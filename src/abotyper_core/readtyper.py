"""Read-level ABO allele typing with per-haplotype consensus (``abotyper-readtyper``).

A complement to the per-position path (pileup/Clair3 -> metrics -> predict-phenotype). Every ONT read of the
exon 2-7 amplicon is one molecule, so instead of summarising reads position by position this scores each read
against each ISBT allele and picks the allele pair that best explains all reads:

1. ISBT alleles (``data/isbt_abo_alleles.tsv``) become expected states at genomic sites on NG_006669.2. cDNA ->
   genomic offsets come per exon from the panel markers, which carry both coordinates; intronic changes
   (``c.374+5G>A``) are anchored on their coding position.
2. Each read's observed state is taken at every site it spans: the base for SNVs, and for indels the net length
   change over the whole homopolymer run (or a small window), so an ONT indel placed anywhere in the run counts.
3. ``log P(read | allele)`` uses a simple error model (SNV 2%, homopolymer indel 20%, other indel 5%), and
   ``log P(reads | a/b) = sum_r log(P(r|a)/2 + P(r|b)/2)`` plus a prior that favours common alleles.
4. Checks on the best pair:
   - allele balance: a second allele on few reads, too few to be sampling noise of a 50:50 heterozygote, is
     called homozygous (contamination or noise);
   - third haplotype: the share of reads that another allele explains much better than either called allele
     (contamination, chimeras, or a wrong pair);
   - haplotype consensus: reads are split by the called pair, a consensus is built per haplotype over the ISBT
     sites and the coding exons, and matched to the ISBT table. A rare allele is kept only if its haplotype's
     consensus carries its defining variants; otherwise the best common allele for that haplotype is used.
     Coding variants not in the ISBT table are reported per haplotype.

Outputs a one-row summary TSV and a per-haplotype TSV. The status column is PASS, REVIEW or NOCALL.
"""

from __future__ import annotations

import argparse
import bisect
import csv
import math
import re
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from abotyper_core.bamreader import read_bam
from abotyper_core.panel import default_panel_path, load_panel

DATA = Path(default_panel_path()).parent
ISBT_TABLE = DATA / "isbt_abo_alleles.tsv"
UNLISTED_TABLE = DATA / "abo_unlisted_haplotypes.tsv"   # recurring haplotypes not fully listed in ISBT

# coding exons of ABO (NM_020469) in cDNA coordinates; exon 1 (c.1-28) is outside the exon 2-7 amplicon
EXON_CDNA = {"Exon 2": (29, 98), "Exon 3": (99, 155), "Exon 4": (156, 203), "Exon 5": (204, 239),
             "Exon 6": (240, 374), "Exon 7": (375, 1065)}
E_SNV, E_HP, E_INDEL = 0.02, 0.20, 0.05
HP_MIN = 3                    # a run of >= 3 identical bases is a homopolymer
WINDOW = 2                    # +-bp around indels outside homopolymers
COMMON = frozenset({"ABO*A1.01", "ABO*A1.02", "ABO*A2.01", "ABO*B.01", "ABO*O.01.01", "ABO*O.01.02",
                    "ABO*O.01.05", "ABO*O.02.01", "ABO*O.02.02"})
RARE_LOG_PRIOR = -3.0         # per rare allele in a pair
MIN_READS = 5                 # fewer informative reads: NOCALL
MAX_CANDIDATES = 40           # alleles considered for the pair (plus the common ones)
REVIEW_READS = 10             # fewer typed reads: REVIEW
MIN_MARGIN = 10.0             # log-likelihood margin over the best pair of another genotype class
MAX_UNEXPLAINED = 0.10        # share of foreign reads (fit neither called allele, not chimeras of the two)
MAX_CHIMERA = 0.12            # chimeric-read share worth flagging (long-amplicon PCR makes some; clean libraries
MIN_CHIMERA_READS = 5         # stay under 0.09), when at least this many reads are chimeric
MIN_AF = 0.20                 # minor allele share below which (and P_het < BALANCE_P) the call is homozygous
SKEW_AF = 0.30                # heterozygote with a minor share below this (and P_het < BALANCE_P): REVIEW, a
BALANCE_P = 1e-3              # contaminant can pose as the second allele (clean heterozygotes stay above 0.34)
CONS_SUPPORT, CONS_MIN_READS = 0.6, 2   # consensus state needs this read share and depth
NOVEL_SUPPORT, NOVEL_MIN_READS = 0.7, 3
# a known (ISBT) site where the haplotype consensus differs from its called allele: reported per haplotype, and a
# REVIEW flag when it is well supported (possible new allele, or a recombinant of known ones)
MISMATCH_SUPPORT, MISMATCH_MIN_READS = 0.8, 5

_CHANGE = re.compile(r"c\.(\d+)([+-]\d+)?(?:_(\d+)([+-]\d+)?)?"
                     r"(?:([ACGT])>([ACGT])|(del|dup)([ACGT]*)|ins([ACGT]+))$")


# --------------------------------------------------------------------------------------------- reference & sites

def read_fasta(path: Path) -> str:
    with open(path) as fh:
        return "".join(line.strip() for line in fh if not line.startswith(">")).upper()


def exon_offsets(panel) -> Dict[str, int]:
    """genomic = cdna + offset, per exon label ('Exon 6'), voted from panel rows with both coordinates."""
    votes: Dict[str, Counter] = defaultdict(Counter)
    for m in panel.variants:
        if m.cdna_pos and m.ng006669_2_pos and m.exon.startswith("Exon"):
            votes[m.exon][m.ng006669_2_pos - m.cdna_pos] += 1
    return {ex: c.most_common(1)[0][0] for ex, c in votes.items()}


def exon_of_cdna(c: int) -> Optional[str]:
    for ex, (lo, hi) in EXON_CDNA.items():
        if lo <= c <= hi:
            return ex
    return None


def genomic_of(c: int, intronic: Optional[str], offsets: Dict[str, int]) -> Optional[int]:
    ex = exon_of_cdna(c)
    if ex is None or ex not in offsets:
        return None
    return c + offsets[ex] + (int(intronic) if intronic else 0)


def coding_map(offsets: Dict[str, int]) -> Dict[int, int]:
    """{genomic position: cDNA position} over the coding exons covered by the offsets."""
    out = {}
    for ex, (lo, hi) in EXON_CDNA.items():
        if ex in offsets:
            for c in range(lo, hi + 1):
                out[c + offsets[ex]] = c
    return out


def run_span(seq: str, g1: int, g2: int) -> Tuple[int, int, Optional[str]]:
    """1-based [g1, g2] extended to the homopolymer run it sits in, if the span is one repeated base."""
    bases = set(seq[g1 - 1:g2])
    if len(bases) != 1:
        return g1, g2, None
    b = bases.pop()
    lo, hi = g1, g2
    while lo > 1 and seq[lo - 2] == b:
        lo -= 1
    while hi < len(seq) and seq[hi] == b:
        hi += 1
    return lo, hi, b if hi - lo + 1 >= HP_MIN else None


@dataclass(frozen=True)
class Site:
    kind: str                 # "snv" or "indel"
    lo: int                   # 1-based genomic span (snv: lo == hi)
    hi: int
    ref: object               # snv: reference base; indel: homopolymer base or None
    err: float

    @property
    def ref_state(self):
        return self.ref if self.kind == "snv" else 0


@dataclass
class AlleleTable:
    names: List[str]
    phenotype: Dict[str, str]
    sites: List[Site]
    expected: List[Dict[int, object]]          # per allele: {site index: expected state}
    skipped: Counter = field(default_factory=Counter)


def load_alleles(path: Path, offsets: Dict[str, int], seq: str) -> AlleleTable:
    site_idx: Dict[tuple, int] = {}
    sites: List[Site] = []
    names, pheno, expected, skipped = [], {}, [], Counter()

    def site(key, s):
        if key not in site_idx:
            site_idx[key] = len(sites)
            sites.append(s)
        return site_idx[key]

    with open(path) as fh:
        rows = list(csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"))
    for r in rows:
        if r.get("obsolete") == "1":
            continue
        state: Dict[int, object] = {}
        for ch in (c.strip() for c in (r.get("changes") or "").split(";")):
            if not ch:
                continue
            m = _CHANGE.match(ch)
            if not m:
                skipped["unparsed"] += 1
                continue
            g1 = genomic_of(int(m.group(1)), m.group(2), offsets)
            g2 = genomic_of(int(m.group(3)), m.group(4), offsets) if m.group(3) else g1
            if g1 is None or g2 is None or not (1 <= g1 <= g2 <= len(seq)):
                skipped["outside amplicon"] += 1
                continue
            if m.group(5):                                           # SNV
                if seq[g1 - 1] != m.group(5):
                    skipped["reference mismatch"] += 1
                    continue
                state[site(("snv", g1), Site("snv", g1, g1, seq[g1 - 1], E_SNV))] = m.group(6)
                continue
            lo, hi, b = run_span(seq, g1, g2)
            delta = len(m.group(9)) if m.group(9) else (-(g2 - g1 + 1) if m.group(7) == "del" else g2 - g1 + 1)
            if b is None:
                lo, hi = lo - WINDOW, hi + WINDOW
            i = site(("indel", lo, hi), Site("indel", lo, hi, b, E_HP if b else E_INDEL))
            state[i] = state.get(i, 0) + delta
        names.append(r["allele"])
        pheno[r["allele"]] = r.get("phenotype", "")
        expected.append(state)
    return AlleleTable(names, pheno, sites, expected, skipped)


# --------------------------------------------------------------------------------------------- read observations

@dataclass
class ReadObs:
    sites: Dict[int, object]          # site index -> observed state
    bases: Dict[int, str]             # genomic position -> base ('-' for deletion), coding/interest positions


def observe_reads(bam: Path, table: AlleleTable, coding: Dict[int, int], ref_name: Optional[str] = None,
                  min_mapq: int = 1) -> Tuple[List[ReadObs], int]:
    """Observed states per primary alignment with MAPQ >= min_mapq; also returns the count of alignments seen."""
    interest = set(coding)
    for s in table.sites:
        interest.update(range(s.lo - 1, s.hi + 1))
    ipos = sorted(interest)
    out, seen = [], 0
    for aln in read_bam(bam, primary_only=True):
        seen += 1
        if aln.mapq < min_mapq or (ref_name and aln.ref != ref_name):
            continue
        base_at: Dict[int, str] = {}
        ins_after: Dict[int, str] = {}
        rp, qp = aln.pos + 1, 0
        for op, n in aln.cigar:
            if op in "M=X":
                for p in ipos[bisect.bisect_left(ipos, rp):bisect.bisect_left(ipos, rp + n)]:
                    base_at[p] = aln.seq[qp + p - rp]
                rp += n
                qp += n
            elif op in "DN":
                for p in ipos[bisect.bisect_left(ipos, rp):bisect.bisect_left(ipos, rp + n)]:
                    base_at[p] = "-"
                rp += n
            elif op == "I":
                if rp - 1 in interest:
                    ins_after[rp - 1] = ins_after.get(rp - 1, "") + aln.seq[qp:qp + n]
                qp += n
            elif op == "S":
                qp += n
        start, end = aln.pos + 1, aln.ref_end
        obs: Dict[int, object] = {}
        for i, s in enumerate(table.sites):
            if s.lo - 1 <= start or s.hi + 1 >= end:
                continue
            if s.kind == "snv":
                obs[i] = base_at.get(s.lo, "-")
            else:
                ins = sum(len(v) if s.ref is None else v.count(s.ref)
                          for p, v in ins_after.items() if s.lo - 1 <= p <= s.hi)
                dels = sum(1 for p in range(s.lo, s.hi + 1) if base_at.get(p) == "-")
                obs[i] = ins - dels
        out.append(ReadObs(obs, {p: b for p, b in base_at.items() if p in coding}))
    return out, seen


# --------------------------------------------------------------------------------------------- likelihoods

def read_likelihoods(reads: Sequence[ReadObs], table: AlleleTable) -> np.ndarray:
    """L[r, a] = log P(read r | allele a)."""
    n_a = len(table.names)
    exp_cols = []
    for i, s in enumerate(table.sites):
        exp_cols.append(np.array([table.expected[a].get(i, s.ref_state) for a in range(n_a)], dtype=object))
    hit = [math.log(1 - s.err) for s in table.sites]
    miss = [math.log(s.err) for s in table.sites]
    cache: Dict[tuple, np.ndarray] = {}
    L = np.zeros((len(reads), n_a))
    for r, ro in enumerate(reads):
        for i, v in ro.sites.items():
            row = cache.get((i, v))
            if row is None:
                row = cache[(i, v)] = np.where(exp_cols[i] == v, hit[i], miss[i])
            L[r] += row
    return L


def allele_class(phenotype: str, name: str = "") -> str:
    p = phenotype or name.replace("ABO*", "")
    if p.startswith("cisAB") or p.startswith("B(A)"):
        return "cisAB" if p.startswith("cisAB") else "B(A)"
    return p[:1] if p[:1] in "ABO" else "?"


def genotype_of(classes: Sequence[str]) -> str:
    simple = all(c in ("A", "B", "O") for c in classes)
    order = {"A": 0, "B": 1, "O": 2}
    cs = sorted(classes, key=lambda c: order.get(c, 3))
    return "".join(cs) if simple else "/".join(cs)


_PHENO_RANK = ["A1", "A2", "A3", "Ax", "Am", "Afinn", "Abantu", "Aweak", "Ael"]


def phenotype_of(p1: str, p2: str) -> str:
    """Expected serological phenotype from two ISBT allele phenotypes (O is recessive)."""
    ps = [p for p in (p1, p2) if p and p != "O"]
    if not ps:
        return "O"
    if any(p.startswith(("cisAB", "B(A)")) for p in ps):
        return next(p for p in ps if p.startswith(("cisAB", "B(A)")))
    a = [p for p in ps if p.startswith("A")]
    b = [p for p in ps if p.startswith("B")]

    def best(group):
        if not group:
            return ""
        def rank(p):
            return next((i for i, k in enumerate(_PHENO_RANK) if p.startswith(k)), len(_PHENO_RANK))
        return min(group, key=lambda p: (rank(p), p))
    pa, pb = best(a), best(b)
    return (pa + pb) if pa and pb else (pa or pb)


# --------------------------------------------------------------------------------------------- haplotypes

@dataclass
class Haplotype:
    allele: str
    reads: int
    consensus: Dict[int, object]       # site index -> consensus state (supported sites only)
    exact: List[str]                   # ISBT alleles matching the consensus at every supported site
    nearest: str
    nearest_mismatches: int
    novel: List[str]
    note: str = ""
    mismatches: List[str] = field(default_factory=list)      # known sites differing from the called allele
    supported_mismatches: int = 0                             # of those, with >= MISMATCH_* support


def site_label(site: "Site", coding: Dict[int, int]) -> str:
    """c. position of a site (g. NG_006669.2 position when it is outside the coding exons)."""
    c = coding.get(site.lo)
    return f"c.{c}" if c is not None else f"g.{site.lo}"


def allele_mismatches(reads: Sequence[ReadObs], idx: Sequence[int], cons: Dict[int, object], allele: int,
                      table: AlleleTable, coding: Dict[int, int]) -> Tuple[List[str], int]:
    """Known sites where the haplotype consensus differs from `allele`, as 'c.220 C>T (12/12 reads)', and how many
    of them reach MISMATCH_SUPPORT and MISMATCH_MIN_READS."""
    out, strong = [], 0
    exp = table.expected[allele]
    for i in sorted(cons, key=lambda i: table.sites[i].lo):
        s = table.sites[i]
        e, v = exp.get(i, s.ref_state), cons[i]
        if e == v:
            continue
        n = sum(1 for r in idx if i in reads[r].sites)
        k = sum(1 for r in idx if reads[r].sites.get(i) == v)
        if s.kind == "snv":
            what = f"{e}>{v}"
        else:
            what = f"indel {int(e):+d}>{int(v):+d} bp" + (f" ({s.ref} run)" if s.ref else "")
        out.append(f"{site_label(s, coding)} {what} ({k}/{n} reads)")
        strong += n >= MISMATCH_MIN_READS and k / n >= MISMATCH_SUPPORT
    return out, strong


def consensus_sites(reads: Sequence[ReadObs], idx: Sequence[int]) -> Dict[int, object]:
    counts: Dict[int, Counter] = defaultdict(Counter)
    for r in idx:
        for i, v in reads[r].sites.items():
            counts[i][v] += 1
    cons = {}
    for i, c in counts.items():
        v, k = c.most_common(1)[0]
        n = sum(c.values())
        if n >= CONS_MIN_READS and k / n >= CONS_SUPPORT:
            cons[i] = v
    return cons


def consensus_mismatches(cons: Dict[int, object], table: AlleleTable) -> Dict[str, int]:
    """{allele: number of consensus sites where the allele expects another state}."""
    return {name: sum(1 for i, v in cons.items() if table.expected[a].get(i, table.sites[i].ref_state) != v)
            for a, name in enumerate(table.names)}


def novel_variants(reads: Sequence[ReadObs], idx: Sequence[int], seq: str, coding: Dict[int, int],
                   table: AlleleTable) -> List[str]:
    """Coding changes carried by most reads of a haplotype that are not ISBT sites (SNVs and non-homopolymer dels)."""
    known = set()
    for s in table.sites:
        known.update(range(s.lo - (1 if s.kind == "indel" else 0), s.hi + 1))
    counts: Dict[int, Counter] = defaultdict(Counter)
    for r in idx:
        for p, b in reads[r].bases.items():
            counts[p][b] += 1
    out = []
    for p in sorted(counts):
        if p in known:
            continue
        b, k = counts[p].most_common(1)[0]
        n = sum(counts[p].values())
        if b == seq[p - 1] or n < NOVEL_MIN_READS or k / n < NOVEL_SUPPORT:
            continue
        c = coding[p]
        if b == "-":
            if run_span(seq, p, p)[2] is None:
                out.append(f"c.{c}del")
        else:
            out.append(f"c.{c}{seq[p - 1]}>{b}")
    return out


# --------------------------------------------------------------------------------------------- typing

@dataclass
class TypingResult:
    summary: Dict[str, object]
    haplotypes: List[Haplotype]


def candidate_alleles(L: np.ndarray, names: Sequence[str]) -> List[int]:
    """Alleles that can be in the best pair: the common ones, and any that fit at least a few reads best."""
    n = L.shape[0]
    best_for = (L >= L.max(axis=1, keepdims=True) - 1e-9).sum(axis=0)
    k = max(2, math.ceil(0.05 * n))                # a real allele carries far more than 5% of the reads
    total = L.sum(axis=0)
    fit = sorted((i for i, c in enumerate(best_for) if c >= k), key=lambda i: (-best_for[i], -total[i]))
    return sorted(set(fit[:MAX_CANDIDATES]) | {i for i, a in enumerate(names) if a in COMMON})


def best_pairs(L: np.ndarray, names: Sequence[str]) -> List[Tuple[float, int, int]]:
    """All pairs of candidate alleles, best first: (log-likelihood + prior, i, j) with i <= j."""
    prior = np.array([0.0 if n in COMMON else RARE_LOG_PRIOR for n in names])
    cand = candidate_alleles(L, names)
    Lc = L[:, cand]
    out = []
    for x, i in enumerate(cand):
        ll = (np.logaddexp(Lc[:, x:x + 1], Lc[:, x:]) + math.log(0.5)).sum(axis=0) + prior[i] + prior[cand[x:]]
        out.extend((float(v), i, cand[x + j]) for j, v in enumerate(ll))
    out.sort(reverse=True)
    return out


def is_chimera(read: ReadObs, table: AlleleTable, a: int, b: int, max_err: int = 1) -> bool:
    """Does the read follow allele a on one side of a breakpoint and allele b on the other (or a-b-a)?

    Along the sites where a and b differ, each observed state is labelled a, b or neither. A chimera has both labels,
    at most two switches between them, and at most `max_err` sites matching neither (sequencing errors)."""
    ea, eb = table.expected[a], table.expected[b]
    labels, err = [], 0
    for i in sorted(read.sites, key=lambda i: table.sites[i].lo):
        s = table.sites[i]
        va, vb, v = ea.get(i, s.ref_state), eb.get(i, s.ref_state), read.sites[i]
        if va == vb:
            err += v != va and s.kind == "snv"          # disagreeing with both at a shared SNV: not a chimera
            continue
        if v == va:
            labels.append("a")
        elif v == vb:
            labels.append("b")
        else:
            err += 1
    switches = sum(1 for x, y in zip(labels, labels[1:]) if x != y)
    return "a" in labels and "b" in labels and switches <= 2 and err <= max_err


def binom_tail(k: int, n: int) -> float:
    """P(X <= k) for X ~ Binomial(n, 0.5), computed in log space."""
    if n == 0:
        return 1.0
    logs = [math.lgamma(n + 1) - math.lgamma(i + 1) - math.lgamma(n - i + 1) - n * math.log(2) for i in range(k + 1)]
    m = max(logs)
    return math.exp(m) * sum(math.exp(x - m) for x in logs)


def load_unlisted(path: Optional[Path] = None) -> Dict[tuple, Tuple[str, str, str]]:
    """{(called allele, frozenset of differing sites): (report_as, phenotype, evidence)} from the table of recurring
    haplotypes that differ from their nearest ISBT allele (see data/abo_unlisted_haplotypes.tsv)."""
    path = UNLISTED_TABLE if path is None else path
    out: Dict[tuple, Tuple[str, str, str]] = {}
    if not path or not Path(path).is_file():
        return out
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t"):
            sites = frozenset(x.strip() for x in r["differs_at"].split(",") if x.strip())
            out[(r["called_allele"].strip(), sites)] = (r["report_as"].strip(), r["phenotype"].strip(),
                                                          r.get("evidence", "").strip())
    return out


_UNLISTED_DEFAULT: Optional[Dict[tuple, Tuple[str, str, str]]] = None


def _default_unlisted():
    global _UNLISTED_DEFAULT
    if _UNLISTED_DEFAULT is None:
        _UNLISTED_DEFAULT = load_unlisted()
    return _UNLISTED_DEFAULT


def type_sample(reads: Sequence[ReadObs], table: AlleleTable, seq: str, coding: Dict[int, int],
                min_af: float = MIN_AF, unlisted: Optional[Dict[tuple, Tuple[str, str, str]]] = None) -> TypingResult:
    names = table.names
    L = read_likelihoods(reads, table)
    # only reads that tell alleles apart are typed: a read covering no distinguishing site fits every allele equally,
    # adds nothing to the call and would dilute the allele balance and unexplained-read shares
    informative = np.ptp(L, axis=1) > 1e-9 if len(reads) else np.zeros(0, dtype=bool)
    reads = [r for r, k in zip(reads, informative) if k]
    L = L[informative]
    if len(reads) < MIN_READS:
        return TypingResult({"reads_typed": len(reads), "genotype": "NOCALL", "status": "NOCALL",
                             "flags": f"only {len(reads)} informative reads"}, [])
    pairs = best_pairs(L, names)
    top = pairs[0]
    def cls(i):
        return allele_class(table.phenotype[names[i]], names[i])
    def gclass(p):
        return genotype_of([cls(p[1]), cls(p[2])])
    ties = [p for p in pairs if top[0] - p[0] < 1e-6]
    alt = next((p for p in pairs if gclass(p) != gclass(top)), None)
    margin = top[0] - alt[0] if alt else float("inf")
    ia, ib = top[1], top[2]
    flags = []

    # allele balance
    to_a = int(np.sum(L[:, ia] > L[:, ib]))
    to_b = int(np.sum(L[:, ib] > L[:, ia]))
    k, n_ab = min(to_a, to_b), to_a + to_b
    minor = k / n_ab if ia != ib and n_ab else 0.0
    if ia != ib and minor < min_af and binom_tail(k, n_ab) < BALANCE_P:
        major = ia if to_a >= to_b else ib
        flags.append(f"minor allele {names[ib if major == ia else ia]} on {minor:.0%} of reads: called homozygous")
        ia = ib = major
    elif ia != ib and minor < SKEW_AF and binom_tail(k, n_ab) < BALANCE_P:
        flags.append(f"allele imbalance: minor allele on {minor:.0%} of reads (a heterozygote gives ~50%): "
                     "possible contamination")

    # third haplotype: reads another allele explains much better than either called allele are PCR chimeras of the
    # two (they switch from one allele to the other along the amplicon) or foreign (contamination)
    gap = L.max(axis=1) - np.maximum(L[:, ia], L[:, ib])
    unexpl = gap > 2 * math.log((1 - E_SNV) / E_SNV)
    chimeric = np.zeros(len(reads), dtype=bool)
    if top[1] != top[2]:
        for r in np.flatnonzero(unexpl):
            chimeric[r] = is_chimera(reads[r], table, top[1], top[2])
    foreign = unexpl & ~chimeric
    third = Counter(names[j] for j in np.argmax(L[foreign], axis=1)) if foreign.any() else Counter()
    unexpl_share = float(foreign.mean())
    chimera_share = float(chimeric.mean())

    # haplotype consensus, rare-allele confirmation, novel variants
    if ia == ib:
        groups = [(ia, list(range(len(reads))))]
    else:
        groups = [(ia, [r for r in range(len(reads)) if L[r, ia] > L[r, ib]]),
                  (ib, [r for r in range(len(reads)) if L[r, ib] > L[r, ia]])]
    haps, final = [], []
    unlisted = _default_unlisted() if unlisted is None else unlisted
    named_pheno: Dict[str, str] = {}
    notes: List[str] = []
    for a, idx in groups:
        cons = consensus_sites(reads, idx)
        mism = consensus_mismatches(cons, table)
        nm = min(mism.values())
        exact = [n for n in names if mism[n] == 0]
        nearest = min(names, key=lambda n: (mism[n], n not in COMMON, n))
        name, note = names[a], ""
        if name not in COMMON:
            commons = [n for n in COMMON if n in mism]
            best_common = min(commons, key=lambda n: (mism[n], -L[idx, names.index(n)].sum()))
            # keep a rare allele only if the haplotype consensus fits it strictly better than every common allele
            if mism[name] >= mism[best_common]:
                note = f"rare allele {name} not confirmed by haplotype consensus: {best_common} used"
                flags.append(note)
                name = best_common
        mm, strong = allele_mismatches(reads, idx, cons, names.index(name), table, coding)
        # a recurring haplotype that differs from its ISBT allele at exactly the listed sites: reported under its own
        # name, with a note instead of a REVIEW flag (any other difference still goes to REVIEW)
        differs = frozenset(m.split(" (")[0] for m in mm)
        hit = unlisted.get((name, differs)) if mm else None
        if hit:
            report_as, pheno, _evidence = hit
            msg = (f"{report_as}: {name} + {', '.join(sorted(differs))}, a recurring haplotype that ISBT does not "
                   "list (see abo_unlisted_haplotypes.tsv)")
            notes.append(msg)
            note = (note + "; " if note else "") + msg
            named_pheno[report_as] = pheno
            name, strong = report_as, 0
        haps.append(Haplotype(name, len(idx), cons, exact, nearest, nm,
                              novel_variants(reads, idx, seq, coding, table), note, mm, strong))
        final.append(name)
    if len(final) == 1:
        final = final * 2
        haps[0].reads = len(reads)
    a1, a2 = final
    pheno_of = {**table.phenotype, **named_pheno}
    p1, p2 = pheno_of.get(a1, ""), pheno_of.get(a2, "")
    genotype = genotype_of([allele_class(p1, a1), allele_class(p2, a2)])
    classes_tied = {gclass(p) for p in ties}
    if len(classes_tied) > 1:
        flags.append("tied pairs of different genotype classes: " + ",".join(sorted(classes_tied)))
    novel = sorted({v for h in haps for v in h.novel})
    off_allele = [f"{h.allele} at {', '.join(m.split(' (')[0] for m in h.mismatches)}"
                  for h in haps if h.supported_mismatches]

    status = "PASS"
    if len(reads) < REVIEW_READS:
        flags.append(f"only {len(reads)} informative reads")
    if margin < MIN_MARGIN:
        flags.append(f"margin {margin:.1f} < {MIN_MARGIN:g}")
    if unexpl_share >= MAX_UNEXPLAINED:
        flags.append(f"{unexpl_share:.0%} of reads fit neither allele (best other: "
                     f"{third.most_common(1)[0][0] if third else '?'}): possible contamination")
    if chimera_share >= MAX_CHIMERA and int(chimeric.sum()) >= MIN_CHIMERA_READS:
        flags.append(f"{chimera_share:.0%} of reads are chimeras of the two alleles")
    if novel:
        flags.append("novel variants: " + ",".join(novel))
    if off_allele:
        flags.append("haplotype differs from its ISBT allele at known site(s): " + "; ".join(off_allele)
                     + " (possible new allele)")
    if flags:
        status = "REVIEW"
    summary = {
        "reads_typed": len(reads), "allele_1": a1, "allele_2": a2, "genotype": genotype,
        "phenotype": phenotype_of(p1, p2),
        "status": status, "loglik": round(top[0], 2), "margin": round(margin, 2) if alt else "",
        "alt_pair": f"{names[alt[1]]}/{names[alt[2]]}" if alt else "", "tied_pairs": len(ties),
        "reads_a_b": f"{to_a}/{to_b}", "minor_allele_share": round(minor, 3),
        "unexplained_share": round(unexpl_share, 3), "chimera_share": round(chimera_share, 3),
        "third_allele": third.most_common(1)[0][0] if third else "",
        "novel_variants": ";".join(novel),
        "haplotype_mismatches": "; ".join(f"{h.allele}: {', '.join(h.mismatches)}" for h in haps if h.mismatches),
        "flags": "; ".join(flags),
        "notes": "; ".join(notes),
    }
    return TypingResult(summary, haps)


SUMMARY_FIELDS = ["sample", "reads_aligned", "reads_typed", "allele_1", "allele_2", "genotype", "phenotype",
                  "status", "loglik", "margin", "alt_pair", "tied_pairs", "reads_a_b", "minor_allele_share",
                  "unexplained_share", "chimera_share", "third_allele", "novel_variants", "haplotype_mismatches",
                  "flags", "notes"]
HAPLOTYPE_FIELDS = ["sample", "haplotype", "allele", "reads", "isbt_exact_matches", "nearest_isbt",
                    "nearest_mismatches", "sites_called", "allele_mismatches", "novel_variants", "note"]


def run(bam: Path, ref: Path, sample: str, panel: Optional[str] = None, isbt: Optional[Path] = None,
        min_af: float = MIN_AF, unlisted: Optional[str] = None) -> TypingResult:
    seq = read_fasta(ref)
    offsets = exon_offsets(load_panel(panel or default_panel_path()))
    table = load_alleles(isbt or ISBT_TABLE, offsets, seq)
    coding = coding_map(offsets)
    reads, seen = observe_reads(bam, table, coding)
    base = {"sample": sample, "reads_aligned": seen, "reads_typed": 0}
    named = None if unlisted is None else ({} if unlisted == "none" else load_unlisted(Path(unlisted)))
    res = type_sample(reads, table, seq, coding, min_af, named)
    res.summary = {k: {**base, **res.summary}.get(k, "") for k in SUMMARY_FIELDS}
    return res


def write_outputs(res: TypingResult, sample: str, out: Path, haplotypes_out: Optional[Path]) -> None:
    fh = sys.stdout if str(out) == "-" else open(out, "w", newline="")
    try:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerow(res.summary)
    finally:
        if fh is not sys.stdout:
            fh.close()
    if haplotypes_out:
        with open(haplotypes_out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=HAPLOTYPE_FIELDS, delimiter="\t")
            w.writeheader()
            for i, h in enumerate(res.haplotypes, 1):
                w.writerow({"sample": sample, "haplotype": i, "allele": h.allele, "reads": h.reads,
                            "isbt_exact_matches": ",".join(h.exact[:20]) + (",..." if len(h.exact) > 20 else ""),
                            "nearest_isbt": h.nearest, "nearest_mismatches": h.nearest_mismatches,
                            "sites_called": len(h.consensus), "allele_mismatches": "; ".join(h.mismatches),
                            "novel_variants": ";".join(h.novel),
                            "note": h.note})


def main(argv: Optional[Sequence[str]] = None) -> int:
    ap = argparse.ArgumentParser(prog="abotyper-readtyper", description=__doc__.split("\n")[0])
    ap.add_argument("--bam", required=True, help="reads aligned to NG_006669.2 (BAM)")
    ap.add_argument("--ref", required=True, help="NG_006669.2 FASTA the reads were aligned to")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True, help="summary TSV ('-' for standard output)")
    ap.add_argument("--haplotypes-out", help="per-haplotype TSV")
    ap.add_argument("--panel", default=None, help="variant panel (default: bundled)")
    ap.add_argument("--isbt", default=None, help="ISBT allele table (default: bundled)")
    ap.add_argument("--min-af", type=float, default=MIN_AF,
                    help=f"minor allele read share below which a skewed het is called homozygous (default {MIN_AF})")
    ap.add_argument("--unlisted", default=None,
                    help="table of recurring haplotypes not listed in ISBT, reported under their own name "
                         "(default: bundled abo_unlisted_haplotypes.tsv; 'none' to report them all as REVIEW)")
    a = ap.parse_args(argv)
    res = run(Path(a.bam), Path(a.ref), a.sample, a.panel, Path(a.isbt) if a.isbt else None, a.min_af,
              a.unlisted)
    write_outputs(res, a.sample, Path(a.out), Path(a.haplotypes_out) if a.haplotypes_out else None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
