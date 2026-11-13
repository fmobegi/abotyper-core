"""Experimental read-level ABO allele typing (prototype, not part of the package API).

Each ONT read of the exon 3-7 amplicon is one molecule, so instead of summarising the reads position by position
(pileup or caller) and then combining per-position calls, this scores every read against every ISBT allele and picks
the allele pair that best explains all reads:

  1. ISBT alleles (data/isbt_abo_alleles.tsv) are turned into expected states at genomic sites on NG_006669.2.
     cDNA -> genomic offsets are derived per exon from the panel markers (which carry both coordinates).
  2. For each read, the observed state is read at every site it spans: the base for SNVs, and for indels the net
     length change over the whole homopolymer run (or a small window), so an ONT indel placed anywhere in the run
     still counts.
  3. log P(read | allele) uses a simple error model (SNV e=0.02, homopolymer indel e=0.20, other indel e=0.05).
  4. log P(reads | pair a/b) = sum_r log(0.5 P(r|a) + 0.5 P(r|b)) (+ a small prior favouring common alleles).
     The best pair, its margin over the best pair with a different genotype class, and tied/indistinguishable pairs
     are reported.

Usage: python readtyper.py --bam S.bam --ref NG_006669v2.fasta [--panel ...] [--isbt ...] --sample S --out S.tsv
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pysam

from abotyper_core.panel import default_panel_path, load_panel

DATA = Path(default_panel_path()).parent
E_SNV, E_HP, E_INDEL = 0.02, 0.20, 0.05
HP_MIN = 3            # a run of >= 3 identical bases is treated as a homopolymer
WINDOW = 2            # +-bp for indels outside homopolymers
COMMON = {"ABO*A1.01", "ABO*A1.02", "ABO*A2.01", "ABO*B.01", "ABO*O.01.01", "ABO*O.01.02", "ABO*O.02.01",
          "ABO*O.02.02", "ABO*O.01.05"}
RARE_LOG_PRIOR = -3.0  # per rare allele in a pair


def read_fasta(path: Path) -> str:
    seq = []
    with open(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq).upper()


def exon_offsets(panel) -> dict:
    """genomic = cdna + offset, per exon label ('Exon 6'), from panel rows with both coordinates."""
    votes = defaultdict(Counter)
    for m in panel.variants:
        if m.cdna_pos and m.ng006669_2_pos and m.exon.startswith("Exon"):
            votes[m.exon][m.ng006669_2_pos - m.cdna_pos] += 1
    return {ex: c.most_common(1)[0][0] for ex, c in votes.items()}


CHANGE = re.compile(r"c\.(\d+)(?:_(\d+))?(?:([ACGT])>([ACGT])|(del|dup)([ACGT]*)|ins([ACGT]+))$")


def run_span(seq: str, g1: int, g2: int):
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


def load_alleles(path: Path, offsets: dict, seq: str):
    """{allele: {site_key: expected}}, site meta {key: (kind, lo, hi, base, err)}, skipped change counts."""
    alleles, sites, skipped, refcheck = {}, {}, Counter(), Counter()
    with open(path) as fh:
        rows = [r for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")]
    for r in rows:
        if r.get("obsolete") == "1":
            continue
        changes = [c.strip() for c in (r["changes"] or "").split(";") if c.strip()]
        exons = [e.strip() for e in (r["exons"] or "").split(";")]
        state = {}
        for i, ch in enumerate(changes):
            ex = exons[i] if i < len(exons) else ""
            m = CHANGE.match(ch)
            if not m or ex not in offsets:
                skipped["unparsed" if not m else f"no offset ({ex or 'none'})"] += 1
                continue
            off = offsets[ex]
            c1 = int(m.group(1))
            c2 = int(m.group(2)) if m.group(2) else c1
            g1, g2 = c1 + off, c2 + off
            if m.group(3):                                   # SNV
                ok = seq[g1 - 1] == m.group(3)
                refcheck["ok" if ok else "mismatch"] += 1
                key = ("snv", g1)
                sites[key] = ("snv", g1, g1, seq[g1 - 1], E_SNV)
                state[key] = m.group(4)
                continue
            if m.group(7):                                   # insertion between c1 and c2
                lo, hi, b = run_span(seq, g1, g2)
                delta = len(m.group(7))
            else:
                lo, hi, b = run_span(seq, g1, g2)
                delta = -(g2 - g1 + 1) if m.group(5) == "del" else (g2 - g1 + 1)
            if b is None:
                lo, hi = lo - WINDOW, hi + WINDOW
            key = ("indel", lo, hi)
            sites[key] = ("indel", lo, hi, b, E_HP if b else E_INDEL)
            state[key] = state.get(key, 0) + delta
        alleles[r["allele"]] = state
    return alleles, sites, skipped, refcheck


def read_observations(bam_path: Path, sites: dict, min_mapq: int = 1, max_reads: int = 0):
    """List of {site_key: observed} per read (primary alignments only)."""
    keys = list(sites)
    obs = []
    with pysam.AlignmentFile(str(bam_path)) as bam:
        for read in bam.fetch(until_eof=True):
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.mapping_quality < min_mapq:
                continue
            seq = read.query_sequence
            base_at, ins_at, last = {}, defaultdict(str), None
            for q, rp in read.get_aligned_pairs():
                if rp is not None:
                    base_at[rp + 1] = seq[q] if q is not None else "-"
                    last = rp + 1
                elif q is not None and last is not None:
                    ins_at[last] += seq[q]
            start, end = read.reference_start + 1, read.reference_end      # 1-based inclusive span
            o = {}
            for k in keys:
                kind, lo, hi, b, _ = sites[k]
                if lo - 1 <= start or hi + 1 >= end:
                    continue
                if kind == "snv":
                    o[k] = base_at.get(lo, "-")
                else:
                    ins = sum(len(s) if b is None else s.count(b) for p, s in ins_at.items() if lo - 1 <= p <= hi)
                    dels = sum(1 for p in range(lo, hi + 1) if base_at.get(p) == "-")
                    o[k] = ins - dels
            obs.append(o)
            if max_reads and len(obs) >= max_reads:
                break
    return obs


def allele_class(name: str) -> str:
    n = name.replace("ABO*", "").upper()
    if n.startswith("CISAB") or n.startswith("B(A)") or n.startswith("BA"):
        return "AB"
    return {"O": "O", "A": "A", "B": "B"}.get(n[:1], "?")


def o_subtype(name: str) -> str:
    m = re.match(r"ABO\*O\.(\d+)", name)
    return f"O{int(m.group(1))}" if m else ""


def genotype_class(a: str, b: str) -> str:
    letters = "".join(sorted((allele_class(a) + allele_class(b)).replace("?", ""), key="ABO".index))
    return letters


def type_reads(obs, alleles, sites, min_af: float = 0.2):
    names = sorted(alleles)
    keys = list(sites)
    kidx = {k: i for i, k in enumerate(keys)}
    ref_state = [sites[k][3] if sites[k][0] == "snv" else 0 for k in keys]
    cols = [np.array([alleles[a].get(k, ref_state[i]) for a in names], dtype=object) for i, k in enumerate(keys)]
    hit = [math.log(1 - sites[k][4]) for k in keys]
    miss = [math.log(sites[k][4]) for k in keys]
    L = np.zeros((len(obs), len(names)))
    for r, o in enumerate(obs):
        for k, v in o.items():
            i = kidx[k]
            L[r] += np.where(cols[i] == v, hit[i], miss[i])
    prior = np.array([0.0 if a in COMMON else RARE_LOG_PRIOR for a in names])
    n = len(names)
    best = []
    for i in range(n):
        pair = np.logaddexp(L[:, i:i + 1], L[:, i:]) + math.log(0.5)
        ll = pair.sum(axis=0) + prior[i] + prior[i:]
        for jj, v in enumerate(ll):
            best.append((float(v), names[i], names[i + jj]))
    best.sort(reverse=True)
    top = best[0]
    ties = [b for b in best if top[0] - b[0] < 1e-6]
    classes = {genotype_class(a, b) for _, a, b in ties}
    top_cls = genotype_class(top[1], top[2])
    alt_cls = next((b for b in best if genotype_class(b[1], b[2]) != top_cls), None)
    margin = top[0] - alt_cls[0] if alt_cls else float("inf")
    ia, ib = names.index(top[1]), names.index(top[2])
    to_a = int(np.sum(L[:, ia] > L[:, ib]))
    to_b = int(np.sum(L[:, ib] > L[:, ia]))
    # third-haplotype check: reads that some other allele explains much better (>= 2 SNV mismatches' worth) than
    # either allele of the call point to contamination, a chimera, or a wrong pair
    gap = L.max(axis=1) - np.maximum(L[:, ia], L[:, ib])
    unexpl = gap > 2 * math.log((1 - E_SNV) / E_SNV)
    third = Counter(names[j] for j in np.argmax(L[unexpl], axis=1)) if unexpl.any() else Counter()
    a, b = top[1], top[2]
    minor = min(to_a, to_b) / max(1, to_a + to_b) if a != b else 0.0
    note = ""
    k, n_ab = min(to_a, to_b), to_a + to_b
    p_het = sum(math.comb(n_ab, i) for i in range(k + 1)) / 2 ** n_ab if n_ab else 1.0   # P(<= k | het, 50:50)
    if a != b and minor < min_af and p_het < 1e-3:
        # allele balance: a second allele on few reads, and too few to be sampling noise of a 50:50 heterozygote,
        # is contamination or noise (at low depth a skewed split is normal, so the het call stands)
        major = a if to_a >= to_b else b
        note = f"minor allele {b if major == a else a} on {minor:.0%} of reads: called homozygous"
        a = b = major
        classes = {genotype_class(a, b)}
    return {
        "best_pair": f"{a}/{b}",
        "genotype": genotype_class(a, b) if len(classes) == 1 else "/".join(sorted(classes)),
        "o_subtypes": "/".join(s for s in (o_subtype(a), o_subtype(b)) if s),
        "minor_allele_share": round(minor, 3),
        "note": note,
        "loglik": round(top[0], 2),
        "margin_vs_other_class": round(margin, 2),
        "alt_class_pair": f"{alt_cls[1]}/{alt_cls[2]}" if alt_cls else "",
        "tied_pairs": len(ties),
        "reads_a_b": f"{to_a}/{to_b}",
        "unexplained_share": round(float(unexpl.mean()), 3),
        "third_allele": third.most_common(1)[0][0] if third else "",
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--bam", required=True)
    ap.add_argument("--ref", required=True)
    ap.add_argument("--panel", default=str(default_panel_path()))
    ap.add_argument("--isbt", default=str(DATA / "isbt_abo_alleles.tsv"))
    ap.add_argument("--sample", default="")
    ap.add_argument("--out", default="-")
    ap.add_argument("--min-af", type=float, default=0.2, help="minimum read share for a second allele (default 0.2)")
    ap.add_argument("--report-sites", action="store_true", help="print the site table to stderr and exit")
    args = ap.parse_args()

    seq = read_fasta(Path(args.ref))
    offsets = exon_offsets(load_panel(args.panel))
    alleles, sites, skipped, refcheck = load_alleles(Path(args.isbt), offsets, seq)
    if args.report_sites:
        print(f"offsets={offsets}\nalleles={len(alleles)} sites={len(sites)} skipped={dict(skipped)} "
              f"snv_refcheck={dict(refcheck)}", file=sys.stderr)
        return 0
    obs = read_observations(Path(args.bam), sites)
    res = {"sample": args.sample, "reads": len(obs)}
    res.update(type_reads(obs, alleles, sites, args.min_af) if obs else {"genotype": "NOCALL"})
    out = sys.stdout if args.out == "-" else open(args.out, "w")
    w = csv.DictWriter(out, fieldnames=list(res), delimiter="\t")
    w.writeheader()
    w.writerow(res)
    return 0


if __name__ == "__main__":
    sys.exit(main())
