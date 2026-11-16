"""Check recurring 'haplotype differs from its ISBT allele' patterns against an independent caller.

usage: python check_recurring_haplotypes.py ALL_READTYPING.tsv AUDIT_PER_LIBRARY.tsv NG_006669v2.fasta [MIN_LIBRARIES]

For every (allele, differing sites) pattern seen in at least MIN_LIBRARIES libraries (default 2): read support of each
site, the sequence context (homopolymer run length), and whether the per-position Clair3 calls of the same libraries
(audit column `vcf_variants`, panel row ids) agree. A site where the haplotype carries a change its allele lacks
("present") should be called by Clair3; a site where it lacks a change its allele has ("absent") should not be,
unless the other allele of the library carries that change (then it cannot be tested from unphased calls).
Prints aggregate numbers only.
"""

import csv
import os
import re
import sys
from collections import defaultdict

from abotyper_core import readtyper as rt
from abotyper_core.panel import default_panel_path, load_panel

rt_path, audit_path, ref = sys.argv[1:4]
min_lib = int(sys.argv[4]) if len(sys.argv) > 4 else 2
rows = {r["sample"]: r for r in csv.DictReader(open(rt_path, encoding="utf-8"), delimiter="\t") if r["sample"] != "sample"}
au = {os.path.basename(r["sample_dir"].rstrip("/\\")): r for r in csv.DictReader(open(audit_path, encoding="utf-8"), delimiter="\t")}
panel = load_panel(default_panel_path())
marker_at = defaultdict(list)                   # cDNA position -> panel row ids
for m in panel.variants:
    pos = getattr(m, "cdna_pos", None) or (m.raw.get("cdna_pos") if hasattr(m, "raw") else None)
    if pos:
        marker_at[int(pos)].append(m.id)
seq = rt.read_fasta(ref)
offsets = rt.exon_offsets(panel)
coding = rt.coding_map(offsets)
table = rt.load_alleles(rt.ISBT_TABLE, offsets, seq)
gpos = {c: g for g, c in coding.items()}


def expects(allele, cpos):
    """Does the ISBT allele carry a change at this cDNA position?"""
    if allele not in table.names:
        return False
    e = table.expected[table.names.index(allele)]
    return any(coding.get(table.sites[i].lo) == cpos for i in e)


def run_len(cpos):
    g = gpos.get(cpos)
    if not g:
        return 0
    b, lo, hi = seq[g - 1], g, g
    while seq[lo - 2] == b:
        lo -= 1
    while seq[hi] == b:
        hi += 1
    return hi - lo + 1


pat = defaultdict(list)
for s, r in rows.items():
    for part in filter(None, (r.get("haplotype_mismatches") or "").split("; ")):
        allele, ms = part.split(": ", 1)
        sites = tuple(sorted(m.split(" (")[0] for m in ms.split(", ")))
        pat[(allele, sites)].append((s, ms))

for (allele, sites), libs in sorted(pat.items(), key=lambda x: -len(x[1])):
    if len(libs) < min_lib:
        continue
    print(f"\n== {allele} + {', '.join(sites)}: {len(libs)} libraries "
          f"(status: {', '.join(f'{k} {v}' for k, v in sorted(__import__('collections').Counter(rows[s]['status'] for s, _ in libs).items()))})")
    others = __import__("collections").Counter(
        rows[s]["allele_2"] if rows[s]["allele_1"] in (allele, "ABO*O.01.12-like") and rows[s]["allele_1"] != rows[s]["allele_2"]
        else rows[s]["allele_1"] for s, _ in libs)
    print(f"   other allele: {dict(others.most_common(5))}")
    for site in sites:
        m = re.match(r"[cg]\.(\d+) (.*)", site)
        cpos, change = int(m.group(1)), m.group(2)
        sup = []
        for s, ms in libs:
            mm = re.search(re.escape(site) + r"(?: \([ACGT] run\))? \((\d+)/(\d+) reads\)", ms)
            if mm:
                sup.append(int(mm.group(1)) / int(mm.group(2)))
        sup.sort()
        if not sup:
            sup = [float("nan")]
        present = not expects(allele, cpos)          # haplotype has a change its allele lacks
        ids = marker_at.get(cpos, [])
        testable = agree = 0
        for s, _ in libs:
            r = rows[s]
            other = r["allele_2"] if r["allele_1"] == allele else r["allele_1"]
            if s not in au or not ids:
                continue
            called = any(i in (au[s]["vcf_variants"] or "").split(";") for i in ids)
            if present:
                testable += 1
                agree += called
            elif not expects(other, cpos):
                testable += 1
                agree += not called
        what = "change present on haplotype" if present else "allele's change absent on haplotype"
        print(f"   {site}: {what}; read support median {sup[len(sup)//2]:.2f} (min {sup[0]:.2f}); "
              f"run length {run_len(cpos)}; panel row {','.join(ids) or '-'}; "
              f"Clair3 agrees in {agree}/{testable} testable libraries")
