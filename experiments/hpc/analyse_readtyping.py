"""Aggregate summary of an all_readtyping.tsv (+ haplotypes) upload. Prints no sample names."""
import csv
import re
import statistics as st
import sys
from collections import Counter

summ, hap = sys.argv[1], sys.argv[2]
rows = [r for r in csv.DictReader(open(summ, encoding="utf-8"), delimiter="\t") if r["sample"] != "sample"]
haps = [r for r in csv.DictReader(open(hap, encoding="utf-8"), delimiter="\t") if r["sample"] != "sample"]
n = len(rows)
print(f"samples: {n}  (unique names {len({r['sample'] for r in rows})}); haplotype rows: {len(haps)} "
      f"from {len({h['sample'] for h in haps})} samples")
st_c = Counter(r["status"] for r in rows)
print("status:", dict(st_c))


def f(x, d=0.0):
    try:
        return float(x)
    except (TypeError, ValueError):
        return d


ra = [int(f(r["reads_aligned"])) for r in rows]
rt = [int(f(r["reads_typed"])) for r in rows]
q = lambda v: (f"min {min(v)} p5 {sorted(v)[len(v)//20]} p25 {sorted(v)[len(v)//4]} median {st.median(v)} "
               f"p75 {sorted(v)[3*len(v)//4]} p95 {sorted(v)[19*len(v)//20]} max {max(v)}") if v else "-"
print("reads aligned:", q(ra))
print("informative reads:", q(rt))
bins = [(0, 5), (5, 10), (10, 20), (20, 50), (50, 200), (200, 1000), (1000, 10**9)]
print("status by informative reads:")
for lo, hi in bins:
    sub = [r for r in rows if lo <= int(f(r["reads_typed"])) < hi]
    if sub:
        c = Counter(r["status"] for r in sub)
        print(f"  {lo:>5}-{'' if hi > 10**8 else hi - 1:<5} n={len(sub):5d}  " + ", ".join(f"{k} {v}" for k, v in c.most_common()))
called = [r for r in rows if r["status"] in ("PASS", "REVIEW")]
print("genotype:", dict(Counter(r["genotype"] for r in called).most_common()))
print("phenotype:", dict(Counter(r["phenotype"] for r in called).most_common()))
common = re.compile(r"ABO\*(A1\.01|A1\.02|A2\.01|B\.01|O\.01\.01|O\.01\.02|O\.02\.01|O\.02\.02)$")
alle = Counter(a for r in called for a in (r["allele_1"], r["allele_2"]))
print("alleles (top 25):", dict(alle.most_common(25)))
rare = Counter(a for a in alle.elements() if not common.match(a))
print(f"rare-allele copies: {sum(rare.values())} in {len({r['sample'] for r in called if not (common.match(r['allele_1']) and common.match(r['allele_2']))})} samples;",
      "top:", dict(rare.most_common(20)))
fl = Counter()
for r in rows:
    for x in filter(None, (r["flags"] or "").split("; ")):
        fl[re.sub(r"ABO\*\S+|c\.\S+|-?[\d.]+%?", "#", x)[:80]] += 1
print("flags:")
for k, v in fl.most_common(20):
    print(f"  {v:6d}  {k}")
het = [r for r in called if r["allele_1"] != r["allele_2"] and int(f(r["reads_typed"])) >= 50]
qs = lambda v, ps=(0.001, 0.01, 0.05, 0.5, 0.95, 0.99): " ".join(f"p{p*100:g}={sorted(v)[min(len(v)-1, int(p*len(v)))]:.3f}" for p in ps) if v else "-"
print(f"minor allele share, hets >=50 reads (n={len(het)}):", qs([f(r['minor_allele_share']) for r in het]))
big = [r for r in called if int(f(r["reads_typed"])) >= 50]
print(f"foreign-read share, >=50 reads (n={len(big)}):", qs([f(r['unexplained_share']) for r in big]))
print(f"chimera share, >=50 reads:", qs([f(r['chimera_share']) for r in big]))
print(f"margin (called):", qs([f(r['margin']) for r in called if r['margin'] not in ('', 'inf')]))
for thr in (0.30, 0.34, 0.38):
    print(f"  hets >=50 reads with minor share < {thr}: {sum(f(r['minor_allele_share']) < thr for r in het)}")
for thr in (0.05, 0.10, 0.15, 0.20):
    print(f"  >=50 reads with foreign share >= {thr}: {sum(f(r['unexplained_share']) >= thr for r in big)}")
for thr in (0.08, 0.12, 0.20):
    print(f"  >=50 reads with chimera share >= {thr}: {sum(f(r['chimera_share']) >= thr for r in big)}")
print("third alleles (foreign >=10%):", dict(Counter(r["third_allele"] for r in rows if f(r["unexplained_share"]) >= 0.10).most_common(10)))
nov = Counter()
for r in rows:
    for v in filter(None, (r["novel_variants"] or "").split(";")):
        nov[v] += 1
print(f"samples with novel variants: {sum(1 for r in rows if r['novel_variants'])}; distinct variants: {len(nov)}")
print("  most frequent:", dict(nov.most_common(15)))
nh = [h for h in haps if h["novel_variants"]]
print("  on haplotype alleles:", dict(Counter(h["allele"] for h in nh).most_common(10)))
mm = Counter(int(f(h["nearest_mismatches"])) for h in haps)
print("haplotype nearest-ISBT mismatches:", dict(sorted(mm.items())))
notes = Counter(re.sub(r"ABO\*\S+", "#", h["note"]) for h in haps if h["note"])
print("haplotype notes:", dict(notes.most_common(5)))
ties = Counter(min(int(f(r["tied_pairs"])), 5) for r in called)
print("tied pairs (capped 5):", dict(sorted(ties.items())))
# replicates by base name (strip run prefixes like R1_ / trailing _barcodeNN)
base = lambda s: re.sub(r"_barcode\d+$", "", s)
grp = {}
for r in called:
    grp.setdefault(base(r["sample"]), []).append(r)
reps = {k: v for k, v in grp.items() if len(v) > 1}
disc = sum(1 for v in reps.values() if len({x["genotype"] for x in v}) > 1)
pdisc = sum(1 for v in reps.values() if len({x["genotype"] for x in v if x['status'] == 'PASS'}) > 1)
print(f"name groups with >1 called library: {len(reps)}; genotype discordant: {disc}; discordant among PASS-only: {pdisc}")
