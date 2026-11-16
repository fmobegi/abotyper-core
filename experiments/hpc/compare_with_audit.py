"""Read-level typing vs the earlier audit (pileup / Clair3 calls) on the same libraries. Aggregates only."""
import csv
import os
import statistics as st
import sys
from collections import Counter

rt_path, audit_path = sys.argv[1], sys.argv[2]
rt = {r["sample"]: r for r in csv.DictReader(open(rt_path, encoding="utf-8"), delimiter="\t") if r["sample"] != "sample"}
au = {}
for r in csv.DictReader(open(audit_path, encoding="utf-8"), delimiter="\t"):
    au[os.path.basename(r["sample_dir"].rstrip("/\\"))] = r
both = [s for s in rt if s in au]
print(f"read-level {len(rt)}, audit {len(au)}, matched by name {len(both)}")


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def cls(g):
    g = (g or "").strip()
    return "".join(sorted(g, key="ABO".index)) if g and set(g) <= set("ABO") and len(g) == 2 else ""


for label, keep in (("NOCALL (0 informative)", lambda r: r["status"] == "NOCALL" and r["reads_typed"] == "0"),
                    ("called", lambda r: r["status"] in ("PASS", "REVIEW"))):
    sub = [s for s in both if keep(rt[s])]
    d = [f(au[s]["pile_min_depth"]) for s in sub if f(au[s]["pile_min_depth"]) is not None]
    vp = Counter(au[s]["vcf_present"] for s in sub)
    pg = Counter(bool(cls(au[s]["pile_geno"])) for s in sub)
    vg = Counter(bool(cls(au[s]["vcf_geno"])) for s in sub)
    print(f"\n{label}: n={len(sub)}")
    if d:
        print(f"  audit pileup min depth at primary sites: median {st.median(d)}, p5 {sorted(d)[len(d)//20]}, "
              f"p95 {sorted(d)[19*len(d)//20]}, zero {sum(x == 0 for x in d)}")
    print(f"  audit vcf_present: {dict(vp)}; pileup genotype called: {dict(pg)}; vcf genotype called: {dict(vg)}")
    print("  aligned reads (read-level):", st.median([int(rt[s]['reads_aligned']) for s in sub]) if sub else "-")

called = [s for s in both if rt[s]["status"] in ("PASS", "REVIEW")]
for src in ("vcf_geno", "pile_geno"):
    pairs = [(cls(au[s][src]), rt[s]["genotype"]) for s in called if cls(au[s][src])]
    agree = sum(a == b for a, b in pairs)
    print(f"\nread-level vs audit {src}: {agree}/{len(pairs)} agree ({100*agree/max(1,len(pairs)):.2f}%)")
    for stt in ("PASS", "REVIEW"):
        p2 = [(cls(au[s][src]), rt[s]["genotype"]) for s in called if cls(au[s][src]) and rt[s]["status"] == stt]
        print(f"  {stt}: {sum(a == b for a, b in p2)}/{len(p2)}")
    print("  disagreements (audit -> read-level):", dict(Counter((a, b) for a, b in pairs if a != b).most_common(12)))
a2 = [(au[s]["vcf_a2"], rt[s]) for s in called if cls(au[s]["vcf_geno"]) == rt[s]["genotype"] and "A" in rt[s]["genotype"]
      and au[s]["vcf_a2"] not in ("", None)]
ok = sum((f(x) or 0) > 0 and any(a.startswith("ABO*A2") for a in (r["allele_1"], r["allele_2"]))
         or (f(x) or 0) == 0 and not any(a.startswith("ABO*A2") for a in (r["allele_1"], r["allele_2"])) for x, r in a2)
print(f"\nA2 agreement with audit vcf_a2 (same class, A carriers): {ok}/{len(a2)}")
