"""Summarise abotyper-readtyper results over a large set, optionally against existing calls.

usage: python summarise_readtyping.py READTYPING_DIR [--calls CALLS.tsv] [--discordant OUT.tsv] > summary.txt

CALLS.tsv: tab-separated with a header; columns `sample` and `genotype` (e.g. AO, or ABO*A1.01/ABO*O.01.01 style
genotypes reduced to the class), optionally `phenotype` and `status`. `sample` must match the BAM name.
The summary on stdout has aggregate numbers only (no sample names), so it can leave the HPC. --discordant writes the
per-sample disagreements for review; keep that file with the data.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
from collections import Counter


def cls(g: str) -> str:
    g = (g or "").strip()
    if "*" in g:                                   # allele pair: take the letters after ABO*
        g = "".join(re.findall(r"ABO\*(A|B|O)", g))
    if not g or not set(g) <= set("ABO"):
        return ""
    return "".join(sorted(g, key="ABO".index))


def q(values, ps=(0.01, 0.05, 0.5, 0.95, 0.99)):
    v = sorted(values)
    return " ".join(f"p{int(p * 100)}={v[min(len(v) - 1, int(p * len(v)))]:.3f}" for p in ps) if v else "-"


def depth_bin(n: int) -> str:
    for lo, hi in ((0, 5), (5, 10), (10, 20), (20, 50), (50, 200), (200, 1000)):
        if n < hi:
            return f"{lo}-{hi - 1}"
    return "1000+"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("dir")
    ap.add_argument("--calls")
    ap.add_argument("--discordant")
    a = ap.parse_args()
    rows = []
    for f in glob.glob(os.path.join(a.dir, "*.ReadTyping.tsv")):
        with open(f) as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    print(f"samples typed: {len(rows)}")
    print("status:", dict(Counter(r["status"] for r in rows)))
    print("genotype:", dict(Counter(r["genotype"] for r in rows).most_common()))
    print("phenotype:", dict(Counter(r["phenotype"] for r in rows).most_common(12)))

    by = Counter((depth_bin(int(r["reads_typed"] or 0)), r["status"]) for r in rows)
    print("\nstatus by informative reads:")
    for b in sorted({k[0] for k in by}, key=lambda s: int(s.split("-")[0].rstrip("+"))):
        print(f"  {b:>9}: " + ", ".join(f"{s} {by[(b, s)]}" for s in ("PASS", "REVIEW", "NOCALL") if by[(b, s)]))

    flags = Counter()
    for r in rows:
        for f in filter(None, (r["flags"] or "").split("; ")):
            flags[re.sub(r"ABO\*\S+|c\.\S+|[\d.]+%?", "#", f)[:70]] += 1
    print("\nflags:")
    for f, n in flags.most_common(15):
        print(f"  {n:6d}  {f}")

    called = [r for r in rows if r["status"] in ("PASS", "REVIEW")]
    het = [r for r in called if r["allele_1"] != r["allele_2"] and int(r["reads_typed"] or 0) >= 50]
    print("\ndistributions (called samples; allele balance for heterozygotes with >=50 informative reads):")
    print("  minor allele share:", q([float(r["minor_allele_share"]) for r in het if r["minor_allele_share"]]))
    print("  foreign-read share:", q([float(r["unexplained_share"]) for r in called if r["unexplained_share"]]))
    print("  chimera share:     ", q([float(r["chimera_share"]) for r in called if r["chimera_share"]]))
    print("  margin:            ", q([float(r["margin"]) for r in called if r["margin"] not in ("", "inf")]))
    print("  rare alleles called:", sum(1 for r in called for x in (r["allele_1"], r["allele_2"])
                                       if not re.match(r"ABO\*(A1\.01|A2\.01|B\.01|O\.01\.0[12]|O\.02\.0[12])$", x)))
    print("  samples with novel variants:", sum(1 for r in called if r["novel_variants"]))

    if a.calls:
        with open(a.calls) as fh:
            ref = {r["sample"]: r for r in csv.DictReader(fh, delimiter="\t")}
        pairs = [(r, ref[r["sample"]]) for r in called if r["sample"] in ref and cls(ref[r["sample"]]["genotype"])]
        print(f"\nagainst existing calls: {len(pairs)} samples in both")
        agree = sum(cls(r["genotype"]) == cls(c["genotype"]) for r, c in pairs)
        print(f"  genotype class agreement: {agree}/{len(pairs)} ({100 * agree / max(1, len(pairs)):.2f}%)")
        for st in ("PASS", "REVIEW"):
            sub = [(r, c) for r, c in pairs if r["status"] == st]
            ag = sum(cls(r["genotype"]) == cls(c["genotype"]) for r, c in sub)
            print(f"  read-level {st}: {ag}/{len(sub)} agree")
        conf = Counter((cls(c["genotype"]), cls(r["genotype"])) for r, c in pairs if cls(r["genotype"]) != cls(c["genotype"]))
        print("  disagreements (existing -> read-level):", dict(conf.most_common(15)))
        if "phenotype" in next(iter(ref.values()), {}):
            ph = [(r, c) for r, c in pairs if c.get("phenotype")]
            sub = Counter((c["phenotype"], r["phenotype"]) for r, c in ph if c["phenotype"] != r["phenotype"])
            print(f"  phenotype agreement: {len(ph) - sum(sub.values())}/{len(ph)}; top differences:",
                  dict(sub.most_common(10)))
        if a.discordant:
            with open(a.discordant, "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["sample", "existing", "read_level", "allele_1", "allele_2", "status", "reads_typed",
                            "flags"])
                for r, c in pairs:
                    if cls(r["genotype"]) != cls(c["genotype"]):
                        w.writerow([r["sample"], c["genotype"], r["genotype"], r["allele_1"], r["allele_2"],
                                    r["status"], r["reads_typed"], r["flags"]])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
