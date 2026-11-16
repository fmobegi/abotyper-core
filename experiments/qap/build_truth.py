"""Build the per-donor truth table for the validation (QAP) libraries.

usage: python build_truth.py AUDIT_PER_LIBRARY.tsv FASTQ_ROOT > truth_audit.tsv

AUDIT_PER_LIBRARY.tsv is the earlier validation audit: one row per library sequenced before (any run), with the
genotype read from that library's variant calls (`vcf_geno`, e.g. AO), the A2 marker (`vcf_a2`, 0/1), the O1 dose
(`vcf_o1`, 0-2), the B marker c.796 dose (`vcf_b796`) and the variants seen (`vcf_variants`). The truth for a donor
is what those earlier VCFs say, not a new call: a donor gets a class only if at least one earlier library has a
VCF-based genotype, and `consistent` is 1 when all its earlier libraries agree.

FASTQ_ROOT holds the new runs (<run>/<donor>_barcodeNN.fastq[.gz]); only their donors are written. Run this where the
restricted data lives and keep the output there (chmod 600); it contains donor identifiers.
"""

import collections
import csv
import glob
import os
import re
import sys


def main() -> int:
    audit, root = sys.argv[1], sys.argv[2]
    donors = {re.match(r"(.+?)_barcode", os.path.basename(f)).group(1)
              for f in glob.glob(os.path.join(root, "*", "*.fastq*"))}
    by = collections.defaultdict(list)
    for r in csv.DictReader(open(audit), delimiter="\t"):
        for d in donors:
            if re.search(re.escape(d) + r"(_|$)", r["sample_dir"]) and r["vcf_geno"]:
                by[d].append(r)
    w = sys.stdout
    w.write("donor\tclass\ta2\to1_dose\tb_dose\tvariants\tn_lib\tconsistent\n")
    for d in sorted(donors):
        rs = by.get(d, [])
        if not rs:
            w.write(f"{d}\t\t\t\t\t\t0\t\n")
            continue
        keys = {(r["vcf_geno"], r["vcf_a2"], r["vcf_o1"], r["vcf_b796"]) for r in rs}
        g, a2, o1, b = sorted(keys)[0]
        w.write(f"{d}\t{g}\t{int(float(a2))}\t{int(float(o1))}\t{int(float(b))}\t{rs[0]['vcf_variants']}\t"
                f"{len(rs)}\t{int(len(keys) == 1)}\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
