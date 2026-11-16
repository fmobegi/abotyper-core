"""Read-by-read evidence at chosen ISBT sites for one library, split by the allele each read is assigned to.

usage: python inspect_sites.py BAM REF SAMPLE_LABEL c.POS [c.POS ...] [--link c.A c.B]

Types the BAM with abotyper-readtyper, then for every requested cDNA site prints, per called allele, how many of its
reads show each state (SNV base or indel length change), and what each allele expects there. --link prints the
joint states of two sites over the reads that cover both (phase), e.g. c.261 and c.796. SAMPLE_LABEL is printed as
given (pass a masked ID); nothing else identifies the sample.
"""

import argparse
import sys
from collections import Counter
from pathlib import Path

import numpy as np

from abotyper_core import readtyper as rt
from abotyper_core.panel import default_panel_path, load_panel


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("bam")
    ap.add_argument("ref")
    ap.add_argument("label")
    ap.add_argument("sites", nargs="+")
    ap.add_argument("--link", nargs=2)
    a = ap.parse_args()
    seq = rt.read_fasta(Path(a.ref))
    offsets = rt.exon_offsets(load_panel(default_panel_path()))
    table = rt.load_alleles(rt.ISBT_TABLE, offsets, seq)
    coding = rt.coding_map(offsets)
    reads, seen = rt.observe_reads(Path(a.bam), table, coding)
    res = rt.type_sample(reads, table, seq, coding)
    s = res.summary
    print(f"== {a.label}: {s['genotype']} {s['allele_1']}/{s['allele_2']} [{s['status']}] aligned {seen}, "
          f"typed {s['reads_typed']}, margin {s['margin']} (next: {s['alt_pair']}), a/b {s['reads_a_b']}, "
          f"foreign {s['unexplained_share']} (third {s['third_allele'] or '-'}), chimeras {s['chimera_share']}")
    if s["flags"]:
        print(f"   flags: {s['flags']}")
    if s.get("notes"):
        print(f"   notes: {s['notes']}")
    for h in res.haplotypes:
        print(f"   haplotype {h.allele}: {h.reads} reads; exact ISBT matches {', '.join(h.exact[:6]) or '-'}; "
              f"mismatches {'; '.join(h.mismatches) or '-'}")

    # informative reads and their assignment, as type_sample does it
    names = table.names
    L = rt.read_likelihoods(reads, table)
    keep = np.ptp(L, axis=1) > 1e-9
    reads = [r for r, k in zip(reads, keep) if k]
    L = L[keep]
    called = [x for x in (s["allele_1"], s["allele_2"]) if x in names]
    idx = {x: names.index(x) for x in called}
    assign = []
    for r in range(len(reads)):
        if len(idx) == 2:
            i, j = idx.values()
            assign.append(called[0] if L[r, i] > L[r, j] else called[1] if L[r, j] > L[r, i] else "tie")
        else:
            assign.append(called[0] if called else "?")

    def site_of(label):
        c = int(label.replace("c.", ""))
        return [i for i, st in enumerate(table.sites) if coding.get(st.lo) == c
                or (st.kind == "indel" and any(coding.get(p) == c for p in range(st.lo, st.hi + 1)))]

    for lab in a.sites:
        for i in site_of(lab):
            st = table.sites[i]
            exp = {x: table.expected[names.index(x)].get(i, st.ref_state) for x in called}
            by = {}
            for r, g in zip(reads, assign):
                if i in r.sites:
                    by.setdefault(g, Counter())[r.sites[i]] += 1
            kind = "SNV" if st.kind == "snv" else "indel (length change)"
            print(f"   {lab} {kind} ref {st.ref_state}; expected " + ", ".join(f"{x}={v}" for x, v in exp.items()))
            for g, c in sorted(by.items()):
                print(f"      reads of {g}: " + ", ".join(f"{k}:{v}" for k, v in c.most_common()))
    if a.link:
        i1, i2 = site_of(a.link[0])[0], site_of(a.link[1])[0]
        joint = Counter((r.sites[i1], r.sites[i2]) for r in reads if i1 in r.sites and i2 in r.sites)
        print(f"   phase {a.link[0]} x {a.link[1]} over reads covering both: " +
              ", ".join(f"{k}:{v}" for k, v in joint.most_common()))
    return 0


if __name__ == "__main__":
    sys.exit(main())
