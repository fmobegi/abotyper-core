"""Score every benchmark batch against the truth encoded in the test sample names (ABO_<truth>_barcodeNN).

usage: python evaluate.py BENCH_DIR > summary.tsv   (also writes BENCH_DIR/calls.tsv with one row per call)
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path


def truth_of(name: str):
    t = re.sub(r"^ABO_|_barcode\d+$", "", name)
    cls = "".join(sorted(re.sub(r"\d", "", t), key="ABO".index))
    return cls, sorted(re.findall(r"O\d", t))


def o_from_extended(ext: str):
    return sorted(f"O{int(n)}" for n in re.findall(r"ABO\*O\.(\d+)", ext or ""))


def pipeline_calls(d: Path):
    """{sequencing_id: (class, o_subtypes, qc_status, extended)} from abotyper-aggregate output."""
    ext = {}
    res = d / "ABO_result.txt"
    if res.is_file():
        rows = list(csv.reader(res.open(), delimiter="\t"))
        hdr = rows[1]
        i_id, i_ext = hdr.index("Sequencing_ID"), hdr.index("ExtendedGenotype")
        for r in rows[2:]:
            if len(r) > i_ext:
                ext[r[i_id]] = r[i_ext]
    out = {}
    fe = d / "final_export.csv"
    if fe.is_file():
        for r in csv.DictReader(fe.open()):
            sid = r["Sample ID"]
            g = (r.get("ABO Geno Type1", "") + r.get("ABO Geno Type2", "")).strip()
            cls = "".join(sorted(g, key=lambda c: "ABO".index(c) if c in "ABO" else 9))
            out[sid] = (cls, o_from_extended(ext.get(sid, "")), r.get("QC Status", ""), ext.get(sid, ""))
    return out


def main() -> int:
    bench = Path(sys.argv[1])
    calls = []
    for b in sorted(bench.glob("d*_r*")):
        m = re.match(r"d(\w+)_r(\d+)", b.name)
        depth, rep = m.group(1), m.group(2)
        samples = sorted(p.stem for p in (b / "bam").glob("*.bam"))
        for s in samples:
            t_cls, t_o = truth_of(s)
            sid = re.sub(r"_barcode\d+$", "", s)
            rt = b / "read" / f"{s}.tsv"
            if rt.is_file():
                r = next(csv.DictReader(rt.open(), delimiter="\t"))
                # prototype (readtyper.py) or packaged abotyper-readtyper output
                pair = r.get("best_pair") or "/".join(x for x in (r.get("allele_1"), r.get("allele_2")) if x)
                o = (sorted(x for x in (r.get("o_subtypes") or "").split("/") if x) if "o_subtypes" in r
                     else o_from_extended(pair))
                call = r["genotype"] if r["genotype"] not in ("NOCALL", "") else ""
                calls.append(dict(depth=depth, rep=rep, method="read", sample=s, truth=t_cls, call=call,
                                  o_truth="/".join(t_o), o_call="/".join(o), qc=r.get("status", ""), extended=pair,
                                  margin=r.get("margin_vs_other_class", r.get("margin", "")),
                                  reads=r.get("reads", r.get("reads_typed", "")),
                                  unexplained=r.get("unexplained_share", ""), note=r.get("note", "")))
            for meth in ("cur", "baq", "q0", "q12", "pipeline"):
                d = b / meth
                if not d.is_dir():
                    continue
                pc = pipeline_calls(d).get(sid)
                cls, o, qc, ext = pc if pc else ("", [], "", "")
                calls.append(dict(depth=depth, rep=rep, method=meth, sample=s, truth=t_cls, call=cls,
                                  o_truth="/".join(t_o), o_call="/".join(o), qc=qc, extended=ext, margin="", reads="",
                                  unexplained="", note=""))
    with open(bench / "calls.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(calls[0]), delimiter="\t")
        w.writeheader()
        w.writerows(calls)

    agg = defaultdict(lambda: defaultdict(int))
    for c in calls:
        k = (c["method"], c["depth"])
        a = agg[k]
        a["n"] += 1
        ok = c["call"] == c["truth"]
        a["class_ok"] += ok
        a["nocall"] += c["call"] in ("", "NOCALL") or "/" in c["call"]
        a["wrong"] += (not ok) and not (c["call"] in ("", "NOCALL") or "/" in c["call"])
        if c["o_truth"]:
            a["o_n"] += 1
            a["o_ok"] += ok and c["o_call"] == c["o_truth"]
        a["review"] += c["qc"].upper().startswith("REVIEW")
        a["wrong_pass"] += (not ok) and c["qc"].upper() == "PASS"

    def dkey(d):
        if d == "full":
            return 10 ** 6
        if d.startswith("mix"):
            return -int(d[3:])          # contamination batches sort after the depth series
        return int(d)
    order = ["read", "pipeline", "cur", "baq", "q0", "q12"]
    print("method\tdepth\tn\tgenotype_correct_%\to_subtype_correct_%\twrong_%\tnocall_%\treview_%\twrong_but_PASS")
    for meth in order:
        for depth in sorted({d for m, d in agg if m == meth}, key=dkey, reverse=True):
            a = agg[(meth, depth)]
            n = a["n"]
            pct = lambda x, d=n: f"{100 * x / d:.0f}" if d else ""
            print(f"{meth}\t{depth}\t{n}\t{pct(a['class_ok'])}\t{pct(a['o_ok'], a['o_n'])}\t{pct(a['wrong'])}\t"
                  f"{pct(a['nocall'])}\t{pct(a['review'])}\t{a['wrong_pass']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
