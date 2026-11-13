"""Score pipeline and read-level calls on a validation set against a per-donor truth table.

usage: python qap_eval.py TRUTH.tsv RESULTS_DIR READ_DIR [OUT_CALLS.tsv]

  TRUTH.tsv    donor, class (AO, ...), a2 (0/1), o1_dose, b_dose, ...   (one row per donor; empty class = unknown)
  RESULTS_DIR  pipeline --outdir (ABO_result.txt)
  READ_DIR     <sample>.ReadTyping.tsv files from abotyper-readtyper

Sample names are <run>_<donor>_barcodeNN. Prints aggregate numbers only (no sample identifiers), so the output can
be shared; OUT_CALLS.tsv has one row per library and stays with the data.
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path


def norm_class(g: str) -> str:
    g = (g or "").strip()
    if not g or g in ("Unknown", "nan", "NOCALL") or not set(g) <= set("ABO"):
        return ""
    return "".join(sorted(g, key="ABO".index))


def dose(ext: str, prefix: str) -> int:
    return sum(1 for a in (ext or "").split("/") if a.strip().startswith(prefix))


def main() -> int:
    truth_path, res_dir, read_dir = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
    out_calls = Path(sys.argv[4]) if len(sys.argv) > 4 else None
    truth = {r["donor"]: r for r in csv.DictReader(truth_path.open(), delimiter="\t")}

    rows = list(csv.reader((res_dir / "ABO_result.txt").open(), delimiter="\t"))
    hdr = rows[1]
    ix = {k: hdr.index(k) for k in ("Barcode", "Sequencing_ID", "Genotype", "ExtendedGenotype", "QCStatus",
                                    "QCFlags")}
    pipe = {}
    for r in rows[2:]:
        sid, bc = r[ix["Sequencing_ID"]], int(r[ix["Barcode"]])
        pipe[(sid, bc)] = {k: r[i] for k, i in ix.items()}

    calls = []
    for f in sorted(read_dir.glob("*.ReadTyping.tsv")):
        rt = next(csv.DictReader(f.open(), delimiter="\t"))
        s = rt["sample"]
        m = re.match(r"(?P<run>R\d+)_(?P<donor>.+)_barcode(?P<bc>\d+)$", s)
        if not m:
            continue
        p = pipe.get((f"{m['run']}_{m['donor']}", int(m["bc"])), {})
        t = truth.get(m["donor"], {})
        rpair = f"{rt['allele_1']}/{rt['allele_2']}"
        calls.append(dict(
            sample=s, run=m["run"], donor=m["donor"], truth=t.get("class", ""), truth_a2=t.get("a2", ""),
            truth_o1=t.get("o1_dose", ""), reads=rt.get("reads_typed", ""),
            pipe=norm_class(p.get("Genotype", "")), pipe_ext=p.get("ExtendedGenotype", ""),
            pipe_qc=p.get("QCStatus", ""), read=norm_class(rt.get("genotype", "")), read_pair=rpair,
            read_status=rt.get("status", ""), read_flags=rt.get("flags", ""),
            unexplained=rt.get("unexplained_share", ""), margin=rt.get("margin", "")))

    if out_calls:
        with out_calls.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(calls[0]), delimiter="\t")
            w.writeheader()
            w.writerows(calls)

    print(f"libraries {len(calls)}; donors {len({c['donor'] for c in calls})}; "
          f"with truth {sum(1 for c in calls if c['truth'])}")
    for meth, ext_key, qc_key in (("pipeline", "pipe_ext", "pipe_qc"), ("read", "read_pair", "read_status")):
        a = Counter()
        for c in calls:
            call = c["pipe" if meth == "pipeline" else "read"]
            qc = c[qc_key].upper()
            a["n"] += 1
            a["nocall"] += not call
            a["review"] += qc.startswith("REVIEW")
            if not c["truth"]:
                continue
            a["t"] += 1
            ok = call == c["truth"]
            a["ok"] += ok
            a["wrong"] += bool(call) and not ok
            a["wrong_pass"] += bool(call) and not ok and qc == "PASS"
            if ok and c["truth_a2"] != "":
                a["a2_n"] += 1
                a["a2_ok"] += (dose(c[ext_key], "ABO*A2") > 0) == (c["truth_a2"] == "1")
            if ok and c["truth_o1"] != "":
                a["o1_n"] += 1
                a["o1_ok"] += dose(c[ext_key], "ABO*O.01") == int(c["truth_o1"])
        pct = lambda x, d: f"{100 * x / d:.1f}%" if d else "-"
        print(f"{meth:9s} genotype {a['ok']}/{a['t']} ({pct(a['ok'], a['t'])})  wrong {a['wrong']} "
              f"(wrong+PASS {a['wrong_pass']})  nocall {a['nocall']}  review {a['review']}/{a['n']}  "
              f"A2 {pct(a['a2_ok'], a['a2_n'])}  O1-dose {pct(a['o1_ok'], a['o1_n'])}")

    both = [c for c in calls if c["pipe"] and c["read"]]
    agree = sum(c["pipe"] == c["read"] for c in both)
    print(f"pipeline vs read genotype agreement {agree}/{len(both)}")
    disagree = Counter((c["truth"] == c["pipe"], c["truth"] == c["read"]) for c in both if c["pipe"] != c["read"]
                       and c["truth"])
    print("  disagreements with truth (pipeline right, read right):", dict(disagree))

    by = defaultdict(list)
    for c in calls:
        by[c["donor"]].append(c)
    for meth in ("pipe", "read"):
        reps = [v for v in by.values() if len(v) > 1 and all(x[meth] for x in v)]
        same = sum(len({x[meth] for x in v}) == 1 for v in reps)
        print(f"replicate concordance ({meth}): {same}/{len(reps)} donors")
    flags = Counter()
    for c in calls:
        for f in (c["read_flags"] or "").split("; "):
            if f:
                flags[re.sub(r"[\d.]+%?|ABO\*\S+", "#", f)[:60]] += 1
    print("read-level flags:", dict(flags.most_common(8)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
