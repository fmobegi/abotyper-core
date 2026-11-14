"""Per-batch detail for the QAP depth series: reads per subsample, genotype class, A1/A2 and read-level PASS calls
(the two full-depth mix-up libraries left out). usage: python depth_detail.py QAP_DIR"""
import csv, glob, os, statistics as st, subprocess, sys
Q = sys.argv[1]
full = {r["sample"]: r for r in csv.DictReader(open(Q + "/qap_calls.tsv"), delimiter="\t")}
mix = {s for s, r in full.items() if r["truth"] and r["pipe"] != r["truth"] and r["pipe"]}
S = "/mnt/storage_c2/bioinformatics/abo_experiments/env/bin/samtools"
def a2(ext): return any(a.strip().startswith("ABO*A2") for a in ext.split("/"))
for b in sorted(glob.glob(Q + "/depth/d*_r*")):
    c = list(csv.DictReader(open(b + "/calls.tsv"), delimiter="\t"))
    n = [int(subprocess.run([S, "view", "-c", f], capture_output=True, text=True).stdout or 0) for f in glob.glob(b + "/bam/*.bam")]
    t = [x for x in c if x["truth"] and x["sample"] not in mix]
    out = [os.path.basename(b), "reads med %s range %s-%s" % (st.median(n), min(n), max(n)), "n=%d" % len(t)]
    for m, ext in (("pipe", "pipe_ext"), ("read", "read_pair")):
        ok = sum(x[m] == x["truth"] for x in t)
        a2t = [x for x in t if x[m] == x["truth"] and x["truth_a2"] == "1"]
        a2n = [x for x in t if x[m] == x["truth"] and x["truth_a2"] == "0" and "A" in x["truth"]]
        out.append(f"{m}: ok {ok}/{len(t)} ({100*ok/len(t):.1f}%) wrong {sum(bool(x[m]) and x[m]!=x['truth'] for x in t)} A2-sens {sum(a2(x[ext]) for x in a2t)}/{len(a2t)} A1-spec {sum(not a2(x[ext]) for x in a2n)}/{len(a2n)}")
    p = [x for x in t if x["read_status"] == "PASS"]
    out.append(f"read PASS {len(p)} correct {sum(x['read'] == x['truth'] for x in p)}")
    print(" | ".join(out))
