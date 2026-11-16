"""Simulate candidate REVIEW rules on existing ReadTyping outputs. Aggregate counts only."""
import csv, glob, re, sys, collections
sys.path.insert(0, "/mnt/storage_c2/bioinformatics/abo_experiments/core_new/src")
from abotyper_core.readtyper import binom_tail
Q = "/mnt/storage_c2/bioinformatics/qap_abo"; X = "/mnt/storage_c2/bioinformatics/abo_experiments"
truth = {}
for b in ["full"] + sorted(glob.glob(Q + "/depth/d*_r*")):
    f = Q + "/qap_calls.tsv" if b == "full" else b + "/calls.tsv"
    for r in csv.DictReader(open(f), delimiter="\t"):
        truth[(b, r["sample"])] = r["truth"]
def rows():
    for f in glob.glob(Q + "/readtyper/*.ReadTyping.tsv"):
        r = next(csv.DictReader(open(f), delimiter="\t")); yield "qap_full", r, truth.get(("full", r["sample"]), "")
    for b in sorted(glob.glob(Q + "/depth/d*_r*")):
        for f in glob.glob(b + "/read/*.ReadTyping.tsv"):
            r = next(csv.DictReader(open(f), delimiter="\t")); yield "qap_depth", r, truth.get((b, r["sample"]), "")
    for b in sorted(glob.glob(X + "/bench/d*_r*")):
        kind = "test_mix" if "mix" in b else "test_depth"
        for f in glob.glob(b + "/read/*.tsv"):
            r = next(csv.DictReader(open(f), delimiter="\t"))
            t = re.sub(r"^ABO_|_barcode\d+$", "", r["sample"]); t = "".join(sorted(re.sub(r"\d", "", t), key="ABO".index))
            yield kind, r, t
SKEW, CH, CHN = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
agg = collections.defaultdict(collections.Counter)
for kind, r, t in rows():
    if r["status"] not in ("PASS", "REVIEW"): continue
    a = agg[kind]; ok = (r["genotype"] == t) if t else None
    het = r["allele_1"] != r["allele_2"]
    to = [int(x) for x in r["reads_a_b"].split("/")]; k, n = min(to), sum(to)
    skew = het and n and k / n < SKEW and binom_tail(k, n) < 1e-3
    nt = int(r["reads_typed"]); ch = float(r["chimera_share"] or 0)
    chim = ch >= CH and round(ch * nt) >= CHN
    new = (skew or chim) and r["status"] == "PASS"
    a["n"] += 1; a["pass"] += r["status"] == "PASS"
    a["wrong_pass"] += r["status"] == "PASS" and ok is False
    a["new_review"] += new
    a["new_review_wrong"] += new and ok is False
    a["new_review_right"] += new and ok is True
    a["skew"] += bool(skew); a["chim"] += bool(chim)
for k, a in agg.items(): print(k, dict(a))
