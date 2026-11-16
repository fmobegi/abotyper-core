"""Reproducibility: retrospective read-level calls (dev2, c.595 rule replayed) vs independent QAP libraries (dev3) of
the same donors. Run on the server; prints aggregates and masked IDs only."""
import csv
import glob
import hashlib
import re
from collections import Counter, defaultdict

Q = "/mnt/storage_c2/bioinformatics/qap_abo"
mask = lambda s: "D" + hashlib.sha1(s.encode()).hexdigest()[:5]
donor = lambda s: re.sub(r"^R\d+_", "", re.sub(r"_barcode\d+.*$", "", s))


def norm(r):
    """(status, genotype, allele pair) with the dev3 naming of O.01.68 + c.595C>T applied to dev2 rows."""
    a = [r["allele_1"], r["allele_2"]]
    mm = r.get("haplotype_mismatches") or ""
    status = r["status"]
    if "ABO*O.01.68: c.595 C>T" in mm and not re.search(r"ABO\*O\.01\.68: c\.595 C>T[^;]*, ", mm):
        a = ["ABO*O.01.12-like" if x == "ABO*O.01.68" else x for x in a]
        other = [f for f in r["flags"].split("; ") if f and not f.startswith("haplotype differs")]
        if status == "REVIEW" and not other:
            status = "PASS"
    return status, r["genotype"], "/".join(sorted(a))


retro = defaultdict(list)
for r in csv.DictReader(open(f"{Q}/retro_dev2_readtyping.tsv"), delimiter="\t"):
    if r["status"] in ("PASS", "REVIEW"):
        retro[donor(r["sample"])].append((r["sample"], *norm(r)))
qap = defaultdict(list)
for f in glob.glob(f"{Q}/readtyper_dev3/*.ReadTyping.tsv"):
    r = next(csv.DictReader(open(f), delimiter="\t"))
    if r["status"] in ("PASS", "REVIEW"):
        qap[donor(r["sample"])].append((r["sample"], r["status"], r["genotype"], "/".join(sorted([r["allele_1"], r["allele_2"]])), r))

both = sorted(set(retro) & set(qap))
print(f"donors with a read-level call in both: {len(both)}  (retrospective libraries {sum(len(retro[d]) for d in both)}, "
      f"QAP libraries {sum(len(qap[d]) for d in both)})")
g_ok = p_ok = 0
dis = []
pass_pairs = Counter()
for d in both:
    rg = {x[2] for x in retro[d]}
    qg = {x[2] for x in qap[d]}
    rp = {x[3] for x in retro[d]}
    qp = {x[3] for x in qap[d]}
    g_ok += rg == qg
    p_ok += rp == qp
    if rg != qg or rp != qp:
        dis.append(d)
    for x in retro[d]:
        for y in qap[d]:
            if x[1] == "PASS" and y[1] == "PASS":
                pass_pairs[(x[2] == y[2], x[3] == y[3])] += 1
print(f"same genotype class in both: {g_ok}/{len(both)}; same ISBT allele pair: {p_ok}/{len(both)}")
print("library pairs where both are PASS -> (same class, same pair):", dict(pass_pairs))
print("\ndonors that differ (masked):")
for d in dis:
    r = "; ".join(f"{x[2]} {x[3]} [{x[1]}]" for x in retro[d])
    q = "; ".join(f"{y[2]} {y[3]} [{y[1]}]" for y in qap[d])
    print(f"  {mask(d)}: retrospective {r} | QAP {q}")

print("\nflagged retrospective donors that were re-sequenced in QAP:")
for d in both:
    flagged = [x for x in retro[d] if x[1] == "REVIEW" or "12-like" in x[3]]
    rows = [r for r in csv.DictReader(open(f"{Q}/retro_dev2_readtyping.tsv"), delimiter="\t") if donor(r["sample"]) == d]
    why = [r["haplotype_mismatches"] or r["novel_variants"] for r in rows if r["haplotype_mismatches"] or r["novel_variants"]]
    if why:
        q = "; ".join(f"{y[2]} {y[3]} [{y[1]}] mm={y[4].get('haplotype_mismatches') or '-'} novel={y[4].get('novel_variants') or '-'} notes={y[4].get('notes') or '-'}" for y in qap[d])
        print(f"  {mask(d)}: retrospective {'; '.join(why)} | QAP {q}")
