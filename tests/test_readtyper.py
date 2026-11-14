"""Read-level typer (abotyper-readtyper): BAM reading, ISBT site model, typing checks and the CLI.

Fixture: tests/readtyper/ABO_AO2_sub.bam is 12 reads of the nf-core/abotyper test sample ABO_AO2_barcode04
(ABO*A1.02/ABO*O.02.01), aligned to NG_006669.2 (tests/readtyper/NG_006669v2.fasta).
"""
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from abotyper_core import readtyper as rt
from abotyper_core.bamreader import read_bam
from abotyper_core.panel import load_panel, default_panel_path

FX = Path(__file__).resolve().parent / "readtyper"
BAM, REF = FX / "ABO_AO2_sub.bam", FX / "NG_006669v2.fasta"

CODON = {a + b + c: aa for (a, b, c), aa in zip(
    ((x, y, z) for x in "TCAG" for y in "TCAG" for z in "TCAG"),
    "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")}


@pytest.fixture(scope="module")
def model():
    seq = rt.read_fasta(REF)
    offsets = rt.exon_offsets(load_panel(default_panel_path()))
    return seq, offsets, rt.load_alleles(rt.ISBT_TABLE, offsets, seq), rt.coding_map(offsets)


def test_bam_reader():
    alns = list(read_bam(BAM))
    assert len(alns) == 12
    a = alns[0]
    assert a.ref == "NG_006669.2" and a.is_primary_mapped and a.mapq > 0
    assert sum(n for op, n in a.cigar if op in "MI=XS") == len(a.seq)
    assert a.ref_end > a.pos + 5000                      # full-length amplicon reads


def test_exon_model_translates(model):
    """The coding exons at the derived offsets give an open reading frame ending in a stop codon."""
    seq, offsets, _, _ = model
    assert set(offsets) >= set(rt.EXON_CDNA)
    cds = "".join(seq[lo + offsets[ex] - 1:hi + offsets[ex]] for ex, (lo, hi) in rt.EXON_CDNA.items())
    cds = cds[2:]                                        # c.29 is the 2nd base of codon 10: start at c.31
    prot = "".join(CODON[cds[i:i + 3]] for i in range(0, len(cds) - 2, 3))
    assert prot.endswith("*") and "*" not in prot[:-1]


def test_allele_table(model):
    _, _, table, _ = model
    assert len(table.names) > 200
    assert table.skipped["reference mismatch"] == 0
    a2 = table.expected[table.names.index("ABO*A2.01")]
    assert any(table.sites[i].kind == "indel" and v == -1 for i, v in a2.items())        # c.1061del
    o1 = table.expected[table.names.index("ABO*O.01.01")]
    assert any(table.sites[i].kind == "indel" and v == -1 for i, v in o1.items())        # c.261del
    assert table.expected[table.names.index("ABO*A1.01")] == {}
    coding = model[3]
    assert any(s.kind == "snv" and s.lo not in coding for s in table.sites)   # splice-site changes are kept


def test_types_fixture(model):
    seq, _, table, coding = model
    reads, seen = rt.observe_reads(BAM, table, coding)
    assert seen == 12 and len(reads) == 12
    res = rt.type_sample(reads, table, seq, coding)
    s = res.summary
    assert {s["allele_1"], s["allele_2"]} == {"ABO*A1.02", "ABO*O.02.01"}
    assert s["genotype"] == "AO" and s["phenotype"] == "A1"
    assert s["unexplained_share"] < 0.1
    assert [h.allele for h in res.haplotypes] and all(h.nearest_mismatches == 0 for h in res.haplotypes)


def _obs(table, allele, n):
    """n synthetic reads carrying exactly `allele` at every site."""
    e = table.expected[table.names.index(allele)]
    return [rt.ReadObs({i: e.get(i, s.ref_state) for i, s in enumerate(table.sites)}, {}) for _ in range(n)]


def test_allele_balance(model):
    seq, _, table, coding = model
    res = rt.type_sample(_obs(table, "ABO*O.01.01", 60) + _obs(table, "ABO*B.01", 3), table, seq, coding)
    assert (res.summary["allele_1"], res.summary["allele_2"]) == ("ABO*O.01.01", "ABO*O.01.01")
    assert "called homozygous" in res.summary["flags"] and res.summary["status"] == "REVIEW"
    # at low depth a skewed split is sampling noise: the heterozygote stands
    res = rt.type_sample(_obs(table, "ABO*O.01.01", 10) + _obs(table, "ABO*B.01", 2), table, seq, coding)
    assert res.summary["genotype"] == "BO"
    # a deep, clearly skewed heterozygote (a contaminant posing as the second allele) stays het but goes to REVIEW
    s = rt.type_sample(_obs(table, "ABO*B.01", 92) + _obs(table, "ABO*A1.01", 23), table, seq, coding).summary
    assert s["genotype"] == "AB" and "allele imbalance" in s["flags"] and s["status"] == "REVIEW"
    # a balanced heterozygote is not flagged
    s = rt.type_sample(_obs(table, "ABO*B.01", 50) + _obs(table, "ABO*A1.01", 40), table, seq, coding).summary
    assert "imbalance" not in s["flags"] and s["status"] == "PASS"


def test_third_haplotype_flag(model):
    seq, _, table, coding = model
    reads = _obs(table, "ABO*O.02.01", 40) + _obs(table, "ABO*O.02.02", 40) + _obs(table, "ABO*A1.01", 12)
    s = rt.type_sample(reads, table, seq, coding).summary
    assert s["unexplained_share"] >= 0.1 and "contamination" in s["flags"] and s["status"] == "REVIEW"


def test_rare_allele_needs_haplotype_support(model):
    """A rare allele is kept when its haplotype carries its defining variants, and replaced when it does not."""
    seq, _, table, coding = model
    rare = "ABO*O.01.68"
    s = rt.type_sample(_obs(table, "ABO*O.01.01", 30) + _obs(table, rare, 30), table, seq, coding).summary
    assert rare in (s["allele_1"], s["allele_2"])
    # hide the defining sites of the rare allele from every read: it can no longer be confirmed
    e_r, e_c = (table.expected[table.names.index(n)] for n in (rare, "ABO*O.01.01"))
    defining = {i for i in set(e_r) | set(e_c) if e_r.get(i) != e_c.get(i)}
    reads = _obs(table, "ABO*O.01.01", 30) + _obs(table, rare, 30)
    for r in reads:
        for i in defining:
            r.sites.pop(i, None)
    s = rt.type_sample(reads, table, seq, coding).summary
    assert rare not in (s["allele_1"], s["allele_2"])


def test_phenotype_rules():
    assert rt.phenotype_of("A1", "A2") == "A1"
    assert rt.phenotype_of("A2", "O") == "A2"
    assert rt.phenotype_of("A1", "B") == "A1B"
    assert rt.phenotype_of("O", "O") == "O"
    assert rt.phenotype_of("cisAB", "O") == "cisAB"
    assert rt.genotype_of(["O", "A"]) == "AO" and rt.genotype_of(["B", "cisAB"]) == "B/cisAB"


def test_cli(tmp_path):
    out, hap = tmp_path / "s.tsv", tmp_path / "s.hap.tsv"
    subprocess.run([sys.executable, "-m", "abotyper_core.readtyper", "--bam", BAM, "--ref", REF,
                    "--sample", "S1", "--out", out, "--haplotypes-out", hap], check=True)
    row = next(csv.DictReader(out.open(), delimiter="\t"))
    assert list(row) == rt.SUMMARY_FIELDS and row["sample"] == "S1" and row["genotype"] == "AO"
    haps = list(csv.DictReader(hap.open(), delimiter="\t"))
    assert len(haps) == 2 and list(haps[0]) == rt.HAPLOTYPE_FIELDS


def test_run_and_nocall(monkeypatch):
    res = rt.run(BAM, REF, "x")
    assert res.summary["reads_typed"] == 12 and np.isfinite(float(res.summary["loglik"]))
    monkeypatch.setattr(rt, "MIN_READS", 13)
    res = rt.run(BAM, REF, "x")
    assert res.summary["status"] == "NOCALL" and res.haplotypes == []


def _rt_row(sample, a1, a2, genotype, unexplained=0.0, third=""):
    row = {k: "" for k in rt.SUMMARY_FIELDS}
    row.update(sample=sample, reads_aligned=90, reads_typed=90, allele_1=a1, allele_2=a2, genotype=genotype,
               status="PASS", unexplained_share=unexplained, third_allele=third)
    return row


def test_aggregate_cross_check(chain, tmp_path):
    """abotyper-aggregate adds a ReadTyping note and sends disagreements and contamination to REVIEW."""
    import shutil
    from conftest import STATE, PANEL, run
    agg = tmp_path / "agg"
    shutil.copytree(STATE["aggregate_dir"] / "per_sample_processing", agg / "per_sample_processing")
    rows = {"ABO_AO2_barcode04": _rt_row("ABO_AO2_barcode04", "ABO*A1.02", "ABO*O.02.01", "AO"),
            "ABO_BB_barcode05": _rt_row("ABO_BB_barcode05", "ABO*A1.01", "ABO*B.01", "AB"),
            "ABO_O2O2_barcode10": _rt_row("ABO_O2O2_barcode10", "ABO*O.02.01", "ABO*O.02.02", "OO", 0.15,
                                          "ABO*A1.01")}
    for s, row in rows.items():
        with open(agg / "per_sample_processing" / s / "combined" / f"{s}.ReadTyping.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=rt.SUMMARY_FIELDS, delimiter="\t")
            w.writeheader()
            w.writerow(row)
    run("abotyper-aggregate", "per_sample_processing", "--panel", PANEL, "--threads", "1", cwd=agg)
    lines = (agg / "ABO_result.txt").read_text().splitlines()
    hdr = lines[1].split("\t")
    res = {r.split("\t")[1]: dict(zip(hdr, r.split("\t"))) for r in lines[2:]}
    assert "ReadTyping" in hdr
    ao2, bb, o2 = res["ABO_AO2"], res["ABO_BB"], res["ABO_O2O2"]
    assert ao2["ReadTyping"].startswith("AO ABO*A1.02/ABO*O.02.01") and "disagree" not in ao2["QCFlags"]
    assert "methods disagree" in bb["QCFlags"] and bb["QCStatus"] == "REVIEW"
    assert "possible contamination" in o2["QCFlags"] and o2["QCStatus"] == "REVIEW"
    assert res["ABO_AA"]["ReadTyping"] == "not run"
    export = (agg / "final_export.csv").read_text().splitlines()[0]
    assert "Read" + " Typing" not in export                    # LIMS export layout is unchanged


def _chimeras(table, a, b, n):
    """n reads that follow allele a before the middle distinguishing site and allele b after it."""
    ea, eb = (table.expected[table.names.index(x)] for x in (a, b))
    diff = sorted((i for i, s in enumerate(table.sites) if ea.get(i, s.ref_state) != eb.get(i, s.ref_state)),
                  key=lambda i: table.sites[i].lo)
    cut = table.sites[diff[len(diff) // 2]].lo
    st = {i: (ea if s.lo < cut else eb).get(i, s.ref_state) for i, s in enumerate(table.sites)}
    return [rt.ReadObs(dict(st), {}) for _ in range(n)]


def test_chimeras_are_not_contamination(model):
    seq, _, table, coding = model
    reads = (_obs(table, "ABO*A1.01", 45) + _obs(table, "ABO*O.01.02", 45)
             + _chimeras(table, "ABO*A1.01", "ABO*O.01.02", 12))
    s = rt.type_sample(reads, table, seq, coding).summary
    assert s["genotype"] == "AO" and s["chimera_share"] > 0.05
    assert s["unexplained_share"] < 0.1 and "contamination" not in s["flags"]
    # an excess of chimeras (well above what PCR makes in clean libraries) goes to REVIEW
    reads = (_obs(table, "ABO*A1.01", 40) + _obs(table, "ABO*O.01.02", 40)
             + _chimeras(table, "ABO*A1.01", "ABO*O.01.02", 20))
    s = rt.type_sample(reads, table, seq, coding).summary
    assert "chimeras" in s["flags"] and s["status"] == "REVIEW"
