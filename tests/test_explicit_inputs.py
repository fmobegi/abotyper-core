"""abotyper-aggregate with explicit inputs (--samples naming each sample's files) and negative controls."""
import subprocess
import sys

from conftest import DATA, PANEL, STATE, run

from abotyper_core.inputs import SAMPLE_FILE_COLUMNS

SUFFIX_TO_COLUMN = {v: k for k, v in SAMPLE_FILE_COLUMNS.items()}


def _samples_table(src, extra_rows=()):
    """A header-format samples.tsv naming every file in the chain's per_sample_processing/<id>/combined/."""
    cols = ["id", "control", "sample", "barcode", *SAMPLE_FILE_COLUMNS]
    rows = ["\t".join(cols)]
    for d in sorted(p for p in src.iterdir() if p.is_dir()):
        sid = d.name
        files = {}
        for f in (d / "combined").iterdir():
            suffix = f.name[len(sid) + 1:]
            if suffix in SUFFIX_TO_COLUMN:
                files[SUFFIX_TO_COLUMN[suffix]] = str(f)
        rows.append("\t".join([sid, "", "", ""] + [files.get(c, "") for c in SAMPLE_FILE_COLUMNS]))
    rows.extend(extra_rows)
    return "\n".join(rows) + "\n"


def test_explicit_inputs_match_folder_layout(chain, tmp_path):
    """Naming the files in --samples gives the same report as the per_sample_processing folder layout."""
    src = STATE["aggregate_dir"] / "per_sample_processing"
    (tmp_path / "samples.tsv").write_text(_samples_table(src), encoding="utf-8")
    run("abotyper-aggregate", "--samples", "samples.tsv", "--panel", PANEL, "--threads", "1", cwd=tmp_path)
    for f in ("ABO_result.txt", "final_export.csv", "ABO_evidence_summary_mqc.tsv"):
        assert (tmp_path / f).read_text(encoding="utf-8") == (STATE["aggregate_dir"] / f).read_text(encoding="utf-8"), f


def test_explicit_inputs_need_files(tmp_path):
    (tmp_path / "samples.tsv").write_text("id\tcontrol\nS1\t\n", encoding="utf-8")
    r = subprocess.run([sys.executable, "-m", "abotyper_core.aggregate", "--samples", "samples.tsv",
                        "--panel", str(PANEL)], cwd=tmp_path, capture_output=True, text=True)
    assert r.returncode == 1 and "names no sample files" in r.stdout


def _empty_sample(tmp_path, sid):
    """The files a sample with no reads gets from the typing steps (empty pileup, zero coverage)."""
    w = tmp_path / sid
    w.mkdir()
    (w / "empty.pileup").write_text("", encoding="utf-8")
    cov = (DATA / "ABO_AA_barcode01" / "ABO_AA_barcode01.coverage.txt").read_text(encoding="utf-8").splitlines()
    hdr, row = cov[0], cov[1].split("\t")
    row[3:] = ["0"] * (len(row) - 3)
    (w / f"{sid}.coverage.txt").write_text(hdr + "\n" + "\t".join(row) + "\n", encoding="utf-8")
    run("abotyper-pileup2metrics", "--pileup", "empty.pileup", "--panel", PANEL, "--sample", sid,
        "--merged-out", f"{sid}.AlignmentStatistics.tsv", "--concordance-out", f"{sid}.concordance.tsv",
        "--status-out", f"{sid}.evidence_status.tsv", cwd=w)
    run("abotyper-predict-phenotype", "-i", f"{sid}.AlignmentStatistics.tsv", "-o", f"{sid}.ABOPhenotype.txt",
        "-c", f"{sid}.coverage.txt", "--panel", PANEL, cwd=w)
    files = {"phenotype": f"{sid}.ABOPhenotype.txt", "evidence_status": f"{sid}.evidence_status.tsv",
             "concordance": f"{sid}.concordance.tsv", "coverage": f"{sid}.coverage.txt"}
    return "\t".join([sid, "negative", "", "95"] + [str(w / files[c]) if c in files else ""
                                                  for c in SAMPLE_FILE_COLUMNS])


def test_clean_negative_control_passes(chain, tmp_path):
    """A declared negative control with no reads is the expected result: PASS, not 'likely amplicon dropout'."""
    src = STATE["aggregate_dir"] / "per_sample_processing"
    ntc = _empty_sample(tmp_path, "NTC")
    agg = tmp_path / "agg"
    agg.mkdir()
    (agg / "samples.tsv").write_text(_samples_table(src, [ntc]), encoding="utf-8")
    run("abotyper-aggregate", "--samples", "samples.tsv", "--panel", PANEL, "--threads", "1", cwd=agg)
    lines = (agg / "ABO_result.txt").read_text(encoding="utf-8").splitlines()
    hdr = lines[1].split("\t")
    row = dict(zip(hdr, next(r.split("\t") for r in lines[2:] if r.split("\t")[1] == "NTC")))
    assert row["QCStatus"] == "PASS", row["QCFlags"]
    assert "dropout" not in row["QCFlags"]
    assert row["DepthCheck"].startswith("Negative control: low depth as expected")
