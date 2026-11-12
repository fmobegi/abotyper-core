"""The ISBT check, offline: a synthetic ISBT page and export stand in for blooddatabase.isbtweb.org."""

import dataclasses
import json
import subprocess
import sys

import pytest

from abotyper_core import isbt
from abotyper_core.panel import load_panel, default_panel_path


def _variant(change, pos):
    return {"hgvs_transcript": f"NM_020469.3:{change}", "sort_position": pos}


def _record(name, pheno, changes, **kw):
    rec = {"id": 1, "system": {"id": 1, "symbol": "ABO"}, "isbt_allele": name, "isbt_phenotype": pheno,
           "isbt_snp": ";".join(changes), "variants": [_variant(c, i) for i, c in enumerate(changes)],
           "exon_intron_change": "<br>".join(kw.pop("exons", ["Exon 7"] * len(changes))),
           "updatedAt": "2026-11-01T00:00:00Z", "approved": True, "obsolete": False, "deleted": False}
    rec.update(kw)
    return rec


def _page(records):
    """An HTML page carrying the records the way the ISBT site does (Next.js flight data)."""
    payload = json.dumps({"alleles": records})
    return f"<html><script>self.__next_f.push([1,{json.dumps(payload)}])</script></html>"


@pytest.fixture
def panel():
    return load_panel(default_panel_path())


@pytest.fixture
def coverage():
    return isbt.load_coverage(isbt.default_coverage_path())


def test_normalise_and_split():
    assert isbt.normalise_change("c.1061delC") == "c.1061del"
    assert isbt.normalise_change("c.804dupG") == "c.804dup"
    assert isbt.split_changes("c.467C>T<br>c.1061delC") == ["c.467C>T", "c.1061del"]
    assert isbt.split_changes("c.539G>C / c.539G>A") == ["c.539G>C", "c.539G>A"]
    assert isbt.split_changes("c.804 region (homopolymer, legacy)") == ["c.804"]


def test_parse_page_uses_structured_variants():
    # the summary field carries a typo; the per-variant HGVS records are right
    rec = _record("ABO*BW.20", "B<sub>weak</sub>", ["c.796C>A", "c.816dup"], isbt_snp="c.815816insG")
    alleles = isbt.parse_alleles(_page([rec, _record("ABO*X.01", "A", ["c.467C>T"], deleted=True)]))
    assert [a.name for a in alleles] == ["ABO*BW.20"]
    assert alleles[0].changes == ["c.796C>A", "c.816dup"]
    assert alleles[0].phenotype == "Bweak"


def test_parse_page_without_alleles_fails():
    with pytest.raises(ValueError, match="layout may have changed"):
        isbt.parse_alleles("<html></html>")


def test_bundled_baseline_matches_panel(panel, coverage):
    """The bundled snapshot is the release the panel was reviewed against: no suggestions."""
    baseline = isbt.load_snapshot(isbt.default_baseline_path())
    assert baseline.version is not None and len(baseline.alleles) > 150
    assert isbt.compare(panel, baseline, coverage, baseline) == []


def test_compare_reports_each_kind(panel, coverage):
    baseline = isbt.load_snapshot(isbt.default_baseline_path())
    live = isbt.IsbtRelease(version=(baseline.version or 0) + 1,
                            alleles=[a for a in baseline.alleles if a.name != "ABO*A2.04"])
    i = next(i for i, a in enumerate(live.alleles) if a.name == "ABO*A2.01")
    live.alleles[i] = dataclasses.replace(live.alleles[i], changes=live.alleles[i].changes + ["c.999A>G"])
    live.alleles.append(isbt.IsbtAllele("ABO*AW.99", "Aweak", ["c.2T>G"], exons="Exon 1"))
    live.alleles.append(isbt.IsbtAllele("ABO*O.01.99", "O", ["c.261del", "c.646T>A"]))

    kinds = {}
    for s in isbt.compare(panel, live, coverage, baseline):
        kinds.setdefault(s.kind, []).append(s)
    assert [s.item for s in kinds["release"]] == [f"ISBT release {live.version}"]
    assert sorted(s.item for s in kinds["new_change"]) == ["c.2T>G", "c.999A>G"]
    assert "assay gap" in next(s for s in kinds["new_change"] if s.item == "c.2T>G").suggestion
    assert sorted(s.item for s in kinds["new_allele"]) == ["ABO*AW.99", "ABO*O.01.99"]
    assert "already handled" in next(s for s in kinds["new_allele"] if s.item == "ABO*O.01.99").suggestion
    assert [s.item for s in kinds["changed_allele"]] == ["ABO*A2.01"]
    assert [s.item for s in kinds["retired_allele"]] == ["ABO*A2.04"]
    assert "panel_change_not_in_isbt" not in kinds


def test_compare_flags_panel_change_no_longer_in_isbt(panel, coverage):
    baseline = isbt.load_snapshot(isbt.default_baseline_path())
    live = isbt.IsbtRelease(version=baseline.version,
                            alleles=[a for a in baseline.alleles if "c.646T>A" not in a.changes])
    flagged = [s for s in isbt.compare(panel, live, coverage, baseline) if s.kind == "panel_change_not_in_isbt"]
    assert "c.646T>A" in [s.item for s in flagged]
    assert "c.467C>T" not in [s.item for s in flagged]


def test_snapshot_round_trip(tmp_path):
    rel = isbt.IsbtRelease(version=7, applied_at="2026-01-02T00:00:00Z",
                           alleles=[isbt.IsbtAllele("ABO*A1.02", "A1", ["c.467C>T"], exons="Exon 7")])
    isbt.write_snapshot(rel, tmp_path / "s.tsv")
    back = isbt.load_snapshot(tmp_path / "s.tsv")
    assert back.version == 7 and back.applied_at == "2026-01-02"
    assert back.alleles[0].changes == ["c.467C>T"]


def test_export_file_is_read(tmp_path):
    export = tmp_path / "abo.csv"
    export.write_text('Phenotype,Allele,DNA Change,Exon/Intron\nA2,ABO*A2.01,"c.467C>T\nc.1061delC","Exon 7\nExon 7"\n'
                      ',not an allele,,\n', encoding="utf-8")
    rel = isbt.fetch_release(str(export))
    assert [(a.name, a.changes) for a in rel.alleles] == [("ABO*A2.01", ["c.467C>T", "c.1061del"])]


def test_cli_offline(tmp_path):
    page = tmp_path / "abo.html"
    page.write_text(_page([_record("ABO*A2.01", "A2", ["c.467C>T", "c.1061delC"]),
                           _record("ABO*AW.99", "Aweak", ["c.2T>G"], exons=["Exon 1"])]), encoding="utf-8")
    report = tmp_path / "report.tsv"
    res = subprocess.run([sys.executable, "-m", "abotyper_core.panel", "--isbt-source", str(page),
                          "--isbt-baseline", "none", "--isbt-report", str(report), "--fail-on-changes"],
                         capture_output=True, text=True)
    assert res.returncode == 3, res.stderr
    rows = report.read_text(encoding="utf-8").splitlines()
    assert rows[0] == "kind\titem\tdetail\tsuggestion"
    assert any(r.startswith("new_change\tc.2T>G\t") for r in rows)
    # panel rows whose changes no allele in this tiny release carries are flagged too
    assert any(r.startswith("panel_change_not_in_isbt\t") for r in rows)
