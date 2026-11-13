"""Unit tests for behaviour fixed during the v1.0.0 review (#71); they guard the refactors that follow."""
import copy
import csv
from types import SimpleNamespace

import pytest
import yaml

from conftest import PANEL

from abotyper_core import panel as abo_panel
from abotyper_core import aggregate as agg
from abotyper_core import clair2haplotypes as c2h


# ---------------------------------------------------------------- panel loading
@pytest.fixture(scope="module")
def panel_yaml():
    return yaml.safe_load(PANEL.read_text(encoding="utf-8"))


def write_yaml(tmp_path, data):
    p = tmp_path / "panel.yaml"
    p.write_text(yaml.safe_dump(data), encoding="utf-8")
    return p


def test_default_panel_loads():
    panel = abo_panel.load_panel(PANEL)
    assert len(panel.variants) == 135
    assert not agg.missing_required_markers(panel)


@pytest.mark.parametrize("key, value", [("ng006669_2_pos", "24281.5"), ("ng006669_2_pos", "abc"), ("alt_threshold_pct", "high")])
def test_bad_values_fail(tmp_path, panel_yaml, key, value):
    data = copy.deepcopy(panel_yaml)
    data["variants"][0][key] = value
    with pytest.raises(abo_panel.PanelError):
        abo_panel.load_panel(write_yaml(tmp_path, data))


def test_integral_float_is_read_as_int(tmp_path):
    rows = list(csv.DictReader(open(PANEL.with_suffix(".csv"), encoding="utf-8")))
    for r in rows:
        if r["id"] == "b_vs_ao_796":
            r["ng006669_2_pos"] = r["ng006669_2_pos"] + ".0"
    p = tmp_path / "panel.csv"
    with open(p, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    m = {v.id: v for v in abo_panel.load_panel(p).variants}
    assert m["b_vs_ao_796"].ng006669_2_pos == 24281


@pytest.mark.parametrize("mutate", [
    lambda v: v[0].__setitem__("alt_treshold_pct", 95),          # misspelt key
    lambda v: v[1].__setitem__("id", v[0]["id"]),                 # duplicate id
    lambda v: v[0].pop("id"),                                     # missing id
])
def test_bad_rows_fail(tmp_path, panel_yaml, mutate):
    data = copy.deepcopy(panel_yaml)
    mutate(data["variants"])
    with pytest.raises(abo_panel.PanelError):
        abo_panel.load_panel(write_yaml(tmp_path, data))


def test_missing_required_role_is_reported(panel_yaml, tmp_path):
    data = copy.deepcopy(panel_yaml)
    data["variants"] = [v for v in data["variants"] if v["id"] != "b_vs_ao_796"]
    panel = abo_panel.load_panel(write_yaml(tmp_path, data))
    assert agg.missing_required_markers(panel) == ["b_796 -> b_vs_ao_796"]


# ---------------------------------------------------------------- sample names and controls
@pytest.fixture(scope="module")
def parser(tmp_path_factory):
    return agg.ABOReportParser(str(tmp_path_factory.mktemp("in")), abo_panel.load_panel(PANEL), qc_report="")


@pytest.mark.parametrize("name, expected", [
    ("S1_barcode07", ("S1", "barcode07", "explicit")),
    ("ABO_AA_barcode01", ("ABO_AA", "barcode01", "explicit")),
    ("Sample-45-44874_barcode25", ("Sample-45-44874", "barcode25", "explicit")),
    ("patient-01", ("patient-01", "barcode00", "default")),
    ("S1.r1", ("S1.r1", "barcode00", "default")),
    ("X_barcode06_rep2", ("X_barcode06_rep2", "barcode00", "default")),
])
def test_sample_names(parser, name, expected):
    assert parser.extract_sample_info(name) == expected


def test_read_samples(tmp_path):
    f = tmp_path / "samples.tsv"
    f.write_text("S1_x\t\tS1\t7\nNTC\tnegative\t\tbarcode95\nS2\n", encoding="utf-8")
    assert agg.read_samples(str(f)) == {
        "S1_x": {"control": "", "sample": "S1", "barcode": "7", "files": {}},
        "NTC": {"control": "negative", "sample": "", "barcode": "barcode95", "files": {}},
        "S2": {"control": "", "sample": "", "barcode": "", "files": {}},
    }
    assert agg.read_controls(str(f)) == {"NTC": "negative"}
    f.write_text("S1\t\tS1\tbc7\n", encoding="utf-8")
    with pytest.raises(ValueError):
        agg.read_samples(str(f))


def test_read_samples_with_header(tmp_path):
    """A header row names the columns (any order) and can give each sample's files."""
    f = tmp_path / "samples.tsv"
    f.write_text("id\tbarcode\tcontrol\tsample\tphenotype\tcoverage\n"
                 "S1_x\t7\t\tS1\tS1_x.ABOPhenotype.txt\tS1_x.coverage.txt\n"
                 "NTC\t95\tnegative\t\t\t\n", encoding="utf-8")
    assert agg.read_samples(str(f)) == {
        "S1_x": {"control": "", "sample": "S1", "barcode": "7",
                 "files": {"phenotype": "S1_x.ABOPhenotype.txt", "coverage": "S1_x.coverage.txt"}},
        "NTC": {"control": "negative", "sample": "", "barcode": "95", "files": {}},
    }
    f.write_text("id\tcontrol\tpileup\nS1\t\tx\n", encoding="utf-8")
    with pytest.raises(ValueError, match="unknown column"):
        agg.read_samples(str(f))


@pytest.mark.parametrize("folder, row, expected", [
    ("whatever", {"sample": "S1", "barcode": "7"}, ("S1", "barcode07", "explicit")),
    ("whatever", {"sample": "", "barcode": "barcode123"}, ("whatever", "barcode123", "explicit")),
    ("P1_barcode05", {"sample": "P1_barcode05", "barcode": ""}, ("P1_barcode05", "barcode00", "default")),
    ("S9_barcode09", {"sample": "", "barcode": ""}, ("S9", "barcode09", "explicit")),  # empty row: folder name
])
def test_sample_names_from_table(tmp_path, folder, row, expected):
    p = agg.ABOReportParser(str(tmp_path), abo_panel.load_panel(PANEL), qc_report="",
                            sample_table={folder: {"control": "", **row}})
    assert p.extract_sample_info(folder) == expected


def test_pileup_bgzipped(tmp_path):
    import gzip
    import shutil
    from conftest import DATA
    from abotyper_core import pileup2metrics as p2m
    src = DATA / "ABO_AA_barcode01" / "ABO_AA_barcode01.pileup.txt"
    gz = tmp_path / "x.mpileup.gz"
    with open(src, "rb") as fi, gzip.open(gz, "wb") as fo:
        shutil.copyfileobj(fi, fo)
    plain = p2m.load_pileup(src)
    packed = p2m.load_pileup(gz)
    assert plain and plain.keys() == packed.keys()
    slots = p2m.PileupColumn.__slots__
    assert all([getattr(plain[k], a) for a in slots] == [getattr(packed[k], a) for a in slots] for k in plain)


def test_read_controls(tmp_path):
    f = tmp_path / "controls.tsv"
    f.write_text("NTC\tnegative\nPOS1\tPositive\nS1\t\nS2\n", encoding="utf-8")
    assert agg.read_controls(str(f)) == {"NTC": "negative", "POS1": "positive"}
    f.write_text("S1\tcontrol\n", encoding="utf-8")
    with pytest.raises(ValueError):
        agg.read_controls(str(f))


# ---------------------------------------------------------------- haplotypes
def rec(pos, ref, alts, gt, ad):
    return c2h.PhasedRecord(pos, ref, alts, gt, "PS1", ad, sum(ad))


def test_haplotype_sides_sized_by_their_own_allele_depth():
    """0|1 and 1|0 at AD 70/30 give the same counts, each side sized by the depth of the allele it carries."""
    panel = abo_panel.load_panel(PANEL)
    pos = next(m.resolved_position() for m in panel.calibrated() if not m.needs_window())
    for gt in [(0, 1), (1, 0)]:
        rows = c2h.build_rows(panel, {pos: rec(pos, "G", ["A"], gt, [70, 30])})
        haps = [r[3] for r in rows]
        assert sorted([haps.count("REF"), len(haps) - haps.count("REF")]) == [30, 70]


def test_window_marker_ignores_a_nearby_snp():
    marker = SimpleNamespace(variant_type="indel", ref_base="REF")
    assert not c2h.record_fits_marker(rec(24288, "G", ["C"], (0, 1), [10, 10]), marker, 24289)
    assert c2h.record_fits_marker(rec(24288, "GG", ["G"], (0, 1), [10, 10]), marker, 24289)


def test_commands_installed():
    """The package installs the nine abotyper-* commands and nothing else (no .py aliases on PATH)."""
    from importlib.metadata import entry_points
    names = {ep.name for ep in entry_points(group="console_scripts") if ep.value.startswith("abotyper_core.")}
    assert names == {"abotyper-panel", "abotyper-panel-to-bed", "abotyper-clair2metrics", "abotyper-clair2haplotypes",
                     "abotyper-pileup2metrics", "abotyper-predict-phenotype", "abotyper-aggregate",
                     "abotyper-calibrate-panel", "abotyper-rename-samples"}
