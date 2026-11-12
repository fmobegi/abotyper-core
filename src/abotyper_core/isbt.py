"""Compare the ABO variant panel with the live ISBT blood group allele database and suggest panel updates.

The ISBT database (https://blooddatabase.isbtweb.org/system/ABO) has no documented allele API: the ABO page embeds
every allele record as Next.js data, and /api/release/latest/info gives the current release. This module reads both,
then sets each ISBT allele and DNA change against what the panel already knows: the panel's own cdna_change values
plus the coverage table (data/abo_isbt_v9_coverage.tsv), which records every ISBT change and how the panel handles it
(on_panel, notation_equivalent, structural_variant, not_diagnostic, pipeline_defined).

Nothing is changed automatically: the result is a list of suggestions for a person to review.
"""

from __future__ import annotations

import csv
import json
import re
import urllib.request
from dataclasses import dataclass, field
from datetime import date
from pathlib import Path
from typing import Iterable, Optional

from . import __version__
from .panel import Panel

ISBT_BASE = "https://blooddatabase.isbtweb.org"
ISBT_ABO_URL = f"{ISBT_BASE}/system/ABO"
ISBT_RELEASE_URL = f"{ISBT_BASE}/api/release/latest/info"
USER_AGENT = f"abotyper-core/{__version__} (+https://github.com/fmobegi/abotyper-core)"


def default_coverage_path() -> Path:
    """The ISBT coverage table shipped with the package."""
    return Path(__file__).resolve().parent / "data" / "abo_isbt_v9_coverage.tsv"


# ---------------------------------------------------------------------------
# Notation
# ---------------------------------------------------------------------------

_DEL_DUP_BASES = re.compile(r"(del|dup)[ACGTN]+$", re.IGNORECASE)


def normalise_change(change: str) -> str:
    """Canonical form of one HGVS c. change, so 'c.1061delC' and 'c.1061del' compare equal."""
    c = change.strip().replace(" ", "")
    c = _DEL_DUP_BASES.sub(lambda m: m.group(1).lower(), c)
    return c


_CHANGE = re.compile(r"c\.[^\s;,/<(]+")


def split_changes(text: str) -> list[str]:
    """The c. changes in a ';', '<br>', ' / ' or newline separated list; free text after a change is ignored
    ('c.804 region (homopolymer)' gives 'c.804')."""
    parts = re.split(r";|<br\s*/?>|\s/\s|,|\n", text or "")
    out = []
    for p in parts:
        m = _CHANGE.match(p.strip())
        if m:
            out.append(normalise_change(m.group(0)))
    return out


def _plain(html: str) -> str:
    return re.sub(r"<[^>]+>", "", html or "").strip()


# ---------------------------------------------------------------------------
# ISBT data
# ---------------------------------------------------------------------------

@dataclass
class IsbtAllele:
    name: str
    phenotype: str
    changes: list[str]
    updated: str = ""          # ISO timestamp of the record's last update
    exons: str = ""
    obsolete: bool = False
    provisional: bool = False
    approved: bool = True

    @property
    def current(self) -> bool:
        return self.approved and not self.obsolete


@dataclass
class IsbtRelease:
    version: Optional[int] = None
    applied_at: str = ""
    alleles: list[IsbtAllele] = field(default_factory=list)


def _http_get(url: str, timeout: float) -> str:
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:  # noqa: S310 (fixed https URL)
        return resp.read().decode("utf-8")


def _next_data(html: str) -> str:
    """Join the Next.js flight-data chunks embedded in a page."""
    chunks = re.findall(r'self\.__next_f\.push\(\[1,("(?:[^"\\]|\\.)*")\]\)', html)
    return "".join(json.loads(c) for c in chunks)


def _json_objects_with(text: str, key: str) -> Iterable[dict]:
    """Yield every top-level JSON object in text that has the given key."""
    decoder = json.JSONDecoder()
    marker = '{"id":'
    i = text.find(marker)
    while i != -1:
        try:
            obj, end = decoder.raw_decode(text, i)
        except ValueError:
            i = text.find(marker, i + 1)
            continue
        if isinstance(obj, dict) and key in obj:
            yield obj
            i = text.find(marker, end)
        else:
            i = text.find(marker, i + 1)


def _record_changes(rec: dict) -> list[str]:
    """An allele's c. changes, from its structured variant records (HGVS against NM_020469.3); the summary
    fields (isbt_snp, nt_change) are only a fallback, as they can carry typos (e.g. 'c.815816insG')."""
    out: list[str] = []
    for v in sorted(rec.get("variants") or [], key=lambda v: v.get("sort_position") or 0):
        hgvs = v.get("hgvs_transcript") or v.get("input") or v.get("freeform_dna_change") or ""
        out += split_changes(hgvs.split(":", 1)[-1])
    return out or split_changes(rec.get("nt_change") or rec.get("isbt_snp") or "")


def parse_alleles(html: str) -> list[IsbtAllele]:
    """The ABO allele records embedded in a saved or fetched ISBT system page."""
    seen: dict[str, IsbtAllele] = {}
    for rec in _json_objects_with(_next_data(html), "isbt_allele"):
        if rec.get("deleted"):
            continue
        name = rec["isbt_allele"]
        seen[name] = IsbtAllele(
            name=name,
            phenotype=_plain(rec.get("isbt_phenotype", "")),
            changes=_record_changes(rec),
            updated=rec.get("updatedAt") or "",
            exons=_plain((rec.get("exon_intron_change") or "").replace("<br>", "; ")),
            obsolete=bool(rec.get("obsolete")),
            provisional=bool(rec.get("provisional")),
            approved=bool(rec.get("approved", True)),
        )
    if not seen:
        raise ValueError("no ISBT allele records found in the page; the site layout may have changed")
    return list(seen.values())


def _cell(row: dict, *names: str) -> str:
    for n in names:
        for k, v in row.items():
            if str(k).strip().lower() == n and v is not None:
                return str(v)
    return ""


def parse_export(path: Path) -> list[IsbtAllele]:
    """The alleles in the site's EXPORT download (.xlsx, or .csv/.tsv): one row per allele, with Allele,
    Phenotype, DNA Change and Exon/Intron columns; multi-change cells are split on newlines, ';' or '<br>'."""
    import pandas as pd

    if path.suffix.lower() in {".xlsx", ".xls"}:
        df = pd.read_excel(path, dtype=str)
    else:
        df = pd.read_csv(path, dtype=str, sep=None, engine="python")
    out: dict[str, IsbtAllele] = {}
    for row in df.fillna("").to_dict("records"):
        name = _cell(row, "allele", "isbt allele", "allele name").strip()
        if not name.startswith("ABO*"):
            continue
        changes = _cell(row, "dna change", "nucleotide change", "isbt_snp").replace("\n", ";")
        exons = re.sub(r"\s*\n\s*", "; ", _cell(row, "exon/intron", "exon/intron change", "exon"))
        out[name] = IsbtAllele(name=name, phenotype=_plain(_cell(row, "phenotype")),
                               changes=split_changes(changes), exons=_plain(exons))
    if not out:
        raise ValueError(f"no ABO alleles found in {path}; expected 'Allele' and 'DNA Change' columns")
    return list(out.values())


def fetch_release(source: Optional[str] = None, timeout: float = 60) -> IsbtRelease:
    """Read the ABO alleles from the ISBT site, a saved page (.html) or the site's EXPORT file (.xlsx/.csv)."""
    if source and not source.startswith(("http://", "https://")):
        path = Path(source)
        if path.suffix.lower() in {".xlsx", ".xls", ".csv", ".tsv"}:
            return IsbtRelease(alleles=parse_export(path))
        return IsbtRelease(alleles=parse_alleles(path.read_text(encoding="utf-8")))
    release = IsbtRelease(alleles=parse_alleles(_http_get(source or ISBT_ABO_URL, timeout)))
    try:
        info = json.loads(_http_get(ISBT_RELEASE_URL, timeout))
        release.version = info.get("releaseVersion")
        release.applied_at = info.get("appliedAt", "")
    except (OSError, ValueError):
        pass  # the release number is informative only
    return release


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

@dataclass
class Suggestion:
    kind: str        # new_change | new_allele | changed_allele | retired_allele | panel_change_not_in_isbt | release
    item: str
    detail: str
    suggestion: str


def load_coverage(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def compare(panel: Panel, release: IsbtRelease, coverage: list[dict],
            baseline: Optional[IsbtRelease] = None) -> list[Suggestion]:
    """Suggestions for bringing the panel and its coverage table in line with an ISBT release.

    DNA changes are checked against the panel and the coverage table. Alleles are checked against a baseline
    snapshot of the ISBT allele list (the one the panel was last reviewed against), when one is given.
    """
    out: list[Suggestion] = []
    meta = panel.meta or {}

    panel_changes: dict[str, str] = {}
    for m in list(panel.variants) + list(panel.structural_variants):
        for c in split_changes(getattr(m, "cdna_change", "")):
            panel_changes.setdefault(c, m.id)
    covered = {normalise_change(r["isbt_change"]): r for r in coverage if r.get("isbt_change")}
    known_changes = set(panel_changes) | set(covered)

    current = [a for a in release.alleles if a.current]
    live_names = {a.name for a in current}
    by_change: dict[str, list[IsbtAllele]] = {}
    for a in current:
        for c in a.changes:
            by_change.setdefault(c, []).append(a)

    seen_release = baseline.version if baseline else meta.get("isbt_release")
    if release.version is not None and str(seen_release or "") != str(release.version):
        out.append(Suggestion(
            "release", f"ISBT release {release.version}",
            f"applied {release.applied_at[:10] or 'unknown'}; the panel was last checked against release "
            f"{seen_release or 'unknown'} ({meta.get('isbt_allele_table_checked', 'date unknown')})",
            "after reviewing the items below, refresh the baseline (--isbt-snapshot) and set "
            f"meta.isbt_release: {release.version}, meta.isbt_allele_table_checked: {date.today().isoformat()}"))

    for c in sorted(set(by_change) - known_changes, key=_change_sort_key):
        alleles = by_change[c]
        names = ";".join(a.name for a in alleles)
        phenos = ", ".join(sorted({a.phenotype for a in alleles if a.phenotype}))
        exon1 = all(_exon_of(a, c) == "Exon 1" for a in alleles)
        hint = ("exon 1 is outside the exon 2-7 amplicon: record it in the coverage table as an assay gap"
                if exon1 else
                "add a panel row (named_marker) or record it in the coverage table with a status")
        out.append(Suggestion("new_change", c, f"{len(alleles)} allele(s): {names} ({phenos})", hint))

    if baseline is not None:
        before = {a.name: a for a in baseline.alleles if a.current}
        for a in current:
            old = before.get(a.name)
            if old is None:
                markers = sorted({panel_changes.get(c) or covered.get(c, {}).get("panel_marker", "")
                                  for c in a.changes} - {""})
                unknown = [c for c in a.changes if c not in known_changes]
                hint = ("its new change(s) are listed under new_change" if unknown else
                        "its changes are already handled (" + ", ".join(markers) + "); "
                        "add it to the coverage table and check the calling can name it")
                out.append(Suggestion("new_allele", a.name,
                                      f"{a.phenotype}; changes {';'.join(a.changes) or 'none listed'}", hint))
            elif sorted(old.changes) != sorted(a.changes) or old.phenotype != a.phenotype:
                out.append(Suggestion(
                    "changed_allele", a.name,
                    f"was {old.phenotype}: {';'.join(old.changes)}; now {a.phenotype}: {';'.join(a.changes)}",
                    "update the coverage table and any panel rows that name this allele"))
        for name in sorted(set(before) - live_names):
            gone = next((a for a in release.alleles if a.name == name), None)
            why = "marked obsolete" if gone and gone.obsolete else "no longer listed"
            out.append(Suggestion("retired_allele", name, f"{why} in ISBT",
                                  "remove it from the coverage table and from any panel associated_alleles"))

    live_changes = set(by_change)
    for c, marker in sorted(panel_changes.items(), key=lambda kv: _change_sort_key(kv[0])):
        if c in live_changes:
            continue
        rows = [r for r in coverage if r.get("panel_marker") == marker or normalise_change(r["isbt_change"]) == c]
        if any(r.get("status") in {"pipeline_defined", "not_diagnostic"} for r in rows):
            continue
        if any(normalise_change(r["isbt_change"]) in live_changes for r in rows):
            continue  # an equivalent notation of the same change is still current
        if any(n.strip() in live_names for r in rows for n in (r.get("alleles") or "").split(";")):
            continue  # e.g. a structural variant ISBT describes in free text
        out.append(Suggestion("panel_change_not_in_isbt", c, f"panel marker {marker}",
                              "no current ISBT allele carries this change; check the notation or retire the row"))
    return out


def _exon_of(allele: IsbtAllele, change: str) -> str:
    """The exon/intron ISBT lists for one of an allele's changes ('' when it cannot tell)."""
    exons = [e.strip() for e in allele.exons.split(";")]
    if change in allele.changes and len(exons) == len(allele.changes):
        return exons[allele.changes.index(change)]
    return exons[0] if len(exons) == 1 else ""


def _change_sort_key(c: str):
    m = re.match(r"c\.(-?\d+)", c)
    return (int(m.group(1)) if m else 10**9, c)


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

def write_report(suggestions: list[Suggestion], path: Path) -> None:
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["kind", "item", "detail", "suggestion"])
        for s in suggestions:
            w.writerow([s.kind, s.item, s.detail, s.suggestion])


def default_baseline_path() -> Path:
    """The ISBT allele list the bundled panel was last reviewed against."""
    return Path(__file__).resolve().parent / "data" / "isbt_abo_alleles.tsv"


def write_snapshot(release: IsbtRelease, path: Path) -> None:
    """Save the allele list (one row per allele); use it as the baseline for the next check."""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write(f"# ISBT ABO alleles; isbt_release={release.version if release.version is not None else ''}; "
                 f"applied={release.applied_at[:10]}; fetched={date.today().isoformat()}; source={ISBT_ABO_URL}\n")
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["allele", "phenotype", "changes", "exons", "updated", "obsolete", "provisional"])
        for a in sorted(release.alleles, key=lambda a: a.name):
            w.writerow([a.name, a.phenotype, ";".join(a.changes), a.exons, a.updated[:10],
                        int(a.obsolete), int(a.provisional)])


def load_snapshot(path: Path) -> IsbtRelease:
    release = IsbtRelease()
    with open(path, encoding="utf-8", newline="") as fh:
        lines = fh.read().splitlines()
    if lines and lines[0].startswith("#"):
        head = dict(kv.strip().split("=", 1) for kv in lines[0][1:].split(";") if "=" in kv)
        release.version = int(head["isbt_release"]) if head.get("isbt_release", "").isdigit() else None
        release.applied_at = head.get("applied", "")
        lines = lines[1:]
    for row in csv.DictReader(lines, delimiter="\t"):
        release.alleles.append(IsbtAllele(
            name=row["allele"], phenotype=row["phenotype"], changes=split_changes(row["changes"]),
            exons=row["exons"], updated=row["updated"], obsolete=row["obsolete"] == "1",
            provisional=row["provisional"] == "1"))
    return release


def format_text(suggestions: list[Suggestion], release: IsbtRelease) -> str:
    head = (f"ISBT ABO: {len([a for a in release.alleles if a.current])} current alleles"
            + (f", release {release.version}" if release.version is not None else ""))
    if not suggestions:
        return head + "\nThe panel and coverage table match this ISBT release; nothing to update."
    lines = [head, f"{len(suggestions)} suggestion(s):"]
    for kind in ["release", "new_change", "new_allele", "changed_allele", "retired_allele",
                 "panel_change_not_in_isbt"]:
        group = [s for s in suggestions if s.kind == kind]
        if not group:
            continue
        lines.append(f"\n[{kind}] ({len(group)})")
        for s in group:
            lines.append(f"  {s.item}: {s.detail}\n    -> {s.suggestion}")
    return "\n".join(lines)
