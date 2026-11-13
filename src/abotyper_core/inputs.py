"""Readers for the per-sample files the aggregation uses besides the phenotype report: samtools stats/flagstat/coverage, the ABO_EVIDENCE status file, NanoPlot NanoStats and the sample table.

Part of abotyper-core (https://github.com/fmobegi/abotyper-core), the ABO typing logic of nf-core/abotyper.
Author: Fredrick M. Mobegi; contributors: Benedict Matern, Mathijs Groeneweg. MIT License (see LICENSE).
"""
from __future__ import annotations

import glob
import gzip
import os
import re
from typing import Dict, Optional


def _first_file(sample_dir: str, patterns) -> Optional[str]:
    """First non-empty file under sample_dir matching any glob pattern (recursive).
    Prefers a combined/ location, then the shallowest path, then alphabetical."""
    found = []
    for pat in patterns:
        found.extend(glob.glob(os.path.join(sample_dir, "**", pat), recursive=True))
    found = [f for f in found if os.path.isfile(f) and os.path.getsize(f) > 0]
    if not found:
        return None

    def rank(path):
        rel = os.path.relpath(path, sample_dir)
        parts = rel.replace("\\", "/").split("/")
        return (0 if "combined" in parts[:-1] else 1, len(parts), rel)

    return sorted(found, key=rank)[0]


def _weighted_median(pairs) -> Optional[float]:
    """Median of values given as (value, count) pairs."""
    pairs = sorted((v, c) for v, c in pairs if c > 0)
    total = sum(c for _, c in pairs)
    if total <= 0:
        return None
    half = total / 2.0
    running = 0
    for i, (v, c) in enumerate(pairs):
        running += c
        if running > half:
            return float(v)
        if running == half:
            nxt = pairs[i + 1][0] if i + 1 < len(pairs) else v
            return (v + nxt) / 2.0
    return float(pairs[-1][0])


def parse_alignment_metrics(sample_dir: str) -> Dict[str, Optional[float]]:
    """Per-sample alignment metrics from samtools stats (preferred) and flagstat (fallback).

    Keys: reads, read_length, insert_size, mappability, median_depth (None when unavailable).
    Never raises. Note samtools stats caps the COV histogram at its last bin
    (default 1000+), so a median depth above the cap is reported as the cap.
    """
    m: Dict[str, Optional[float]] = {
        "reads": None, "read_length": None, "insert_size": None,
        "mappability": None, "median_depth": None,
    }
    mapped = None
    try:
        stats_file = _first_file(sample_dir, ["*.stats"])
        if stats_file:
            cov = []
            with open(stats_file, "r", encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if line.startswith("SN\t"):
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) < 3:
                            continue
                        key, val = parts[1].rstrip(":").strip(), parts[2].strip()
                        try:
                            if key == "raw total sequences":
                                m["reads"] = float(val)
                            elif key == "reads mapped":
                                mapped = float(val)
                            elif key == "average length":
                                m["read_length"] = float(val)
                            elif key == "insert size average":
                                m["insert_size"] = float(val)
                        except ValueError:
                            pass
                    elif line.startswith("COV\t"):
                        parts = line.rstrip("\n").split("\t")
                        if len(parts) >= 4:
                            try:
                                cov.append((int(parts[2]), int(parts[3])))
                            except ValueError:
                                pass
            m["median_depth"] = _weighted_median(cov)
    except OSError:
        pass

    try:
        if m["reads"] is None or mapped is None:
            flag_file = _first_file(sample_dir, ["*.flagstat"])
            if flag_file:
                with open(flag_file, "r", encoding="utf-8", errors="replace") as fh:
                    for line in fh:
                        t = re.match(r"^(\d+) \+ \d+ in total", line)
                        if t and m["reads"] is None:
                            m["reads"] = float(t.group(1))
                        g = re.match(r"^(\d+) \+ \d+ mapped \(", line)
                        if g and mapped is None:
                            mapped = float(g.group(1))
    except OSError:
        pass

    m["mapped"] = mapped
    if m["reads"] and mapped is not None:
        m["mappability"] = 100.0 * mapped / m["reads"]
    if m["insert_size"] is not None and m["insert_size"] <= 0:
        m["insert_size"] = None   # 0 = unpaired/long reads: not meaningful
    return m


def parse_coverage_file(sample_dir: str) -> Dict[str, Optional[float]]:
    """ABO-specific figures from a samtools coverage table (*.coverage.txt, one row per contig).

    abo_reads     reads mapped to the ABO reference (numreads)
    covered_depth mean depth over the bases that have any coverage (meandepth * length / covbases)
    Both None when the file is absent. In an ABO-HLA duplex run the HLA reads do not map to the
    ABO reference, so these describe the ABO amplicon rather than the whole library.
    """
    out: Dict[str, Optional[float]] = {"abo_reads": None, "covered_depth": None}
    try:
        f = _first_file(sample_dir, ["*.coverage.txt"])
        if not f:
            return out
        with open(f, "r", encoding="utf-8", errors="replace") as fh:
            rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip() and not ln.startswith("#")]
        reads = depth_x_len = cov = 0.0
        for r in rows:
            if len(r) < 7:
                continue
            start, end, numreads, covbases, meandepth = (float(r[1]), float(r[2]), float(r[3]),
                                                         float(r[4]), float(r[6]))
            reads += numreads
            cov += covbases
            depth_x_len += meandepth * (end - start + 1)
        out["abo_reads"] = reads
        if cov > 0:
            out["covered_depth"] = depth_x_len / cov
    except (OSError, ValueError):
        pass
    return out


def parse_evidence_status(sample_dir: str) -> Dict[str, object]:
    """Two-line evidence summary written by the pipeline's ABO_EVIDENCE step (*.evidence_status.tsv).

    evidence_source  CLAIR3+BAM | CLAIR3_ONLY | BAM_ONLY | NONE ('' when no file was found)
    filled_primary   primary typing positions taken from the BAM pileup (Clair3 had no evidence)
    discordant       panel positions where Clair3 and the BAM pileup disagree
    primary_names    the primary typing positions among them
    """
    out: Dict[str, object] = {"evidence_source": "", "filled_primary": None, "bam_only": None,
                              "discordant": None, "discordant_primary": None, "primary_names": "",
                              "bam_max_depth": None}
    try:
        f = _first_file(sample_dir, ["*.evidence_status.tsv"])
        if not f:
            return out
        with open(f, "r", encoding="utf-8", errors="replace") as fh:
            lines = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
        if len(lines) >= 2:
            row = dict(zip(lines[0], lines[1]))
            out["evidence_source"] = row.get("evidence_source", "")
            out["filled_primary"] = int(row.get("filled_primary_from_bam", "0"))
            out["bam_only"] = int(row.get("bam_only", "0"))
            if row.get("bam_max_depth", "") != "":
                out["bam_max_depth"] = int(row["bam_max_depth"])
            out["discordant"] = int(row.get("discordant", "0"))
            out["discordant_primary"] = int(row.get("discordant_primary", "0"))
            out["primary_names"] = row.get("discordant_primary_names", "").replace(";", ", ")
    except (OSError, ValueError):
        pass
    return out


def parse_nanostats(path: str) -> Dict[str, Optional[float]]:
    """Read count and mean read length from a NanoPlot NanoStats.txt (None when absent)."""
    out: Dict[str, Optional[float]] = {"reads": None, "read_length": None}
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                key, _, val = line.partition(":")
                val = val.strip().replace(",", "")
                try:
                    if key.strip() == "Number of reads":
                        out["reads"] = float(val)
                    elif key.strip() == "Mean read length":
                        out["read_length"] = float(val)
                except ValueError:
                    pass
    except OSError:
        pass
    return out


def find_nanostats(sample_dir: str, nanoplot_dir: Optional[str], sample_name: str,
                   barcode: str) -> Optional[str]:
    """NanoStats.txt inside the sample folder, else under nanoplot_dir in a path that
    contains the sample name (preferred) or the barcode token."""
    found = _first_file(sample_dir, ["*NanoStats.txt"])
    if found or not nanoplot_dir or not os.path.isdir(nanoplot_dir):
        return found
    cands = [f for f in glob.glob(os.path.join(nanoplot_dir, "**", "*NanoStats.txt"), recursive=True)
             if os.path.getsize(f) > 0]
    for token in (sample_name, barcode):
        hits = sorted(f for f in cands if token and token.lower() in
                      os.path.relpath(f, nanoplot_dir).lower())
        if token == barcode and len(hits) > 1:
            continue   # barcode alone is ambiguous across runs
        if hits:
            return hits[0]
    return None


def open_text(path):
    """Open a text file for reading, gzip/bgzip-compressed or not (detected from the magic bytes)."""
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "r", encoding="utf-8")


# per-sample files that abotyper-aggregate reads, by --samples column -> the file suffix it looks for
SAMPLE_FILE_COLUMNS = {
    "phenotype": "ABOPhenotype.txt",          # abotyper-predict-phenotype
    "haplotypes": "Haplotypes.tsv",           # abotyper-clair2haplotypes
    "evidence_status": "evidence_status.tsv", # abotyper-pileup2metrics --status-out
    "concordance": "concordance.tsv",         # abotyper-pileup2metrics --concordance-out
    "coverage": "coverage.txt",               # samtools coverage
    "stats": "stats",                         # samtools stats
    "flagstat": "flagstat",                   # samtools flagstat
    "nanostats": "NanoStats.txt",             # NanoPlot
}
SAMPLE_TABLE_COLUMNS = ["id", "control", "sample", "barcode", *SAMPLE_FILE_COLUMNS]


def read_samples(path: Optional[str]) -> Dict[str, Dict[str, object]]:
    """The sample table from the samplesheet: {sample id: {"control", "sample", "barcode", "files"}}.

    Tab-separated, one row per sample id. Two layouts are accepted:

    - with a header row starting with ``id``: named columns from SAMPLE_TABLE_COLUMNS, in any order. The file
      columns (phenotype, haplotypes, ...) give each sample's files directly, so no sample folders are needed;
      empty cells are files the sample does not have.
    - without a header: id, control type and, optionally, the sample name and barcode to report.

    Control is 'negative', 'positive' or empty. When name and barcode are both empty, they are taken from the
    sample id (see ABOReportParser.extract_sample_info).
    """
    samples: Dict[str, Dict[str, object]] = {}
    if not path:
        return samples
    with open(path, encoding="utf-8") as fh:
        lines = [(n, ln.rstrip("\r\n")) for n, ln in enumerate(fh, start=1) if ln.strip()]
    header = None
    if lines and lines[0][1].split("\t")[0].strip().lower() == "id":
        header = [h.strip().lower() for h in lines[0][1].split("\t")]
        unknown = [h for h in header if h not in SAMPLE_TABLE_COLUMNS]
        if unknown:
            raise ValueError(f"{path} line 1: unknown column(s) {', '.join(unknown)}; "
                             f"expected {', '.join(SAMPLE_TABLE_COLUMNS)}")
        lines = lines[1:]
    for n, line in lines:
        parts = [p.strip() for p in line.split("\t")]
        if header:
            row = dict(zip(header, parts + [""] * (len(header) - len(parts))))
        else:
            row = dict(zip(["id", "control", "sample", "barcode"], parts + ["", "", ""]))
        sample_id, kind = row.get("id", ""), row.get("control", "").lower()
        name, barcode = row.get("sample", ""), row.get("barcode", "")
        if not sample_id:
            raise ValueError(f"{path} line {n}: empty sample id")
        if kind not in ("", "negative", "positive"):
            raise ValueError(f"{path} line {n}: control must be empty, 'negative' or 'positive', got '{kind}'")
        if barcode and not re.fullmatch(r"(barcode)?\d+", barcode, re.IGNORECASE):
            raise ValueError(f"{path} line {n}: barcode must be a number or barcode<number>, got '{barcode}'")
        files = {c: row[c] for c in SAMPLE_FILE_COLUMNS if row.get(c)}
        samples[sample_id] = {"control": kind, "sample": name, "barcode": barcode, "files": files}
    return samples


def stage_sample_files(samples: Dict[str, Dict[str, object]], root: str) -> int:
    """Lay out the files named in a sample table as <root>/<id>/combined/<id>.<suffix> (symlinks, else copies), the
    folder layout abotyper-aggregate reads. Returns the number of samples staged."""
    import shutil
    n = 0
    for sid, row in samples.items():
        files = row.get("files") or {}
        if not files:
            continue
        d = os.path.join(root, sid, "combined")
        os.makedirs(d, exist_ok=True)
        for col, src in files.items():
            if not os.path.isfile(src):
                raise ValueError(f"sample {sid}: {col} file '{src}' does not exist")
            dst = os.path.join(d, f"{sid}.{SAMPLE_FILE_COLUMNS[col]}")
            try:
                os.symlink(os.path.abspath(src), dst)
            except OSError:
                shutil.copyfile(src, dst)
        n += 1
    return n


def read_controls(path: Optional[str]) -> Dict[str, str]:
    """Controls declared in the samplesheet: {sample id: 'negative' | 'positive'} (see read_samples)."""
    return {sid: row["control"] for sid, row in read_samples(path).items() if row["control"]}
