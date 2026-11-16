"""Two lines of evidence in one report: the read-level typer (A, abotyper-readtyper) and the per-position caller
(B, Clair3 + BAM pileup -> predict-phenotype).

A is the primary call when it makes one. B is the secondary line: it may have run on every sample (`always`), only on
the samples A could not settle (`review`), or not at all (`off`). Per sample:

    A called, B absent              -> A's call and status
    A called, B agrees              -> A's call (its full ISBT names), A's status; B's own QC kept in PositionQC
    A called, B disagrees / no call -> A's call, REVIEW ("methods disagree")
    A no-call, B called             -> B's call, REVIEW ("read-level no-call")
    A no-call, B absent             -> Unknown, REVIEW (a clean negative control passes)
    A not run                       -> B's call, unchanged

Controls are checked last, on the final call: a negative control passes when it has no more reads than --min-depth
and is flagged as possible contamination otherwise; a positive control with an expected genotype (samplesheet column
`expected`) must match it. After all samples are merged, a sample with foreign reads (a third haplotype) or a skewed
allele balance is matched against the alleles called in the other samples of the run: a match names a possible
source of the contamination. Novel coding variants of every haplotype are written to ABO_novel_variants.tsv and a
read-level table to the MultiQC report.
"""

from __future__ import annotations

import csv
import glob
import os
import re
from collections import Counter
from typing import Dict, List, Optional

import pandas as pd

READ_TYPING_DIRNAME = "read_typing"        # per_sample_processing/<id>/read_typing/<id>.ReadTyping.tsv
READ_TYPING_SUFFIX = ".ReadTyping.tsv"
READ_TYPING_HAP_SUFFIX = ".ReadTypingHaplotypes.tsv"
READ_TYPING_MAX_UNEXPLAINED = 0.10         # foreign-read share that sends a sample to REVIEW
EVIDENCE_COLUMNS = ["ReadTyping", "PrimaryMethod", "MethodAgreement", "Margin", "AlleleBalance", "PositionQC"]
NOVEL_VARIANTS_FILE = "ABO_novel_variants.tsv"
READ_TYPING_MQC_FILE = "ABO_readtyping_mqc.tsv"


def _blank(v) -> bool:
    return v is None or (isinstance(v, float) and pd.isna(v)) or str(v).strip() in ("", "nan", "None")


def _float(v, default=0.0) -> float:
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def genotype_class(text: str) -> str:
    """Genotype class (AO, BB, ...) from a class or an ISBT allele pair; '' when it has none."""
    t = str(text or "").strip()
    if "*" in t:
        t = "".join(re.findall(r"ABO\*(A|B|O)", t))
    if not t or not set(t) <= set("ABO") or len(t) != 2:
        return ""
    return "".join(sorted(t, key="ABO".index))


def phenotype_class(genotype: str, fallback: str = "") -> str:
    """Serological class reported in the Phenotype column (A, B, AB, O) from a genotype class."""
    g = genotype_class(genotype)
    if not g:
        return fallback or "Unknown"
    letters = "".join(c for c in "AB" if c in g)
    return letters or "O"


class ReconcileMixin:
    """Mixed into ABOReportParser: finds read-level results, reconciles them with the per-position call, and writes
    the run-level read-level outputs. Uses self.min_depth, self.sample_table, self.sample_metrics, self.qc_records."""

    # ------------------------------------------------------------------ files
    def _read_typing_paths(self, sample_dir: str):
        """(summary, haplotypes) of a sample: read_typing/ first, combined/ for older layouts."""
        for sub in (READ_TYPING_DIRNAME, "combined"):
            s = glob.glob(os.path.join(sample_dir, sub, f"*{READ_TYPING_SUFFIX}"))
            if s:
                h = glob.glob(os.path.join(sample_dir, sub, f"*{READ_TYPING_HAP_SUFFIX}"))
                return s[0], (h[0] if h else None)
        return None, None

    def _detect_read_typing(self) -> bool:
        return any(glob.glob(os.path.join(str(self.input_dir), "*", sub, f"*{READ_TYPING_SUFFIX}"))
                   for sub in (READ_TYPING_DIRNAME, "combined"))

    @staticmethod
    def _load_read_typing(path: Optional[str]):
        if not path:
            return None, ""
        try:
            with open(path, newline="") as fh:
                return next(csv.DictReader(fh, delimiter="\t")), ""
        except (OSError, StopIteration, csv.Error) as e:
            return None, f"unreadable ({e})"

    @staticmethod
    def _load_haplotypes(path: Optional[str]) -> List[dict]:
        if not path:
            return []
        try:
            with open(path, newline="") as fh:
                return list(csv.DictReader(fh, delimiter="\t"))
        except (OSError, csv.Error):
            return []

    # ------------------------------------------------------------------ per sample
    def _row_from_read_typing(self, result_df, rt: dict):
        """Fill the Result/Notes cells of a sample the per-position caller did not type."""
        geno = rt.get("genotype", "")
        called = geno and geno != "NOCALL"
        pair = f"{rt.get('allele_1', '')}/{rt.get('allele_2', '')}" if called else ""
        result_df.at[0, ("Result", "Genotype")] = genotype_class(geno) or (geno if called else "Unknown")
        result_df.at[0, ("Result", "Phenotype")] = phenotype_class(geno, rt.get("phenotype", "")) if called \
            else "Unknown"
        result_df.at[0, ("Result", "ExtendedGenotype")] = pair
        pheno = rt.get("phenotype", "")
        result_df.at[0, ("Notes", "ASubtype")] = pheno if pheno.startswith("A") and called else ""
        result_df.at[0, ("Notes", "ReadReliability")] = f"{rt.get('reads_typed', 0)} informative reads (read-level)"
        result_df.at[0, ("Notes", "DepthCheck")] = (f"read-level: {rt.get('reads_aligned', 0)} reads aligned, "
                                                    f"{rt.get('reads_typed', 0)} informative")
        result_df.at[0, ("Notes", "PhaseConfidence")] = "read-level (each read is a haplotype)"
        result_df.at[0, ("Notes", "BwSubtype")] = ""
        result_df.at[0, ("Notes", "A2Trace")] = "n/a (read-level)"
        result_df.at[0, ("Notes", "EvidenceBasis")] = "read-level only"
        result_df.at[0, ("Notes", "QCStatus")] = rt.get("status", "") if called else "REVIEW"
        result_df.at[0, ("Notes", "QCFlags")] = rt.get("flags", "")

    def _reconcile(self, result_df, sample_dir: str, sample_id: str, key, has_position_call: bool):
        """Combine read-level (A) and per-position (B) evidence for one sample; see the module docstring."""
        rt_path, hap_path = self._read_typing_paths(sample_dir)
        rt, err = self._load_read_typing(rt_path)
        flags: List[str] = []
        put = lambda col, v: result_df.__setitem__(("Notes", col), v)   # noqa: E731

        b_geno = str(result_df.at[0, ("Result", "Genotype")]) if has_position_call else ""
        b_ext = str(result_df.at[0, ("Result", "ExtendedGenotype")]) if has_position_call else ""
        b_status = str(result_df.at[0, ("Notes", "QCStatus")]) if has_position_call else ""
        b_flags = str(result_df.at[0, ("Notes", "QCFlags")] or "") if has_position_call else ""
        b_called = has_position_call and b_geno not in ("", "nan", "Unknown")
        put("PositionQC", (f"{b_status}" + (f": {b_flags}" if b_flags and b_flags != "nan" else ""))
            if has_position_call else "not run")

        if rt is None:
            put("ReadTyping", err or "not run")
            put("PrimaryMethod", "position" if b_called else "none")
            put("MethodAgreement", "n/a")
            put("Margin", "")
            put("AlleleBalance", "")
            if not has_position_call:
                result_df.at[0, ("Result", "Genotype")] = "Unknown"
                result_df.at[0, ("Result", "Phenotype")] = "Unknown"
                put("QCStatus", "REVIEW")
                put("QCFlags", "no typing: neither read-level nor position-based results")
            return None

        geno, status = rt.get("genotype", ""), rt.get("status", "")
        a_called = bool(geno) and geno != "NOCALL"
        pair = f"{rt.get('allele_1', '')}/{rt.get('allele_2', '')}" if a_called else ""
        put("ReadTyping", (f"{geno} {pair} [{status}; {rt.get('reads_typed', '')} reads]"
                           + (f" -- note: {rt['notes']}" if rt.get("notes") else "")) if a_called
            else f"NOCALL [{rt.get('flags', '')}]")
        put("Margin", rt.get("margin", ""))
        put("AlleleBalance", f"{rt.get('reads_a_b', '')} (minor {rt.get('minor_allele_share', '')})"
            if a_called and rt.get("allele_1") != rt.get("allele_2") else "")
        a_flags = [f for f in str(rt.get("flags") or "").split("; ") if f]

        if a_called:
            agree = b_called and genotype_class(geno) == genotype_class(b_geno) and \
                self._typings_compatible(b_ext, pair)
            if not has_position_call:
                self._row_from_read_typing(result_df, rt)
                put("MethodAgreement", "position caller not run")
            else:
                b_pheno = str(result_df.at[0, ("Result", "Phenotype")])
                self._row_from_read_typing(result_df, rt)
                if agree:
                    result_df.at[0, ("Result", "Phenotype")] = b_pheno   # same class; keep B's annotated wording
                    put("EvidenceBasis", "read-level + position (agree)")
                    put("MethodAgreement", "agree")
                elif b_called:
                    put("EvidenceBasis", "read-level (position caller disagrees)")
                    put("MethodAgreement", f"disagree: position {b_geno} {b_ext}")
                    flags.append(f"methods disagree: read-level {genotype_class(geno) or geno} {pair}, "
                                 f"position {b_geno} {b_ext}")
                else:
                    put("EvidenceBasis", "read-level (position caller: no call)")
                    put("MethodAgreement", "position caller no call")
                    flags.append("position caller made no call")
            put("PrimaryMethod", "read-level")
            flags = a_flags + flags
            if any(a.startswith("ABO*AEL") for a in (rt.get("allele_1", ""), rt.get("allele_2", ""))):
                flags.append("Ael allele (c.804 homopolymer indel): very weak or silent A expression, so serology "
                             "may read O; ONT reads are error-prone in this 7-G run, confirm the indel")
            put("QCStatus", "REVIEW" if flags or status != "PASS" else "PASS")
            put("QCFlags", "; ".join(flags))
        else:
            put("AlleleBalance", "")
            if b_called:
                put("PrimaryMethod", "position")
                put("MethodAgreement", "read-level no-call")
                old = b_flags if b_flags != "nan" else ""
                put("QCFlags", "; ".join(x for x in (old, f"read-level typing made no call ({rt.get('flags', '')})")
                                         if x))
                put("QCStatus", "REVIEW")
            else:
                if not has_position_call:
                    self._row_from_read_typing(result_df, rt)
                result_df.at[0, ("Result", "Genotype")] = "Unknown"
                result_df.at[0, ("Result", "Phenotype")] = "Unknown"
                put("PrimaryMethod", "none")
                put("MethodAgreement", "no call by either method" if has_position_call else "read-level no-call")
                put("QCStatus", "REVIEW")
                put("QCFlags", f"no call: read-level ({rt.get('flags', '')})"
                    + ("" if has_position_call else ", position caller not run"))

        m = self.sample_metrics.setdefault(key, {})
        m["read_typing"] = {k: rt.get(k, "") for k in ("genotype", "allele_1", "allele_2", "status", "reads_aligned",
                                                       "reads_typed", "margin", "minor_allele_share", "reads_a_b",
                                                       "unexplained_share", "chimera_share", "third_allele",
                                                       "novel_variants")}
        m["read_typing_haplotypes"] = self._load_haplotypes(hap_path)
        if not has_position_call:
            self.qc_records.append({"sample": sample_id, "genotype": str(result_df.at[0, ("Result", "Genotype")]),
                                    "ext": pair, "a2trace": "n/a", "evidence": "read-level only", "s261": None,
                                    "reads_uniform": False})
        return rt

    def _apply_controls(self, result_df, rt: Optional[dict], has_position_call: bool, sample_folder: str):
        """Negative and positive controls, judged on the final call (B's own control logic stays for B-only rows)."""
        kind = self._current_control
        sid = str(result_df.at[0, ("", "Sequencing_ID")])
        if not kind and getattr(self, "negative_control_re", None) and self.negative_control_re.search(sid):
            kind = "negative"
        if not kind and getattr(self, "positive_control_re", None) and self.positive_control_re.search(sid):
            kind = "positive"
        if not kind:
            return
        put = lambda col, v: result_df.__setitem__(("Notes", col), v)   # noqa: E731
        if kind == "negative":
            if rt is None and has_position_call:
                return                      # B-only row: the per-position caller's own control logic stands
            reads = int(_float(rt.get("reads_aligned"), 0)) if rt else 0
            b_clean = getattr(self, "_current_negative_clean", False) if has_position_call else True
            if reads <= self.min_depth and b_clean:
                put("DepthCheck", f"Negative control: no reads over the ABO region beyond {self.min_depth}, as expected"
                    + (f" ({reads} read-level reads)" if rt else ""))
                put("QCStatus", "PASS")
                put("QCFlags", "")
            elif reads > self.min_depth:
                old = str(result_df.at[0, ("Notes", "QCFlags")] or "")
                msg = f"negative control has reads (read-level {reads} > {self.min_depth}): possible contamination"
                if msg not in old:
                    put("QCFlags", "; ".join(x for x in (old if old != "nan" else "", msg) if x))
                put("QCStatus", "REVIEW")
            return
        expected = (self.sample_table.get(sample_folder, {}) or {}).get("expected", "")
        if not expected:
            return
        got = str(result_df.at[0, ("Result", "Genotype")])
        old = str(result_df.at[0, ("Notes", "QCFlags")] or "")
        old = "" if old == "nan" else old
        if genotype_class(expected) and genotype_class(expected) == genotype_class(got):
            put("DepthCheck", f"Positive control: expected {expected}, typed {got}")
        else:
            put("QCFlags", "; ".join(x for x in (old, f"positive control: expected {expected}, typed {got}") if x))
            put("QCStatus", "REVIEW")

    # ------------------------------------------------------------------ run level
    def _apply_cross_sample_checks(self, final_df):
        """Name possible sources of contamination: another sample of the run whose called allele is the foreign
        (third) allele, or the minor allele of a skewed heterozygote."""
        self.contamination_sources = []
        rows = []
        for i in final_df.index:
            try:
                key = (str(final_df.at[i, ("", "Sequencing_ID")]), int(final_df.at[i, ("", "Barcode")]))
            except (TypeError, ValueError):
                continue
            rows.append((i, key, self.sample_metrics.get(key, {}).get("read_typing") or {}))
        carriers: Dict[str, List[str]] = {}
        for i, key, _ in rows:
            ext = str(final_df.at[i, ("Result", "ExtendedGenotype")])
            for a in set(ext.split("/")):
                if a.startswith("ABO*"):
                    carriers.setdefault(a, []).append(f"{key[0]} (barcode{key[1]:02d})")
        for i, key, rt in rows:
            suspects = []
            if _float(rt.get("unexplained_share")) >= READ_TYPING_MAX_UNEXPLAINED and rt.get("third_allele"):
                suspects.append(rt["third_allele"])
            flags = str(final_df.at[i, ("Notes", "QCFlags")] or "")
            if "allele imbalance" in flags and rt.get("allele_1") and rt.get("allele_1") != rt.get("allele_2"):
                to = [int(x) for x in str(rt.get("reads_a_b", "0/0")).split("/") if x.isdigit()] or [0, 0]
                suspects.append(rt["allele_1"] if len(to) == 2 and to[0] < to[1] else rt["allele_2"])
            me = f"{key[0]} (barcode{key[1]:02d})"
            for allele in suspects:
                src = [s for s in carriers.get(allele, []) if s != me]
                if not src:
                    continue
                shown = ", ".join(src[:5]) + (f" and {len(src) - 5} more" if len(src) > 5 else "")
                msg = f"possible contamination source in this run ({allele}): {shown}"
                final_df.at[i, ("Notes", "QCFlags")] = (flags + "; " if flags and flags != "nan" else "") + msg
                final_df.at[i, ("Notes", "QCStatus")] = "REVIEW"
                flags = str(final_df.at[i, ("Notes", "QCFlags")])
                self.contamination_sources.append((me, allele, shown))

    def write_novel_variants(self, final_df, path: str = NOVEL_VARIANTS_FILE):
        """Run-level list of coding variants not in the ISBT table, one row per sample and haplotype."""
        cols = ["Sample", "Barcode", "Haplotype", "Allele", "Reads", "Novel variants", "Nearest ISBT allele",
                "Mismatches to nearest", "QC status"]
        out = []
        for _, r in final_df.iterrows():
            try:
                key = (str(r["Sequencing_ID"]), int(r["Barcode"]))
            except (TypeError, ValueError):
                continue
            for h in self.sample_metrics.get(key, {}).get("read_typing_haplotypes") or []:
                if h.get("novel_variants"):
                    out.append([key[0], f"barcode{key[1]:02d}", h.get("haplotype", ""), h.get("allele", ""),
                                h.get("reads", ""), h.get("novel_variants", ""), h.get("nearest_isbt", ""),
                                h.get("nearest_mismatches", ""), r.get("QCStatus", "")])
        try:
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write("\t".join(cols) + "\n")
                for row in out:
                    fh.write("\t".join(str(x) for x in row) + "\n")
            print(f"Novel variants written to {path} ({len(out)} haplotypes)")
        except OSError as e:
            self.output_errors.append(f"{path}: {e}")
        self.novel_variant_rows = len(out)

    def write_read_typing_mqc(self, final_df, path: str = READ_TYPING_MQC_FILE):
        cols = ["Sample", "Barcode", "Read-level call", "Read-level status", "Informative reads", "Margin",
                "Minor allele share", "Foreign-read share", "Chimera share", "Primary method", "Method agreement",
                "Final QC"]
        rows = []
        for _, r in final_df.iterrows():
            try:
                key = (str(r["Sequencing_ID"]), int(r["Barcode"]))
            except (TypeError, ValueError):
                continue
            rt = self.sample_metrics.get(key, {}).get("read_typing")
            if not rt:
                continue
            call = f"{rt['allele_1']}/{rt['allele_2']}" if rt.get("genotype") not in ("", "NOCALL") else "NOCALL"
            rows.append([key[0], f"barcode{key[1]:02d}", call, rt.get("status", ""), rt.get("reads_typed", ""),
                         rt.get("margin", ""), rt.get("minor_allele_share", ""), rt.get("unexplained_share", ""),
                         rt.get("chimera_share", ""), r.get("PrimaryMethod", ""), r.get("MethodAgreement", ""),
                         r.get("QCStatus", "")])
        if not rows:
            return
        try:
            with open(path, "w", encoding="utf-8", newline="\n") as fh:
                fh.write("# id: 'abo_read_typing'\n# section_name: 'ABO read-level typing'\n"
                         "# description: 'Read-level typing (abotyper-readtyper): the allele pair that best explains "
                         "the reads, its likelihood margin, allele balance, foreign and chimeric reads, and agreement "
                         "with the position-based call.'\n# plot_type: 'table'\n")
                fh.write("\t".join(cols) + "\n")
                for row in rows:
                    fh.write("\t".join("" if x is None else str(x) for x in row) + "\n")
            print(f"Read-level MultiQC table written to {path}")
        except OSError as e:
            self.output_errors.append(f"{path}: {e}")

    def read_typing_qc_lines(self, final_df) -> List[str]:
        """Run-QC section for the read-level typing."""
        lines = ["Read-level typing (abotyper-readtyper):"]
        st = Counter()
        for (sid, bc), m in self.sample_metrics.items():
            rt = m.get("read_typing")
            if rt:
                st[rt.get("status") or "?"] += 1
        lines.append("  read-level status: " + (", ".join(f"{k}={v}" for k, v in st.most_common()) or "none"))
        for col, label in (("PrimaryMethod", "primary method"), ("MethodAgreement", "method agreement")):
            if col in final_df.columns:
                vals = Counter(re.sub(r":.*", "", str(v)) for v in final_df[col])
                lines.append(f"  {label}: " + ", ".join(f"{k}={v}" for k, v in vals.most_common()))
        lines.append(f"  haplotypes with novel coding variants: {getattr(self, 'novel_variant_rows', 0)} "
                     f"(see {NOVEL_VARIANTS_FILE})")
        srcs = getattr(self, "contamination_sources", [])
        if srcs:
            lines.append("  possible contamination sources within the run:")
            for me, allele, shown in srcs:
                lines.append(f"    {me}: {allele} also called in {shown}")
        return lines
