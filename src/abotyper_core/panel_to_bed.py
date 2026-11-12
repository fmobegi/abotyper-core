"""BED file of the diagnostic positions of the ABO variant panel, for Clair3's --bed_fn and bcftools view -R
(command abotyper-panel-to-bed, also installed as panel_to_bed.py).

Every calibrated panel row becomes one BED interval anchored on its
resolved_position() (amplicon_pos, falling back to ng006669_2_pos / legacy
coordinates -- see abo_panel.VariantMarker.resolved_position). Plain SNP
rows get an exact 1bp interval since their position never shifts.

Indel, homopolymer, and multi-offset intronic-splice rows get padded
(--pad, default 15bp) on both sides. This matters because Clair3 (and
variant callers generally) anchor an indel candidate wherever the
event's left-flank base sits, not at the panel's nominal cDNA-derived
position -- empirically confirmed in this pipeline for both o1_marker
(c.261delG: called at genomic pos-1) and a1_a2_1061del (c.1061delC: same
offset). A 1bp-exact BED interval for these rows excludes the position
Clair3 actually calls at, which silently reproduces the exact detection
failure --gvcf/padding is meant to fix. Plain SNP rows are intentionally
left unpadded: several panel positions sit only 2-5bp apart, and padding
every row would let a neighbouring marker's own evidence bleed into an
unrelated marker's window.

Row order in the output is sorted by start coordinate (not panel
declaration order), which is what Clair3/bcftools expect from a BED file.

Part of abotyper-core (https://github.com/fmobegi/abotyper-core), the ABO typing logic of nf-core/abotyper.
Author: Fredrick M. Mobegi; contributors: Benedict Matern, Mathijs Groeneweg. MIT License (see LICENSE).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from .panel import load_panel, VariantMarker, default_panel_path

def marker_to_interval(marker: VariantMarker, pad: int) -> tuple[int, int]:
    """Return a 0-based half-open (start, end) BED interval for one marker."""
    pos = marker.resolved_position()
    if pos is None:
        raise ValueError(f"{marker.id} has no resolved position (not calibrated)")
    if marker.needs_window():
        lo, hi = marker.search_bounds(pad)
        start, end = max(0, lo - 1), hi
    else:
        start, end = pos - 1, pos
    return start, end


def build_bed_rows(panel, chrom: str, pad: int) -> list[tuple[int, int, str]]:
    rows = []
    for marker in panel.calibrated():
        start, end = marker_to_interval(marker, pad)
        rows.append((start, end, marker.id))
    rows.sort(key=lambda r: (r[0], r[1]))
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Generate a BED file of ABO panel diagnostic positions "
                    "(for Clair3 --bed_fn / bcftools view -R).",
    )
    ap.add_argument("--panel", default=str(default_panel_path()),
                     help="Path to the variant panel (.yaml/.json/.csv/.tsv); default: the bundled panel")
    ap.add_argument("--output", required=True, help="Output BED path")
    ap.add_argument("--chrom", default="NG_006669.2",
                     help="Reference contig name to use in the BED (default: %(default)s, "
                          "matching the pipeline's single combined ABO reference)")
    ap.add_argument("--pad", type=int, default=15,
                     help="Bases of padding applied to each side of indel/homopolymer/"
                          "multi-offset rows (default: %(default)s)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args()

    panel = load_panel(args.panel)
    rows = build_bed_rows(panel, args.chrom, args.pad)

    out_path = Path(args.output)
    with open(out_path, "w", encoding="utf-8", newline="\n") as f:
        for start, end, name in rows:
            f.write(f"{args.chrom}\t{start}\t{end}\t{name}\n")

    n_padded = sum(1 for m in panel.calibrated() if m.needs_window())
    print(f"Wrote {len(rows)} BED intervals to {out_path} "
          f"({n_padded} padded ±{args.pad}bp, {len(rows) - n_padded} exact 1bp SNP rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
