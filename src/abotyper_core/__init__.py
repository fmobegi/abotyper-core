"""abotyper-core: ABO blood group typing logic from long-read amplicon data.

The panel-driven typing used by nf-core/abotyper: loading and validating the ABO variant panel, turning Clair3 and
BAM-pileup evidence into per-position calls, resolving phenotype, genotype and ISBT allele names, and writing the
result tables. Each step is also a command-line tool (see the README).

Research use only; not validated for clinical diagnosis.
"""

__version__ = "0.1.0"

from .panel import Panel, PanelError, VariantMarker, default_panel_path, load_panel  # noqa: E402

__all__ = ["Panel", "PanelError", "VariantMarker", "default_panel_path", "load_panel", "__version__"]
