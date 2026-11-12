# Changelog

## 0.3.0 (2026-11-13)

- One version for the whole package. The modules no longer carry their own version numbers; `--version` of
  every command, the About sheet and the run QC report (`ABO_run_qc.txt`) show the abotyper-core version, under
  the label `abotyper-core` (formerly `aggregate_abo_reports.py: <script version>`).
- Every module has the same header: a docstring that names its command, ending with the package, author and
  licence lines. The per-module `__author__`/`__credits__`/`__version__` fields are removed.
- The bundled panel's `isbt_coverage_table` field points at the copy in the package.
- The `.py` command aliases (`abo_panel.py`, `panel_to_bed.py`, `clair2metrics.py`, `clair2haplotypes.py`,
  `pileup2metrics.py`, `predict_abo_phenotype.py`, `aggregate_abo_reports.py`) are removed; the package installs
  only the nine `abotyper-*` commands. nf-core/abotyper calls these commands.

## 0.2.0 (2026-11-12)

- `pileup2metrics` reads bgzipped pileups (`*.mpileup.gz`, as written by the nf-core `samtools/mpileup` module)
  as well as plain ones.
- `aggregate` (2.6.0) takes a sample table, `--samples` (the former `--controls`, which still works): sample
  folder name, control type and, optionally, the sample name and barcode to report. The pipeline passes them
  from the samplesheet, so they no longer have to be parsed from the folder name. Without them the folder name
  is parsed as before.
- The bundled panel's comments name `abotyper-calibrate-panel`.

## 0.1.0 (2026-11-12)

First release. The typing logic of nf-core/abotyper (refactor branch, `bin/` scripts) as an installable package.

- Modules: `panel`, `naming`, `inputs`, `phase`, `calling`, `report`, `aggregate`, `panel_to_bed`,
  `clair2metrics`, `clair2haplotypes`, `pileup2metrics`, `predict_phenotype`, `calibrate_panel`,
  `rename_samples`.
- Commands `abotyper-*` for each step, plus the pipeline's script names (`abo_panel.py`,
  `predict_abo_phenotype.py`, ...) as aliases.
- The ABO variant panel (v2.3.0, YAML and CSV) and the ISBT coverage table ship with the package; `--panel`
  defaults to the bundled panel.
- `abotyper-panel --check-isbt` compares the panel with the ISBT blood group allele database
  (blooddatabase.isbtweb.org) and lists suggested updates: new DNA changes, new, changed or retired alleles, and
  panel rows no current allele carries. It also reads the site's EXPORT file (`--isbt-source abo.xlsx`). The
  baseline it diffs alleles against is the bundled snapshot of ISBT release 21 (2026-11-01).
- A monthly workflow runs the check and opens an issue when it has suggestions.
- Golden tests: per-sample and aggregated outputs are byte-identical to the pipeline's `bin/` scripts.
