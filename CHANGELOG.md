# Changelog

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
