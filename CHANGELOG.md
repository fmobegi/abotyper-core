# Changelog

## 0.4.1 (2026-11-16)

- Three more recurring haplotypes in `abo_unlisted_haplotypes.tsv`, each seen in 4-6 libraries of the
  retrospective set with at least 90% read support and confirmed by Clair3 in every library:
  `ABO*O.01.83[+c.189C>T;+c.220C>T]`, `ABO*O.01.36[-c.829G>A]` and `ABO*O.01.24[+c.220C>T]` (names say how the
  haplotype differs from the ISBT allele; they are not ISBT names). They are reported with a note instead of
  REVIEW.
- `abotyper-aggregate`: a haplotype named from that table agrees with a per-position call of the ISBT allele it was
  called as (for example `ABO*O.01.24` and `ABO*O.01.24[+c.220C>T]`, `ABO*O.01.68` and `ABO*O.01.12-like`), so it no
  longer reads as "methods disagree".

## 0.4.0 (2026-11-16)

Development versions 0.4.0.dev0 to 0.4.0.dev3 were published to TestPyPI for testing on nf-core/abotyper's
read-level typing branch.

- New command `abotyper-readtyper`: read-level typing. Every read of the exon 2-7 amplicon is scored against every
  ISBT ABO allele (error rates per SNV, homopolymer indel and other indel) and the allele pair that best explains all
  reads is called, with its likelihood margin to the best pair of another genotype class. Needs only the BAM and
  the NG_006669.2 FASTA (a pure-Python BAM reader; no Clair3, pysam or samtools).
  - At least 5 informative reads, else NOCALL; REVIEW below 10.
  - Allele balance: a heterozygote whose minor allele is under 20% of reads (and unlikely under 50:50,
    P < 10^-3) is called homozygous; under 30% it stays heterozygous and goes to REVIEW (a contaminant can pose as
    the second allele).
  - Third haplotype: reads that another allele explains much better than either called allele are split into
    chimeras of the two (PCR) and foreign reads; 10% or more foreign reads, or 12% and 5 reads or more chimeras,
    go to REVIEW.
  - Haplotype consensus: the reads of each allele are pooled. A rare allele is kept only if the consensus fits it
    better than every common allele; coding variants not in the ISBT table are reported on their haplotype; known
    sites where the consensus differs from its called allele are reported (`haplotype_mismatches`, and per
    haplotype `allele_mismatches`), and send the sample to REVIEW ("possible new allele") when backed by at least
    5 reads and 80% of them.
  - Recurring haplotypes that ISBT does not list are named from a bundled table (`abo_unlisted_haplotypes.tsv`,
    `--unlisted`) and reported with a note (`notes`) instead of a REVIEW flag. First entry: ABO*O.01.68 + c.595C>T,
    reported as `ABO*O.01.12-like` (104 of 2,776 typed libraries of a retrospective set, all confirmed by Clair3).
- `abotyper-aggregate` combines the read-level call with the per-position call: the read-level call is primary,
  a disagreement goes to REVIEW, and a sample the read-level typer cannot call keeps the per-position call with
  REVIEW. A sample with only read-level typing is reported on its own. New evidence columns (`ReadTyping`,
  `PrimaryMethod`, `MethodAgreement`, `Margin`, `AlleleBalance`, `PositionQC`), a positive control's expected
  genotype (`expected` column of `--samples`), possible contamination sources named from the other samples of the
  run, and two run-level files: `ABO_novel_variants.tsv` and `ABO_readtyping_mqc.tsv` (MultiQC). Read-level files
  are read from `<id>/read_typing/` or named in `--samples` (`readtyping`, `readtyping_haplotypes`).

## 0.3.1 (2026-11-13)

- `abotyper-aggregate` takes explicit inputs: with a header row, `--samples` can name each sample's files
  (columns `phenotype`, `haplotypes`, `evidence_status`, `concordance`, `coverage`, `stats`, `flagstat`,
  `nanostats`), and the input directory is then optional. The report is the same as with the
  `per_sample_processing/<id>/combined/` folder layout, which still works. nf-core/abotyper uses this, so its
  aggregation module no longer builds the folder layout.
- A declared negative control with no reads over the ABO region now passes QC. Before, it was sent to REVIEW as
  "No reads over ABO region (likely amplicon dropout)", although the Typing sheet said low depth was expected.

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
