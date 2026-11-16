# abotyper-core

A Python package to assign ABO genotypes and alleles from variant calls using the ISBT database panel.

The ABO blood group typing logic from [nf-core/abotyper](https://github.com/nf-core/abotyper), packaged as a Python
library and a set of command-line tools.

It loads and validates the ABO variant panel, turns Clair3 calls and BAM pileups into per-position evidence,
resolves phenotype, genotype and ISBT allele names, and writes the per-sample and aggregated result tables.

> **Research use only.** Not validated for clinical diagnosis or transfusion decisions.

## Install

```bash
pip install abotyper-core                     # PyPI
conda install -c conda-forge -c bioconda abotyper-core   # Bioconda
```

Containers: `quay.io/biocontainers/abotyper-core` (BioContainers, built from the Bioconda package).

Python 3.11 or newer. Dependencies: pandas, numpy, pyyaml, xlsxwriter and openpyxl.

## Commands

| Command | What it does |
| --- | --- |
| `abotyper-panel` | Validate and summarise a variant panel; `--check-isbt` compares it with the ISBT allele table |
| `abotyper-panel-to-bed` | BED of panel positions for Clair3 `--bed_fn` and `bcftools view -R` |
| `abotyper-clair2metrics` | Per-position metrics from Clair3 VCFs |
| `abotyper-clair2haplotypes` | Phased haplotypes from Clair3 VCFs |
| `abotyper-pileup2metrics` | Per-position metrics from BAM pileups (plain or bgzipped `samtools mpileup` output) |
| `abotyper-predict-phenotype` | Phenotype, genotype and ISBT allele calls per sample |
| `abotyper-aggregate` | Aggregate per-sample results into the final report (xlsx, csv) |
| `abotyper-calibrate-panel` | Fill amplicon coordinates into a panel from a reference |
| `abotyper-rename-samples` | Rename samples from a mapping sheet |
| `abotyper-readtyper` | Read-level allele-pair typing with haplotype consensus, allele-balance and contamination checks (complements the per-position calls) |

Each command takes `--help`. Commands that read the panel default to the bundled
`abotyper_core/data/abo_variant_panel.yaml`; pass `--panel` to use another one.

### Read-level typing

```bash
abotyper-readtyper --bam S1.bam --ref NG_006669v2.fasta --sample S1 \
    --out S1.ReadTyping.tsv --haplotypes-out S1.ReadTypingHaplotypes.tsv
```

Each long read spans the exon 2-7 amplicon, so each read is a haplotype. `abotyper-readtyper` scores every read
against every ISBT ABO allele and calls the pair that best explains all of them. It needs only the BAM (aligned to
NG_006669.2) and that FASTA.

`S1.ReadTyping.tsv` has one row: the two ISBT alleles, genotype and phenotype, status (`PASS`, `REVIEW` or
`NOCALL`), the likelihood margin and next-best pair, allele balance (`reads_a_b`, `minor_allele_share`), the share
of foreign reads and of chimeras, novel coding variants, `haplotype_mismatches` (known sites where a haplotype
differs from its allele, with read support), `flags` (the reasons for REVIEW) and `notes` (information that does not
change the status). `S1.ReadTypingHaplotypes.tsv` has one row per haplotype: allele, reads, the ISBT alleles its
consensus matches exactly, the nearest allele and mismatches, `allele_mismatches` and novel variants.

A haplotype that differs from its ISBT allele at known sites goes to REVIEW as a possible new allele, unless it is
in the bundled table of recurring haplotypes that ISBT does not list (`abo_unlisted_haplotypes.tsv`): those are
reported under the table's name with a note. `--unlisted <file>` uses another table, `--unlisted none` none.

`abotyper-aggregate` picks up `<id>/read_typing/<id>.ReadTyping.tsv` (or the `readtyping` column of `--samples`) and
makes the read-level call the primary call, with the per-position call as the second line of evidence.

### Checking the panel against ISBT

```bash
abotyper-panel --check-isbt                          # fetch the live ISBT ABO allele table
abotyper-panel --isbt-source ABO_allele_table.xlsx   # or the file the site's EXPORT button downloads
```

The check reads the ABO alleles from the [ISBT blood group allele database](https://blooddatabase.isbtweb.org/system/ABO)
and compares them with the panel and its coverage table (`abo_isbt_v9_coverage.tsv`, which records how the panel
handles every ISBT DNA change). It lists suggestions; it never edits the panel:

| Kind | Meaning |
| --- | --- |
| `new_change` | a DNA change in a current allele that neither the panel nor the coverage table knows (exon 1 changes are flagged as assay gaps) |
| `new_allele` | an allele not in the baseline snapshot, and whether its changes are already handled |
| `changed_allele` | an allele whose changes or phenotype differ from the baseline |
| `retired_allele` | a baseline allele that is gone or marked obsolete |
| `panel_change_not_in_isbt` | a panel row whose change no current allele carries |
| `release` | the ISBT release differs from the one the panel was last checked against |

Alleles are compared with a baseline snapshot, by default the bundled `isbt_abo_alleles.tsv` (ISBT release 21,
2026-11-01). After reviewing and updating the panel, save a new baseline with `--isbt-snapshot`. Use
`--isbt-report` to write the suggestions as a TSV and `--fail-on-changes` (exit status 3) in scheduled jobs; the
repository's monthly `ISBT check` workflow opens an issue when there is something to review.

The ISBT site has no documented allele API: the check reads the allele records the ABO page embeds, and fails
with a clear message if the page layout changes. `--isbt-source` with an export file always works.

## Library use

```python
from abotyper_core import load_panel, default_panel_path

panel = load_panel(default_panel_path())
for marker in panel.calibrated():
    print(marker.id, marker.resolved_position())
```

## Development

```bash
pip install -e ".[test]"
ruff check src tests
pytest
```

The golden tests run every step on ten fixture samples and compare the outputs byte for byte (xlsx cell by cell)
with `tests/golden/`. After an intended output change, regenerate them with `UPDATE_GOLDEN=1 pytest` and review the
diff. `tests/requirements-pinned.txt` lists the versions the golden files were made with.

### Branches and releases

- Open pull requests against `dev`, the default branch. CI (ruff and pytest on Python 3.11 and 3.13) must pass
  before a merge.
- `main` holds released code only. A release is a pull request from `dev` into `main`, then a GitHub release
  tagged `vX.Y.Z` from `main`, which publishes to PyPI. Urgent fixes can go from a `hotfix/*` branch.

## Citation

If you use this package, please cite [nf-core/abotyper](https://github.com/nf-core/abotyper) and its publication:
[PubMed 40564908](https://pubmed.ncbi.nlm.nih.gov/40564908/).

## Licence

MIT. See [LICENSE](LICENSE).
