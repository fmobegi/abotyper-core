# Supplementary: in-depth testing and stress testing of ABO typing

Draft for the manuscript supplement (November 2026). It covers the methods, the scripts and the aggregate numbers.
Per-library calls of the validation data stay on the analysis server with the restricted data. No sample or donor
identifiers appear here.

Marked **[pending]**: results still running when this was written.

## S1. Methods compared

| Arm | Description |
| --- | --- |
| Pipeline | nf-core/abotyper typing path: Clair3 v2.0.0 (ONT model `r1041_e82_400bps_sup_v500`, region NG_006669.2:17950–24700, phasing on, gVCF), `abotyper-clair2metrics`/`clair2haplotypes`, BAM pileup (`samtools mpileup -a -B -d 0 -Q 7 -q 1` at the 135 panel positions), `abotyper-pileup2metrics` → `abotyper-predict-phenotype` → `abotyper-aggregate` |
| Pileup only (`cur`) | Same as the pipeline, without Clair3 |
| `baq`, `q0`, `q12` | Pileup only, with BAQ on (no `-B`), or base-quality filter `-Q 0` / `-Q 12` |
| Read-level typer | `abotyper-readtyper`. Each primary read is scored against every ISBT ABO allele: per-site error rates of 0.02 for SNVs, 0.20 for homopolymer indels and 0.05 for other indels. Reads that do not discriminate between alleles are set aside. The allele pair with the highest summed log-likelihood over all reads is called (candidates: alleles best for ≥5% of reads, top 40, plus the common alleles). The call is made homozygous when the minor allele is <20% of reads *and* P(≤k minor reads \| 50:50 heterozygote) < 10⁻³. Reported QC: likelihood margin to the next pair, tied pairs, allele balance, foreign-read share (reads that another allele explains by ≥2 SNVs, i.e. a third haplotype), chimera share (reads that switch between the two called alleles, at most 2 switches), and novel coding variants per haplotype from a read-backed consensus. A rare allele is kept only if it fits strictly better than every common allele. |

**Software:**
- abotyper-core 0.3.0 for the pipeline steps;
- readtyper from abotyper-core branch `experiment/read-level-typing`, commit 0c1cfec;
- samtools 1.22.1;
- Clair3 v2.0.0;
- Nextflow 26.10.6.

**Truth genotype class:** the unordered pair of A, B and O. Further truth fields:
- **A2:** presence of an A2 allele;
- **O1 dose:** the number of ABO\*O.01 alleles.

## S2. Test samples: depth series

**Data:** the 10 samples of the pipeline `test` profile (nf-core/test-datasets, branch `abotyper`).
- One library of each common genotype: AA, AB, AO1, AO2, BB, BO1, BO2, O1O1, O1O2, O2O2.
- 67–198 reads per library, nearly all spanning the whole exon 2–7 amplicon (median 6.7 kb).
- Truth comes from the sample names.

**Design:**
- Each BAM was downsampled with `samtools view --subsample` to 40, 20, 12, 8 and 5 reads, with 5 seeds per depth.
- Clair3 (the pipeline arm) ran on seeds 1–3; every other arm ran on all 5.
- Scripts: `bench_batch.sh`, `run_all.sh`, `evaluate.py`.

**Genotype class correct (%)**, n = 10 samples × seeds:

| Depth (reads) | Read-level | Pipeline | Pileup `-Q 7` | BAQ | `-Q 0` | `-Q 12` |
| ------------- | ---------- | -------- | ------------- | --- | ------ | ------- |
| Full (67–198) | 100 | 100 | 100 | 100 | 100 | 100 |
| 40 | 98 | 100 | 100 | 100 | 94 | 98 |
| 20 | 96 | 97 | 80 | 80 | 80 | 92 |
| 12 | 96 | 83 | 76 | 72 | 72 | 84 |
| 8 | 94 | 70 | 68 | 66 | 68 | 74 |
| 5 | 74 | 43 | 64 | 58 | 62 | 62 |

**O1/O2 subtype correct (%, samples carrying an O):**

| Depth (reads) | Read-level | Pipeline |
| --- | --- | --- |
| 20 | 94 | 95 |
| 12 | 97 | 81 |
| 8 | 94 | 52 |

**Pipeline wrong or no-call (%), and REVIEW:**

| Depth (reads) | Wrong | No-call | REVIEW | Wrong with PASS |
| --- | --- | --- | --- | --- |
| 40 | 0 | 0 | 40 | 0 |
| 20 | 3 | 0 | 100 | 0 |
| 12 | 7 | 10 | 100 | 0 |
| 8 | 23 | 7 | 100 | 0 |
| 5 | 27 | 30 | 100 | 0 |

**Notes on these numbers:**
- At low depth, every wrong pipeline call was flagged REVIEW; none passed QC.
- Two seeds (8 and 12 reads, seed 2, O2O2) got 0 reads from samtools subsampling. They are counted as no-calls in every arm.
- The read-level column above comes from the prototype (`readtyper.py`). The series is being re-typed with the packaged `abotyper-readtyper` (`rerun_read_pkg.sh` → `bench/summary_pkg.tsv`). **[pending]**

## S3. Test samples: contamination (stress) series

**Design:**
- Each test sample was spiked with reads from a sample of a different genotype, making up 5, 10 or 20% of the final reads (`mix_all.sh`, one replicate, seed 7).
- The truth stays the main sample's genotype.

**Results.** Each cell is genotype correct (%) / calls that were wrong but PASS:

| Spike | Read-level | Pipeline | Pileup only |
| ----- | ---------- | -------- | ----------- |
| 5% | 90 / – | 90 / 1 | 100 / 0 |
| 10% | 90 / – | 80 / 1 | 100 / 0 |
| 20% | 80 / – | 80 / 1 | 50 / 0 |

**What happened in the wrong pipeline PASS:**
- **Setup:** the O2O2 sample already carries about 12% reads that look like A1. With 5% more foreign reads, the pipeline called ABO\*A1.01/O.02.02 with QC PASS, phenotype A instead of O.
- **Read-level typer:** it made the same call. Unlike the pipeline, it flagged the sample: 14–29% of reads fit neither called allele.
- **Review rule:** sending samples above 10% of such reads to REVIEW catches every wrong O2O2 call.

The packaged rerun (S2) will add the read-level wrong-PASS column. **[pending]**

## S4. Processing settings (from S2)

- **BAQ:** never better than `-B`. At full depth it changes both A2 samples from ABO\*A2.01 to ABO\*A1.02, because it down-weights the c.1061delC homopolymer evidence. `-B` is kept.
- **`-Q 12` vs `-Q 7`:** `-Q 12` beats the current `-Q 7` below 40 reads: 92 vs 80% at 20 reads, 84 vs 76% at 12.
- **`-Q 0`:** no better than `-Q 7`.

## S5. Validation runs (QAP): full depth

**Data:**
- 129 ONT ABO amplicon libraries from 5 sequencing runs (12, 20, 24, 22 and 51 libraries).
- 87 donors: 42 sequenced in two runs, 45 in one.
- Typed reads per library: median 2,045 (IQR 551–5,057; range 0–17,280); 4 libraries under 20 reads.
- Pipeline: nf-core/abotyper, `-profile conda`, `map-ont`, defaults otherwise.
- Read-level typer: run on the pipeline's BAMs.
- Scripts: `qap/qap_setup.sh`, `qap/qap_read.sh`, `qap_eval.py`.

**Truth:**
- **Source:** the variant calls (VCFs) of each donor's earlier libraries in the validation audit (`qap/build_truth.py`). This is per donor: genotype class, A2, O1 dose and B dose.
- **Coverage:** 82 donors have a truth genotype (AO 34, OO 23, AB 13, BO 7, AA 4, BB 1); 33 carry A2.
- **Consistency:** for 80 of the 82 donors, all earlier libraries agree; for 2 they don't. 5 donors have no earlier VCF-based call.
- **Libraries scored:** 123 of the 129 have truth.

**Results:**

| | Pipeline | Read-level |
| --- | --- | --- |
| Genotype class correct | 121/123 (98.4%) | 121/123 (98.4%) |
| Wrong (all with PASS) | 2 | 2 |
| No-call | 2 | 3 |
| REVIEW | 8/129 | 3/129 |
| A2 correct (among correct class) | 97.5% | 96.7% |
| O1 dose correct (among correct class) | 100% | 100% |
| Donors whose 2 libraries agree | 40/42 | 40/42 |

**Notes on these results:**
- **Agreement:** the two methods gave the same genotype class on all 126 libraries both of them called.
- **The two wrong calls:**
  - Both are AB donors with one library typed A2.01/O.01.01 with QC PASS, the same call by both methods (one of these libraries has 7,088 reads).
  - Each donor's other library, in a different run, was AB.
  - The audit truth is also inconsistent for both donors.
  - This pattern points to sample or barcode mix-ups at the bench, not to typing errors. Without these two libraries, both methods are 121/121.
- **No-calls:** two libraries had 0 reads. One had 3 informative reads; the pipeline called it with REVIEW.
- **Novel variant:** the read-level typer flagged c.468G>A on an O.01.02 haplotype. It was seen in the same donor in two runs and is not in the ISBT table (codon 156, likely synonymous).
- **Replicate check:** the pipeline's replicate-concordance check did not fire on the mix-ups. In this test the library names carried the run prefix, so the check did not see them as replicates. With the same sample ID across runs it flags exactly this case.

## S6. Validation runs (QAP): depth series

**Design:**
- Each QAP library was downsampled to 8, 12 and 20 reads, with seeds 1 and 2.
- Only reads overlapping the typed region NG_006669.2:18064–24544 were counted and kept. Many QAP reads align partly outside it, so subsampling by total mapped reads under-samples the amplicon.
- Each subsample was typed by the pipeline arm (Clair3 + pileup) and by the packaged read-level typer.
- Scripts: `qap/qap_depth_batch.sh`, `qap/qap_depth_all.sh`, scored with `qap_eval.py`.

**Results:** **[pending]**

## S7. Alleles the exon 2–7 assay cannot resolve

| ISBT allele | Defining changes | Reported as |
| --- | --- | --- |
| ABO\*AW.13 | c.2T>C (exon 1) | ABO\*A1.01 |
| ABO\*AW.16 | c.1A>G (exon 1), c.467C>T, c.1061del | ABO\*A2.01 |
| ABO\*O.16 | c.106G>T, c.188G>A, c.204-98\_\*3252del | the other allele may appear homozygous |

See the pipeline's `docs/intended_use.md`.

## Caveats

- **Test samples:** one library per genotype, from one run. The S2–S4 numbers show direction; they are not performance claims.
- **QAP truth:** it comes from earlier libraries of the same donors, typed with the earlier pipeline's VCFs. It is not an independent serological or Sanger reference.
- **Software versions:** the read-level typer is a development version (not released). The pipeline steps used abotyper-core 0.3.0. Version 0.3.1 changes only the aggregation inputs and the negative-control QC, not typing.
