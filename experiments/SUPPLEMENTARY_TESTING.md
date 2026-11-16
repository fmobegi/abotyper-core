# Supplementary: in-depth testing and stress testing of ABO typing

Draft for the manuscript supplement (November 2026). It covers the methods, the scripts and the aggregate numbers.
Per-library calls of the validation data stay on the analysis server with the restricted data. No sample or donor
identifiers appear here.


## S1. Methods compared

| Arm | Description |
| --- | --- |
| Pipeline | nf-core/abotyper typing path: Clair3 v2.0.0 (ONT model `r1041_e82_400bps_sup_v500`, region NG_006669.2:17950–24700, phasing on, gVCF), `abotyper-clair2metrics`/`clair2haplotypes`, BAM pileup (`samtools mpileup -a -B -d 0 -Q 7 -q 1` at the 135 panel positions), `abotyper-pileup2metrics` → `abotyper-predict-phenotype` → `abotyper-aggregate` |
| Pileup only (`cur`) | Same as the pipeline, without Clair3 |
| `baq`, `q0`, `q12` | Pileup only, with BAQ on (no `-B`), or base-quality filter `-Q 0` / `-Q 12` |
| Read-level typer | `abotyper-readtyper`. Each primary read is scored against every ISBT ABO allele: per-site error rates of 0.02 for SNVs, 0.20 for homopolymer indels and 0.05 for other indels. Reads that do not discriminate between alleles are set aside. The allele pair with the highest summed log-likelihood over all reads is called (candidates: alleles best for ≥5% of reads, top 40, plus the common alleles). The call is made homozygous when the minor allele is <20% of reads *and* P(≤k minor reads \| 50:50 heterozygote) < 10⁻³. Reported QC: likelihood margin to the next pair, tied pairs, allele balance, foreign-read share (reads that another allele explains by ≥2 SNVs, i.e. a third haplotype), chimera share (reads that switch between the two called alleles, at most 2 switches), and novel coding variants per haplotype from a read-backed consensus. A rare allele is kept only if it fits strictly better than every common allele. |

**Software:**
- abotyper-core 0.3.0 for the pipeline steps;
- readtyper from abotyper-core 0.4.0.dev0 (branch `experiment/read-level-typing`; development version, not released);
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
- The read-level column above comes from the prototype (`readtyper.py`).

**Packaged `abotyper-readtyper`** (abotyper-core 0.4.0.dev0, with the review rules in S3). Same BAMs, re-typed with `rerun_read_pkg.sh` and `qap/retype_all.sh`:

| Depth (reads) | Correct % | Wrong % | No-call % | REVIEW % | Wrong with PASS |
| --- | --- | --- | --- | --- | --- |
| Full | 100 | 0 | 0 | 20 | 0 |
| 40 | 98 | 2 | 0 | 8 | 0 |
| 20 | 96 | 4 | 0 | 26 | 0 |
| 12 | 94 | 0 | 6 | 46 | 0 |
| 8 | 80 | 0 | 20 | 66 | 0 |
| 5 | 38 | 0 | 62 | 36 | 0 |

The packaged typer calls nothing on fewer than 5 informative reads. Where the prototype guessed at 5–8 reads, it gives a no-call: no wrong calls at 12, 8 or 5 reads, and no wrong call with PASS at any depth.

## S3. Test samples: contamination (stress) series

**Design:**
- Each test sample was spiked with reads from a sample of a different genotype, making up 5, 10 or 20% of the final reads (`mix_all.sh`, one replicate, seed 7).
- The truth stays the main sample's genotype.

**Results.** Each cell is genotype correct (%) / calls that were wrong but PASS:

| Spike | Read-level | Pipeline | Pileup only |
| ----- | ---------- | -------- | ----------- |
| 5% | 90 / 0 | 90 / 1 | 100 / 0 |
| 10% | 90 / 0 | 80 / 1 | 100 / 0 |
| 20% | 80 / 0 | 80 / 1 | 50 / 0 |

The read-level column is the packaged typer with the review rules below. It sends 40%, 40% and 90% of the spiked samples to REVIEW.

**What happened in the wrong pipeline PASS:**
- **Setup:** the O2O2 sample already carries about 12% reads that look like A1. With 5% more foreign reads, the pipeline called ABO\*A1.01/O.02.02 with QC PASS, phenotype A instead of O.
- **Read-level typer:** it made the same call. Unlike the pipeline, it flagged the sample: 14–29% of reads fit neither called allele.
- **Review rule:** sending samples above 10% of such reads to REVIEW catches every wrong O2O2 call.

**Review rules added after this test:**
- **What went wrong:** in the packaged typer, the contaminant could pose as the second allele. O2O2 + 10% or 20% A1 was called A1/O2, and BB + 20% AA was called A1/B. Both passed: the foreign-read share stays low because the contaminant is "explained" by the called allele.
- **Rule 1, allele imbalance:** a heterozygote whose minor allele is under 30% of reads, when P(≤k \| 50:50) < 10⁻³, goes to REVIEW.
- **Rule 2, excess chimeras:** a chimera share of at least 12%, with at least 5 chimeric reads, goes to REVIEW.
- **Why these cut-offs:** in the 129 QAP libraries at full depth, heterozygotes stay at a minor share of 0.34 or more, and the chimera share stays under 0.09.
- **Simulation before the change (`qap/sim_rules.py`):** the rules sent all 3 wrong contamination calls to REVIEW and changed no other call: 126 QAP full-depth calls, 716 QAP depth-series calls and 219 test-sample depth calls.
- **Re-typing after the change:** confirmed the same.

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

**How it was scored:**
- 121 libraries with truth. The two suspected mix-up libraries (S5) are left out, because they are wrong at every depth for both methods.
- The pipeline arm used abotyper-core 0.3.0; the read-level arm used the packaged readtyper.
- Scripts: `qap_eval.py` and `qap/depth_detail.py`.

**Genotype class, per seed.** Each cell is correct / wrong / no-call:

| Target reads | Seed | Reads in subsample, median (range) | Pipeline | Read-level | Read-level PASS (all correct) |
| --- | --- | --- | --- | --- | --- |
| 20 | 1 | 21 (3–30) | 117 / 2 / 2 | 120 / 0 / 1 | 109 |
| 20 | 2 | 19 (3–35) | 120 / 1 / 0 | 120 / 0 / 1 | 101 |
| 12 | 1 | 13 (3–20) | 112 / 7 / 2 | 119 / 1 / 1 | 62 |
| 12 | 2 | 10.5 (3–26) | 113 / 5 / 3 | 113 / 5 / 3 | 54 |
| 8 | 1 | 8 (2–15) | 99 / 15 / 7 | 111 / 2 / 8 | 19 |
| 8 | 2 | 7 (1–16) | 98 / 17 / 6 | 95 / 2 / 24 | 20 |

**QC behaviour:**
- **Pipeline:** every library is REVIEW at every reduced depth (all 128), so it never makes a wrong call with PASS.
- **Read-level typer:** it passes 85–90% of libraries at about 20 reads, 45–51% at 12 and 16% at 8. Every one of those PASS calls has the right genotype class. Every wrong call it made was REVIEW.
- **Read-level no-calls:** its no-calls are libraries with fewer than 5 informative reads (minimum not met). This is why seed 2 at 8 reads, which drew fewer reads (median 7), has 24 no-calls.

**A subgroup (A1/A2), among libraries whose class was right:**
- **Read-level typer:** found A2 in 54/55 A2 carriers at 20 reads, 53/55 and 49/50 at 12, and 49/51 and 44/46 at 8. It called A2 in 2 of about 22 non-A2 A carriers at each depth. At full depth, A2 was right 96.7% of the time.
- **Pipeline:** at these depths it withholds the subgroup. It reports `ABO*A` with phenotype A (or AB) and REVIEW; at full depth A2 was right 97.5% of the time.
- **Result:** with 8–20 reads, only the read-level typer resolves A1/A2.

**Allele-pair stability:**
- **Read-level typer:** 369 of its 372 PASS calls at reduced depth give the same ISBT allele pair as its full-depth call on the same library.
- **Pipeline:** its extended genotype matches its own full-depth one in 19–26% of libraries, because it reports `ABO*A`/`ABO*O.01` without the full allele names at low depth.

**Summary:**
- **Genotype class:** both methods stay above 95% down to about 20 reads, and the read-level typer is ahead at 8–12 reads.
- **Confidence:** below full depth, the pipeline sends every sample to REVIEW. The read-level typer makes most calls at about 20 reads with PASS: 85–90% of libraries, all right, with A1/A2 and the full allele pair resolved.

## S6b. Large retrospective set: read-level typer vs the per-position audit

**Data:** 6,241 libraries, the full retrospective set. They were typed earlier with the per-position path (pileup and Clair3; the validation audit). The packaged read-level typer (0.4.0.dev0) was then run on the same BAMs.

**Libraries neither method can type:**
- 3,348 libraries (54%) have reads aligned (median about 4,900) but none over the typed sites.
- For the read-level typer, that means 0 informative reads.
- In the audit, these same libraries also had 0 pileup depth at every primary position and no genotype from the pileup or Clair3.
- So these BAMs cannot be typed by either method. This is a property of these BAMs, not a typer failure.

**The rest (2,893 libraries):**

| | n |
| --- | --- |
| Read-level PASS | 2,587 |
| Read-level REVIEW | 189 |
| Read-level no-call (1–4 informative reads) | 117 |

**Agreement with the per-position calls:**

| Compared with | Agree | of which read-level PASS | of which read-level REVIEW |
| --- | --- | --- | --- |
| Clair3 genotype (audit `vcf_geno`) | 2,528/2,536 (99.68%) | 2,490/2,494 (99.84%) | 38/42 |
| Pileup genotype (audit `pile_geno`) | 2,531/2,538 (99.72%) | 2,491/2,495 (99.84%) | 40/43 |

**A2:** among A carriers where both methods give the same class, the read-level typer and Clair3 agree on A2 status for 1,099/1,112 (98.8%).

**Disagreements: 8 libraries with a Clair3 genotype.**
- 4 are read-level REVIEW:
  - a tie between AO and OO;
  - two skewed BB/AB calls with 12–15% minor allele;
  - one call on 8 reads.
- 4 are read-level PASS, all with high margins:
  - 3 call a B(A) or cisAB allele where the per-position path reported B: ABO\*BA.03 twice, ABO\*cisAB.06 once. The panel has no rules for these rare alleles, so the read-level call may be the more specific one.
  - 1 calls A1/O where the per-position path called BO.
- All 8 are listed for review with the restricted data.

**Calls only the read-level typer made:**
- 238 libraries had no per-position genotype but a read-level call: 92 PASS and 146 REVIEW.
- These libraries are thin: median 16 informative reads, median pileup depth 12 at the primary sites.

**Rare alleles and novel variants:**
- The typer called a non-common ISBT allele on 672 haplotypes in 633 libraries. Most are O.01 subtypes: O.01.68, O.01.26, O.01.75, O.01.67.
- 26 libraries carry a coding variant not in the ISBT table: 19 distinct variants.
- Of all haplotypes, 4,783 of 4,937 (96.9%) match an ISBT allele exactly at every site they cover.

**Review cut-offs checked at this scale** (heterozygotes or libraries with at least 50 informative reads):

| Measure | n | Result | Cut-off |
| --- | --- | --- | --- |
| Minor-allele share | 1,932 heterozygotes | 0.1st percentile 0.321, none below 0.30 | 0.30 (allele-imbalance flag) |
| Chimera share | 2,491 libraries | 99th percentile 0.076, none at or above 0.12 | 0.12 (excess-chimera flag) |
| Foreign-read share | 2,491 libraries | 4 at or above 0.10 | 0.10 (contamination flag) |

The allele-imbalance and excess-chimera rules (S3) add no REVIEW anywhere in this set.

**Review of the 8 disagreements** (per-library evidence from the audit and the read-level haplotypes):
- **3 libraries (2 donors), read-level PASS:** the read-level typer names a rare allele where the per-position path reported B.
  - ABO\*BA.03 (twice, one donor, two runs): B.01 without c.703G>A. The audit's own variant list for both libraries has every B marker except c.703.
  - ABO\*cisAB.06: B.01 without c.526C>G. The audit lists every B marker except c.526.
  - Both haplotypes match their ISBT alleles exactly at 134 sites, with balanced allele splits and margins of 665–1,152.
  - The per-position panel has no rules for these alleles, so it fell back to B.
- **1 library, read-level PASS in 0.4.0.dev0:** A1/O versus the per-position BO, a phenotype A versus B difference.
  - Phased reads put c.261del on the same haplotype as the B markers: ABO\*O.01.24, an O allele on a B background.
  - That haplotype also carries c.220C>T (172/190 reads), which O.01.24 does not list.
  - From 0.4.0.dev2 this goes to REVIEW as a possible new allele.
- **4 libraries, read-level REVIEW:**
  - a tie between O.09.01 and AW.31.02-05, which ISBT defines by the same four exon 7 changes;
  - two skewed BB/AB libraries with 11–13% foreign reads;
  - one library of 15 reads.

**Haplotypes that differ from their ISBT allele (0.4.0.dev2):**
- **Effect on calls:** none of the 6,241 calls change between 0.4.0.dev0 and dev2. 134 PASS calls go to REVIEW, all for this reason; 111 of them have at least 200 informative reads.
- **One recurring pattern:** 104 of the 134 are the same haplotype, ABO\*O.01.68 + c.595C>T.
  - Read support: median 97%, lowest 87%, at a median of 315 reads per haplotype.
  - Sequence context: no homopolymer.
  - Clair3 confirmation: the per-position path's Clair3 calls have c.595C>T in all 104 libraries, and in only 1 of the 119 O.01.68 libraries without the mismatch.
  - Backgrounds: O.01.01 (51), A1.01 (20), O.01.02 (13), B.01 (6), A2.01 (5) and others.
- **Interpretation:** c.595C>T defines ABO\*O.01.12, which ISBT lists without the exon 3–4 changes (c.106, c.188, c.189). Many O alleles were first described from exons 6–7 only. The haplotype is most likely O.01.12 carrying exon 3–4 changes that its definition does not list. It is present in about 3.7% of typed libraries.
- **How 0.4.0.dev3 reports it:** as ABO\*O.01.12-like, from a bundled table of recurring haplotypes that ISBT does not list (`abo_unlisted_haplotypes.tsv`). It carries a note and keeps PASS; any other difference still goes to REVIEW. Replayed on the dev2 results, this gives:

| Version | PASS | REVIEW | No-call |
| --- | --- | --- | --- |
| 0.4.0.dev0 | 2,587 | 189 | 3,465 |
| 0.4.0.dev2 | 2,453 | 323 | 3,465 |
| 0.4.0.dev3 (replayed) | 2,555 | 221 | 3,465 |

- **The remaining "possible new allele" REVIEWs:** 32, for example c.189C>T and c.220C>T on O.01.83 or O.16, and c.829A>G on O.01.36.

**Overlap with the QAP set (S5):**
- 103 of the 104 retrospective libraries from QAP donors are the QAP libraries themselves, with identical read counts. The QAP set is therefore not an independent check of this set.
- Re-typed with 0.4.0.dev3, every donor gets the same genotype class as before.
- The one independent library pair of the same donor differs at the allele level:
  - O.01.01/O.01.02 on 606 reads, against O.01.01/O.01.75 on 9,932 reads.
  - O.01.75 is O.01.02 + c.542G>A.
  - This points to a sample mix-up or a site missed in one library. It is listed for review.

**Read-by-read review of the 8 disagreements and the c.542 pair** (`hpc/inspect_sites.py`, 0.4.0.dev3; read counts are reads assigned to each allele):

| Case | Deciding evidence | Verdict |
| --- | --- | --- |
| cisAB.06 (1 library) | B-type haplotype (437 reads): c.526 C on 433 (B.01 needs G); c.657T, c.703A, c.796A, c.803C, no c.261del | read-level correct: cisAB.06, not B.01 |
| B(A).03 (1 donor, 2 runs) | B-type haplotype: c.703 G on 172/173 and 191/194 (B.01 needs A); c.526G, c.657T, c.796A, c.803C present | read-level correct in both runs: B(A).03, not B.01 |
| A1/O vs BO (1 library) | c.261del and c.796A on the same reads: 167 reads del+A, 163 neither, 34 mixed (chimeras); O haplotype also c.220T on 172/190 | per-position BO is wrong: A1 + an O.01.24-like allele with c.220T (phenotype A, not B); the c.220T allele is a candidate new allele |
| O.09.01 vs AW.31.02-05 (1 library) | second haplotype (1,974 reads): no c.261del, c.646A, c.681A, c.771T, c.829A, c.467C | unresolvable from exons 2–7: the two alleles have the same definition there; REVIEW is right; serology (O vs weak A) needed |
| BB/AB (2 libraries) | thousands of reads aligned but only 36 and 52 reach the typed sites; skewed split, 11–13% foreign reads | not typable from these libraries; re-test |
| 15-read library | 8 informative reads | not typable; re-test |
| c.542 pair (1 donor) | retrospective library: c.542 G on 91/99 reads of the O.01.02-type haplotype; QAP library: c.542 A on 4,310/4,520 | the two libraries come from different individuals (O.01.02 vs O.01.75): a sample identity problem, not a typing error |

So of the 4 confident disagreements, all 4 favour the read-level call:
- 3 because the per-position panel has no rules for B(A) and cisAB alleles;
- 1 because only whole-read phasing shows that c.261del sits on the B-background allele.

None of the 9 cases is a read-level typing error.

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
