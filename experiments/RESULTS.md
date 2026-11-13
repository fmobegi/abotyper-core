# Processing and typing experiments (November 2026)

Exploratory work, not part of the package. Goal: see whether the BAM → pileup/Clair3 → per-position evidence →
allele matching path could be more accurate or more sensitive, and try a different way of assigning alleles.

## Data and set-up

- The 10 samples of the pipeline `test` profile (AA, AB, AO1, AO2, BB, BO1, BO2, O1O1, O1O2, O2O2). The truth is
  taken from the sample names; the BAMs come from a completed `-profile test,conda` run, so no pipeline step was
  re-run. 67–198 reads per sample, nearly all spanning the whole exon 2–7 amplicon (median 6.7 kb).
- **Depth series:** each BAM was downsampled to 40, 20, 12, 8 and 5 reads (`samtools view --subsample`), with five
  seeds per depth. Clair3 ran on three of the five.
- **Contamination series:** each sample was spiked with reads from a sample of another genotype, making up 5, 10 or
  20% of the final reads.
- **Methods compared:**

| arm | what it is |
|---|---|
| `pipeline` | Clair3 (pipeline arguments) + clair2metrics/haplotypes + current pileup → predict-phenotype → aggregate |
| `cur` | pileup only, current settings `-a -B -d 0 -Q 7 -q 1` |
| `baq` | as `cur`, BAQ on (no `-B`) |
| `q0` / `q12` | as `cur`, base-quality filter `-Q 0` / `-Q 12` |
| `read` | `readtyper.py`: read-level allele-pair likelihood against the ISBT table (below) |

Scripts:

- `bench_batch.sh`, `run_all.sh`: the depth series.
- `mix_all.sh`: the contamination series.
- `rerun_read.sh`: re-runs only the read-level typer.
- `evaluate.py`: writes `calls.tsv`, one row per call, and a summary.

## Results

Genotype class correct (%). n = 10 samples × reps (50 per depth; 30 for `pipeline`, which is the Clair3 subset).

| depth (reads) | read | pipeline | cur | baq | q0 | q12 |
|---|---|---|---|---|---|---|
| full (67–198) | 100 | 100 | 100 | 100 | 100 | 100 |
| 40 | 98 | 100 | 100 | 100 | 94 | 98 |
| 20 | 96 | 97 | 80 | 80 | 80 | 92 |
| 12 | 96 | 83 | 76 | 72 | 72 | 84 |
| 8 | 94 | 70 | 68 | 66 | 68 | 74 |
| 5 | 74 | 43 | 64 | 58 | 62 | 62 |

O1/O2 subtype correct (% of samples carrying an O) follows the same pattern:

| method | 20 reads | 12 reads | 8 reads |
|---|---|---|---|
| read | 94 | 97 | 94 |
| pipeline | 95 | 81 | 52 |

Contamination, one replicate each. Each cell is genotype correct / calls that were wrong but PASS:

| spike | read | pipeline | cur |
|---|---|---|---|
| 5% | 90 / – | 90 / **1** | 100 / 0 |
| 10% | 90 / – | 80 / **1** | 100 / 0 |
| 20% | 80 / – | 80 / **1** | 50 / 0 |

### What the numbers say

1. **At full depth every method is right.** On these 10 samples, the processing settings only matter once depth
   drops or the sample is impure.
2. **`-B` (no BAQ) is the right setting.** BAQ is never better. At full depth it changes the extended call in both
   A2 samples, ABO\*A2.01 → ABO\*A1.02, because BAQ down-weights the c.1061delC homopolymer evidence. Keep `-B`.
3. **`-Q 12` beats the current `-Q 7` below 40 reads** (20 reads: 92 vs 80%; 12 reads: 84 vs 76%). `-Q 0` is no
   better than `-Q 7`. This is worth testing on real data before changing it: on R10/SUP the base qualities around
   homopolymers are what the `-Q` filter acts on.
4. **Most low-depth failures are no-calls and hom/het confusions, not wrong alleles.** Below about 20 reads, the
   pipeline and pileup arms often call AO as AA or BO as "Unknown". They are always REVIEW, never a wrong PASS. The
   current QC gating does its job.
5. **Read-level typing is much more sensitive.** It stays at 94–96% correct down to 8 reads, where the pipeline is at
   70%. It keeps the O1/O2 subtype at 94% where the pipeline drops to 52%. Each long read is a whole haplotype, so
   phase comes for free. Ten reads that agree across all 129 ISBT sites are stronger evidence than ten independent
   per-position calls.
6. **Contamination can produce a wrong PASS in the pipeline.** The O2O2 test sample already carries about 12% reads
   that look like A1 (reference sequence at every O2-defining site). At full depth both methods still call O2O2, but
   adding as little as 5% foreign reads tips the pipeline to ABO\*A1.01/O.02.02 with **QC PASS** (phenotype A instead
   of O). The read-level typer made the same wrong call. Unlike the pipeline, it flags the sample: 14–29% of reads fit
   neither called allele (`unexplained_share`, below).

### Read-level typer failures, and what was added for them

- **A third haplotype was present** (the O2O2 cases above).
  - Added `unexplained_share` and `third_allele`: the share of reads that some other ISBT allele explains at least two
    SNV mismatches better than either called allele, and which allele that is.
  - With a review threshold of 10%, all wrong O2O2 calls are flagged.
  - It also flags the undownsampled O2O2 (12%). That sample is worth a look.
- **Low-level minor alleles.**
  - The O1O1 sample has 4/90 reads that match ABO\*O.01.68.
  - Added an allele-balance rule: the call becomes homozygous when the minor allele is under 20% of reads *and*
    P(≤ k minor reads | 50:50 heterozygote) < 10⁻³.
  - Without the binomial test, the rule wrongly made a 2-of-12 heterozygote homozygous.
  - With it, the rule wrongly made only one heterozygote homozygous: AO1 at 20 reads, a 3/22 split.
- **Rare alleles one SNV from a common pair.**
  - At 8–12 reads, two calls of ABO\*B.01/cisAB.0x were made where the truth was B.01/A1.01 or B.01/B.01.
  - A flat log-prior of −3 for rare alleles is not enough. Allele-frequency priors would fix this, for example from
    the 50-donor validation set or the published frequencies. So would requiring a minimum number of reads that
    carry the rare allele's defining variant.
- **Subsampling artefact.** At 8 and 12 reads with seed 2, samtools returned 0 reads for O2O2. This counts as a
  no-call in every arm; it is not a typing error.

## Recommendations

### Processing (can go into a minor release; validate on the 50-donor set first)

1. Keep `-B` and `-a -d 0`. They are justified by the A2 result above.
2. Try `-Q 12` (or `-Q 10`) in place of `-Q 7`, and validate on the donor set. It is a one-line change in
   `modules.config`, and it helped at every depth below 40 reads in this test.
3. Add a contamination / third-haplotype check to QC. Even without read-level typing, it can be done from the
   phased Clair3 VCF or from read-backed haplotype counts: count reads whose haplotype matches neither called allele
   at the diagnostic sites. Any sample above about 10% goes to REVIEW. This would have turned the wrong PASS above into
   REVIEW.
4. Report allele balance per heterozygous call (reads per allele) next to the existing depth QC, and flag skewed
   splits that cannot be explained by sampling.
5. Make the minimum-depth gate depth-aware. Below about 20 reads, pileup-based calls fall to about 80%. Either raise
   the REVIEW threshold or use the read-level result there (below).

### New ways to assign alleles and genotypes (later releases)

1. **Read-level allele-pair likelihood** (this prototype).
   - Each read is scored against every ISBT allele, and the pair that best explains the reads is chosen. The output
     includes a likelihood margin, tied pairs, allele balance and the unexplained-read share.
   - Run it first as a cross-check column in the report. "Agrees / disagrees with the pileup call" is a cheap and
     strong QC signal.
   - Before it could be a primary caller, it needs:
     - allele-frequency priors;
     - a calibrated error model, fitted per chemistry from the homopolymer indel and SNV error rates;
     - handling of the 7 ISBT changes the parser skips, and of exon 1.
2. **Haplotype consensus.**
   - Cluster reads into two haplotypes on the diagnostic sites (or with WhatsHap/longphase on the Clair3 VCF), then
     build a per-haplotype consensus with medaka, or spoa plus a polish.
   - Look up each consensus in the ISBT table by full-sequence identity. This also reports novel variants on the
     right haplotype, which per-position matching can only list.
3. **EM over the allele table with a contamination component.** Fit the pair plus an "other" component with its own
   weight. This turns contamination into an estimated fraction instead of a threshold, and allows a mixed or
   chimeric sample to be reported as such.
4. **`bcftools mpileup -a AD,ADF,ADR,DP` + `bcftools call`** in place of the custom pileup parser for the per-position
   evidence. The format is standard, strand counts are free (strand bias is a common artefact flag), and the VCF can
   be phased with WhatsHap.
5. **Calibrated confidence.** Report a genotype posterior or a likelihood margin, and set the PASS/REVIEW cut-offs
   from the validation set, so that thresholds stop being per-position heuristics.

### Caveats

- 10 samples, one sequencing run, one replicate per contamination level. The numbers show direction, not
  performance claims.
- Truth for the test samples is the genotype class and O subtype from their names. Extended alleles were compared
  with the pipeline's own full-depth calls.
- Nothing here used the restricted validation data. Calibrating the read-level model and the `-Q` change needs the
  50-donor set, on the server where it lives.
