#!/usr/bin/env bash
# Run abotyper-readtyper on many BAMs as a Slurm job array: one array task per chunk of BAMs.
#
# 1. List the BAMs (reads aligned to NG_006669.2, e.g. per_sample_processing/*/alignment/*.bam from nf-core/abotyper
#    or the earlier pipeline), one path per line:
#        find /path/to/runs -path '*alignment/*.bam' > bams.txt
# 2. Split into chunks and submit (here 100 BAMs per task):
#        split -l 100 -d -a 3 bams.txt chunk_ && ls chunk_* > chunks.txt
#        sbatch --array=1-$(wc -l < chunks.txt)%50 readtyper_array.slurm
# 3. Summarise: python summarise_readtyping.py OUTDIR [--calls existing_calls.tsv] > summary.txt
#
# Each BAM takes about 5-20 s and well under 1 GB. Output: OUTDIR/<sample>.ReadTyping.tsv and
# OUTDIR/<sample>.ReadTypingHaplotypes.tsv. The sample name is the BAM file name without .bam.
#SBATCH --job-name=readtyper
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=02:00:00
#SBATCH --output=logs/readtyper_%A_%a.log
set -euo pipefail

: "${REF:?set REF to the NG_006669.2 FASTA (assets/refs/NG_006669v2.fasta in nf-core/abotyper)}"
: "${OUTDIR:=readtyping}"
# activate the environment that has abotyper-core (with readtyper) installed, e.g.
#   module load miniforge && conda activate abotyper-readtyper
mkdir -p "$OUTDIR" logs
chunk=$(sed -n "${SLURM_ARRAY_TASK_ID}p" chunks.txt)
while read -r bam; do
    s=$(basename "$bam" .bam)
    [ -s "$OUTDIR/$s.ReadTyping.tsv" ] && continue          # resumable: skip samples already typed
    abotyper-readtyper --bam "$bam" --ref "$REF" --sample "$s" \
        --out "$OUTDIR/$s.ReadTyping.tsv" --haplotypes-out "$OUTDIR/$s.ReadTypingHaplotypes.tsv" \
        || echo "FAILED $bam" >&2
done < "$chunk"
