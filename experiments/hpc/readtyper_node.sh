#!/usr/bin/env bash
# Run abotyper-readtyper on many BAMs directly on one node, in parallel (no scheduler).
#
# usage: bash readtyper_node.sh BAM_LIST REF OUTDIR [JOBS]
#   BAM_LIST  one BAM path per line (reads aligned to NG_006669.2), e.g.
#               find /path/to/runs -path '*alignment/*.bam' > bams.txt
#   REF       NG_006669.2 FASTA (assets/refs/NG_006669v2.fasta in nf-core/abotyper)
#   OUTDIR    <sample>.ReadTyping.tsv and <sample>.ReadTypingHaplotypes.tsv per BAM (sample = BAM name without .bam)
#   JOBS      parallel processes (default: all cores; each needs < 1 GB)
#
# Resumable: samples with an output already are skipped, so it can be re-run after an interruption. Run it inside
# tmux/screen (or with nohup) so it survives a dropped connection. Progress: ls OUTDIR/*.ReadTyping.tsv | wc -l.
# Then: python summarise_readtyping.py OUTDIR [--calls existing_calls.tsv] > summary.txt
set -uo pipefail
BAM_LIST=$1; REF=$2; OUTDIR=$3; JOBS=${4:-$(nproc)}
mkdir -p "$OUTDIR"
export REF OUTDIR
one() {
    bam=$1; s=$(basename "$bam" .bam)
    [ -s "$OUTDIR/$s.ReadTyping.tsv" ] && return 0
    abotyper-readtyper --bam "$bam" --ref "$REF" --sample "$s" \
        --out "$OUTDIR/$s.ReadTyping.tsv.tmp" --haplotypes-out "$OUTDIR/$s.ReadTypingHaplotypes.tsv" \
        2> "$OUTDIR/$s.err" \
        && mv "$OUTDIR/$s.ReadTyping.tsv.tmp" "$OUTDIR/$s.ReadTyping.tsv" && rm -f "$OUTDIR/$s.err" \
        || echo "FAILED $bam (see $OUTDIR/$s.err)"
}
export -f one
start=$(date +%s)
grep -v '^\s*$' "$BAM_LIST" | xargs -P "$JOBS" -I{} bash -c 'one "$@"' _ {} | tee "$OUTDIR/failed.txt"
n=$(ls "$OUTDIR"/*.ReadTyping.tsv 2>/dev/null | wc -l)
echo "typed $n of $(grep -cv '^\s*$' "$BAM_LIST") BAMs in $(( $(date +%s) - start )) s; failures: $(grep -c FAILED "$OUTDIR/failed.txt")"
