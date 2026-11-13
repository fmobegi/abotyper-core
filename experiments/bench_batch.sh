#!/usr/bin/env bash
# One benchmark batch: downsample every sample BAM to DEPTH reads (seed REP), then type it with
#   read      - experiments/readtyper.py (read-level allele-pair likelihood)
#   cur       - BAM pileup only, current pipeline settings (-a -B -d 0 -Q 7 -q 1) -> pileup2metrics -> predict
#   baq       - as cur, BAQ on (no -B)
#   q0        - as cur, no base-quality filter (-Q 0)
#   q12       - as cur, -Q 12
#   pipeline  - Clair3 (pipeline args) + clair2metrics/haplotypes + cur pileup evidence -> predict (only if CLAIR=1)
# and aggregate each method's batch with abotyper-aggregate.
#
# usage: bench_batch.sh DEPTH REP   (DEPTH = number of reads, or "full")
# env:   ROOT (experiment dir), BAMS (glob of input BAMs), REF, CLAIR_ENV, CLAIR_MODEL, CLAIR_BED, CLAIR=0|1
set -euo pipefail
DEPTH=$1; REP=$2
B="$ROOT/bench/d${DEPTH}_r${REP}"
E="$ROOT/env/bin"
PY="$E/python"
mkdir -p "$B/bam" "$B/read"
"$E/abotyper-panel-to-bed" --output "$B/panel.bed" >/dev/null 2>&1

declare -A MP=( [cur]="-a -B -d 0 -Q 7 -q 1" [baq]="-a -d 0 -Q 7 -q 1" [q0]="-a -B -d 0 -Q 0 -q 1" [q12]="-a -B -d 0 -Q 12 -q 1" )

typeit() {   # method sample pileup [clair_metrics haplotypes]
    local m=$1 s=$2 pile=$3 clair=${4:-} hap=${5:-}
    local d="$B/$m/per_sample_processing/$s/combined"
    mkdir -p "$d"
    "$E/abotyper-pileup2metrics" --pileup "$pile" ${clair:+--clair-metrics "$clair"} --sample "$s" \
        --merged-out "$d/$s.AlignmentStatistics.tsv" --concordance-out "$d/$s.concordance.tsv" \
        --status-out "$d/$s.evidence_status.tsv" 2>/dev/null
    "$E/abotyper-predict-phenotype" -i "$d/$s.AlignmentStatistics.tsv" -o "$d/$s.ABOPhenotype.txt" \
        -c "$B/bam/$s.coverage.txt" >/dev/null 2>&1
    cp "$B/bam/$s.coverage.txt" "$B/bam/$s.stats" "$B/bam/$s.flagstat" "$d/"
    [ -n "$hap" ] && cp "$hap" "$d/$s.Haplotypes.tsv"
    return 0
}

for src in $BAMS; do
    s=$(basename "$src" .bam)
    bam="$B/bam/$s.bam"
    if [ "$DEPTH" = full ] || [[ "$DEPTH" == mix* ]]; then
        cp "$src" "$bam"
    else
        n=$("$E/samtools" view -c -F 0x904 -q 1 "$src")
        frac=$(awk -v d="$DEPTH" -v n="$n" 'BEGIN{f=d/n; if(f>0.999) f=0.999; printf "%.4f", f}')
        "$E/samtools" view -b -F 0x904 -q 1 --subsample "$frac" --subsample-seed "$REP" -o "$bam" "$src"
    fi
    "$E/samtools" index "$bam"
    "$E/samtools" coverage "$bam" > "$B/bam/$s.coverage.txt"
    "$E/samtools" stats "$bam" > "$B/bam/$s.stats"
    "$E/samtools" flagstat "$bam" > "$B/bam/$s.flagstat"

    "$PY" "$ROOT/core/experiments/readtyper.py" --bam "$bam" --ref "$REF" --sample "$s" --out "$B/read/$s.tsv"

    for m in cur baq q0 q12; do
        "$E/samtools" mpileup ${MP[$m]} -l "$B/panel.bed" -f "$REF" "$bam" 2>/dev/null > "$B/bam/$s.$m.pileup"
        typeit "$m" "$s" "$B/bam/$s.$m.pileup"
    done

    if [ "${CLAIR:-0}" = 1 ]; then
        c="$B/clair3/$s"; mkdir -p "$c"
        ( export PATH="$CLAIR_ENV/bin:$PATH" CONDA_PREFIX="$CLAIR_ENV"
          run_clair3.sh --bam_fn="$bam" --ref_fn="$REF" --threads=2 --output="$c" --platform=ont \
              --model_path="$CLAIR_MODEL" --bed_fn="$CLAIR_BED" --include_all_ctgs --enable_phasing \
              --enable_long_indel --sample_name="$s" --gvcf > "$c/clair3.log" 2>&1 ) || true
        if [ -s "$c/merge_output.gvcf.gz" ]; then
            "$E/abotyper-clair2metrics" -i "$c/merge_output.gvcf.gz" -o "$c/$s.AlignmentStatistics.tsv" 2>/dev/null
            hap=""
            if [ -s "$c/phased_merge_output.vcf.gz" ]; then
                "$E/abotyper-clair2haplotypes" -i "$c/phased_merge_output.vcf.gz" -o "$c/$s.Haplotypes.tsv" 2>/dev/null && hap="$c/$s.Haplotypes.tsv"
            fi
            typeit pipeline "$s" "$B/bam/$s.cur.pileup" "$c/$s.AlignmentStatistics.tsv" "$hap"
        else
            typeit pipeline "$s" "$B/bam/$s.cur.pileup"      # Clair3 failed: the pipeline falls back to the pileup
        fi
    fi
done

for m in cur baq q0 q12 pipeline; do
    [ -d "$B/$m" ] || continue
    ( cd "$B/$m" && "$E/abotyper-aggregate" per_sample_processing > aggregate.log 2>&1 || true )
done
echo "done $DEPTH $REP"
