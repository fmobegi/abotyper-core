#!/usr/bin/env bash
# Contamination test: spike each sample with reads of a sample of a different genotype (5, 10, 20% of the final
# reads), then type the mixed BAMs with every method (bench_batch.sh, depth "full"). Truth stays the main sample.
# usage: source benv.sh; bash mix_all.sh
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
E="$ROOT/env/bin"
SRC=$(dirname "$(echo $BAMS | cut -d' ' -f1)")/../..
declare -A CONTAM=( [ABO_AA_barcode01]=ABO_BB_barcode05 [ABO_AB_barcode02]=ABO_O1O1_barcode08
    [ABO_AO1_barcode03]=ABO_BB_barcode05 [ABO_AO2_barcode04]=ABO_O1O2_barcode09 [ABO_BB_barcode05]=ABO_AA_barcode01
    [ABO_BO1_barcode06]=ABO_AO2_barcode04 [ABO_BO2_barcode07]=ABO_AA_barcode01 [ABO_O1O1_barcode08]=ABO_AB_barcode02
    [ABO_O1O2_barcode09]=ABO_BB_barcode05 [ABO_O2O2_barcode10]=ABO_AA_barcode01 )
bam_of() { echo "$SRC/$1/alignment/$1.bam"; }
for pct in 5 10 20; do
    out="$ROOT/mix/f$pct"; mkdir -p "$out"
    for main in "${!CONTAM[@]}"; do
        c=${CONTAM[$main]}
        n=$("$E/samtools" view -c -F 0x904 -q 1 "$(bam_of "$main")")
        nc=$("$E/samtools" view -c -F 0x904 -q 1 "$(bam_of "$c")")
        frac=$(awk -v p="$pct" -v n="$n" -v nc="$nc" 'BEGIN{f=(p/(100-p))*n/nc; if(f>0.999) f=0.999; printf "%.4f", f}')
        "$E/samtools" view -b -F 0x904 -q 1 "$(bam_of "$main")" -o "$out/$main.main.bam"
        "$E/samtools" view -b -F 0x904 -q 1 --subsample "$frac" --subsample-seed 7 "$(bam_of "$c")" -o "$out/$main.contam.bam"
        "$E/samtools" merge -f -o - "$out/$main.main.bam" "$out/$main.contam.bam" | "$E/samtools" sort -o "$out/$main.bam"
        rm "$out/$main.main.bam" "$out/$main.contam.bam"
    done
done
for pct in 5 10 20; do
    echo "mix$pct"
done | xargs -P 3 -I{} bash -c 'p={}; BAMS="$ROOT/mix/f${p#mix}/*.bam" CLAIR=1 bash '"$here"'/bench_batch.sh "$p" 1 >/dev/null 2>"$ROOT/bench/err_d${p}_r1.log"; echo "done $p"'
