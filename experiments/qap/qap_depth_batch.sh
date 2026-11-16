#!/usr/bin/env bash
# Depth benchmark on the QAP libraries: downsample each BAM to DEPTH reads (seed REP), then type with the pipeline
# path (Clair3 + pileup -> pileup2metrics -> predict -> aggregate) and abotyper-readtyper. Restricted data: stays in qap_abo/.
set -uo pipefail
DEPTH=$1; REP=$2
X=/mnt/storage_c2/bioinformatics/abo_experiments; Q=/mnt/storage_c2/bioinformatics/qap_abo
source $X/benv.sh
B=$Q/depth/d${DEPTH}_r${REP}; E=$X/env/bin
mkdir -p $B/bam $B/read
$E/abotyper-panel-to-bed --output $B/panel.bed >/dev/null 2>&1
for src in $Q/results/per_sample_processing/*/alignment/*.bam; do
  s=$(basename $src .bam); bam=$B/bam/$s.bam
  R=NG_006669.2:18064-24544; n=$($E/samtools view -c -F 0x904 -q 1 $src $R); [ "$n" -eq 0 ] && continue
  frac=$(awk -v d=$DEPTH -v n=$n 'BEGIN{f=d/n; if(f>0.999) f=0.999; printf "%.4f", f}')
  $E/samtools view -b -F 0x904 -q 1 --subsample $frac --subsample-seed $REP -o $bam $src $R && $E/samtools index $bam
  $E/samtools coverage $bam > $B/bam/$s.coverage.txt; $E/samtools stats $bam > $B/bam/$s.stats; $E/samtools flagstat $bam > $B/bam/$s.flagstat
  PYTHONPATH=$X/core_new/src $E/python -m abotyper_core.readtyper --bam $bam --ref $REF --sample $s --out $B/read/$s.ReadTyping.tsv >/dev/null 2>&1
  $E/samtools mpileup -a -B -d 0 -Q 7 -q 1 -l $B/panel.bed -f $REF $bam 2>/dev/null > $B/bam/$s.pileup
  c=$B/clair3/$s; mkdir -p $c
  ( export PATH="$CLAIR_ENV/bin:$PATH" CONDA_PREFIX="$CLAIR_ENV"
    run_clair3.sh --bam_fn=$bam --ref_fn=$REF --threads=2 --output=$c --platform=ont --model_path=$CLAIR_MODEL \
      --bed_fn=$CLAIR_BED --include_all_ctgs --enable_phasing --enable_long_indel --sample_name=$s --gvcf > $c/clair3.log 2>&1 ) || true
  d=$B/pipeline/per_sample_processing/$s/combined; mkdir -p $d; cm=""; 
  if [ -s $c/merge_output.gvcf.gz ]; then
    $E/abotyper-clair2metrics -i $c/merge_output.gvcf.gz -o $c/$s.AlignmentStatistics.tsv 2>/dev/null && cm=$c/$s.AlignmentStatistics.tsv
    [ -s $c/phased_merge_output.vcf.gz ] && $E/abotyper-clair2haplotypes -i $c/phased_merge_output.vcf.gz -o $d/$s.Haplotypes.tsv 2>/dev/null
  fi
  $E/abotyper-pileup2metrics --pileup $B/bam/$s.pileup ${cm:+--clair-metrics $cm} --sample $s --merged-out $d/$s.AlignmentStatistics.tsv \
     --concordance-out $d/$s.concordance.tsv --status-out $d/$s.evidence_status.tsv 2>/dev/null
  $E/abotyper-predict-phenotype -i $d/$s.AlignmentStatistics.tsv -o $d/$s.ABOPhenotype.txt -c $B/bam/$s.coverage.txt >/dev/null 2>&1
  cp $B/bam/$s.coverage.txt $B/bam/$s.stats $B/bam/$s.flagstat $d/
  rm -rf $c/tmp $c/log
done
( cd $B/pipeline && $E/abotyper-aggregate per_sample_processing > aggregate.log 2>&1 )
echo "done d$DEPTH r$REP"
