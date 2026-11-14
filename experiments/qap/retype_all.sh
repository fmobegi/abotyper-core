#!/usr/bin/env bash
# Re-type every dataset with the current readtyper (after a rule change), then re-score. Restricted data stays here.
X=/mnt/storage_c2/bioinformatics/abo_experiments; Q=/mnt/storage_c2/bioinformatics/qap_abo
source $X/benv.sh
{
  for b in $Q/results/per_sample_processing/*/alignment/*.bam; do s=$(basename $b .bam); echo "$b $Q/readtyper/$s.ReadTyping.tsv $s"; done
  for d in $Q/depth/d*_r*; do for b in $d/bam/*.bam; do s=$(basename $b .bam); echo "$b $d/read/$s.ReadTyping.tsv $s"; done; done
  for d in $X/bench/d*_r*; do for b in $d/bam/*.bam; do s=$(basename $b .bam); echo "$b $d/read/$s.tsv $s"; done; done
} | xargs -P 5 -L 1 bash -c 'PYTHONPATH='$X'/core_new/src '$X'/env/bin/python -m abotyper_core.readtyper --bam "$0" --ref "$REF" --sample "$2" --out "$1" >/dev/null 2>&1'
$X/env/bin/python $X/core_new/experiments/qap_eval.py $Q/truth_audit.tsv $Q/results $Q/readtyper $Q/qap_calls.tsv > $Q/retype_full.txt
for b in $Q/depth/d*_r*; do echo "== $(basename $b)"; $X/env/bin/python $X/core_new/experiments/qap_eval.py $Q/truth_audit.tsv $b/pipeline $b/read $b/calls.tsv; done > $Q/depth/summary.txt 2>&1
$X/env/bin/python $Q/depth_detail.py $Q > $Q/depth/detail.txt
$X/env/bin/python $X/core_new/experiments/evaluate.py $X/bench > $X/bench/summary_pkg.tsv
echo RETYPE_DONE > $Q/retype.done
