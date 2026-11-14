#!/usr/bin/env bash
Q=/mnt/storage_c2/bioinformatics/qap_abo; X=/mnt/storage_c2/bioinformatics/abo_experiments
printf '%s\n' "8 1" "12 1" "20 1" "8 2" "12 2" "20 2" | xargs -P 2 -L 1 bash -c 'bash '$Q'/qap_depth_batch.sh $0 $1 2>'$Q'/depth/err_d$0_r$1.log'
for b in $Q/depth/d*_r*; do echo "== $(basename $b)"; $X/env/bin/python $X/core_new/experiments/qap_eval.py $Q/truth_audit.tsv $b/pipeline $b/read $b/calls.tsv; done > $Q/depth/summary.txt 2>&1
echo ALL_DONE >> $Q/depth/summary.txt
