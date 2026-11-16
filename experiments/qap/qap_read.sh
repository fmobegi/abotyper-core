#!/usr/bin/env bash
# Run abotyper-readtyper (development copy) on every QAP BAM, then score both methods. Restricted data: all output
# stays under qap_abo/; the evaluation prints aggregate numbers only.
set -uo pipefail
X=/mnt/storage_c2/bioinformatics/abo_experiments
Q=/mnt/storage_c2/bioinformatics/qap_abo
source $X/benv.sh
mkdir -p $Q/readtyper
ls $Q/results/per_sample_processing/*/alignment/*.bam | xargs -P 4 -I{} bash -c 'b={}; s=$(basename $b .bam); PYTHONPATH='$X'/core_new/src '$X'/env/bin/python -m abotyper_core.readtyper --bam $b --ref $REF --sample $s --out '$Q'/readtyper/$s.ReadTyping.tsv --haplotypes-out '$Q'/readtyper/$s.ReadTypingHaplotypes.tsv'
$X/env/bin/python $X/core_new/experiments/qap_eval.py $Q/truth_audit.tsv $Q/results $Q/readtyper $Q/qap_calls.tsv
