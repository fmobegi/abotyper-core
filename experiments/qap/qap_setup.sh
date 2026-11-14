#!/usr/bin/env bash
# Build the samplesheet for the validation (QAP) runs and the server resource config. Restricted data: everything
# stays under $Q on the server.
# usage: bash qap_setup.sh RUN_DIR_1 RUN_DIR_2 ...   (one folder of barcoded FASTQs per sequencing run, in run order;
#        they are renamed R1, R2, ... and each library becomes <run>_<donor>_barcodeNN)
set -euo pipefail
Q=/mnt/storage_c2/bioinformatics/qap_abo
cd "$Q/fastq"
i=1
for d in "$@"; do [ -d "$d" ] && mv "$d" "R$i"; i=$((i + 1)); done
echo "sample,fastq_1,fastq_2,control" > "$Q/samplesheet.csv"
for f in R*/*.fastq.gz; do
    r=${f%%/*}; b=$(basename "$f" .fastq.gz)
    echo "${r}_${b},$Q/fastq/$f,," >> "$Q/samplesheet.csv"
done
wc -l "$Q/samplesheet.csv"
cat > "$Q/server.config" <<'CFG'
process { resourceLimits = [ cpus: 6, memory: 12.GB, time: 24.h ] }
CFG
# then: nextflow run nf-core/abotyper -r <revision> -profile conda -c $Q/server.config \
#         --input $Q/samplesheet.csv --outdir $Q/results
