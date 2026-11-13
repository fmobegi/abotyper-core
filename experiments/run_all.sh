#!/usr/bin/env bash
# Run the downsampling benchmark: depths x replicates, 3 batches at a time (Clair3 uses 2 threads each).
# Clair3 (the full pipeline path) runs for replicates 1-3 only; the pileup and read-level methods for all 5.
# usage: source benv.sh; bash run_all.sh
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
jobs=()
for d in 40 20 12 8 5; do
    for r in 1 2 3 4 5; do
        jobs+=("$d $r $([ "$r" -le 3 ] && echo 1 || echo 0)")
    done
done
printf '%s\n' "${jobs[@]}" | xargs -P 3 -L 1 bash -c 'CLAIR=$2 bash '"$here"'/bench_batch.sh $0 $1 > /dev/null 2>"$ROOT/bench/err_d$0_r$1.log"; echo "done d$0 r$1"'
"$ROOT/env/bin/python" "$here/evaluate.py" "$ROOT/bench" > "$ROOT/bench/summary.tsv"
cat "$ROOT/bench/summary.tsv"
