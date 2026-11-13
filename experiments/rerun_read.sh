#!/usr/bin/env bash
# Re-run only the read-level typer on every benchmark batch (after a change to readtyper.py), then re-evaluate.
# usage: source benv.sh; bash rerun_read.sh
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
for b in "$ROOT"/bench/d*_r*; do
    for bam in "$b"/bam/*.bam; do
        s=$(basename "$bam" .bam)
        echo "$bam $b/read/$s.tsv $s"
    done
done | xargs -P 6 -L 1 bash -c '"$ROOT/env/bin/python" '"$here"'/readtyper.py --bam "$0" --ref "$REF" --sample "$2" --out "$1"'
"$ROOT/env/bin/python" "$here/evaluate.py" "$ROOT/bench" > "$ROOT/bench/summary.tsv"
cat "$ROOT/bench/summary.tsv"
