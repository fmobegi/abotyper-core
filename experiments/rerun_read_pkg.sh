#!/usr/bin/env bash
# Re-type every test-sample benchmark batch (depth and contamination series) with the packaged abotyper-readtyper
# instead of the prototype readtyper.py. The prototype's outputs are kept in read_proto/; the packaged ones replace
# read/, and evaluate.py writes bench/summary_pkg.tsv.
# usage: source benv.sh; SRC=<abotyper-core checkout>/src bash rerun_read_pkg.sh
set -uo pipefail
here=$(cd "$(dirname "$0")" && pwd)
: "${SRC:=$here/../src}"
for b in "$ROOT"/bench/d*_r*; do
    [ -d "$b/read_proto" ] || mv "$b/read" "$b/read_proto"
    mkdir -p "$b/read"
    for bam in "$b"/bam/*.bam; do
        s=$(basename "$bam" .bam)
        echo "$bam $b/read/$s.tsv $s"
    done
done | xargs -P 4 -L 1 bash -c 'PYTHONPATH='"$SRC"' "$ROOT/env/bin/python" -m abotyper_core.readtyper --bam "$0" --ref "$REF" --sample "$2" --out "$1" >/dev/null 2>&1'
"$ROOT/env/bin/python" "$here/evaluate.py" "$ROOT/bench" > "$ROOT/bench/summary_pkg.tsv"
cat "$ROOT/bench/summary_pkg.tsv"
