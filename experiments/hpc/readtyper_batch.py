"""Type many BAMs with abotyper-readtyper in parallel, on any OS (Windows, Linux, macOS).

usage: python readtyper_batch.py --ref NG_006669v2.fasta --out OUTDIR (--bams BAM_LIST.txt | --root DIR) [--jobs N]

  --bams   text file, one BAM path per line
  --root   folder searched recursively for *.bam (e.g. a drive or a folder of pipeline runs)
  --jobs   parallel worker processes (default: all cores; each needs well under 1 GB)

Writes OUTDIR/<sample>.ReadTyping.tsv and <sample>.ReadTypingHaplotypes.tsv per BAM (sample = file name without
.bam) and OUTDIR/failed.tsv. Resumable: samples with an output already are skipped, so re-run the same command after
an interruption. The allele table is loaded once per worker, not once per BAM. Needs abotyper-core >= 0.4.0.dev0.
Then: python summarise_readtyping.py OUTDIR [--calls existing_calls.tsv] > summary.txt
"""

from __future__ import annotations

import argparse
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from abotyper_core import readtyper as rt

_W = {}


def _init(ref: str) -> None:
    seq = rt.read_fasta(Path(ref))
    offsets = rt.exon_offsets(rt.load_panel(rt.default_panel_path()))
    _W.update(seq=seq, table=rt.load_alleles(rt.ISBT_TABLE, offsets, seq), coding=rt.coding_map(offsets))


def _type(bam: str, out: str) -> tuple[str, str]:
    s = Path(bam).name[:-4] if bam.lower().endswith(".bam") else Path(bam).name
    try:
        reads, seen = rt.observe_reads(Path(bam), _W["table"], _W["coding"])
        res = rt.type_sample(reads, _W["table"], _W["seq"], _W["coding"])
        base = {"sample": s, "reads_aligned": seen, "reads_typed": 0}
        res.summary = {k: {**base, **res.summary}.get(k, "") for k in rt.SUMMARY_FIELDS}
        tmp = Path(out) / f"{s}.ReadTyping.tsv.tmp"
        rt.write_outputs(res, s, tmp, Path(out) / f"{s}.ReadTypingHaplotypes.tsv")
        os.replace(tmp, Path(out) / f"{s}.ReadTyping.tsv")     # only complete outputs count as done
        return s, ""
    except Exception as e:                                       # noqa: BLE001 - report and carry on
        return s, f"{type(e).__name__}: {e}".replace("\t", " ").replace("\n", " ")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", required=True)
    ap.add_argument("--out", required=True)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--bams")
    src.add_argument("--root")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    a = ap.parse_args()

    if a.bams:
        bams = [ln.strip() for ln in open(a.bams) if ln.strip()]
    else:
        bams = sorted(str(p) for p in Path(a.root).rglob("*.bam"))
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    names = {}
    for b in bams:
        names.setdefault(Path(b).name[:-4], []).append(b)
    dupes = {k: v for k, v in names.items() if len(v) > 1}
    if dupes:
        print(f"warning: {len(dupes)} sample names occur in more than one BAM; only the first of each is typed "
              f"(e.g. {next(iter(dupes))})", file=sys.stderr)
    todo = [v[0] for k, v in names.items() if not (out / f"{k}.ReadTyping.tsv").is_file()]
    print(f"{len(bams)} BAMs, {len(names) - len(todo)} already typed, {len(todo)} to type with {a.jobs} workers",
          file=sys.stderr)

    t0, done, failed = time.time(), 0, []
    with ProcessPoolExecutor(max_workers=a.jobs, initializer=_init, initargs=(a.ref,)) as ex:
        futs = [ex.submit(_type, b, str(out)) for b in todo]
        for f in as_completed(futs):
            s, err = f.result()
            done += 1
            if err:
                failed.append((s, err))
            if done % 100 == 0 or done == len(todo):
                rate = done / max(1e-9, time.time() - t0)
                print(f"  {done}/{len(todo)} ({rate:.1f}/s, ~{(len(todo) - done) / max(rate, 1e-9) / 60:.0f} min left)"
                      f"{f', {len(failed)} failed' if failed else ''}", file=sys.stderr)
    with open(out / "failed.tsv", "w") as fh:
        fh.write("sample\terror\n")
        fh.writelines(f"{s}\t{e}\n" for s, e in failed)
    print(f"done: {done - len(failed)} typed, {len(failed)} failed (see {out / 'failed.tsv'})", file=sys.stderr)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.exit(130)
    except Exception:
        traceback.print_exc()
        sys.exit(1)
