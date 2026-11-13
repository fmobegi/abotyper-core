"""Minimal pure-Python BAM reader (no pysam/samtools needed).

BAM is BGZF, which is plain multi-member gzip, so the standard library can decompress it; the record layout is
fixed by the SAM/BAM specification (https://samtools.github.io/hts-specs/SAMv1.pdf, section 4.2). Only what the
read-level typer needs is decoded: reference id, position, flag, MAPQ, CIGAR and sequence. Amplicon BAMs are small
(hundreds to a few thousand reads), so streaming the whole file is fine.
"""

from __future__ import annotations

import gzip
import struct
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterator, List, Tuple, Union

import numpy as np

_SEQ = "=ACMGRSVTWYHKDBN"
_CIGAR_OPS = "MIDNSHP=X"
# flags
UNMAPPED, SECONDARY, SUPPLEMENTARY = 0x4, 0x100, 0x800


@dataclass(frozen=True)
class Alignment:
    ref: str
    pos: int                          # 0-based leftmost reference position
    flag: int
    mapq: int
    cigar: Tuple[Tuple[str, int], ...]
    seq: str

    @property
    def is_primary_mapped(self) -> bool:
        return not self.flag & (UNMAPPED | SECONDARY | SUPPLEMENTARY)

    @property
    def ref_end(self) -> int:
        """0-based exclusive reference end."""
        return self.pos + sum(n for op, n in self.cigar if op in "MDN=X")


def _read_exact(fh: BinaryIO, n: int) -> bytes:
    b = fh.read(n)
    if len(b) != n:
        raise EOFError
    return b


_SEQ_LUT = np.frombuffer(_SEQ.encode(), dtype=np.uint8)


def _decode_seq(packed: bytes, l_seq: int) -> str:
    """4-bit packed BAM sequence -> string (two bases per byte, high nibble first)."""
    a = np.frombuffer(packed, dtype=np.uint8)
    out = np.empty(2 * len(a), dtype=np.uint8)
    out[0::2] = _SEQ_LUT[a >> 4]
    out[1::2] = _SEQ_LUT[a & 15]
    return out[:l_seq].tobytes().decode()


def read_bam(path: Union[str, Path], primary_only: bool = False) -> Iterator[Alignment]:
    """Yield the alignments in a BAM file (with primary_only, skip unmapped, secondary and supplementary records
    without decoding them)."""
    skip = UNMAPPED | SECONDARY | SUPPLEMENTARY if primary_only else 0
    with gzip.open(path, "rb") as fh:
        if _read_exact(fh, 4) != b"BAM\1":
            raise ValueError(f"{path}: not a BAM file")
        (l_text,) = struct.unpack("<i", _read_exact(fh, 4))
        fh.read(l_text)
        (n_ref,) = struct.unpack("<i", _read_exact(fh, 4))
        refs: List[str] = []
        for _ in range(n_ref):
            (l_name,) = struct.unpack("<i", _read_exact(fh, 4))
            refs.append(_read_exact(fh, l_name)[:-1].decode())
            fh.read(4)
        while True:
            try:
                (block,) = struct.unpack("<i", _read_exact(fh, 4))
            except EOFError:
                return
            rec = _read_exact(fh, block)
            ref_id, pos, l_name, mapq, _bin, n_cig, flag, l_seq = struct.unpack_from("<iiBBHHHi", rec, 0)
            if flag & skip:
                continue
            off = 32 + l_name
            cig = struct.unpack_from(f"<{n_cig}I", rec, off)
            off += 4 * n_cig
            seq = _decode_seq(rec[off:off + (l_seq + 1) // 2], l_seq)
            yield Alignment(
                ref=refs[ref_id] if 0 <= ref_id < len(refs) else "*",
                pos=pos, flag=flag, mapq=mapq,
                cigar=tuple((_CIGAR_OPS[c & 15], c >> 4) for c in cig),
                seq=seq,
            )
