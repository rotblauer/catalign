"""Sequence I/O – FASTA reading and writing (plain and gzipped)."""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from pathlib import Path
from typing import Generator, Iterable, Tuple, Union


@dataclass
class Sequence:
    """A named biological sequence."""

    name: str
    seq: str


def read_fasta(filepath: Union[str, Path]) -> Generator[Tuple[str, str], None, None]:
    """Yield (name, sequence) tuples from a FASTA file.

    Supports plain-text and gzip-compressed files (.gz).
    """
    filepath = Path(filepath)
    opener = gzip.open if filepath.suffix == ".gz" else open
    mode = "rt"

    name: str | None = None
    parts: list[str] = []

    with opener(filepath, mode) as fh:  # type: ignore[arg-type]
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            yield name, "".join(parts)


def write_fasta(
    filepath: Union[str, Path],
    sequences: Iterable[Union[Tuple[str, str], Sequence]],
    line_width: int = 80,
) -> None:
    """Write sequences to a FASTA file.

    *sequences* can be an iterable of ``(name, seq)`` tuples or
    ``Sequence`` objects.  Supports gzip when *filepath* ends with ``.gz``.
    """
    filepath = Path(filepath)
    opener = gzip.open if filepath.suffix == ".gz" else open
    mode = "wt"

    with opener(filepath, mode) as fh:  # type: ignore[arg-type]
        for item in sequences:
            if isinstance(item, Sequence):
                name, seq = item.name, item.seq
            else:
                name, seq = item
            fh.write(f">{name}\n")
            if seq:
                for i in range(0, len(seq), line_width):
                    fh.write(seq[i : i + line_width] + "\n")
            else:
                fh.write("\n")
