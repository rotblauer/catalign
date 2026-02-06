"""Sequence sketching module for long-range attraction (coarse-grained matching)."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple

import numpy as np


def _hash_kmer(kmer: str) -> int:
    """Simple polynomial rolling hash for a k-mer string."""
    h = 0
    for ch in kmer:
        h = (h * 31 + ord(ch)) & 0xFFFFFFFFFFFFFFFF
    return h


@dataclass
class Sketch:
    """A collection of (hash, position, kmer) entries representing a sequence sketch."""
    entries: List[Tuple[int, int, str]] = field(default_factory=list)

    def hashes(self) -> set:
        return {h for h, _, _ in self.entries}

    def hash_to_positions(self) -> dict:
        """Map hash -> list of (position, kmer)."""
        d: dict[int, List[Tuple[int, str]]] = {}
        for h, pos, kmer in self.entries:
            d.setdefault(h, []).append((pos, kmer))
        return d


def kmer_sketch(sequence: str, k: int = 15) -> Sketch:
    """Extract all k-mers with their positions from *sequence*."""
    seq = sequence.upper()
    entries = []
    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if "N" not in kmer:
            entries.append((_hash_kmer(kmer), i, kmer))
    return Sketch(entries=entries)


def minimizer_sketch(sequence: str, k: int = 15, w: int = 50) -> Sketch:
    """Extract minimizer sketch (one representative k-mer per window)."""
    seq = sequence.upper()
    if len(seq) < k:
        return Sketch()

    # Pre-compute all k-mer hashes
    n_kmers = len(seq) - k + 1
    hashes = []
    for i in range(n_kmers):
        kmer = seq[i : i + k]
        if "N" in kmer:
            hashes.append((0xFFFFFFFFFFFFFFFF, i, kmer))
        else:
            hashes.append((_hash_kmer(kmer), i, kmer))

    seen: set[int] = set()
    entries = []
    for win_start in range(max(n_kmers - w + 1, 1)):
        win_end = min(win_start + w, n_kmers)
        best = min(hashes[win_start:win_end], key=lambda x: x[0])
        h, pos, kmer = best
        if h == 0xFFFFFFFFFFFFFFFF:
            continue
        if pos not in seen:
            seen.add(pos)
            entries.append((h, pos, kmer))

    entries.sort(key=lambda x: x[1])
    return Sketch(entries=entries)


def find_anchors(
    sketch_a: Sketch, sketch_b: Sketch
) -> List[Tuple[int, int, str]]:
    """Find matching k-mers between two sketches.

    Returns list of (position_a, position_b, kmer) anchor triples,
    sorted by position_a.
    """
    b_map = sketch_b.hash_to_positions()
    anchors = []
    for h, pos_a, kmer_a in sketch_a.entries:
        if h in b_map:
            for pos_b, kmer_b in b_map[h]:
                if kmer_a == kmer_b:
                    anchors.append((pos_a, pos_b, kmer_a))
    anchors.sort(key=lambda x: (x[0], x[1]))
    return anchors
