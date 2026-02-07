"""Metrics tiling for multi-scale visualization.

Computes alignment metrics at multiple resolution levels for efficient
genome-wide visualization.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional
import numpy as np

from catalign.viewer.cali_format import TileMetrics, DEFAULT_TILE_SIZES


class TileLevel(Enum):
    """Pre-defined tile levels for visualization."""
    BASE = 1
    FINE = 1_000
    MEDIUM = 10_000
    COARSE = 100_000
    CHROMOSOME = 1_000_000


@dataclass
class AlignmentPosition:
    """A single position in an alignment."""
    query_pos: Optional[int]
    target_pos: Optional[int]
    operation: str  # M, X, I, D
    query_base: str = ""
    target_base: str = ""
    energy: float = 0.0


@dataclass
class TileAccumulator:
    """Accumulator for computing tile metrics incrementally."""
    start: int
    end: int
    matches: int = 0
    mismatches: int = 0
    insertions: int = 0
    deletions: int = 0
    total_energy: float = 0.0
    covered_positions: int = 0

    def add_position(self, pos: AlignmentPosition) -> None:
        """Add a position to the accumulator."""
        if pos.operation == "M":
            self.matches += 1
            self.covered_positions += 1
        elif pos.operation == "X":
            self.mismatches += 1
            self.covered_positions += 1
        elif pos.operation == "I":
            self.insertions += 1
        elif pos.operation == "D":
            self.deletions += 1
            self.covered_positions += 1
        self.total_energy += pos.energy

    def to_metrics(self) -> TileMetrics:
        """Convert to TileMetrics."""
        tile_size = self.end - self.start
        aligned = self.matches + self.mismatches
        gap_total = self.insertions + self.deletions

        return TileMetrics(
            start=self.start,
            end=self.end,
            identity=self.matches / max(aligned, 1),
            coverage=self.covered_positions / max(tile_size, 1),
            energy=self.total_energy,
            gap_rate=gap_total / max(self.covered_positions + self.insertions, 1),
            quality_score=self._compute_quality_score(),
            match_count=self.matches,
            mismatch_count=self.mismatches,
            insertion_count=self.insertions,
            deletion_count=self.deletions,
        )

    def _compute_quality_score(self) -> float:
        """Compute a quality score (0-100)."""
        if self.covered_positions == 0:
            return 0.0
        aligned = self.matches + self.mismatches
        identity = self.matches / max(aligned, 1)
        gap_rate = (self.insertions + self.deletions) / max(self.covered_positions + self.insertions, 1)
        return max(0, min(100, identity * 100 * (1 - gap_rate * 0.5)))


class MetricsTiler:
    """Computes multi-scale metrics tiles from alignment data."""

    def __init__(
        self,
        chromosome_length: int,
        tile_sizes: Optional[List[int]] = None,
    ):
        self.chromosome_length = chromosome_length
        self.tile_sizes = tile_sizes or list(DEFAULT_TILE_SIZES)

        # Initialize accumulators for each tile size
        self._accumulators: Dict[int, List[TileAccumulator]] = {}
        for ts in self.tile_sizes:
            num_tiles = (chromosome_length + ts - 1) // ts
            self._accumulators[ts] = [
                TileAccumulator(start=i * ts, end=min((i + 1) * ts, chromosome_length))
                for i in range(num_tiles)
            ]

    def add_position(self, target_pos: int, pos: AlignmentPosition) -> None:
        """Add a position to all relevant tiles."""
        for ts, accumulators in self._accumulators.items():
            tile_idx = target_pos // ts
            if 0 <= tile_idx < len(accumulators):
                accumulators[tile_idx].add_position(pos)

    def get_tiles(self, tile_size: int) -> List[TileMetrics]:
        """Get all tiles at a specific resolution."""
        if tile_size not in self._accumulators:
            raise ValueError(f"Unknown tile size: {tile_size}")
        return [acc.to_metrics() for acc in self._accumulators[tile_size]]

    def get_all_tiles(self) -> Dict[int, List[TileMetrics]]:
        """Get tiles at all resolution levels."""
        return {ts: self.get_tiles(ts) for ts in self.tile_sizes}

    def to_numpy(self, tile_size: int) -> Dict[str, np.ndarray]:
        """Get tiles as numpy arrays for efficient processing."""
        tiles = self.get_tiles(tile_size)
        return {
            "start": np.array([t.start for t in tiles], dtype=np.int64),
            "end": np.array([t.end for t in tiles], dtype=np.int64),
            "identity": np.array([t.identity for t in tiles], dtype=np.float32),
            "coverage": np.array([t.coverage for t in tiles], dtype=np.float32),
            "energy": np.array([t.energy for t in tiles], dtype=np.float32),
            "gap_rate": np.array([t.gap_rate for t in tiles], dtype=np.float32),
            "quality_score": np.array([t.quality_score for t in tiles], dtype=np.float32),
            "matches": np.array([t.match_count for t in tiles], dtype=np.int32),
            "mismatches": np.array([t.mismatch_count for t in tiles], dtype=np.int32),
            "insertions": np.array([t.insertion_count for t in tiles], dtype=np.int32),
            "deletions": np.array([t.deletion_count for t in tiles], dtype=np.int32),
        }


def compute_tiles_from_alignment(
    alignment,
    query_seq: str,
    target_seq: str,
    target_length: int,
    tile_sizes: Optional[List[int]] = None,
) -> Dict[int, List[TileMetrics]]:
    """Compute multi-scale tiles from a Catalign alignment."""
    tiler = MetricsTiler(target_length, tile_sizes)

    for qpos, tpos, op in alignment.aligned_pairs:
        if tpos is not None:
            qb = query_seq[qpos] if qpos is not None and qpos < len(query_seq) else ""
            tb = target_seq[tpos] if tpos < len(target_seq) else ""
            pos = AlignmentPosition(
                query_pos=qpos,
                target_pos=tpos,
                operation=op,
                query_base=qb,
                target_base=tb,
            )
            tiler.add_position(tpos, pos)

    return tiler.get_all_tiles()
