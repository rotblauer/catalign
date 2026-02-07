"""CALI file format - Catalign Index for multi-scale alignment metrics.

Binary format specification:
- Efficient random access via tiled storage
- Multiple resolution levels (1kb, 10kb, 100kb, 1Mb)
- Pre-computed metrics for fast visualization
"""

from __future__ import annotations

import struct
import mmap
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, BinaryIO
import numpy as np

# Magic number for .cali files
CALI_MAGIC = b"CALI"
CALI_VERSION = 1

# Tile sizes in base pairs
DEFAULT_TILE_SIZES = [1_000, 10_000, 100_000, 1_000_000]


@dataclass
class TileMetrics:
    """Metrics for a single tile."""

    start: int
    end: int
    identity: float
    coverage: float
    energy: float
    gap_rate: float
    quality_score: float
    match_count: int
    mismatch_count: int
    insertion_count: int
    deletion_count: int

    def to_bytes(self) -> bytes:
        """Serialize to bytes."""
        return struct.pack(
            "<IIffffffff IIII",
            self.start,
            self.end,
            self.identity,
            self.coverage,
            self.energy,
            self.gap_rate,
            self.quality_score,
            0.0,  # reserved
            0.0,  # reserved
            0.0,  # reserved
            self.match_count,
            self.mismatch_count,
            self.insertion_count,
            self.deletion_count,
        )

    @classmethod
    def from_bytes(cls, data: bytes) -> "TileMetrics":
        """Deserialize from bytes."""
        values = struct.unpack("<IIffffffff IIII", data)
        return cls(
            start=values[0],
            end=values[1],
            identity=values[2],
            coverage=values[3],
            energy=values[4],
            gap_rate=values[5],
            quality_score=values[6],
            match_count=values[10],
            mismatch_count=values[11],
            insertion_count=values[12],
            deletion_count=values[13],
        )

    @staticmethod
    def byte_size() -> int:
        """Size in bytes."""
        return 56  # 2*4 + 8*4 + 4*4


@dataclass
class ChromosomeInfo:
    """Information about a chromosome in the CALI file."""

    name: str
    length: int
    data_offset: int = 0
    data_length: int = 0
    tile_counts: Dict[int, int] = field(default_factory=dict)  # tile_size -> count

    def to_bytes(self) -> bytes:
        """Serialize to bytes."""
        name_bytes = self.name.encode("utf-8")[:31].ljust(32, b"\x00")

        # Encode tile counts (up to 8 levels)
        tile_data = b""
        for tile_size in DEFAULT_TILE_SIZES:
            count = self.tile_counts.get(tile_size, 0)
            tile_data += struct.pack("<II", tile_size, count)

        return name_bytes + struct.pack(
            "<QQQ",
            self.length,
            self.data_offset,
            self.data_length,
        ) + tile_data

    @classmethod
    def from_bytes(cls, data: bytes) -> "ChromosomeInfo":
        """Deserialize from bytes."""
        name = data[:32].rstrip(b"\x00").decode("utf-8")
        length, data_offset, data_length = struct.unpack("<QQQ", data[32:56])

        # Parse tile counts
        tile_counts = {}
        offset = 56
        for _ in range(4):
            if offset + 8 <= len(data):
                tile_size, count = struct.unpack("<II", data[offset:offset + 8])
                if tile_size > 0:
                    tile_counts[tile_size] = count
                offset += 8

        return cls(
            name=name,
            length=length,
            data_offset=data_offset,
            data_length=data_length,
            tile_counts=tile_counts,
        )

    @staticmethod
    def byte_size() -> int:
        """Size in bytes."""
        return 32 + 24 + 32  # name + offsets + tile counts


@dataclass
class CaliHeader:
    """Header for CALI file."""

    version: int = CALI_VERSION
    flags: int = 0
    num_chromosomes: int = 0
    num_samples: int = 1
    tile_sizes: List[int] = field(default_factory=lambda: DEFAULT_TILE_SIZES.copy())
    reference_name: str = ""
    sample_name: str = ""

    def to_bytes(self) -> bytes:
        """Serialize header to bytes."""
        ref_bytes = self.reference_name.encode("utf-8")[:63].ljust(64, b"\x00")
        sample_bytes = self.sample_name.encode("utf-8")[:63].ljust(64, b"\x00")

        # Encode tile sizes (up to 8)
        tile_bytes = b""
        for ts in self.tile_sizes[:8]:
            tile_bytes += struct.pack("<I", ts)
        tile_bytes = tile_bytes.ljust(32, b"\x00")

        return (
            CALI_MAGIC +
            struct.pack("<IIII", self.version, self.flags, self.num_chromosomes, self.num_samples) +
            tile_bytes +
            ref_bytes +
            sample_bytes +
            b"\x00" * 56  # Reserved
        )

    @classmethod
    def from_bytes(cls, data: bytes) -> "CaliHeader":
        """Deserialize header from bytes."""
        if data[:4] != CALI_MAGIC:
            raise ValueError("Invalid CALI file magic number")

        version, flags, num_chroms, num_samples = struct.unpack("<IIII", data[4:20])

        # Parse tile sizes
        tile_sizes = []
        for i in range(8):
            ts = struct.unpack("<I", data[20 + i*4:24 + i*4])[0]
            if ts > 0:
                tile_sizes.append(ts)

        ref_name = data[52:116].rstrip(b"\x00").decode("utf-8")
        sample_name = data[116:180].rstrip(b"\x00").decode("utf-8")

        return cls(
            version=version,
            flags=flags,
            num_chromosomes=num_chroms,
            num_samples=num_samples,
            tile_sizes=tile_sizes,
            reference_name=ref_name,
            sample_name=sample_name,
        )

    @staticmethod
    def byte_size() -> int:
        """Size in bytes."""
        return 256


class CaliWriter:
    """Writer for CALI files."""

    def __init__(
        self,
        filepath: Path | str,
        reference_name: str = "",
        sample_name: str = "",
        tile_sizes: Optional[List[int]] = None,
    ):
        self.filepath = Path(filepath)
        self.tile_sizes = tile_sizes or DEFAULT_TILE_SIZES.copy()

        self.header = CaliHeader(
            tile_sizes=self.tile_sizes,
            reference_name=reference_name,
            sample_name=sample_name,
        )

        self.chromosomes: List[ChromosomeInfo] = []
        self._chromosome_tiles: Dict[str, Dict[int, List[TileMetrics]]] = {}
        self._file: Optional[BinaryIO] = None

    def __enter__(self) -> "CaliWriter":
        self._file = open(self.filepath, "wb")
        # Write placeholder header
        self._file.write(b"\x00" * CaliHeader.byte_size())
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._file:
            self._finalize()
            self._file.close()

    def add_chromosome(self, name: str, length: int) -> None:
        """Add a chromosome to the index."""
        chrom = ChromosomeInfo(name=name, length=length)
        self.chromosomes.append(chrom)
        self._chromosome_tiles[name] = {ts: [] for ts in self.tile_sizes}

    def add_tile(self, chrom: str, tile_size: int, metrics: TileMetrics) -> None:
        """Add a tile's metrics."""
        if chrom not in self._chromosome_tiles:
            raise ValueError(f"Unknown chromosome: {chrom}")
        if tile_size not in self._chromosome_tiles[chrom]:
            raise ValueError(f"Unknown tile size: {tile_size}")

        self._chromosome_tiles[chrom][tile_size].append(metrics)

    def add_tiles_from_array(
        self,
        chrom: str,
        tile_size: int,
        starts: np.ndarray,
        ends: np.ndarray,
        identities: np.ndarray,
        coverages: np.ndarray,
        energies: np.ndarray,
        gap_rates: np.ndarray,
        quality_scores: np.ndarray,
    ) -> None:
        """Add multiple tiles from numpy arrays (efficient bulk add)."""
        for i in range(len(starts)):
            metrics = TileMetrics(
                start=int(starts[i]),
                end=int(ends[i]),
                identity=float(identities[i]),
                coverage=float(coverages[i]),
                energy=float(energies[i]),
                gap_rate=float(gap_rates[i]),
                quality_score=float(quality_scores[i]),
                match_count=0,
                mismatch_count=0,
                insertion_count=0,
                deletion_count=0,
            )
            self.add_tile(chrom, tile_size, metrics)

    def _finalize(self) -> None:
        """Write all data and finalize the file."""
        if not self._file:
            return

        # First pass: calculate chromosome data positions
        # Chromosome table will be right after the header
        chrom_table_start = CaliHeader.byte_size()
        chrom_table_size = len(self.chromosomes) * ChromosomeInfo.byte_size()
        data_start = chrom_table_start + chrom_table_size

        current_offset = data_start
        for chrom in self.chromosomes:
            chrom.data_offset = current_offset
            # Calculate size: 4 bytes per tile count + tile data
            data_size = 0
            for tile_size in self.tile_sizes:
                tiles = self._chromosome_tiles[chrom.name][tile_size]
                chrom.tile_counts[tile_size] = len(tiles)
                data_size += 4 + len(tiles) * TileMetrics.byte_size()
            chrom.data_length = data_size
            current_offset += data_size

        # Write chromosome table
        self._file.seek(chrom_table_start)
        for chrom in self.chromosomes:
            self._file.write(chrom.to_bytes())

        # Write chromosome data
        for chrom in self.chromosomes:
            self._file.seek(chrom.data_offset)
            for tile_size in self.tile_sizes:
                tiles = self._chromosome_tiles[chrom.name][tile_size]
                # Write tile count for this level
                self._file.write(struct.pack("<I", len(tiles)))
                # Write tiles
                for tile in tiles:
                    self._file.write(tile.to_bytes())


        # Update and write header
        self.header.num_chromosomes = len(self.chromosomes)
        self._file.seek(0)
        self._file.write(self.header.to_bytes())


class CaliFile:
    """Reader for CALI files with memory-mapped access."""

    def __init__(self, filepath: Path | str):
        self.filepath = Path(filepath)
        self._file: Optional[BinaryIO] = None
        self._mmap: Optional[mmap.mmap] = None
        self.header: Optional[CaliHeader] = None
        self.chromosomes: Dict[str, ChromosomeInfo] = {}

        self._open()

    def _open(self) -> None:
        """Open the file and read header."""
        self._file = open(self.filepath, "rb")
        self._mmap = mmap.mmap(self._file.fileno(), 0, access=mmap.ACCESS_READ)

        # Read header
        header_data = self._mmap[:CaliHeader.byte_size()]
        self.header = CaliHeader.from_bytes(header_data)

        # Read chromosome table (at end of file)
        # For now, scan through the file to find chromosomes
        self._read_chromosome_table()

    def _read_chromosome_table(self) -> None:
        """Read the chromosome table (stored right after header)."""
        offset = CaliHeader.byte_size()

        for _ in range(self.header.num_chromosomes):
            # Read chromosome info
            chrom_data = self._mmap[offset:offset + ChromosomeInfo.byte_size()]
            chrom = ChromosomeInfo.from_bytes(chrom_data)
            self.chromosomes[chrom.name] = chrom
            # Move to next chromosome entry (sequential in table)
            offset += ChromosomeInfo.byte_size()

    def close(self) -> None:
        """Close the file."""
        if self._mmap:
            self._mmap.close()
        if self._file:
            self._file.close()

    def __enter__(self) -> "CaliFile":
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def get_chromosomes(self) -> List[str]:
        """Get list of chromosome names."""
        return list(self.chromosomes.keys())

    def get_chromosome_length(self, chrom: str) -> int:
        """Get length of a chromosome."""
        return self.chromosomes[chrom].length

    def get_tiles(
        self,
        chrom: str,
        start: int,
        end: int,
        tile_size: int,
    ) -> List[TileMetrics]:
        """Get tiles for a region at a specific resolution."""
        if chrom not in self.chromosomes:
            raise ValueError(f"Unknown chromosome: {chrom}")

        chrom_info = self.chromosomes[chrom]

        # Find the appropriate tile level
        if tile_size not in chrom_info.tile_counts:
            # Find closest available tile size
            available = list(chrom_info.tile_counts.keys())
            if not available:
                return []
            tile_size = min(available, key=lambda x: abs(x - tile_size))

        # Read tiles (simplified - in real implementation would use index)
        tiles = []
        offset = chrom_info.data_offset

        for ts in self.header.tile_sizes:
            count_data = self._mmap[offset:offset + 4]
            count = struct.unpack("<I", count_data)[0]
            offset += 4

            if ts == tile_size:
                # Read relevant tiles
                for _ in range(count):
                    tile_data = self._mmap[offset:offset + TileMetrics.byte_size()]
                    tile = TileMetrics.from_bytes(tile_data)

                    # Check if tile overlaps region
                    if tile.end >= start and tile.start <= end:
                        tiles.append(tile)

                    offset += TileMetrics.byte_size()
            else:
                # Skip this level
                offset += count * TileMetrics.byte_size()

        return tiles

    def get_region_summary(
        self,
        chrom: str,
        start: int,
        end: int,
    ) -> Dict:
        """Get summary metrics for a region (auto-select resolution)."""
        region_size = end - start

        # Select appropriate tile size
        for tile_size in reversed(self.header.tile_sizes):
            if region_size / tile_size >= 10:  # At least 10 tiles
                break

        tiles = self.get_tiles(chrom, start, end, tile_size)

        if not tiles:
            return {}

        # Aggregate metrics
        identities = [t.identity for t in tiles if t.coverage > 0]
        coverages = [t.coverage for t in tiles]
        energies = [t.energy for t in tiles if t.coverage > 0]

        return {
            "chrom": chrom,
            "start": start,
            "end": end,
            "tile_size": tile_size,
            "num_tiles": len(tiles),
            "mean_identity": np.mean(identities) if identities else 0,
            "mean_coverage": np.mean(coverages) if coverages else 0,
            "mean_energy": np.mean(energies) if energies else 0,
            "min_identity": min(identities) if identities else 0,
            "max_identity": max(identities) if identities else 0,
        }

    def to_numpy(
        self,
        chrom: str,
        tile_size: int,
    ) -> Dict[str, np.ndarray]:
        """Get all tiles for a chromosome as numpy arrays."""
        chrom_info = self.chromosomes[chrom]
        tiles = self.get_tiles(chrom, 0, chrom_info.length, tile_size)

        return {
            "start": np.array([t.start for t in tiles]),
            "end": np.array([t.end for t in tiles]),
            "identity": np.array([t.identity for t in tiles]),
            "coverage": np.array([t.coverage for t in tiles]),
            "energy": np.array([t.energy for t in tiles]),
            "gap_rate": np.array([t.gap_rate for t in tiles]),
            "quality_score": np.array([t.quality_score for t in tiles]),
            "matches": np.array([t.match_count for t in tiles]),
            "mismatches": np.array([t.mismatch_count for t in tiles]),
            "insertions": np.array([t.insertion_count for t in tiles]),
            "deletions": np.array([t.deletion_count for t in tiles]),
        }
