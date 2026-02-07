"""BAM/CRAM to CALI converter."""

from __future__ import annotations
from pathlib import Path
from typing import Optional, List

try:
    import pysam
    HAS_PYSAM = True
except ImportError:
    HAS_PYSAM = False

from catalign.viewer.cali_format import CaliWriter, DEFAULT_TILE_SIZES
from catalign.viewer.metrics_tiler import MetricsTiler, AlignmentPosition


def bam_to_cali(
    bam_path: Path | str,
    output_path: Path | str,
    reference_path: Optional[Path | str] = None,
    tile_sizes: Optional[List[int]] = None,
    sample_name: str = "",
    min_mapq: int = 0,
    progress_callback=None,
) -> None:
    """Convert BAM file to CALI format.

    Parameters
    ----------
    bam_path : Path
        Input BAM file
    output_path : Path
        Output CALI file
    reference_path : Path, optional
        Reference FASTA (for sequence access)
    tile_sizes : list of int, optional
        Tile sizes to compute (default: [1000, 10000, 100000, 1000000])
    sample_name : str, optional
        Sample name for metadata
    min_mapq : int
        Minimum mapping quality filter
    progress_callback : callable, optional
        Progress callback function
    """
    if not HAS_PYSAM:
        raise ImportError("pysam is required for BAM conversion. Install with: pip install pysam")

    bam_path = Path(bam_path)
    output_path = Path(output_path)
    tile_sizes = tile_sizes or list(DEFAULT_TILE_SIZES)

    bam = pysam.AlignmentFile(str(bam_path), "rb")
    ref_name = reference_path.name if reference_path else bam_path.stem

    with CaliWriter(output_path, reference_name=ref_name, sample_name=sample_name or bam_path.stem, tile_sizes=tile_sizes) as writer:
        for chrom, length in zip(bam.references, bam.lengths):
            if progress_callback:
                progress_callback(f"Processing {chrom}...")

            writer.add_chromosome(chrom, length)
            tiler = MetricsTiler(length, tile_sizes)

            for read in bam.fetch(chrom):
                if read.is_unmapped or read.mapping_quality < min_mapq:
                    continue
                _process_read(read, tiler)

            for ts in tile_sizes:
                for tile in tiler.get_tiles(ts):
                    writer.add_tile(chrom, ts, tile)

    bam.close()


def cram_to_cali(
    cram_path: Path | str,
    reference_path: Path | str,
    output_path: Path | str,
    tile_sizes: Optional[List[int]] = None,
    sample_name: str = "",
    min_mapq: int = 0,
    progress_callback=None,
) -> None:
    """Convert CRAM file to CALI format.

    Parameters
    ----------
    cram_path : Path
        Input CRAM file
    reference_path : Path
        Reference FASTA (required for CRAM)
    output_path : Path
        Output CALI file
    tile_sizes : list of int, optional
        Tile sizes to compute
    sample_name : str, optional
        Sample name for metadata
    min_mapq : int
        Minimum mapping quality filter
    progress_callback : callable, optional
        Progress callback function
    """
    if not HAS_PYSAM:
        raise ImportError("pysam is required for CRAM conversion. Install with: pip install pysam")

    cram_path = Path(cram_path)
    reference_path = Path(reference_path)
    output_path = Path(output_path)
    tile_sizes = tile_sizes or list(DEFAULT_TILE_SIZES)

    cram = pysam.AlignmentFile(str(cram_path), "rc", reference_filename=str(reference_path))

    with CaliWriter(output_path, reference_name=reference_path.name, sample_name=sample_name or cram_path.stem, tile_sizes=tile_sizes) as writer:
        for chrom, length in zip(cram.references, cram.lengths):
            if progress_callback:
                progress_callback(f"Processing {chrom}...")

            writer.add_chromosome(chrom, length)
            tiler = MetricsTiler(length, tile_sizes)

            for read in cram.fetch(chrom):
                if read.is_unmapped or read.mapping_quality < min_mapq:
                    continue
                _process_read(read, tiler)

            for ts in tile_sizes:
                for tile in tiler.get_tiles(ts):
                    writer.add_tile(chrom, ts, tile)

    cram.close()


def _process_read(read, tiler: MetricsTiler) -> None:
    """Process a single read and add positions to tiler."""
    if read.cigartuples is None:
        return

    ref_pos = read.reference_start
    query_pos = 0
    query_seq = read.query_sequence or ""

    for op, length in read.cigartuples:
        if op == 0:  # M - match/mismatch
            for i in range(length):
                pos = AlignmentPosition(
                    query_pos=query_pos + i,
                    target_pos=ref_pos + i,
                    operation="M",
                    query_base=query_seq[query_pos + i] if query_pos + i < len(query_seq) else "N",
                )
                tiler.add_position(ref_pos + i, pos)
            ref_pos += length
            query_pos += length
        elif op == 1:  # I - insertion
            for i in range(length):
                pos = AlignmentPosition(query_pos=query_pos + i, target_pos=None, operation="I")
                if ref_pos > 0:
                    tiler.add_position(ref_pos - 1, pos)
            query_pos += length
        elif op == 2:  # D - deletion
            for i in range(length):
                pos = AlignmentPosition(query_pos=None, target_pos=ref_pos + i, operation="D")
                tiler.add_position(ref_pos + i, pos)
            ref_pos += length
        elif op == 3:  # N - skip
            ref_pos += length
        elif op == 4:  # S - soft clip
            query_pos += length
        elif op == 5:  # H - hard clip
            pass
        elif op == 7:  # = - sequence match
            for i in range(length):
                pos = AlignmentPosition(query_pos=query_pos + i, target_pos=ref_pos + i, operation="M")
                tiler.add_position(ref_pos + i, pos)
            ref_pos += length
            query_pos += length
        elif op == 8:  # X - sequence mismatch
            for i in range(length):
                pos = AlignmentPosition(query_pos=query_pos + i, target_pos=ref_pos + i, operation="X")
                tiler.add_position(ref_pos + i, pos)
            ref_pos += length
            query_pos += length
