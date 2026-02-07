"""Tests for the viewer module."""

import pytest
import tempfile
from pathlib import Path
import numpy as np


class TestCaliFormat:
    """Test the CALI file format."""

    def test_tile_metrics_serialization(self):
        """Test TileMetrics serialization."""
        from catalign.viewer.cali_format import TileMetrics

        tile = TileMetrics(
            start=1000,
            end=2000,
            identity=0.95,
            coverage=0.80,
            energy=-100.0,
            gap_rate=0.02,
            quality_score=92.5,
            match_count=950,
            mismatch_count=50,
            insertion_count=10,
            deletion_count=10,
        )

        # Serialize and deserialize
        data = tile.to_bytes()
        assert len(data) == TileMetrics.byte_size()

        restored = TileMetrics.from_bytes(data)
        assert restored.start == tile.start
        assert restored.end == tile.end
        assert abs(restored.identity - tile.identity) < 0.001
        assert abs(restored.coverage - tile.coverage) < 0.001

    def test_cali_writer_reader(self):
        """Test writing and reading CALI files."""
        from catalign.viewer.cali_format import CaliWriter, CaliFile, TileMetrics

        with tempfile.NamedTemporaryFile(suffix=".cali", delete=False) as f:
            test_path = Path(f.name)

        try:
            # Write file
            with CaliWriter(test_path, reference_name="test_ref", sample_name="test_sample") as writer:
                writer.add_chromosome("chr1", 10000)

                for i in range(10):
                    tile = TileMetrics(
                        start=i * 1000,
                        end=(i + 1) * 1000,
                        identity=0.90 + i * 0.01,
                        coverage=0.80,
                        energy=-100.0,
                        gap_rate=0.01,
                        quality_score=90.0 + i,
                        match_count=900 + i * 10,
                        mismatch_count=100 - i * 10,
                        insertion_count=5,
                        deletion_count=5,
                    )
                    writer.add_tile("chr1", 1000, tile)

            # Read file
            cali = CaliFile(test_path)

            assert cali.header.reference_name == "test_ref"
            assert cali.header.sample_name == "test_sample"
            assert len(cali.chromosomes) == 1
            assert "chr1" in cali.chromosomes

            cali.close()

        finally:
            test_path.unlink(missing_ok=True)

    def test_chromosome_info(self):
        """Test ChromosomeInfo serialization."""
        from catalign.viewer.cali_format import ChromosomeInfo

        chrom = ChromosomeInfo(
            name="chr1",
            length=248956422,
            data_offset=1024,
            data_length=102400,
        )

        data = chrom.to_bytes()
        assert len(data) == ChromosomeInfo.byte_size()

        restored = ChromosomeInfo.from_bytes(data)
        assert restored.name == chrom.name
        assert restored.length == chrom.length


class TestMetricsTiler:
    """Test the metrics tiler."""

    def test_tiler_basic(self):
        """Test basic tiling functionality."""
        from catalign.viewer.metrics_tiler import MetricsTiler, AlignmentPosition

        tiler = MetricsTiler(10000, tile_sizes=[1000, 5000])

        # Add some positions
        for i in range(100):
            pos = AlignmentPosition(
                query_pos=i,
                target_pos=i,
                operation="M",
            )
            tiler.add_position(i, pos)

        # Get tiles
        tiles_1k = tiler.get_tiles(1000)
        assert len(tiles_1k) == 10  # 10000 / 1000

        tiles_5k = tiler.get_tiles(5000)
        assert len(tiles_5k) == 2  # 10000 / 5000

        # Check first tile
        first_tile = tiles_1k[0]
        assert first_tile.start == 0
        assert first_tile.end == 1000
        assert first_tile.match_count == 100  # All 100 positions in first tile

    def test_tiler_with_gaps(self):
        """Test tiling with gaps."""
        from catalign.viewer.metrics_tiler import MetricsTiler, AlignmentPosition

        tiler = MetricsTiler(5000, tile_sizes=[1000])

        # Add matches
        for i in range(50):
            pos = AlignmentPosition(query_pos=i, target_pos=i, operation="M")
            tiler.add_position(i, pos)

        # Add mismatches
        for i in range(50, 80):
            pos = AlignmentPosition(query_pos=i, target_pos=i, operation="X")
            tiler.add_position(i, pos)

        # Add insertions
        for i in range(5):
            pos = AlignmentPosition(query_pos=80+i, target_pos=None, operation="I")
            tiler.add_position(79, pos)

        tiles = tiler.get_tiles(1000)
        first_tile = tiles[0]

        assert first_tile.match_count == 50
        assert first_tile.mismatch_count == 30
        assert first_tile.insertion_count == 5

    def test_to_numpy(self):
        """Test numpy array conversion."""
        from catalign.viewer.metrics_tiler import MetricsTiler, AlignmentPosition

        tiler = MetricsTiler(5000, tile_sizes=[1000])

        for i in range(100):
            pos = AlignmentPosition(query_pos=i, target_pos=i, operation="M")
            tiler.add_position(i, pos)

        arrays = tiler.to_numpy(1000)

        assert "start" in arrays
        assert "end" in arrays
        assert "identity" in arrays
        assert len(arrays["start"]) == 5


class TestBamConverter:
    """Test BAM converter (basic import test)."""

    def test_import(self):
        """Test module imports."""
        # Import should work even without pysam
        from catalign.viewer import bam_converter

        # Functions should exist
        assert hasattr(bam_converter, 'bam_to_cali')
        assert hasattr(bam_converter, 'cram_to_cali')
