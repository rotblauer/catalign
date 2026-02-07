"""Catalign Viewer - Multi-scale genome alignment visualization.

This module provides:
- CaliFile: Reader/writer for .cali multi-scale metrics files
- CaliWriter: Generate .cali files from alignments
- Integration with BAM/CRAM via pysam
- Tile-based metric computation for efficient visualization
"""


def __getattr__(name):
    """Lazy imports to avoid circular dependencies."""
    if name in ("CaliFile", "CaliWriter", "CaliHeader", "TileMetrics"):
        from catalign.viewer.cali_format import CaliFile, CaliWriter, CaliHeader, TileMetrics
        return locals()[name]
    elif name in ("MetricsTiler", "TileLevel", "AlignmentPosition"):
        from catalign.viewer.metrics_tiler import MetricsTiler, TileLevel, AlignmentPosition
        return locals()[name]
    elif name in ("bam_to_cali", "cram_to_cali"):
        from catalign.viewer.bam_converter import bam_to_cali, cram_to_cali
        return locals()[name]
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "CaliFile",
    "CaliWriter",
    "CaliHeader",
    "TileMetrics",
    "MetricsTiler",
    "TileLevel",
    "AlignmentPosition",
    "bam_to_cali",
    "cram_to_cali",
]
