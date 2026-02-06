"""Catalign: Where sequences find each other like cats find the warmest spot."""

__version__ = "0.1.0"

from catalign.energy import EnergyModel
from catalign.align import CatalignAligner, Alignment
from catalign.quality import evaluate_quality, AlignmentQuality
from catalign.io import read_fasta, write_fasta, Sequence
from catalign.metrics import compute_metrics, AlignmentMetrics

__all__ = [
    "CatalignAligner",
    "Alignment",
    "EnergyModel",
    "evaluate_quality",
    "AlignmentQuality",
    "read_fasta",
    "write_fasta",
    "Sequence",
    "compute_metrics",
    "AlignmentMetrics",
]
