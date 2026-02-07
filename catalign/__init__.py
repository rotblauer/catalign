"""
Catalign: Where sequences find each other like cats find the warmest spot.

⚠️  PROTOTYPE - NOT FOR PRODUCTION USE  ⚠️

This project is an experimental prototype developed using AI-assisted "vibe coding".
The methods and algorithms have NOT been validated for correctness or biological accuracy.
Results should NOT be used for clinical, research, or production purposes.

If you need a production aligner, use established tools like minimap2, BWA-MEM, or LAST.
"""

import warnings as _warnings

__version__ = "0.1.0-alpha"

# Show warning on first import
_warnings.warn(
    "catalign is a PROTOTYPE and should NOT be used for production. "
    "Results have not been validated for correctness. "
    "See README.md for details.",
    UserWarning,
    stacklevel=2
)

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
