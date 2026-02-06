"""Catalign Visualization Suite.

Interactive visualization tools for alignment analysis:
- Dot plots
- Energy heatmaps
- CIGAR visualization
- Quality dashboards
- Alignment comparison
"""

from catalign.viz.dotplot import create_dotplot, DotPlot
from catalign.viz.heatmap import (
    create_energy_heatmap,
    create_identity_heatmap,
    create_quality_heatmap,
)
from catalign.viz.cigar_viz import visualize_cigar, CigarVisualization
from catalign.viz.alignment_view import AlignmentView, create_alignment_view
from catalign.viz.dashboard import launch_dashboard

__all__ = [
    "create_dotplot",
    "DotPlot",
    "create_energy_heatmap",
    "create_identity_heatmap",
    "create_quality_heatmap",
    "visualize_cigar",
    "CigarVisualization",
    "AlignmentView",
    "create_alignment_view",
    "launch_dashboard",
]
