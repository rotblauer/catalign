"""Heatmap visualizations for alignment energy and quality metrics."""

from __future__ import annotations
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from catalign.align import Alignment
    from catalign.quality import AlignmentQuality

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def _check_deps():
    if not HAS_PLOTLY:
        raise ImportError("plotly is required")
    if not HAS_NUMPY:
        raise ImportError("numpy is required")


def create_energy_heatmap(
    alignment: "Alignment", query_seq: str, target_seq: str,
    window_size: int = 50, step: int = 10, width: int = 900, height: int = 400,
) -> "go.Figure":
    """Create a heatmap showing alignment energy across the alignment."""
    _check_deps()
    from catalign.energy import EnergyModel
    em = EnergyModel()

    pairs = alignment.aligned_pairs
    if not pairs:
        fig = go.Figure()
        fig.update_layout(title="No alignment data")
        return fig

    energies = []
    for qpos, tpos, op in pairs:
        if op == "M":
            energies.append(em.match_energy)
        elif op == "X":
            qb = query_seq[qpos] if qpos is not None and qpos < len(query_seq) else "N"
            tb = target_seq[tpos] if tpos is not None and tpos < len(target_seq) else "N"
            energies.append(em.compute_base_energy(qb, tb))
        else:
            energies.append(em.gap_open_energy)

    n = len(energies)
    positions, windowed = [], []
    for i in range(0, n - window_size + 1, step):
        positions.append(i + window_size // 2)
        windowed.append(sum(energies[i:i + window_size]))

    energy_matrix = np.array(windowed).reshape(1, -1)

    fig = go.Figure(data=go.Heatmap(
        z=energy_matrix, x=positions, y=["Energy"],
        colorscale="RdBu_r", colorbar=dict(title="Energy"),
    ))
    fig.update_layout(
        title=f"Energy Landscape: {alignment.query_name} vs {alignment.target_name}",
        xaxis_title="Alignment Position", width=width, height=height,
    )
    return fig


def create_identity_heatmap(
    alignment: "Alignment", window_size: int = 50, step: int = 10,
    width: int = 900, height: int = 400,
) -> "go.Figure":
    """Create a heatmap showing local identity across the alignment."""
    _check_deps()

    pairs = alignment.aligned_pairs
    if not pairs:
        fig = go.Figure()
        fig.update_layout(title="No alignment data")
        return fig

    is_match = [1 if op == "M" else 0 for _, _, op in pairs]
    n = len(is_match)
    positions, identities = [], []

    for i in range(0, n - window_size + 1, step):
        positions.append(i + window_size // 2)
        identities.append(sum(is_match[i:i + window_size]) / window_size)

    identity_matrix = np.array(identities).reshape(1, -1)

    fig = go.Figure(data=go.Heatmap(
        z=identity_matrix, x=positions, y=["Identity"],
        colorscale="Viridis", zmin=0, zmax=1,
        colorbar=dict(title="Identity", tickformat=".0%"),
    ))
    fig.update_layout(
        title=f"Local Identity: {alignment.query_name} vs {alignment.target_name}",
        xaxis_title="Alignment Position", width=width, height=height,
    )
    return fig


def create_quality_heatmap(
    quality: "AlignmentQuality", width: int = 900, height: int = 500,
) -> "go.Figure":
    """Create a multi-track heatmap showing quality metrics."""
    _check_deps()

    if not quality.base_qualities:
        fig = go.Figure()
        fig.update_layout(title="No quality data")
        return fig

    n = len(quality.base_qualities)
    positions = list(range(n))

    op_map = {"M": 2, "X": 1, "I": 0, "D": 0}
    operations = [op_map.get(bq.operation, 0) for bq in quality.base_qualities]
    energies = [bq.energy for bq in quality.base_qualities]

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.1,
        subplot_titles=("Operation Type", "Per-Base Energy"))

    op_matrix = np.array(operations).reshape(1, -1)
    fig.add_trace(go.Heatmap(
        z=op_matrix, x=positions, y=["Op"],
        colorscale=[[0, "red"], [0.5, "yellow"], [1, "green"]], showscale=False,
    ), row=1, col=1)

    energy_matrix = np.array(energies).reshape(1, -1)
    fig.add_trace(go.Heatmap(
        z=energy_matrix, x=positions, y=["Energy"],
        colorscale="RdBu_r", colorbar=dict(title="Energy", x=1.02),
    ), row=2, col=1)

    fig.update_layout(title="Multi-Scale Quality Assessment", width=width, height=height)
    fig.update_xaxes(title_text="Alignment Position", row=2, col=1)

    return fig


def create_block_quality_chart(
    quality: "AlignmentQuality", width: int = 900, height: int = 400,
) -> "go.Figure":
    """Create a bar chart showing quality metrics for each alignment block."""
    _check_deps()

    blocks = quality.block_qualities
    if not blocks:
        fig = go.Figure()
        fig.update_layout(title="No block data")
        return fig

    block_labels = [f"Block {i+1}" for i in range(len(blocks))]
    identities = [b.identity * 100 for b in blocks]
    energies = [b.energy for b in blocks]

    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.15,
        subplot_titles=("Block Identity", "Block Energy"))

    fig.add_trace(go.Bar(x=block_labels, y=identities, name="Identity %",
        marker_color="steelblue"), row=1, col=1)
    fig.add_trace(go.Bar(x=block_labels, y=energies, name="Energy",
        marker_color=["green" if e < 0 else "red" for e in energies]), row=2, col=1)

    fig.update_layout(title="Block-Level Quality Metrics", width=width, height=height)
    fig.update_yaxes(title_text="Identity (%)", row=1, col=1)
    fig.update_yaxes(title_text="Energy", row=2, col=1)

    return fig
