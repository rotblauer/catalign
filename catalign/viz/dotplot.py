"""Dot plot visualization for sequence alignments."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from catalign.align import Alignment

try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


@dataclass
class DotPlot:
    """Container for dot plot data and visualization."""

    query_name: str
    target_name: str
    query_length: int
    target_length: int
    matches: List[Tuple[int, int]] = field(default_factory=list)
    mismatches: List[Tuple[int, int]] = field(default_factory=list)
    k: int = 1

    def to_figure(
        self,
        title: Optional[str] = None,
        show_mismatches: bool = False,
        color_matches: str = "blue",
        color_mismatches: str = "red",
        marker_size: int = 2,
        width: int = 800,
        height: int = 800,
    ) -> "go.Figure":
        """Create a Plotly figure from the dot plot data."""
        if not HAS_PLOTLY:
            raise ImportError("plotly is required for visualization")

        fig = go.Figure()

        if self.matches:
            match_x, match_y = zip(*self.matches)
            fig.add_trace(go.Scattergl(
                x=match_x, y=match_y, mode="markers",
                marker=dict(size=marker_size, color=color_matches),
                name="Matches",
                hovertemplate="Query: %{x}<br>Target: %{y}<extra></extra>",
            ))

        if show_mismatches and self.mismatches:
            mm_x, mm_y = zip(*self.mismatches)
            fig.add_trace(go.Scattergl(
                x=mm_x, y=mm_y, mode="markers",
                marker=dict(size=marker_size, color=color_mismatches),
                name="Mismatches",
            ))

        fig.update_layout(
            title=title or f"Dot Plot: {self.query_name} vs {self.target_name}",
            xaxis_title=f"Query: {self.query_name} ({self.query_length} bp)",
            yaxis_title=f"Target: {self.target_name} ({self.target_length} bp)",
            width=width, height=height,
        )
        fig.update_xaxes(range=[0, self.query_length])
        fig.update_yaxes(range=[0, self.target_length], scaleanchor="x", scaleratio=1)

        return fig

    def to_html(self, filepath: str, **kwargs) -> None:
        """Save dot plot as interactive HTML file."""
        fig = self.to_figure(**kwargs)
        fig.write_html(filepath)


def create_dotplot(
    query: str, target: str,
    query_name: str = "query", target_name: str = "target",
    k: int = 11, step: int = 1,
) -> DotPlot:
    """Create a dot plot from two sequences using k-mer matching."""
    if not HAS_NUMPY:
        raise ImportError("numpy is required")

    query, target = query.upper(), target.upper()

    target_kmers: dict[str, List[int]] = {}
    for i in range(0, len(target) - k + 1, step):
        kmer = target[i:i + k]
        if "N" not in kmer:
            target_kmers.setdefault(kmer, []).append(i)

    matches: List[Tuple[int, int]] = []
    for i in range(0, len(query) - k + 1, step):
        kmer = query[i:i + k]
        if kmer in target_kmers:
            for j in target_kmers[kmer]:
                matches.append((i, j))

    return DotPlot(
        query_name=query_name, target_name=target_name,
        query_length=len(query), target_length=len(target),
        matches=matches, k=k,
    )


def create_dotplot_from_alignment(
    alignment: "Alignment", query_seq: str, target_seq: str,
    include_background: bool = True, background_k: int = 11,
) -> DotPlot:
    """Create a dot plot highlighting the alignment path."""
    if include_background:
        dp = create_dotplot(query_seq, target_seq,
            query_name=alignment.query_name, target_name=alignment.target_name,
            k=background_k)
    else:
        dp = DotPlot(
            query_name=alignment.query_name, target_name=alignment.target_name,
            query_length=len(query_seq), target_length=len(target_seq),
        )

    for qpos, tpos, op in alignment.aligned_pairs:
        if qpos is not None and tpos is not None:
            if op == "M":
                dp.matches.append((qpos, tpos))
            elif op == "X":
                dp.mismatches.append((qpos, tpos))

    return dp
