"""CIGAR string visualization."""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, TYPE_CHECKING
import re

if TYPE_CHECKING:
    from catalign.align import Alignment

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


CIGAR_COLORS = {
    "M": "#4CAF50", "X": "#FF9800", "I": "#2196F3", "D": "#F44336",
    "=": "#4CAF50", "S": "#9E9E9E", "H": "#607D8B", "N": "#9C27B0",
}


@dataclass
class CigarBlock:
    """A single CIGAR operation block."""
    operation: str
    length: int
    query_start: int
    query_end: Optional[int]
    target_start: int
    target_end: Optional[int]


@dataclass
class CigarVisualization:
    """Container for CIGAR visualization data."""
    cigar_string: str
    blocks: List[CigarBlock] = field(default_factory=list)
    query_length: int = 0
    target_length: int = 0
    query_name: str = "query"
    target_name: str = "target"

    def to_figure(self, width: int = 1000, height: int = 300, show_legend: bool = True) -> "go.Figure":
        """Create a Plotly figure visualizing the CIGAR string."""
        if not HAS_PLOTLY:
            raise ImportError("plotly is required")

        fig = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.1,
            subplot_titles=("Query Alignment", "Target Alignment"), row_heights=[0.5, 0.5])

        for block in self.blocks:
            color = CIGAR_COLORS.get(block.operation, "#CCCCCC")
            if block.query_end is not None:
                fig.add_trace(go.Bar(
                    x=[block.query_end - block.query_start], y=["Query"],
                    orientation="h", base=block.query_start,
                    marker=dict(color=color), name=block.operation,
                    showlegend=False,
                ), row=1, col=1)
            if block.target_end is not None:
                fig.add_trace(go.Bar(
                    x=[block.target_end - block.target_start], y=["Target"],
                    orientation="h", base=block.target_start,
                    marker=dict(color=color), name=block.operation,
                    showlegend=False,
                ), row=2, col=1)

        fig.update_layout(
            title=f"CIGAR: {self.cigar_string[:50]}{'...' if len(self.cigar_string) > 50 else ''}",
            barmode="stack", width=width, height=height, showlegend=show_legend,
        )
        fig.update_xaxes(title_text="Position (bp)", row=2, col=1)
        return fig

    def to_text(self) -> str:
        """Generate a text-based representation."""
        lines = [f"CIGAR: {self.cigar_string}", ""]
        op_counts = {}
        for block in self.blocks:
            op_counts[block.operation] = op_counts.get(block.operation, 0) + block.length
        lines.append("Operation Summary:")
        for op, count in sorted(op_counts.items()):
            lines.append(f"  {op}: {count} bp")
        return "\n".join(lines)


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    """Parse a CIGAR string into (length, operation) tuples."""
    return [(int(m.group(1)), m.group(2)) for m in re.finditer(r"(\d+)([MIDNSHP=X])", cigar)]


def visualize_cigar(
    cigar: str, query_name: str = "query", target_name: str = "target",
    query_start: int = 0, target_start: int = 0,
) -> CigarVisualization:
    """Create a visualization from a CIGAR string."""
    ops = parse_cigar(cigar)
    blocks = []
    q_pos, t_pos = query_start, target_start

    for length, op in ops:
        block = CigarBlock(op, length, q_pos, q_pos, t_pos, t_pos)

        if op in ("M", "=", "X"):
            block.query_end = q_pos + length
            block.target_end = t_pos + length
            q_pos += length
            t_pos += length
        elif op in ("I", "S"):
            block.query_end = q_pos + length
            block.target_end = None
            q_pos += length
        elif op in ("D", "N"):
            block.query_end = None
            block.target_end = t_pos + length
            t_pos += length

        blocks.append(block)

    return CigarVisualization(
        cigar_string=cigar, blocks=blocks,
        query_length=q_pos - query_start, target_length=t_pos - target_start,
        query_name=query_name, target_name=target_name,
    )


def visualize_alignment_cigar(alignment: "Alignment") -> CigarVisualization:
    """Create a visualization from an Alignment object."""
    return visualize_cigar(
        alignment.cigar, alignment.query_name, alignment.target_name,
        alignment.query_start, alignment.target_start,
    )
