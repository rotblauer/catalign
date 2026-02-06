"""Alignment view - detailed side-by-side alignment visualization."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from catalign.align import Alignment

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


@dataclass
class AlignmentView:
    """Container for detailed alignment visualization."""

    query_name: str
    target_name: str
    query_aligned: str
    target_aligned: str
    match_line: str
    query_seq: str = ""
    target_seq: str = ""

    def to_text(self, line_width: int = 60, show_positions: bool = True) -> str:
        """Generate text representation of the alignment."""
        lines = []
        lines.append(f"Query:  {self.query_name}")
        lines.append(f"Target: {self.target_name}")
        lines.append("")

        q_pos = 0
        t_pos = 0

        for i in range(0, len(self.query_aligned), line_width):
            chunk_q = self.query_aligned[i:i + line_width]
            chunk_t = self.target_aligned[i:i + line_width]
            chunk_m = self.match_line[i:i + line_width]

            q_consumed = sum(1 for c in chunk_q if c != "-")
            t_consumed = sum(1 for c in chunk_t if c != "-")

            q_start = q_pos + 1
            t_start = t_pos + 1
            q_end = q_pos + q_consumed
            t_end = t_pos + t_consumed

            if show_positions:
                lines.append(f"Query  {q_start:>6}  {chunk_q}  {q_end}")
                lines.append(f"              {chunk_m}")
                lines.append(f"Target {t_start:>6}  {chunk_t}  {t_end}")
            else:
                lines.append(f"Query:  {chunk_q}")
                lines.append(f"        {chunk_m}")
                lines.append(f"Target: {chunk_t}")

            lines.append("")
            q_pos = q_end
            t_pos = t_end

        return "\n".join(lines)

    def to_html(self, line_width: int = 80) -> str:
        """Generate HTML representation with color-coded alignment."""
        html_parts = [
            "<style>",
            ".alignment { font-family: monospace; white-space: pre; }",
            ".match { background-color: #c8e6c9; }",
            ".mismatch { background-color: #ffcdd2; }",
            ".gap { background-color: #e3f2fd; }",
            "</style>",
            "<div class='alignment'>",
        ]

        for i in range(0, len(self.query_aligned), line_width):
            chunk_q = self.query_aligned[i:i + line_width]
            chunk_t = self.target_aligned[i:i + line_width]
            chunk_m = self.match_line[i:i + line_width]

            q_html = "Query:  "
            for j, (qc, mc) in enumerate(zip(chunk_q, chunk_m)):
                if qc == "-":
                    q_html += f"<span class='gap'>{qc}</span>"
                elif mc == "|":
                    q_html += f"<span class='match'>{qc}</span>"
                else:
                    q_html += f"<span class='mismatch'>{qc}</span>"
            html_parts.append(q_html)

            html_parts.append(f"        {chunk_m}")

            t_html = "Target: "
            for j, (tc, mc) in enumerate(zip(chunk_t, chunk_m)):
                if tc == "-":
                    t_html += f"<span class='gap'>{tc}</span>"
                elif mc == "|":
                    t_html += f"<span class='match'>{tc}</span>"
                else:
                    t_html += f"<span class='mismatch'>{tc}</span>"
            html_parts.append(t_html)
            html_parts.append("")

        html_parts.append("</div>")
        return "\n".join(html_parts)

    def to_figure(
        self,
        start: int = 0,
        end: Optional[int] = None,
        width: int = 1000,
        height: int = 400,
    ) -> "go.Figure":
        """Create an interactive Plotly figure for the alignment."""
        if not HAS_PLOTLY:
            raise ImportError("plotly is required for visualization")

        if end is None:
            end = len(self.query_aligned)

        q_slice = self.query_aligned[start:end]
        t_slice = self.target_aligned[start:end]
        m_slice = self.match_line[start:end]

        colors_q = []
        colors_t = []
        for qc, tc, mc in zip(q_slice, t_slice, m_slice):
            if qc == "-":
                colors_q.append("#e3f2fd")
            elif mc == "|":
                colors_q.append("#c8e6c9")
            else:
                colors_q.append("#ffcdd2")

            if tc == "-":
                colors_t.append("#e3f2fd")
            elif mc == "|":
                colors_t.append("#c8e6c9")
            else:
                colors_t.append("#ffcdd2")

        fig = make_subplots(
            rows=2, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.05,
            subplot_titles=("Query", "Target"),
        )

        fig.add_trace(go.Heatmap(
            z=[[1] * len(q_slice)],
            x=list(range(start, end)),
            y=[self.query_name],
            text=[[c for c in q_slice]],
            texttemplate="%{text}",
            colorscale=[[0, "white"], [1, "white"]],
            showscale=False,
            hovertemplate="Pos: %{x}<br>Base: %{text}<extra></extra>",
        ), row=1, col=1)

        for i, (qc, col) in enumerate(zip(q_slice, colors_q)):
            fig.add_shape(
                type="rect",
                x0=start + i - 0.5, x1=start + i + 0.5,
                y0=-0.5, y1=0.5,
                fillcolor=col,
                line=dict(width=0),
                row=1, col=1,
            )

        fig.add_trace(go.Heatmap(
            z=[[1] * len(t_slice)],
            x=list(range(start, end)),
            y=[self.target_name],
            text=[[c for c in t_slice]],
            texttemplate="%{text}",
            colorscale=[[0, "white"], [1, "white"]],
            showscale=False,
            hovertemplate="Pos: %{x}<br>Base: %{text}<extra></extra>",
        ), row=2, col=1)

        for i, (tc, col) in enumerate(zip(t_slice, colors_t)):
            fig.add_shape(
                type="rect",
                x0=start + i - 0.5, x1=start + i + 0.5,
                y0=-0.5, y1=0.5,
                fillcolor=col,
                line=dict(width=0),
                row=2, col=1,
            )

        fig.update_layout(
            title="Alignment Detail View",
            width=width,
            height=height,
        )

        fig.update_xaxes(title_text="Alignment Position", row=2, col=1)

        return fig


def create_alignment_view(
    alignment: "Alignment",
    query_seq: str,
    target_seq: str,
) -> AlignmentView:
    """Create an alignment view from an Alignment object."""
    query_aligned = []
    target_aligned = []
    match_line = []

    for qpos, tpos, op in alignment.aligned_pairs:
        if op == "M":
            qb = query_seq[qpos] if qpos is not None and qpos < len(query_seq) else "N"
            tb = target_seq[tpos] if tpos is not None and tpos < len(target_seq) else "N"
            query_aligned.append(qb)
            target_aligned.append(tb)
            match_line.append("|")
        elif op == "X":
            qb = query_seq[qpos] if qpos is not None and qpos < len(query_seq) else "N"
            tb = target_seq[tpos] if tpos is not None and tpos < len(target_seq) else "N"
            query_aligned.append(qb)
            target_aligned.append(tb)
            match_line.append(".")
        elif op == "I":
            qb = query_seq[qpos] if qpos is not None and qpos < len(query_seq) else "N"
            query_aligned.append(qb)
            target_aligned.append("-")
            match_line.append(" ")
        elif op == "D":
            tb = target_seq[tpos] if tpos is not None and tpos < len(target_seq) else "N"
            query_aligned.append("-")
            target_aligned.append(tb)
            match_line.append(" ")

    return AlignmentView(
        query_name=alignment.query_name,
        target_name=alignment.target_name,
        query_aligned="".join(query_aligned),
        target_aligned="".join(target_aligned),
        match_line="".join(match_line),
        query_seq=query_seq,
        target_seq=target_seq,
    )


def compare_alignments(
    alignment1: "Alignment",
    alignment2: "Alignment",
    query_seq: str,
    target_seq: str,
    width: int = 1000,
    height: int = 600,
) -> "go.Figure":
    """Create a side-by-side comparison of two alignments."""
    if not HAS_PLOTLY:
        raise ImportError("plotly is required for visualization")

    view1 = create_alignment_view(alignment1, query_seq, target_seq)
    view2 = create_alignment_view(alignment2, query_seq, target_seq)

    fig = make_subplots(
        rows=4, cols=1,
        shared_xaxes=False,
        vertical_spacing=0.08,
        subplot_titles=(
            f"Alignment 1: Query ({alignment1.query_name})",
            f"Alignment 1: Target ({alignment1.target_name})",
            f"Alignment 2: Query ({alignment2.query_name})",
            f"Alignment 2: Target ({alignment2.target_name})",
        ),
    )

    _add_alignment_track(fig, view1.query_aligned, row=1)
    _add_alignment_track(fig, view1.target_aligned, row=2)
    _add_alignment_track(fig, view2.query_aligned, row=3)
    _add_alignment_track(fig, view2.target_aligned, row=4)

    fig.update_layout(
        title="Alignment Comparison",
        width=width,
        height=height,
        showlegend=False,
    )

    return fig


def _add_alignment_track(fig: "go.Figure", sequence: str, row: int) -> None:
    """Helper to add an alignment track to a figure."""
    fig.add_trace(go.Heatmap(
        z=[[ord(c) for c in sequence]],
        x=list(range(len(sequence))),
        y=[""],
        colorscale="Viridis",
        showscale=False,
        text=[[c for c in sequence]],
        texttemplate="%{text}",
        hovertemplate="Pos: %{x}<br>Base: %{text}<extra></extra>",
    ), row=row, col=1)
