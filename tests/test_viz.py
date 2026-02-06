"""Tests for visualization components."""

import pytest
from unittest.mock import MagicMock

from catalign import CatalignAligner, evaluate_quality
from catalign.metrics import compute_metrics


# Check if visualization dependencies are available
try:
    import plotly
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

try:
    import streamlit
    HAS_STREAMLIT = True
except ImportError:
    HAS_STREAMLIT = False


pytestmark = pytest.mark.viz


class TestDotPlot:
    """Tests for dot plot visualization."""

    @pytest.mark.skipif(not HAS_PLOTLY, reason="plotly not installed")
    def test_create_dotplot_basic(self):
        """Test basic dot plot creation."""
        from catalign.viz.dotplot import create_dotplot

        query = "ACGTACGTACGT"
        target = "ACGTACGTACGT"

        dp = create_dotplot(query, target, k=4)

        assert dp.query_length == len(query)
        assert dp.target_length == len(target)
        assert len(dp.matches) > 0

    @pytest.mark.skipif(not HAS_PLOTLY, reason="plotly not installed")
    def test_dotplot_to_figure(self):
        """Test dot plot figure generation."""
        from catalign.viz.dotplot import create_dotplot

        query = "ACGTACGTACGTACGT"
        target = "ACGTACGTACGTACGT"

        dp = create_dotplot(query, target, k=4)
        fig = dp.to_figure()

        assert fig is not None
        assert hasattr(fig, 'data')

    @pytest.mark.skipif(not HAS_PLOTLY, reason="plotly not installed")
    def test_dotplot_from_alignment(self):
        """Test dot plot creation from alignment."""
        from catalign.viz.dotplot import create_dotplot_from_alignment

        query = "ACGTACGTACGT"
        target = "ACGTACGTACGT"

        aligner = CatalignAligner()
        aln = aligner.align(query, target)

        dp = create_dotplot_from_alignment(aln, query, target, background_k=4)

        assert len(dp.matches) > 0


class TestHeatmaps:
    """Tests for heatmap visualizations."""

    @pytest.mark.skipif(not HAS_PLOTLY, reason="plotly not installed")
    def test_energy_heatmap(self):
        """Test energy heatmap creation."""
        from catalign.viz.heatmap import create_energy_heatmap

        query = "A" * 100
        target = "A" * 100

        aligner = CatalignAligner()
        aln = aligner.align(query, target)

        fig = create_energy_heatmap(aln, query, target, window_size=10, step=5)

        assert fig is not None

    @pytest.mark.skipif(not HAS_PLOTLY, reason="plotly not installed")
    def test_identity_heatmap(self):
        """Test identity heatmap creation."""
        from catalign.viz.heatmap import create_identity_heatmap

        query = "A" * 100
        target = "A" * 100

        aligner = CatalignAligner()
        aln = aligner.align(query, target)

        fig = create_identity_heatmap(aln, window_size=10, step=5)

        assert fig is not None

    @pytest.mark.skipif(not HAS_PLOTLY, reason="plotly not installed")
    def test_quality_heatmap(self):
        """Test quality heatmap creation."""
        from catalign.viz.heatmap import create_quality_heatmap

        query = "ACGT" * 25
        target = "ACGT" * 25

        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        quality = evaluate_quality(aln, query, target)

        fig = create_quality_heatmap(quality)

        assert fig is not None


class TestCigarVisualization:
    """Tests for CIGAR visualization."""

    @pytest.mark.skipif(not HAS_PLOTLY, reason="plotly not installed")
    def test_visualize_cigar(self):
        """Test CIGAR string visualization."""
        from catalign.viz.cigar_viz import visualize_cigar

        cigar = "10M2I5M3D8M"

        viz = visualize_cigar(cigar)

        assert viz.cigar_string == cigar
        assert len(viz.blocks) == 5  # 5 operations: M, I, M, D, M

    @pytest.mark.skipif(not HAS_PLOTLY, reason="plotly not installed")
    def test_cigar_to_figure(self):
        """Test CIGAR figure generation."""
        from catalign.viz.cigar_viz import visualize_cigar

        cigar = "10M2I5M"

        viz = visualize_cigar(cigar)
        fig = viz.to_figure()

        assert fig is not None

    def test_cigar_to_text(self):
        """Test CIGAR text representation."""
        from catalign.viz.cigar_viz import visualize_cigar

        cigar = "10M2I5M"

        viz = visualize_cigar(cigar)
        text = viz.to_text()

        assert "CIGAR:" in text
        assert "M:" in text
        assert "I:" in text


class TestAlignmentView:
    """Tests for alignment view visualization."""

    def test_create_alignment_view(self):
        """Test alignment view creation."""
        from catalign.viz.alignment_view import create_alignment_view

        query = "ACGTACGT"
        target = "ACGTACGT"

        aligner = CatalignAligner()
        aln = aligner.align(query, target)

        view = create_alignment_view(aln, query, target)

        assert len(view.query_aligned) > 0
        assert len(view.target_aligned) > 0
        assert len(view.match_line) > 0

    def test_alignment_view_to_text(self):
        """Test alignment view text output."""
        from catalign.viz.alignment_view import create_alignment_view

        query = "ACGTACGT"
        target = "ACGTACGT"

        aligner = CatalignAligner()
        aln = aligner.align(query, target)

        view = create_alignment_view(aln, query, target)
        text = view.to_text()

        assert "Query" in text
        assert "Target" in text

    def test_alignment_view_to_html(self):
        """Test alignment view HTML output."""
        from catalign.viz.alignment_view import create_alignment_view

        query = "ACGTACGT"
        target = "ACGTACGT"

        aligner = CatalignAligner()
        aln = aligner.align(query, target)

        view = create_alignment_view(aln, query, target)
        html = view.to_html()

        assert "<div" in html
        assert "alignment" in html


class TestDashboard:
    """Tests for dashboard functionality."""

    def test_dashboard_import(self):
        """Test that dashboard module can be imported."""
        from catalign.viz import dashboard

        assert hasattr(dashboard, 'launch_dashboard')

    @pytest.mark.skipif(not HAS_STREAMLIT, reason="streamlit not installed")
    def test_dashboard_content_function_exists(self):
        """Test dashboard content function."""
        from catalign.viz.dashboard import create_dashboard_content

        # Just verify the function exists and is callable
        assert callable(create_dashboard_content)
