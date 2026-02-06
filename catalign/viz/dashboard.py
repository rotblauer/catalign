"""Interactive Streamlit dashboard for alignment analysis."""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional, List, TYPE_CHECKING

# Check for streamlit availability
try:
    import streamlit as st
    HAS_STREAMLIT = True
except ImportError:
    HAS_STREAMLIT = False

try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False


def launch_dashboard(port: int = 8501, host: str = "localhost") -> None:
    """Launch the Catalign visualization dashboard.

    This starts a Streamlit server with the alignment analysis dashboard.

    Parameters
    ----------
    port : int
        Port to run the dashboard on
    host : str
        Host address
    """
    if not HAS_STREAMLIT:
        raise ImportError(
            "streamlit is required for the dashboard. "
            "Install with: pip install catalign[viz]"
        )

    import subprocess
    dashboard_path = Path(__file__).parent / "dashboard_app.py"

    subprocess.run([
        sys.executable, "-m", "streamlit", "run",
        str(dashboard_path),
        "--server.port", str(port),
        "--server.address", host,
    ])


def create_dashboard_content():
    """Create the main dashboard content (for use within Streamlit)."""
    if not HAS_STREAMLIT:
        raise ImportError("streamlit is required")

    st.set_page_config(
        page_title="Catalign Dashboard",
        page_icon="🐱",
        layout="wide",
    )

    st.title("🐱 Catalign Alignment Dashboard")
    st.markdown("*Where sequences find each other like cats find the warmest spot*")

    # Sidebar for inputs
    with st.sidebar:
        st.header("Input Sequences")

        input_method = st.radio(
            "Input method",
            ["Paste sequences", "Upload FASTA files"],
        )

        if input_method == "Paste sequences":
            query_seq = st.text_area(
                "Query sequence",
                height=100,
                placeholder="Paste DNA sequence here...",
            )
            target_seq = st.text_area(
                "Target sequence",
                height=100,
                placeholder="Paste DNA sequence here...",
            )
            query_name = "query"
            target_name = "target"
        else:
            query_file = st.file_uploader("Query FASTA", type=["fa", "fasta", "fna"])
            target_file = st.file_uploader("Target FASTA", type=["fa", "fasta", "fna"])
            query_seq, query_name = _parse_uploaded_fasta(query_file)
            target_seq, target_name = _parse_uploaded_fasta(target_file)

        st.header("Alignment Parameters")
        k = st.slider("K-mer size", min_value=5, max_value=31, value=15, step=2)
        w = st.slider("Window size", min_value=10, max_value=200, value=50, step=10)

        st.header("Energy Model")
        match_e = st.number_input("Match energy", value=-2.0, step=0.5)
        mismatch_e = st.number_input("Mismatch energy", value=3.0, step=0.5)
        gap_open = st.number_input("Gap open energy", value=5.0, step=0.5)
        gap_extend = st.number_input("Gap extend energy", value=1.0, step=0.5)

        run_alignment = st.button("🚀 Run Alignment", type="primary")

    # Main content area
    if run_alignment and query_seq and target_seq:
        _run_alignment_and_display(
            query_seq, target_seq,
            query_name, target_name,
            k, w,
            match_e, mismatch_e, gap_open, gap_extend,
        )
    elif not (query_seq and target_seq):
        st.info("👈 Enter sequences in the sidebar and click 'Run Alignment' to begin")

        # Show example
        with st.expander("📖 Quick Start Guide"):
            st.markdown("""
            ### How to use this dashboard
            
            1. **Input sequences**: Either paste DNA sequences directly or upload FASTA files
            2. **Adjust parameters**: Fine-tune k-mer size, window size, and energy model
            3. **Run alignment**: Click the button to compute the alignment
            4. **Explore results**: View alignment metrics, visualizations, and quality assessments
            
            ### Example sequences
            
            Try these test sequences:
            
            **Query:**
            ```
            ACGTACGTACGTACGTACGTACGTACGTACGT
            ```
            
            **Target:**
            ```
            ACGTACGTACATACGTACGTACGTACGTACGT
            ```
            
            This pair has a single SNP (G→A) at position 11.
            """)


def _parse_uploaded_fasta(uploaded_file) -> tuple:
    """Parse an uploaded FASTA file."""
    if uploaded_file is None:
        return "", ""

    content = uploaded_file.read().decode("utf-8")
    lines = content.strip().split("\n")

    name = ""
    seq_parts = []

    for line in lines:
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
        else:
            seq_parts.append(line)

    return "".join(seq_parts), name


def _run_alignment_and_display(
    query_seq: str,
    target_seq: str,
    query_name: str,
    target_name: str,
    k: int,
    w: int,
    match_e: float,
    mismatch_e: float,
    gap_open: float,
    gap_extend: float,
):
    """Run alignment and display results."""
    if not HAS_STREAMLIT:
        return

    from catalign import CatalignAligner, EnergyModel, evaluate_quality
    from catalign.viz.dotplot import create_dotplot, create_dotplot_from_alignment
    from catalign.viz.heatmap import (
        create_energy_heatmap,
        create_identity_heatmap,
        create_quality_heatmap,
        create_block_quality_chart,
    )
    from catalign.viz.cigar_viz import visualize_alignment_cigar
    from catalign.viz.alignment_view import create_alignment_view

    # Clean sequences
    query_seq = "".join(c.upper() for c in query_seq if c.upper() in "ACGTN")
    target_seq = "".join(c.upper() for c in target_seq if c.upper() in "ACGTN")

    # Progress indicator
    with st.spinner("Computing alignment..."):
        # Create aligner and run
        em = EnergyModel(
            match_energy=match_e,
            mismatch_energy=mismatch_e,
            gap_open_energy=gap_open,
            gap_extend_energy=gap_extend,
        )
        aligner = CatalignAligner(energy_model=em, k=k, w=w)

        alignment = aligner.align(query_seq, target_seq, query_name, target_name)
        quality = evaluate_quality(alignment, query_seq, target_seq, em)

    st.success("✅ Alignment complete!")

    # Summary metrics
    st.header("📊 Alignment Summary")

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Identity", f"{quality.overall_identity:.1%}")
    with col2:
        st.metric("Quality Score", f"{quality.quality_score:.1f}/100")
    with col3:
        st.metric("Total Energy", f"{quality.total_energy:.2f}")
    with col4:
        aln_length = len(alignment.aligned_pairs)
        st.metric("Alignment Length", f"{aln_length:,} bp")

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Matches", f"{quality.total_matches:,}")
    with col2:
        st.metric("Mismatches", f"{quality.total_mismatches:,}")
    with col3:
        st.metric("Insertions", f"{quality.total_insertions:,}")
    with col4:
        st.metric("Deletions", f"{quality.total_deletions:,}")

    # CIGAR string
    st.subheader("CIGAR String")
    cigar = alignment.cigar
    if len(cigar) > 200:
        st.code(cigar[:200] + f"... ({len(cigar)} total characters)")
    else:
        st.code(cigar)

    # Tabs for different visualizations
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "🎯 Dot Plot",
        "📈 Energy Landscape",
        "📊 Quality Metrics",
        "🔤 CIGAR View",
        "📝 Text Alignment",
    ])

    with tab1:
        st.subheader("Dot Plot")

        dot_col1, dot_col2 = st.columns([3, 1])
        with dot_col2:
            dot_k = st.slider("K-mer size for dots", 5, 21, 11, key="dotplot_k")
            show_mismatch = st.checkbox("Show mismatches", value=False)

        with dot_col1:
            try:
                dp = create_dotplot_from_alignment(
                    alignment, query_seq, target_seq,
                    include_background=True,
                    background_k=dot_k,
                )
                fig = dp.to_figure(show_mismatches=show_mismatch)
                st.plotly_chart(fig, use_container_width=True)
            except Exception as e:
                st.error(f"Could not create dot plot: {e}")

    with tab2:
        st.subheader("Energy & Identity Landscape")

        land_col1, land_col2 = st.columns([3, 1])
        with land_col2:
            window = st.slider("Window size", 10, 200, 50, key="energy_window")
            step = st.slider("Step size", 1, 50, 10, key="energy_step")

        with land_col1:
            try:
                energy_fig = create_energy_heatmap(
                    alignment, query_seq, target_seq,
                    window_size=window, step=step,
                )
                st.plotly_chart(energy_fig, use_container_width=True)

                identity_fig = create_identity_heatmap(
                    alignment, window_size=window, step=step,
                )
                st.plotly_chart(identity_fig, use_container_width=True)
            except Exception as e:
                st.error(f"Could not create landscape: {e}")

    with tab3:
        st.subheader("Multi-Scale Quality")

        try:
            qual_fig = create_quality_heatmap(quality)
            st.plotly_chart(qual_fig, use_container_width=True)

            if quality.block_qualities:
                block_fig = create_block_quality_chart(quality)
                st.plotly_chart(block_fig, use_container_width=True)

            # Region quality details
            if quality.region_quality:
                st.subheader("Region Quality")
                rq = quality.region_quality

                rcol1, rcol2, rcol3 = st.columns(3)
                with rcol1:
                    st.metric("Query Coverage", f"{rq.coverage_query:.1%}")
                with rcol2:
                    st.metric("Target Coverage", f"{rq.coverage_target:.1%}")
                with rcol3:
                    st.metric("Concordance", f"{rq.concordance:.1%}")

                st.write(f"**Aligned blocks:** {rq.num_blocks}")
                st.write(f"**Total aligned bases:** {rq.total_aligned_bases:,}")
        except Exception as e:
            st.error(f"Could not create quality visualization: {e}")

    with tab4:
        st.subheader("CIGAR Visualization")

        try:
            cigar_viz = visualize_alignment_cigar(alignment)
            fig = cigar_viz.to_figure()
            st.plotly_chart(fig, use_container_width=True)

            st.text(cigar_viz.to_text())
        except Exception as e:
            st.error(f"Could not create CIGAR visualization: {e}")

    with tab5:
        st.subheader("Text Alignment")

        try:
            view = create_alignment_view(alignment, query_seq, target_seq)

            line_width = st.slider("Line width", 40, 120, 60, key="text_width")

            # Show text alignment
            text_aln = view.to_text(line_width=line_width)
            st.code(text_aln, language=None)

            # Download button
            st.download_button(
                "📥 Download alignment",
                text_aln,
                file_name="alignment.txt",
                mime="text/plain",
            )
        except Exception as e:
            st.error(f"Could not create text alignment: {e}")


# Entry point for running as script
if __name__ == "__main__":
    if HAS_STREAMLIT:
        create_dashboard_content()
    else:
        print("Streamlit is required. Install with: pip install streamlit")
