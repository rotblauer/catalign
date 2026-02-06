#!/usr/bin/env python3
"""
Quick Mitochondrial Genome Demo

Downloads human mitochondrial reference and a sample variant,
then runs full Catalign analysis with visualizations.

The mitochondrial genome is only ~16.5kb so this runs quickly.

Usage:
    python examples/mito_demo.py
"""

import sys
import urllib.request
from pathlib import Path

# Add parent to path for development
sys.path.insert(0, str(Path(__file__).parent.parent))


# NCBI URLs for mitochondrial sequences
# NC_012920.1 - Human mitochondrial reference (rCRS)
MITO_REF_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text"

# KC911424.1 - A human mitochondrial variant sample
MITO_SAMPLE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=KC911424.1&rettype=fasta&retmode=text"


def download_fasta(url: str, dest: Path) -> None:
    """Download a FASTA file from URL."""
    print(f"Downloading: {dest.name}")
    try:
        with urllib.request.urlopen(url, timeout=30) as response:
            content = response.read().decode('utf-8')
            dest.write_text(content)

            # Count sequence length
            lines = content.strip().split('\n')
            seq = ''.join(l for l in lines if not l.startswith('>'))
            print(f"  Downloaded {len(seq):,} bp")
    except Exception as e:
        print(f"  Error: {e}")
        raise


def main():
    print("="*60)
    print("🧬 Catalign Mitochondrial Genome Demo")
    print("="*60)
    print()

    # Setup output directory
    output_dir = Path("mito_demo_output")
    data_dir = output_dir / "data"
    data_dir.mkdir(parents=True, exist_ok=True)

    # Download sequences
    ref_path = data_dir / "mito_reference.fa"
    sample_path = data_dir / "mito_sample.fa"

    if not ref_path.exists():
        print("\nDownloading mitochondrial reference (rCRS)...")
        download_fasta(MITO_REF_URL, ref_path)
    else:
        print(f"\nUsing cached reference: {ref_path}")

    if not sample_path.exists():
        print("Downloading sample mitochondrial sequence...")
        download_fasta(MITO_SAMPLE_URL, sample_path)
    else:
        print(f"Using cached sample: {sample_path}")

    # Import catalign
    print("\nLoading Catalign...")
    from catalign import CatalignAligner, EnergyModel, evaluate_quality, read_fasta
    from catalign.metrics import compute_metrics

    # Load sequences
    print("\nLoading sequences...")
    ref_records = list(read_fasta(ref_path))
    sample_records = list(read_fasta(sample_path))

    ref_name, ref_seq = ref_records[0]
    sample_name, sample_seq = sample_records[0]

    print(f"  Reference: {ref_name[:50]}... ({len(ref_seq):,} bp)")
    print(f"  Sample: {sample_name[:50]}... ({len(sample_seq):,} bp)")

    # Configure aligner
    print("\nConfiguring aligner...")
    energy_model = EnergyModel(
        match_energy=-2.0,
        mismatch_energy=3.0,
        gap_open_energy=5.0,
        gap_extend_energy=1.0,
        transition_energy=2.0,
    )

    aligner = CatalignAligner(energy_model=energy_model, k=11, w=20)

    # Run alignment
    print("\nRunning alignment...")
    import time
    start = time.perf_counter()

    alignment = aligner.align(sample_seq, ref_seq, "sample", "reference")

    elapsed = time.perf_counter() - start
    print(f"  Completed in {elapsed:.2f} seconds")

    # Evaluate quality
    print("\nEvaluating quality...")
    quality = evaluate_quality(alignment, sample_seq, ref_seq, energy_model)
    metrics = compute_metrics(alignment, sample_seq, ref_seq, quality)

    # Print results
    print("\n" + "="*60)
    print("ALIGNMENT RESULTS")
    print("="*60)
    print(f"  Alignment length: {metrics.alignment_length:,} bp")
    print(f"  Identity:         {metrics.identity*100:.2f}%")
    print(f"  Matches:          {metrics.matches:,}")
    print(f"  Mismatches:       {metrics.mismatches:,}")
    print(f"  Insertions:       {metrics.insertions:,}")
    print(f"  Deletions:        {metrics.deletions:,}")
    print(f"  Quality Score:    {metrics.quality_score:.1f}/100")
    print(f"  Energy:           {metrics.total_energy:.2f}")

    # Generate visualizations
    print("\n" + "="*60)
    print("GENERATING VISUALIZATIONS")
    print("="*60)

    try:
        from catalign.viz import (
            create_dotplot,
            create_energy_heatmap,
            create_identity_heatmap,
            create_quality_heatmap,
            visualize_cigar,
            create_alignment_view,
        )
        from catalign.viz.dotplot import create_dotplot_from_alignment

        viz_dir = output_dir / "visualizations"
        viz_dir.mkdir(exist_ok=True)

        # 1. Dot Plot
        print("\n1. Creating Dot Plot...")
        dp = create_dotplot_from_alignment(alignment, sample_seq, ref_seq, background_k=11)
        fig = dp.to_figure(title="Mitochondrial Genome Alignment - Dot Plot", marker_size=2)
        fig.write_html(str(viz_dir / "dotplot.html"))
        print(f"   Saved: {viz_dir / 'dotplot.html'}")

        # 2. Energy Landscape
        print("\n2. Creating Energy Landscape...")
        energy_fig = create_energy_heatmap(alignment, sample_seq, ref_seq, window_size=100, step=20)
        energy_fig.update_layout(title="Energy Landscape - Mitochondrial Alignment")
        energy_fig.write_html(str(viz_dir / "energy_landscape.html"))
        print(f"   Saved: {viz_dir / 'energy_landscape.html'}")

        # 3. Identity Heatmap
        print("\n3. Creating Identity Heatmap...")
        id_fig = create_identity_heatmap(alignment, window_size=100, step=20)
        id_fig.update_layout(title="Local Identity - Mitochondrial Alignment")
        id_fig.write_html(str(viz_dir / "identity_heatmap.html"))
        print(f"   Saved: {viz_dir / 'identity_heatmap.html'}")

        # 4. Quality Heatmap
        print("\n4. Creating Quality Assessment...")
        qual_fig = create_quality_heatmap(quality)
        qual_fig.update_layout(title="Quality Assessment - Mitochondrial Alignment")
        qual_fig.write_html(str(viz_dir / "quality_heatmap.html"))
        print(f"   Saved: {viz_dir / 'quality_heatmap.html'}")

        # 5. CIGAR Visualization
        print("\n5. Creating CIGAR Visualization...")
        cigar_viz = visualize_cigar(alignment.cigar)
        cigar_fig = cigar_viz.to_figure()
        cigar_fig.update_layout(title="CIGAR String Visualization")
        cigar_fig.write_html(str(viz_dir / "cigar_view.html"))
        print(f"   Saved: {viz_dir / 'cigar_view.html'}")

        # 6. Text Alignment
        print("\n6. Creating Text Alignment View...")
        view = create_alignment_view(alignment, sample_seq, ref_seq)
        text_aln = view.to_text(line_width=80)
        (viz_dir / "alignment.txt").write_text(text_aln)
        print(f"   Saved: {viz_dir / 'alignment.txt'}")

        # Save HTML version
        html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Mitochondrial Genome Alignment</title>
    <style>
        body {{ font-family: monospace; padding: 20px; }}
        .match {{ background-color: #c8e6c9; }}
        .mismatch {{ background-color: #ffcdd2; }}
        .gap {{ background-color: #e3f2fd; }}
        pre {{ white-space: pre-wrap; line-height: 1.5; }}
    </style>
</head>
<body>
    <h1>Mitochondrial Genome Alignment</h1>
    <p><strong>Identity:</strong> {metrics.identity*100:.2f}%</p>
    <p><strong>Matches:</strong> {metrics.matches:,} | <strong>Mismatches:</strong> {metrics.mismatches:,}</p>
    <hr>
    <pre>{text_aln[:50000]}</pre>
</body>
</html>"""
        (viz_dir / "alignment.html").write_text(html_content)
        print(f"   Saved: {viz_dir / 'alignment.html'}")

        print("\n" + "="*60)
        print("✅ All visualizations generated!")
        print("="*60)
        print(f"\nOpen the visualizations in your browser:")
        print(f"  {viz_dir.absolute()}")
        print()
        print("Recommended viewing order:")
        print("  1. dotplot.html        - Overview of sequence similarity")
        print("  2. identity_heatmap.html - Where are the differences?")
        print("  3. energy_landscape.html - Alignment quality profile")
        print("  4. quality_heatmap.html  - Detailed quality metrics")
        print("  5. cigar_view.html     - Alignment operations")
        print("  6. alignment.html      - Base-by-base alignment")

    except ImportError as e:
        print(f"\n⚠️  Visualization dependencies not installed: {e}")
        print("Install with: pip install catalign[viz]")
        print("\nAlignment completed successfully - metrics above are valid.")

    print("\n" + "="*60)
    print("Demo complete!")
    print("="*60)


if __name__ == "__main__":
    main()
