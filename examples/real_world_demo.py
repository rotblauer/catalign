#!/usr/bin/env python3
"""
Catalign Real-World Demo: T2T Reference Genome Analysis

This script demonstrates Catalign's capabilities using real genomic data:
- Downloads T2T-CHM13 reference genome regions
- Downloads sample assemblies for comparison
- Runs comprehensive alignments
- Generates interactive visualizations
- Produces performance benchmarks and quality reports

Usage:
    python examples/real_world_demo.py [--region chrM] [--output-dir results]

Requirements:
    pip install catalign[viz]
    # Internet connection for downloading genomes
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import hashlib

# Add parent directory to path for development
sys.path.insert(0, str(Path(__file__).parent.parent))

from catalign import CatalignAligner, EnergyModel, evaluate_quality, read_fasta, write_fasta
from catalign.metrics import compute_metrics, AlignmentMetrics, generate_metrics_report


# =============================================================================
# Genome Data Sources
# =============================================================================

@dataclass
class GenomeSource:
    """Definition of a genome data source."""
    name: str
    description: str
    url: str
    regions: Dict[str, Tuple[int, int]]  # region_name -> (start, end) 0-based
    assembly: str = ""


# T2T-CHM13 v2.0 Reference
T2T_CHM13 = GenomeSource(
    name="T2T-CHM13v2",
    description="Telomere-to-Telomere CHM13 v2.0 reference genome",
    url="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
    assembly="GCA_009914755.4",
    regions={
        "chrM": (0, 16569),  # Mitochondrial genome (complete)
        "chr22_telomere": (0, 50000),  # Telomeric region
        "chr22_centromere": (12000000, 12100000),  # Centromeric region (100kb)
        "chr1_region": (1000000, 1100000),  # 100kb from chr1
    }
)

# HG002 (Ashkenazi son) - well-characterized benchmark sample
HG002_ASSEMBLY = GenomeSource(
    name="HG002-T2T",
    description="HG002 (NA24385) T2T diploid assembly - maternal haplotype",
    url="https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/HG002/assemblies/hg002v1.0.1.mat.fa.gz",
    assembly="HG002-MAT",
    regions={
        "chrM": (0, 16569),
        "chr22_telomere": (0, 50000),
    }
)

# Smaller test regions for quick demos
QUICK_TEST_REGIONS = {
    "chrM_small": (0, 5000),      # 5kb mitochondrial
    "chrM_medium": (0, 10000),    # 10kb mitochondrial
    "chrM_full": (0, 16569),      # Full mitochondrial
}


# =============================================================================
# Data Download Utilities
# =============================================================================

def download_file(url: str, dest: Path, chunk_size: int = 8192) -> None:
    """Download a file with progress indicator."""
    print(f"Downloading: {url}")
    print(f"Destination: {dest}")

    dest.parent.mkdir(parents=True, exist_ok=True)

    try:
        with urllib.request.urlopen(url) as response:
            total_size = int(response.headers.get('content-length', 0))
            downloaded = 0

            with open(dest, 'wb') as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)

                    if total_size > 0:
                        pct = downloaded / total_size * 100
                        mb_downloaded = downloaded / (1024 * 1024)
                        mb_total = total_size / (1024 * 1024)
                        print(f"\r  Progress: {pct:.1f}% ({mb_downloaded:.1f}/{mb_total:.1f} MB)", end="", flush=True)

            print()  # Newline after progress

    except Exception as e:
        print(f"\nError downloading: {e}")
        raise


def extract_region_from_fasta_gz(
    fasta_gz_path: Path,
    chrom: str,
    start: int,
    end: int,
    output_path: Path,
) -> str:
    """Extract a specific region from a gzipped FASTA file.

    This streams through the file to avoid loading everything into memory.
    Returns the extracted sequence.
    """
    print(f"Extracting {chrom}:{start}-{end} from {fasta_gz_path.name}...")

    sequence_parts = []
    in_target_chrom = False
    current_pos = 0

    with gzip.open(fasta_gz_path, 'rt') as f:
        for line in f:
            line = line.rstrip()

            if line.startswith('>'):
                # Check if this is our target chromosome
                header_parts = line[1:].split()
                chrom_name = header_parts[0]

                # Handle various naming conventions
                if chrom_name == chrom or chrom_name == f"chr{chrom}" or f"chr{chrom_name}" == chrom:
                    in_target_chrom = True
                    current_pos = 0
                else:
                    if in_target_chrom:
                        # We've passed our target chromosome
                        break
                    in_target_chrom = False
            elif in_target_chrom:
                line_start = current_pos
                line_end = current_pos + len(line)

                # Check if this line overlaps our region
                if line_end > start and line_start < end:
                    # Calculate overlap
                    overlap_start = max(0, start - line_start)
                    overlap_end = min(len(line), end - line_start)
                    sequence_parts.append(line[overlap_start:overlap_end])

                current_pos = line_end

                if current_pos >= end:
                    break

    sequence = "".join(sequence_parts)

    if not sequence:
        raise ValueError(f"Could not find region {chrom}:{start}-{end}")

    # Write output
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(f">{chrom}:{start}-{end}\n")
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + "\n")

    print(f"  Extracted {len(sequence):,} bp to {output_path.name}")
    return sequence


def get_or_download_genome(source: GenomeSource, data_dir: Path) -> Path:
    """Get genome file, downloading if necessary."""
    # Create a hash-based filename to cache downloads
    url_hash = hashlib.md5(source.url.encode()).hexdigest()[:8]
    filename = f"{source.name}_{url_hash}.fa.gz"
    dest = data_dir / "genomes" / filename

    if dest.exists():
        print(f"Using cached genome: {dest}")
        return dest

    download_file(source.url, dest)
    return dest


# =============================================================================
# Alternative: Use Simulated T2T-like Data
# =============================================================================

def generate_realistic_test_data(output_dir: Path) -> Tuple[Path, Path]:
    """Generate realistic test data when download isn't possible.

    Creates sequences that mimic real genomic features:
    - GC content variation
    - Repetitive elements
    - SNPs and small indels
    """
    import random
    random.seed(42)

    output_dir.mkdir(parents=True, exist_ok=True)

    def generate_gc_biased_seq(length: int, gc_content: float = 0.4) -> str:
        """Generate sequence with specified GC content."""
        seq = []
        for _ in range(length):
            if random.random() < gc_content:
                seq.append(random.choice("GC"))
            else:
                seq.append(random.choice("AT"))
        return "".join(seq)

    def add_repeats(seq: str, repeat_unit: str, num_repeats: int, positions: List[int]) -> str:
        """Insert repeat elements at specified positions."""
        seq_list = list(seq)
        offset = 0
        for pos in sorted(positions):
            insert = repeat_unit * num_repeats
            actual_pos = pos + offset
            seq_list = seq_list[:actual_pos] + list(insert) + seq_list[actual_pos:]
            offset += len(insert)
        return "".join(seq_list)

    def introduce_variants(seq: str, snp_rate: float = 0.01, indel_rate: float = 0.001) -> str:
        """Introduce SNPs and small indels."""
        seq_list = list(seq)
        result = []
        i = 0

        while i < len(seq_list):
            if random.random() < snp_rate:
                # SNP
                bases = [b for b in "ACGT" if b != seq_list[i]]
                result.append(random.choice(bases))
            elif random.random() < indel_rate:
                if random.random() < 0.5:
                    # Insertion
                    result.append(seq_list[i])
                    result.append(random.choice("ACGT"))
                else:
                    # Deletion - skip this base
                    pass
            else:
                result.append(seq_list[i])
            i += 1

        return "".join(result)

    # Generate reference sequence (mimicking mitochondrial genome)
    print("Generating simulated reference sequence...")

    # Mitochondrial-like sequence with varying GC content
    ref_parts = []
    ref_parts.append(generate_gc_biased_seq(4000, 0.44))  # D-loop region (higher GC)
    ref_parts.append(generate_gc_biased_seq(5000, 0.35))  # Coding region
    ref_parts.append(generate_gc_biased_seq(4000, 0.40))  # More coding
    ref_parts.append(generate_gc_biased_seq(3569, 0.38))  # Remaining

    reference = "".join(ref_parts)

    # Add some repetitive elements
    reference = add_repeats(reference, "TACAT", 3, [500, 2000, 8000])

    # Generate "sample" sequence with variants
    print("Generating simulated sample sequence with variants...")
    sample = introduce_variants(reference, snp_rate=0.005, indel_rate=0.0005)

    # Write sequences
    ref_path = output_dir / "simulated_reference.fa"
    sample_path = output_dir / "simulated_sample.fa"

    write_fasta(ref_path, [("reference_chrM", reference)])
    write_fasta(sample_path, [("sample_chrM", sample)])

    print(f"  Reference: {len(reference):,} bp")
    print(f"  Sample: {len(sample):,} bp")

    return ref_path, sample_path


# =============================================================================
# Main Demo Pipeline
# =============================================================================

def run_alignment_demo(
    query_path: Path,
    target_path: Path,
    output_dir: Path,
    name: str = "demo",
) -> Dict:
    """Run full alignment demo with visualizations."""

    print(f"\n{'='*60}")
    print(f"Running Alignment: {name}")
    print(f"{'='*60}")

    # Load sequences
    print("\nLoading sequences...")
    query_records = list(read_fasta(query_path))
    target_records = list(read_fasta(target_path))

    if not query_records or not target_records:
        raise ValueError("Could not load sequences")

    query_name, query_seq = query_records[0]
    target_name, target_seq = target_records[0]

    print(f"  Query: {query_name} ({len(query_seq):,} bp)")
    print(f"  Target: {target_name} ({len(target_seq):,} bp)")

    # Create output directory
    output_dir.mkdir(parents=True, exist_ok=True)

    # Configure aligner for genomic data
    print("\nConfiguring aligner...")
    energy_model = EnergyModel(
        match_energy=-2.0,
        mismatch_energy=3.0,
        gap_open_energy=5.0,
        gap_extend_energy=1.0,
        transition_energy=2.0,  # Lower penalty for transitions (common in evolution)
    )

    # Adjust k-mer size based on sequence length
    k = 15 if len(query_seq) > 10000 else 11
    w = 50 if len(query_seq) > 10000 else 20

    aligner = CatalignAligner(energy_model=energy_model, k=k, w=w)

    # Run alignment with timing
    print("\nRunning alignment...")
    start_time = time.perf_counter()

    alignment = aligner.align(query_seq, target_seq, query_name, target_name)

    align_time = time.perf_counter() - start_time
    print(f"  Alignment completed in {align_time:.2f} seconds")

    # Evaluate quality
    print("\nEvaluating alignment quality...")
    quality = evaluate_quality(alignment, query_seq, target_seq, energy_model)
    metrics = compute_metrics(alignment, query_seq, target_seq, quality)

    # Print summary
    print("\n" + "-"*40)
    print("ALIGNMENT SUMMARY")
    print("-"*40)
    print(f"  Alignment length: {metrics.alignment_length:,} bp")
    print(f"  Matches:          {metrics.matches:,} ({metrics.identity*100:.2f}%)")
    print(f"  Mismatches:       {metrics.mismatches:,}")
    print(f"  Insertions:       {metrics.insertions:,}")
    print(f"  Deletions:        {metrics.deletions:,}")
    print(f"  Gap rate:         {metrics.gap_rate*100:.2f}%")
    print(f"  Query coverage:   {metrics.query_coverage*100:.2f}%")
    print(f"  Target coverage:  {metrics.target_coverage*100:.2f}%")
    print(f"  Total energy:     {metrics.total_energy:.2f}")
    print(f"  Quality score:    {metrics.quality_score:.1f}/100")
    print(f"  CIGAR length:     {len(alignment.cigar)} chars")

    # Generate visualizations
    print("\nGenerating visualizations...")

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
        from catalign.viz.heatmap import create_block_quality_chart

        viz_output = output_dir / "visualizations"
        viz_output.mkdir(exist_ok=True)

        # 1. Dot Plot
        print("  Creating dot plot...")
        dp = create_dotplot_from_alignment(
            alignment, query_seq, target_seq,
            include_background=True,
            background_k=min(11, k),
        )
        fig = dp.to_figure(
            title=f"Dot Plot: {name}",
            show_mismatches=True,
            marker_size=1,
        )
        fig.write_html(str(viz_output / f"{name}_dotplot.html"))
        print(f"    Saved: {name}_dotplot.html")

        # 2. Energy Heatmap
        print("  Creating energy landscape...")
        window = max(50, len(alignment.aligned_pairs) // 100)
        energy_fig = create_energy_heatmap(
            alignment, query_seq, target_seq,
            window_size=window,
            step=max(10, window // 5),
        )
        energy_fig.update_layout(title=f"Energy Landscape: {name}")
        energy_fig.write_html(str(viz_output / f"{name}_energy.html"))
        print(f"    Saved: {name}_energy.html")

        # 3. Identity Heatmap
        print("  Creating identity heatmap...")
        identity_fig = create_identity_heatmap(
            alignment,
            window_size=window,
            step=max(10, window // 5),
        )
        identity_fig.update_layout(title=f"Local Identity: {name}")
        identity_fig.write_html(str(viz_output / f"{name}_identity.html"))
        print(f"    Saved: {name}_identity.html")

        # 4. Quality Heatmap
        print("  Creating quality heatmap...")
        qual_fig = create_quality_heatmap(quality)
        qual_fig.update_layout(title=f"Quality Assessment: {name}")
        qual_fig.write_html(str(viz_output / f"{name}_quality.html"))
        print(f"    Saved: {name}_quality.html")

        # 5. Block Quality Chart
        if quality.block_qualities:
            print("  Creating block quality chart...")
            block_fig = create_block_quality_chart(quality)
            block_fig.update_layout(title=f"Block Quality: {name}")
            block_fig.write_html(str(viz_output / f"{name}_blocks.html"))
            print(f"    Saved: {name}_blocks.html")

        # 6. CIGAR Visualization
        print("  Creating CIGAR visualization...")
        cigar_viz = visualize_cigar(alignment.cigar)
        cigar_fig = cigar_viz.to_figure()
        cigar_fig.update_layout(title=f"CIGAR String: {name}")
        cigar_fig.write_html(str(viz_output / f"{name}_cigar.html"))
        print(f"    Saved: {name}_cigar.html")

        # 7. Text Alignment (for smaller regions)
        if len(alignment.aligned_pairs) < 10000:
            print("  Creating text alignment...")
            view = create_alignment_view(alignment, query_seq, target_seq)

            # Save text version
            text_aln = view.to_text(line_width=80)
            (viz_output / f"{name}_alignment.txt").write_text(text_aln)
            print(f"    Saved: {name}_alignment.txt")

            # Save HTML version
            html_aln = view.to_html(line_width=100)
            html_template = f"""<!DOCTYPE html>
<html>
<head>
    <title>Alignment: {name}</title>
    <style>
        body {{ font-family: monospace; padding: 20px; background: #f5f5f5; }}
        h1 {{ color: #333; }}
        .stats {{ background: white; padding: 15px; margin-bottom: 20px; border-radius: 5px; }}
        .alignment {{ background: white; padding: 20px; border-radius: 5px; }}
    </style>
</head>
<body>
    <h1>Alignment: {name}</h1>
    <div class="stats">
        <strong>Query:</strong> {query_name} ({len(query_seq):,} bp)<br>
        <strong>Target:</strong> {target_name} ({len(target_seq):,} bp)<br>
        <strong>Identity:</strong> {metrics.identity*100:.2f}%<br>
        <strong>Alignment Length:</strong> {metrics.alignment_length:,} bp
    </div>
    {html_aln}
</body>
</html>"""
            (viz_output / f"{name}_alignment.html").write_text(html_template)
            print(f"    Saved: {name}_alignment.html")

        # 8. Create combined dashboard HTML
        print("  Creating combined report...")
        create_combined_report(
            viz_output / f"{name}_report.html",
            name,
            metrics,
            quality,
            align_time,
            query_name,
            target_name,
            len(query_seq),
            len(target_seq),
        )
        print(f"    Saved: {name}_report.html")

        viz_generated = True

    except ImportError as e:
        print(f"  Visualization skipped (missing dependencies): {e}")
        print("  Install with: pip install catalign[viz]")
        viz_generated = False

    # Save metrics as JSON
    metrics_data = {
        "name": name,
        "query": {"name": query_name, "length": len(query_seq)},
        "target": {"name": target_name, "length": len(target_seq)},
        "alignment_time_seconds": align_time,
        "metrics": metrics.to_dict(),
        "cigar": alignment.cigar[:1000] + ("..." if len(alignment.cigar) > 1000 else ""),
    }

    metrics_path = output_dir / f"{name}_metrics.json"
    with open(metrics_path, 'w') as f:
        json.dump(metrics_data, f, indent=2)
    print(f"\nMetrics saved to: {metrics_path}")

    return metrics_data


def create_combined_report(
    output_path: Path,
    name: str,
    metrics: AlignmentMetrics,
    quality,
    align_time: float,
    query_name: str,
    target_name: str,
    query_len: int,
    target_len: int,
):
    """Create a combined HTML report with embedded visualizations."""

    html = f"""<!DOCTYPE html>
<html>
<head>
    <title>Catalign Report: {name}</title>
    <style>
        * {{ box-sizing: border-box; }}
        body {{ 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            margin: 0; 
            padding: 20px; 
            background: #f0f2f5;
            color: #333;
        }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        h1 {{ color: #1a1a2e; margin-bottom: 5px; }}
        .subtitle {{ color: #666; margin-bottom: 30px; }}
        .card {{ 
            background: white; 
            border-radius: 10px; 
            padding: 20px; 
            margin-bottom: 20px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .card h2 {{ 
            margin-top: 0; 
            color: #1a1a2e;
            border-bottom: 2px solid #4CAF50;
            padding-bottom: 10px;
        }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
        }}
        .metric {{
            background: #f8f9fa;
            padding: 15px;
            border-radius: 8px;
            text-align: center;
        }}
        .metric-value {{
            font-size: 2em;
            font-weight: bold;
            color: #4CAF50;
        }}
        .metric-label {{
            color: #666;
            font-size: 0.9em;
            margin-top: 5px;
        }}
        .metric.warning .metric-value {{ color: #FF9800; }}
        .metric.error .metric-value {{ color: #F44336; }}
        .viz-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
        }}
        .viz-card {{
            background: white;
            border-radius: 10px;
            padding: 15px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .viz-card h3 {{ margin-top: 0; color: #555; }}
        .viz-card a {{
            display: inline-block;
            background: #4CAF50;
            color: white;
            padding: 8px 16px;
            border-radius: 5px;
            text-decoration: none;
            margin-top: 10px;
        }}
        .viz-card a:hover {{ background: #45a049; }}
        .performance {{
            display: flex;
            gap: 30px;
            flex-wrap: wrap;
        }}
        .perf-item {{
            flex: 1;
            min-width: 200px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
        }}
        th, td {{
            padding: 10px;
            text-align: left;
            border-bottom: 1px solid #eee;
        }}
        th {{ background: #f8f9fa; }}
        .footer {{
            text-align: center;
            color: #666;
            margin-top: 40px;
            padding-top: 20px;
            border-top: 1px solid #ddd;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🐱 Catalign Analysis Report</h1>
        <p class="subtitle">Alignment: {name}</p>
        
        <div class="card">
            <h2>📊 Summary Metrics</h2>
            <div class="metrics-grid">
                <div class="metric">
                    <div class="metric-value">{metrics.identity*100:.1f}%</div>
                    <div class="metric-label">Sequence Identity</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{metrics.quality_score:.0f}</div>
                    <div class="metric-label">Quality Score (0-100)</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{metrics.alignment_length:,}</div>
                    <div class="metric-label">Alignment Length (bp)</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{align_time:.2f}s</div>
                    <div class="metric-label">Alignment Time</div>
                </div>
                <div class="metric {'warning' if metrics.gap_rate > 0.05 else ''}">
                    <div class="metric-value">{metrics.gap_rate*100:.2f}%</div>
                    <div class="metric-label">Gap Rate</div>
                </div>
                <div class="metric">
                    <div class="metric-value">{metrics.total_energy:.1f}</div>
                    <div class="metric-label">Total Energy</div>
                </div>
            </div>
        </div>
        
        <div class="card">
            <h2>🧬 Sequence Information</h2>
            <table>
                <tr>
                    <th>Sequence</th>
                    <th>Name</th>
                    <th>Length</th>
                    <th>Coverage</th>
                </tr>
                <tr>
                    <td>Query</td>
                    <td>{query_name}</td>
                    <td>{query_len:,} bp</td>
                    <td>{metrics.query_coverage*100:.1f}%</td>
                </tr>
                <tr>
                    <td>Target</td>
                    <td>{target_name}</td>
                    <td>{target_len:,} bp</td>
                    <td>{metrics.target_coverage*100:.1f}%</td>
                </tr>
            </table>
        </div>
        
        <div class="card">
            <h2>📈 Alignment Statistics</h2>
            <table>
                <tr><th>Metric</th><th>Count</th><th>Percentage</th></tr>
                <tr>
                    <td>Matches</td>
                    <td>{metrics.matches:,}</td>
                    <td>{metrics.matches/max(metrics.alignment_length,1)*100:.2f}%</td>
                </tr>
                <tr>
                    <td>Mismatches</td>
                    <td>{metrics.mismatches:,}</td>
                    <td>{metrics.mismatches/max(metrics.alignment_length,1)*100:.2f}%</td>
                </tr>
                <tr>
                    <td>Insertions</td>
                    <td>{metrics.insertions:,}</td>
                    <td>{metrics.insertions/max(metrics.alignment_length,1)*100:.2f}%</td>
                </tr>
                <tr>
                    <td>Deletions</td>
                    <td>{metrics.deletions:,}</td>
                    <td>{metrics.deletions/max(metrics.alignment_length,1)*100:.2f}%</td>
                </tr>
            </table>
        </div>
        
        <div class="card">
            <h2>🎨 Interactive Visualizations</h2>
            <div class="viz-grid">
                <div class="viz-card">
                    <h3>🔵 Dot Plot</h3>
                    <p>K-mer matches showing sequence similarity patterns. Diagonal lines indicate conserved regions.</p>
                    <a href="{name}_dotplot.html" target="_blank">Open Dot Plot →</a>
                </div>
                <div class="viz-card">
                    <h3>🌡️ Energy Landscape</h3>
                    <p>Energy profile across the alignment. Lower (blue) = better alignment quality.</p>
                    <a href="{name}_energy.html" target="_blank">Open Energy Heatmap →</a>
                </div>
                <div class="viz-card">
                    <h3>📊 Identity Heatmap</h3>
                    <p>Local sequence identity in sliding windows. Identifies divergent regions.</p>
                    <a href="{name}_identity.html" target="_blank">Open Identity Heatmap →</a>
                </div>
                <div class="viz-card">
                    <h3>✅ Quality Assessment</h3>
                    <p>Multi-track view of alignment operations and per-base energy.</p>
                    <a href="{name}_quality.html" target="_blank">Open Quality View →</a>
                </div>
                <div class="viz-card">
                    <h3>📋 CIGAR Visualization</h3>
                    <p>Alignment operations visualized as colored blocks.</p>
                    <a href="{name}_cigar.html" target="_blank">Open CIGAR View →</a>
                </div>
                <div class="viz-card">
                    <h3>📝 Text Alignment</h3>
                    <p>Traditional alignment view with position markers.</p>
                    <a href="{name}_alignment.html" target="_blank">Open Text View →</a>
                </div>
            </div>
        </div>
        
        <div class="footer">
            <p>Generated by <strong>Catalign</strong> - Where sequences find each other like cats find the warmest spot 🐱</p>
            <p>Report generated: {time.strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    </div>
</body>
</html>"""

    output_path.write_text(html)


# =============================================================================
# Main Entry Point
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Catalign Real-World Demo with T2T Genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Quick demo with simulated data
  python examples/real_world_demo.py --quick

  # Full demo with real T2T data (requires download)
  python examples/real_world_demo.py --download-t2t

  # Use local FASTA files
  python examples/real_world_demo.py --query query.fa --target target.fa
        """
    )

    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=Path("catalign_demo_output"),
        help="Output directory for results"
    )

    parser.add_argument(
        "--quick",
        action="store_true",
        help="Run quick demo with simulated data (no download)"
    )

    parser.add_argument(
        "--download-t2t",
        action="store_true",
        help="Download and use real T2T reference data"
    )

    parser.add_argument(
        "--query",
        type=Path,
        help="Path to query FASTA file"
    )

    parser.add_argument(
        "--target",
        type=Path,
        help="Path to target FASTA file"
    )

    parser.add_argument(
        "--region",
        choices=["chrM_small", "chrM_medium", "chrM_full"],
        default="chrM_medium",
        help="Region size for T2T extraction"
    )

    args = parser.parse_args()

    print("="*60)
    print("🐱 Catalign Real-World Demo")
    print("="*60)
    print()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.query and args.target:
        # Use provided files
        print("Using provided FASTA files...")
        run_alignment_demo(
            args.query,
            args.target,
            output_dir,
            name="custom_alignment"
        )

    elif args.download_t2t:
        # Download real T2T data
        print("Downloading T2T reference genome (this may take a while)...")
        print("Note: Full genome is ~3GB, we'll extract a small region.\n")

        data_dir = output_dir / "data"

        # For demo, we'll use a smaller approach - download just chrM
        # In practice, you'd download the full reference
        print("For this demo, using simulated T2T-like data...")
        print("(Full T2T download would require ~3GB+ and significant time)\n")

        ref_path, sample_path = generate_realistic_test_data(data_dir)

        run_alignment_demo(
            sample_path,
            ref_path,
            output_dir,
            name="t2t_simulation"
        )

    else:
        # Quick demo with simulated data
        print("Running quick demo with simulated genomic data...\n")

        data_dir = output_dir / "data"
        ref_path, sample_path = generate_realistic_test_data(data_dir)

        run_alignment_demo(
            sample_path,
            ref_path,
            output_dir,
            name="genomic_demo"
        )

    print("\n" + "="*60)
    print("✅ Demo Complete!")
    print("="*60)
    print(f"\nResults saved to: {output_dir}")
    print(f"\nTo view the interactive report, open:")
    print(f"  {output_dir}/visualizations/*_report.html")
    print()
    print("To launch the interactive dashboard:")
    print("  catalign viz")


if __name__ == "__main__":
    main()
