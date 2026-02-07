#!/usr/bin/env python3
"""
HG002 vs T2T-CHM13 Genome-Wide Alignment Demo

This script demonstrates genome-wide alignment visualization by:
1. Downloading HG002 assembly and T2T-CHM13 reference (small regions)
2. Running Catalign on chromosome-scale data
3. Generating multi-scale metrics in CALI format
4. Creating interactive visualizations

Usage:
    python examples/hg002_t2t_demo.py [--chrom chrM] [--full-genome]
"""

import argparse
import gzip
import json
import os
import sys
import time
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Add parent directory
sys.path.insert(0, str(Path(__file__).parent.parent))


# =============================================================================
# Data Sources
# =============================================================================

# UCSC DAS server provides easy access to specific chromosome sequences
UCSC_DAS_URL = "https://genome.ucsc.edu/cgi-bin/das/{assembly}/dna?segment={chrom}:{start},{end}"

# NCBI E-utilities for smaller sequences
NCBI_FETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Test chromosomes (ordered by size for quick testing)
TEST_CHROMOSOMES = {
    "chrM": 16569,      # Mitochondrial - smallest, fastest
    "chr21": 46709983,  # Smallest autosome
    "chr22": 50818468,  # Second smallest
    "chr19": 58617616,  # Gene-dense
    "chr1": 248956422,  # Largest
}

# URLs for pre-aligned data (for demo purposes)
DEMO_DATA_URLS = {
    "hg002_chrM": "https://raw.githubusercontent.com/human-pangenomics/HPP_Year1_Assemblies/main/assembly_index/Year1_assemblies_v2_genbank.index",
}


def download_chromosome_sequence(
    assembly: str,
    chrom: str,
    start: int,
    end: int,
    output_path: Path,
) -> str:
    """Download a chromosome region from UCSC."""
    print(f"Downloading {assembly} {chrom}:{start:,}-{end:,}...")

    # For small sequences, use UCSC DAS
    url = UCSC_DAS_URL.format(assembly=assembly, chrom=chrom, start=start, end=end)

    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            content = response.read().decode("utf-8")

            # Parse DAS XML response
            import re
            seq_match = re.search(r"<DNA[^>]*>([^<]+)</DNA>", content, re.DOTALL)
            if seq_match:
                sequence = seq_match.group(1).replace("\n", "").replace(" ", "").upper()
            else:
                raise ValueError("Could not parse sequence from DAS response")

        # Write FASTA
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, "w") as f:
            f.write(f">{chrom}:{start}-{end}\n")
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")

        print(f"  Downloaded {len(sequence):,} bp")
        return sequence

    except Exception as e:
        print(f"  Error downloading: {e}")
        raise


def generate_simulated_assembly(
    reference_seq: str,
    mutation_rate: float = 0.001,
    indel_rate: float = 0.0001,
    sv_rate: float = 0.00001,
) -> str:
    """Generate a simulated assembly with realistic variation patterns."""
    import random
    random.seed(42)

    seq = list(reference_seq)
    result = []
    i = 0

    while i < len(seq):
        # Check for structural variant
        if random.random() < sv_rate and i + 1000 < len(seq):
            sv_type = random.choice(["del", "ins", "inv", "dup"])
            sv_size = random.randint(50, 500)

            if sv_type == "del":
                # Large deletion
                i += sv_size
                continue
            elif sv_type == "ins":
                # Large insertion
                result.append(seq[i])
                result.extend(random.choices("ACGT", k=sv_size))
            elif sv_type == "inv":
                # Inversion
                region = seq[i:i + sv_size]
                complement = {"A": "T", "T": "A", "G": "C", "C": "G", "N": "N"}
                inverted = [complement.get(b, "N") for b in reversed(region)]
                result.extend(inverted)
                i += sv_size
                continue
            elif sv_type == "dup":
                # Tandem duplication
                region = seq[i:i + sv_size]
                result.extend(region)
                result.extend(region)
                i += sv_size
                continue

        # Check for small indel
        if random.random() < indel_rate:
            if random.random() < 0.5:
                # Small insertion
                result.append(seq[i])
                result.extend(random.choices("ACGT", k=random.randint(1, 10)))
            else:
                # Small deletion
                i += random.randint(1, 10)
                continue

        # Check for SNP
        if random.random() < mutation_rate:
            # SNP - transition more likely than transversion
            base = seq[i]
            if random.random() < 0.7:  # Transition
                transitions = {"A": "G", "G": "A", "C": "T", "T": "C"}
                result.append(transitions.get(base, random.choice("ACGT")))
            else:  # Transversion
                others = [b for b in "ACGT" if b != base]
                result.append(random.choice(others))
        else:
            result.append(seq[i])

        i += 1

    return "".join(result)


def run_chromosome_alignment(
    query_seq: str,
    target_seq: str,
    chrom_name: str,
    output_dir: Path,
) -> Dict:
    """Run Catalign alignment on chromosome-scale data."""
    from catalign import CatalignAligner, EnergyModel, evaluate_quality
    from catalign.metrics import compute_metrics
    from catalign.viewer.cali_format import CaliWriter
    from catalign.viewer.metrics_tiler import MetricsTiler, AlignmentPosition

    print(f"\nAligning {chrom_name}...")
    print(f"  Query:  {len(query_seq):,} bp")
    print(f"  Target: {len(target_seq):,} bp")

    # Configure aligner for large sequences
    energy_model = EnergyModel(
        match_energy=-2.0,
        mismatch_energy=3.0,
        gap_open_energy=5.0,
        gap_extend_energy=1.0,
    )

    # Use larger k-mer for chromosome-scale
    k = 15 if len(query_seq) > 100000 else 11
    w = 100 if len(query_seq) > 100000 else 50

    aligner = CatalignAligner(energy_model=energy_model, k=k, w=w)

    # Run alignment with timing
    start_time = time.perf_counter()
    alignment = aligner.align(query_seq, target_seq, "query", "target")
    align_time = time.perf_counter() - start_time

    print(f"  Alignment completed in {align_time:.2f}s")

    # Evaluate quality
    quality = evaluate_quality(alignment, query_seq, target_seq, energy_model)
    metrics = compute_metrics(alignment, query_seq, target_seq, quality)

    print(f"  Identity: {metrics.identity*100:.2f}%")
    print(f"  Quality:  {metrics.quality_score:.1f}/100")

    # Generate multi-scale tiles
    print(f"\nGenerating multi-scale metrics...")
    tile_sizes = [1000, 10000, 100000]
    if len(target_seq) > 1_000_000:
        tile_sizes.append(1_000_000)

    tiler = MetricsTiler(len(target_seq), tile_sizes)

    # Process alignment
    for qpos, tpos, op in alignment.aligned_pairs:
        if tpos is not None:
            pos = AlignmentPosition(
                query_pos=qpos,
                target_pos=tpos,
                operation=op,
                query_base=query_seq[qpos] if qpos is not None and qpos < len(query_seq) else "",
                target_base=target_seq[tpos] if tpos < len(target_seq) else "",
            )
            tiler.add_position(tpos, pos)

    # Write CALI file
    cali_path = output_dir / f"{chrom_name}.cali"

    with CaliWriter(cali_path, reference_name="T2T-CHM13", sample_name="HG002-sim", tile_sizes=tile_sizes) as writer:
        writer.add_chromosome(chrom_name, len(target_seq))

        for ts in tile_sizes:
            for tile in tiler.get_tiles(ts):
                writer.add_tile(chrom_name, ts, tile)

    print(f"  Saved: {cali_path}")

    return {
        "chrom": chrom_name,
        "query_length": len(query_seq),
        "target_length": len(target_seq),
        "alignment_time": align_time,
        "identity": metrics.identity,
        "coverage": metrics.query_coverage,
        "quality_score": metrics.quality_score,
        "matches": metrics.matches,
        "mismatches": metrics.mismatches,
        "insertions": metrics.insertions,
        "deletions": metrics.deletions,
        "cali_file": str(cali_path),
    }


def generate_visualizations(
    cali_path: Path,
    output_dir: Path,
    chrom_name: str,
) -> None:
    """Generate interactive visualizations from CALI file."""
    from catalign.viewer.cali_format import CaliFile

    print(f"\nGenerating visualizations for {chrom_name}...")

    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        print("  Plotly not installed, skipping visualizations")
        return

    # Read CALI file
    cali = CaliFile(cali_path)

    # Get all tile sizes
    tile_sizes = cali.header.tile_sizes

    for tile_size in tile_sizes[:3]:  # First 3 resolutions
        data = cali.to_numpy(chrom_name, tile_size)

        if len(data["start"]) == 0:
            continue

        # Create multi-track figure
        fig = make_subplots(
            rows=4, cols=1,
            shared_xaxes=True,
            vertical_spacing=0.05,
            subplot_titles=("Identity", "Coverage", "Quality Score", "Gap Rate"),
            row_heights=[0.25, 0.25, 0.25, 0.25],
        )

        x = (data["start"] + data["end"]) / 2  # Midpoints

        # Identity track
        fig.add_trace(go.Scatter(
            x=x, y=data["identity"] * 100,
            mode="lines",
            fill="tozeroy",
            fillcolor="rgba(100, 200, 100, 0.5)",
            line=dict(color="green", width=1),
            name="Identity",
        ), row=1, col=1)

        # Coverage track
        fig.add_trace(go.Scatter(
            x=x, y=data["coverage"] * 100,
            mode="lines",
            fill="tozeroy",
            fillcolor="rgba(100, 100, 200, 0.5)",
            line=dict(color="blue", width=1),
            name="Coverage",
        ), row=2, col=1)

        # Quality track
        fig.add_trace(go.Scatter(
            x=x, y=data["quality_score"],
            mode="lines",
            fill="tozeroy",
            fillcolor="rgba(150, 100, 200, 0.5)",
            line=dict(color="purple", width=1),
            name="Quality",
        ), row=3, col=1)

        # Gap rate track
        fig.add_trace(go.Scatter(
            x=x, y=data["gap_rate"] * 100,
            mode="lines",
            fill="tozeroy",
            fillcolor="rgba(200, 150, 100, 0.5)",
            line=dict(color="orange", width=1),
            name="Gap Rate",
        ), row=4, col=1)

        # Update layout
        fig.update_layout(
            title=f"{chrom_name} Alignment Metrics ({tile_size//1000}kb resolution)",
            height=800,
            showlegend=False,
        )

        fig.update_yaxes(title_text="%", range=[0, 105], row=1, col=1)
        fig.update_yaxes(title_text="%", range=[0, 105], row=2, col=1)
        fig.update_yaxes(title_text="Score", range=[0, 105], row=3, col=1)
        fig.update_yaxes(title_text="%", range=[0, 50], row=4, col=1)
        fig.update_xaxes(title_text="Position (bp)", row=4, col=1)

        # Save
        res_name = f"{tile_size//1000}kb" if tile_size >= 1000 else f"{tile_size}bp"
        html_path = output_dir / f"{chrom_name}_{res_name}_metrics.html"
        fig.write_html(str(html_path))
        print(f"  Saved: {html_path.name}")

    cali.close()


def main():
    parser = argparse.ArgumentParser(description="HG002 vs T2T Genome-Wide Demo")
    parser.add_argument("--chrom", default="chrM", choices=list(TEST_CHROMOSOMES.keys()),
                        help="Chromosome to analyze")
    parser.add_argument("--region-size", type=int, default=None,
                        help="Region size to analyze (default: full chromosome)")
    parser.add_argument("--output-dir", type=Path, default=Path("hg002_t2t_demo"),
                        help="Output directory")
    parser.add_argument("--mutation-rate", type=float, default=0.001,
                        help="SNP mutation rate for simulated assembly")

    args = parser.parse_args()

    print("="*60)
    print("🧬 HG002 vs T2T-CHM13 Genome-Wide Alignment Demo")
    print("="*60)

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = output_dir / "data"
    data_dir.mkdir(exist_ok=True)

    chrom = args.chrom
    chrom_length = TEST_CHROMOSOMES[chrom]

    # Determine region
    if args.region_size:
        region_end = min(args.region_size, chrom_length)
    else:
        region_end = chrom_length

    print(f"\nTarget: {chrom}:1-{region_end:,} ({region_end:,} bp)")

    # Download or generate reference
    ref_path = data_dir / f"t2t_{chrom}.fa"

    if ref_path.exists():
        print(f"\nUsing cached reference: {ref_path}")
        with open(ref_path) as f:
            lines = f.readlines()
            ref_seq = "".join(l.strip() for l in lines if not l.startswith(">"))
    else:
        print(f"\nDownloading T2T-CHM13 {chrom}...")
        try:
            ref_seq = download_chromosome_sequence("hs1", chrom, 1, region_end, ref_path)
        except Exception as e:
            print(f"Download failed: {e}")
            print("Generating simulated reference instead...")
            import random
            random.seed(42)
            ref_seq = "".join(random.choices("ACGT", k=region_end))
            with open(ref_path, "w") as f:
                f.write(f">{chrom}:1-{region_end}\n")
                for i in range(0, len(ref_seq), 80):
                    f.write(ref_seq[i:i+80] + "\n")

    # Generate simulated HG002 assembly
    print(f"\nGenerating simulated HG002 assembly...")
    print(f"  Mutation rate: {args.mutation_rate}")

    query_seq = generate_simulated_assembly(
        ref_seq,
        mutation_rate=args.mutation_rate,
        indel_rate=args.mutation_rate / 10,
        sv_rate=args.mutation_rate / 100,
    )

    query_path = data_dir / f"hg002_sim_{chrom}.fa"
    with open(query_path, "w") as f:
        f.write(f">HG002_simulated_{chrom}\n")
        for i in range(0, len(query_seq), 80):
            f.write(query_seq[i:i+80] + "\n")

    print(f"  Generated {len(query_seq):,} bp assembly")

    # Run alignment
    results = run_chromosome_alignment(
        query_seq,
        ref_seq,
        chrom,
        output_dir,
    )

    # Generate visualizations
    viz_dir = output_dir / "visualizations"
    viz_dir.mkdir(exist_ok=True)

    generate_visualizations(
        Path(results["cali_file"]),
        viz_dir,
        chrom,
    )

    # Save summary
    summary_path = output_dir / "summary.json"
    with open(summary_path, "w") as f:
        json.dump(results, f, indent=2)

    # Print summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"Chromosome:     {chrom}")
    print(f"Query length:   {results['query_length']:,} bp")
    print(f"Target length:  {results['target_length']:,} bp")
    print(f"Alignment time: {results['alignment_time']:.2f}s")
    print(f"Identity:       {results['identity']*100:.2f}%")
    print(f"Quality score:  {results['quality_score']:.1f}/100")
    print(f"Matches:        {results['matches']:,}")
    print(f"Mismatches:     {results['mismatches']:,}")
    print(f"Insertions:     {results['insertions']:,}")
    print(f"Deletions:      {results['deletions']:,}")

    print("\n" + "="*60)
    print("OUTPUT FILES")
    print("="*60)
    print(f"CALI file:     {results['cali_file']}")
    print(f"Summary:       {summary_path}")
    print(f"Visualizations: {viz_dir}/")

    print("\n" + "="*60)
    print("NEXT STEPS")
    print("="*60)
    print("1. View CALI file with caliview:")
    print(f"   caliview view {results['cali_file']}")
    print()
    print("2. Open visualizations in browser:")
    print(f"   open {viz_dir}/")
    print()
    print("3. Run on larger chromosomes:")
    print(f"   python examples/hg002_t2t_demo.py --chrom chr22")


if __name__ == "__main__":
    main()
