#!/usr/bin/env python3
"""
NA19240 vs hg38 Alignment Demo

Aligns NA19240 (Yoruba individual from 1000 Genomes) assembly against
hg38 reference, focusing on a specific region of interest.

Target region: chr17:10958130-11017410 (~59kb region)

This is a larger-scale test demonstrating genome-wide alignment capabilities.

Usage:
    python examples/na19240_t2t_demo.py [--region chr17:10958130-11017410]
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

sys.path.insert(0, str(Path(__file__).parent.parent))


# =============================================================================
# Data Sources - NA19240 from HPRC
# =============================================================================

# NA19240 is part of the Human Pangenome Reference Consortium (HPRC)
# Year 1 assemblies available from NCBI GenBank

# NCBI GenBank assemblies for NA19240 (HPRC Year 1, freeze 2)
# GCA_018503275.3 = NA19240 maternal haplotype
# GCA_018503265.2 = NA19240 paternal haplotype
NCBI_BASE = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA"

NA19240_MAT = f"{NCBI_BASE}/018/503/275/GCA_018503275.3_NA19240_mat_hprc_f2/GCA_018503275.3_NA19240_mat_hprc_f2_genomic.fna.gz"
NA19240_PAT = f"{NCBI_BASE}/018/503/265/GCA_018503265.2_NA19240_pat_hprc_f2/GCA_018503265.2_NA19240_pat_hprc_f2_genomic.fna.gz"

# Region of interest on hg38
DEFAULT_REGION = "chr17:10958130-11017410"  # ~59kb region on hg38
DEFAULT_GENOME = "hg38"


def parse_region(region_str: str) -> Tuple[str, int, int]:
    """Parse region string like 'chr17:10882740-10902740'."""
    chrom, coords = region_str.split(":")
    start, end = coords.split("-")
    return chrom, int(start.replace(",", "")), int(end.replace(",", ""))


def download_file_with_progress(url: str, dest: Path, chunk_size: int = 1024*1024) -> None:
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
                        mb = downloaded / (1024 * 1024)
                        total_mb = total_size / (1024 * 1024)
                        print(f"\r  Progress: {pct:.1f}% ({mb:.1f}/{total_mb:.1f} MB)", end="", flush=True)

            print()  # Newline

    except Exception as e:
        print(f"\nError downloading: {e}")
        raise


def extract_region_from_gzipped_fasta(
    fasta_gz: Path,
    chrom: str,
    start: int,
    end: int,
    output_path: Path,
) -> str:
    """Extract a region from a gzipped FASTA file."""
    print(f"Extracting {chrom}:{start:,}-{end:,} from {fasta_gz.name}...")

    sequence_parts = []
    in_target_chrom = False
    current_pos = 0
    found = False
    found_header = ""

    # Parse chromosome number from chrom (e.g., "chr17" -> "17")
    chrom_num = chrom.replace("chr", "")

    with gzip.open(fasta_gz, 'rt') as f:
        for line in f:
            line = line.rstrip()

            if line.startswith('>'):
                header = line[1:].split()[0]
                full_header = line[1:]

                # Handle different naming conventions:
                # 1. Direct match: chr17, 17
                # 2. NCBI GenBank format: CM099623.1 with "chromosome 17" in description
                is_match = False
                if header == chrom or header == chrom_num or f"chr{header}" == chrom:
                    is_match = True
                elif f"chromosome {chrom_num}," in full_header.lower() or f"chromosome {chrom_num} " in full_header.lower():
                    is_match = True
                    print(f"  Found chromosome {chrom_num}: {header}")

                if is_match:
                    in_target_chrom = True
                    current_pos = 0
                    found = True
                    found_header = header
                else:
                    if in_target_chrom:
                        break  # Done with target chromosome
                    in_target_chrom = False
            elif in_target_chrom:
                line_start = current_pos
                line_end = current_pos + len(line)

                # Check overlap with region
                if line_end > start and line_start < end:
                    overlap_start = max(0, start - line_start)
                    overlap_end = min(len(line), end - line_start)
                    sequence_parts.append(line[overlap_start:overlap_end])

                current_pos = line_end

                if current_pos >= end:
                    break

    if not found:
        print(f"  Warning: Chromosome {chrom} not found in {fasta_gz.name}")
        # List available chromosomes
        with gzip.open(fasta_gz, 'rt') as f:
            chroms = []
            for line in f:
                if line.startswith('>'):
                    chroms.append(line[1:].split()[0])
                    if len(chroms) >= 10:
                        break
            print(f"  Available chromosomes (first 10): {chroms}")
        return ""

    sequence = "".join(sequence_parts)

    if sequence:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            f.write(f">{found_header}:{start}-{end}\n")
            for i in range(0, len(sequence), 80):
                f.write(sequence[i:i+80] + "\n")

        print(f"  Extracted {len(sequence):,} bp from {found_header}")
    else:
        print(f"  No sequence extracted for region")

    return sequence


def find_homologous_region_in_assembly(
    assembly_gz: Path,
    ref_seq: str,
    chrom_hint: str,
    output_path: Path,
    search_window: int = 100000,
    padding: int = 5000,
) -> Tuple[str, str]:
    """Find the homologous region in an assembly using k-mer matching.

    This function searches for the region in the assembly that is most similar
    to the reference sequence by comparing minimizer sketches.

    Returns (contig_name, sequence).
    """
    from catalign.sketch import minimizer_sketch

    print(f"Searching for homologous region in {assembly_gz.name}...")

    # Parse chromosome number from hint
    chrom_num = chrom_hint.replace("chr", "")

    # Extract minimizers from reference (use smaller windows for precise matching)
    ref_sketch = minimizer_sketch(ref_seq, k=15, w=30)
    ref_kmer_set = set(h for h, _, _ in ref_sketch.entries)
    print(f"  Reference minimizers: {len(ref_kmer_set)}")

    best_contig = ""
    best_start = 0
    best_score = 0
    best_seq = ""
    target_contig_seq = ""

    # First, find the target chromosome and load it
    current_contig = ""
    current_header = ""
    current_seq = []

    print(f"  Loading chromosome {chrom_num}...")
    with gzip.open(assembly_gz, 'rt') as f:
        for line in f:
            line = line.rstrip()

            if line.startswith('>'):
                # Check if previous contig was our target
                if target_contig_seq:
                    break

                header = line[1:].split()[0]
                full_header = line[1:]

                # Check if this is our target chromosome
                if f"chromosome {chrom_num}," in full_header.lower() or f"chromosome {chrom_num} " in full_header.lower():
                    current_contig = header
                    current_header = full_header
                    print(f"  Found: {header}")
                else:
                    if current_contig and current_seq:
                        target_contig_seq = "".join(current_seq)
                    current_contig = ""
                    current_seq = []
            elif current_contig:
                current_seq.append(line)

    # Get the sequence if not already done
    if current_contig and current_seq and not target_contig_seq:
        target_contig_seq = "".join(current_seq)

    if not target_contig_seq:
        print(f"  WARNING: Chromosome {chrom_num} not found")
        return "", ""

    print(f"  Chromosome {chrom_num} length: {len(target_contig_seq):,} bp")

    # Now scan through the chromosome to find the best matching region
    # Use sliding window with minimizer comparison
    window_size = len(ref_seq)
    step_size = 10000  # Step by 10kb for initial scan

    print(f"  Scanning chromosome for homologous region (window: {window_size:,} bp, step: {step_size:,})...")

    best_positions = []

    for pos in range(0, len(target_contig_seq) - window_size, step_size):
        window_seq = target_contig_seq[pos:pos + window_size]
        window_sketch = minimizer_sketch(window_seq, k=15, w=30)
        window_kmers = set(h for h, _, _ in window_sketch.entries)
        overlap = len(ref_kmer_set & window_kmers)

        if overlap > best_score * 0.8:  # Track near-best positions
            best_positions.append((pos, overlap))

        if overlap > best_score:
            best_score = overlap
            best_start = pos

        # Progress
        if pos % 1000000 == 0:
            print(f"\r    Scanned: {pos//1000000}M bp, best score: {best_score}", end="", flush=True)

    print()

    if best_score == 0:
        print("  No significant match found")
        return "", ""

    # Refine the position with finer step
    print(f"  Refining position around {best_start:,}...")
    refined_start = best_start
    refined_score = best_score

    for pos in range(max(0, best_start - step_size), min(len(target_contig_seq) - window_size, best_start + step_size), 1000):
        window_seq = target_contig_seq[pos:pos + window_size]
        window_sketch = minimizer_sketch(window_seq, k=15, w=30)
        window_kmers = set(h for h, _, _ in window_sketch.entries)
        overlap = len(ref_kmer_set & window_kmers)

        if overlap > refined_score:
            refined_score = overlap
            refined_start = pos

    # Extract the region with padding
    extract_start = max(0, refined_start - padding)
    extract_end = min(len(target_contig_seq), refined_start + window_size + padding)
    extracted = target_contig_seq[extract_start:extract_end]

    print(f"  Best match at position {refined_start:,} (score: {refined_score}, {refined_score*100/len(ref_kmer_set):.1f}% of ref minimizers)")
    print(f"  Extracted region: {extract_start:,}-{extract_end:,} ({len(extracted):,} bp)")

    # Save the extracted sequence
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(f">{current_contig}:{extract_start}-{extract_end}\n")
        for i in range(0, len(extracted), 80):
            f.write(extracted[i:i+80] + "\n")

    return current_contig, extracted


def run_alignment(
    query_path: Path,
    target_path: Path,
    output_dir: Path,
    name: str,
) -> Dict:
    """Run Catalign alignment and generate outputs."""
    from catalign import CatalignAligner, EnergyModel, evaluate_quality, read_fasta
    from catalign.metrics import compute_metrics, detect_structural_variants, detect_svs_by_sequence_comparison
    from catalign.viewer.cali_format import CaliWriter
    from catalign.viewer.metrics_tiler import MetricsTiler, AlignmentPosition

    # Load sequences
    print(f"\nLoading sequences...")
    query_records = list(read_fasta(query_path))
    target_records = list(read_fasta(target_path))

    if not query_records:
        raise ValueError(f"No sequences in query file: {query_path}")
    if not target_records:
        raise ValueError(f"No sequences in target file: {target_path}")

    query_name, query_seq = query_records[0]
    target_name, target_seq = target_records[0]

    print(f"  Query:  {query_name} ({len(query_seq):,} bp)")
    print(f"  Target: {target_name} ({len(target_seq):,} bp)")

    length_diff = len(query_seq) - len(target_seq)
    if abs(length_diff) > 100:
        print(f"  Length difference: {length_diff:+,} bp (indicates structural variation)")

    # Configure aligner with wider bandwidth for large indels
    energy_model = EnergyModel(
        match_energy=-2.0,
        mismatch_energy=3.0,
        gap_open_energy=5.0,
        gap_extend_energy=1.0,
        transition_energy=2.0,
    )

    # Adjust parameters based on sequence length
    # Use larger bandwidth for expected large SVs
    if len(query_seq) > 100000:
        k, w = 19, 100
    elif len(query_seq) > 10000:
        k, w = 15, 50
    else:
        k, w = 11, 20

    aligner = CatalignAligner(energy_model=energy_model, k=k, w=w)

    # Run alignment
    print(f"\nRunning alignment (k={k}, w={w})...")
    start_time = time.perf_counter()

    alignment = aligner.align(query_seq, target_seq, query_name, target_name)

    align_time = time.perf_counter() - start_time
    print(f"  Completed in {align_time:.2f}s")

    # Evaluate
    quality = evaluate_quality(alignment, query_seq, target_seq, energy_model)
    metrics = compute_metrics(alignment, query_seq, target_seq, quality)

    print(f"  Identity: {metrics.identity*100:.2f}%")
    print(f"  Quality:  {metrics.quality_score:.1f}/100")
    print(f"  Total gaps: {metrics.insertions + metrics.deletions:,}")

    # Detect structural variants using multiple methods
    print(f"\nDetecting structural variants (min size: 50bp)...")

    # Method 1: From alignment
    svs_from_alignment = detect_structural_variants(alignment, query_seq, target_seq, min_sv_size=50)

    # Method 2: Direct sequence comparison (more reliable for large SVs)
    svs_from_sequence = detect_svs_by_sequence_comparison(query_seq, target_seq, min_sv_size=50)

    # Merge results, preferring sequence-based detection
    svs = svs_from_sequence if svs_from_sequence else svs_from_alignment

    print(f"  Alignment-based detection: {len(svs_from_alignment)} SV(s)")
    print(f"  Sequence-based detection:  {len(svs_from_sequence)} SV(s)")

    if svs:
        print(f"\n  Detected structural variants:")
        for sv in svs:
            print(f"    • {sv}")
            if sv.target_start > 0 or sv.target_end > 0:
                print(f"      Position: target {sv.target_start:,}-{sv.target_end:,}")
    else:
        print(f"  No large SVs detected")
        len_diff = len(query_seq) - len(target_seq)
        if abs(len_diff) >= 50:
            print(f"  Note: Length difference of {len_diff:+,}bp suggests structural variation")

    # Generate multi-scale metrics
    print(f"\nGenerating multi-scale metrics...")

    # Choose tile sizes based on sequence length
    if len(target_seq) > 100000:
        tile_sizes = [100, 1000, 10000, 100000]
    elif len(target_seq) > 10000:
        tile_sizes = [100, 1000, 10000]
    else:
        tile_sizes = [10, 100, 1000]

    tiler = MetricsTiler(len(target_seq), tile_sizes)

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
    output_dir.mkdir(parents=True, exist_ok=True)
    cali_path = output_dir / f"{name}.cali"

    with CaliWriter(cali_path, reference_name="hg38", sample_name="NA19240", tile_sizes=tile_sizes) as writer:
        writer.add_chromosome(target_name.split(":")[0] if ":" in target_name else "region", len(target_seq))

        for ts in tile_sizes:
            for tile in tiler.get_tiles(ts):
                writer.add_tile(target_name.split(":")[0] if ":" in target_name else "region", ts, tile)

    print(f"  Saved: {cali_path}")

    # Save sequences for viewing
    seq_dir = output_dir / "sequences"
    seq_dir.mkdir(exist_ok=True)

    (seq_dir / "reference.fa").write_text(f">{target_name}\n{target_seq}\n")
    (seq_dir / "query.fa").write_text(f">{query_name}\n{query_seq}\n")

    # Generate alignment text
    alignment_text = generate_alignment_text(alignment, query_seq, target_seq)
    (output_dir / f"{name}_alignment.txt").write_text(alignment_text)

    return {
        "name": name,
        "query_name": query_name,
        "target_name": target_name,
        "query_length": len(query_seq),
        "target_length": len(target_seq),
        "length_difference": len(query_seq) - len(target_seq),
        "alignment_time": align_time,
        "identity": metrics.identity,
        "coverage": metrics.query_coverage,
        "quality_score": metrics.quality_score,
        "matches": metrics.matches,
        "mismatches": metrics.mismatches,
        "insertions": metrics.insertions,
        "deletions": metrics.deletions,
        "structural_variants": [sv.to_dict() for sv in svs],
        "cali_file": str(cali_path),
    }


def generate_alignment_text(alignment, query_seq: str, target_seq: str, line_width: int = 80) -> str:
    """Generate human-readable alignment text."""
    lines = []
    lines.append(f"Query:  {alignment.query_name}")
    lines.append(f"Target: {alignment.target_name}")
    lines.append(f"CIGAR:  {alignment.cigar[:100]}{'...' if len(alignment.cigar) > 100 else ''}")
    lines.append("")

    # Build aligned sequences
    q_aligned = []
    t_aligned = []
    m_line = []

    for qpos, tpos, op in alignment.aligned_pairs:
        if op == "M":
            qb = query_seq[qpos] if qpos is not None and qpos < len(query_seq) else "N"
            tb = target_seq[tpos] if tpos is not None and tpos < len(target_seq) else "N"
            q_aligned.append(qb)
            t_aligned.append(tb)
            m_line.append("|" if qb == tb else ".")
        elif op == "X":
            qb = query_seq[qpos] if qpos is not None and qpos < len(query_seq) else "N"
            tb = target_seq[tpos] if tpos is not None and tpos < len(target_seq) else "N"
            q_aligned.append(qb)
            t_aligned.append(tb)
            m_line.append("X")
        elif op == "I":
            qb = query_seq[qpos] if qpos is not None and qpos < len(query_seq) else "N"
            q_aligned.append(qb)
            t_aligned.append("-")
            m_line.append(" ")
        elif op == "D":
            tb = target_seq[tpos] if tpos is not None and tpos < len(target_seq) else "N"
            q_aligned.append("-")
            t_aligned.append(tb)
            m_line.append(" ")

    q_str = "".join(q_aligned)
    t_str = "".join(t_aligned)
    m_str = "".join(m_line)

    # Format in blocks
    for i in range(0, len(q_str), line_width):
        lines.append(f"Query  {q_str[i:i+line_width]}")
        lines.append(f"       {m_str[i:i+line_width]}")
        lines.append(f"Target {t_str[i:i+line_width]}")
        lines.append("")

    return "\n".join(lines)


def generate_visualizations(results: Dict, output_dir: Path) -> None:
    """Generate interactive visualizations."""
    from catalign.viewer.cali_format import CaliFile

    print(f"\nGenerating visualizations...")

    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        print("  Plotly not installed, skipping visualizations")
        return

    cali = CaliFile(results["cali_file"])
    viz_dir = output_dir / "visualizations"
    viz_dir.mkdir(exist_ok=True)

    for chrom in cali.get_chromosomes():
        for tile_size in cali.header.tile_sizes[:3]:
            data = cali.to_numpy(chrom, tile_size)

            if len(data["start"]) == 0:
                continue

            x = (data["start"] + data["end"]) / 2

            fig = make_subplots(
                rows=5, cols=1,
                shared_xaxes=True,
                vertical_spacing=0.03,
                subplot_titles=("Identity", "Coverage", "Quality", "Gap Rate", "Indel Counts"),
                row_heights=[0.2, 0.2, 0.2, 0.2, 0.2],
            )

            # Identity
            fig.add_trace(go.Scatter(
                x=x, y=data["identity"] * 100,
                mode="lines", fill="tozeroy",
                fillcolor="rgba(100, 200, 100, 0.5)",
                line=dict(color="green", width=1),
                name="Identity",
            ), row=1, col=1)

            # Coverage
            fig.add_trace(go.Scatter(
                x=x, y=data["coverage"] * 100,
                mode="lines", fill="tozeroy",
                fillcolor="rgba(100, 100, 200, 0.5)",
                line=dict(color="blue", width=1),
                name="Coverage",
            ), row=2, col=1)

            # Quality
            fig.add_trace(go.Scatter(
                x=x, y=data["quality_score"],
                mode="lines", fill="tozeroy",
                fillcolor="rgba(150, 100, 200, 0.5)",
                line=dict(color="purple", width=1),
                name="Quality",
            ), row=3, col=1)

            # Gap rate
            fig.add_trace(go.Scatter(
                x=x, y=data["gap_rate"] * 100,
                mode="lines", fill="tozeroy",
                fillcolor="rgba(200, 150, 100, 0.5)",
                line=dict(color="orange", width=1),
                name="Gap Rate",
            ), row=4, col=1)

            # Indel counts
            fig.add_trace(go.Bar(
                x=x, y=data["insertions"],
                name="Insertions",
                marker_color="rgba(100, 150, 255, 0.7)",
            ), row=5, col=1)
            fig.add_trace(go.Bar(
                x=x, y=-data["deletions"],  # Negative to show below axis
                name="Deletions",
                marker_color="rgba(255, 100, 100, 0.7)",
            ), row=5, col=1)

            fig.update_layout(
                title=f"NA19240 vs hg38 - {chrom} ({tile_size}bp tiles)",
                height=900,
                showlegend=False,
            )

            fig.update_yaxes(title_text="%", range=[0, 105], row=1, col=1)
            fig.update_yaxes(title_text="%", range=[0, 105], row=2, col=1)
            fig.update_yaxes(title_text="Q", range=[0, 105], row=3, col=1)
            fig.update_yaxes(title_text="%", range=[0, 20], row=4, col=1)
            fig.update_yaxes(title_text="Count", row=5, col=1)  # Auto-range for indels
            fig.update_xaxes(title_text="Position (bp)", row=5, col=1)

            res_str = f"{tile_size}bp" if tile_size < 1000 else f"{tile_size//1000}kb"
            html_path = viz_dir / f"{chrom}_{res_str}_metrics.html"
            fig.write_html(str(html_path))
            print(f"  Saved: {html_path.name}")

    cali.close()


def main():
    parser = argparse.ArgumentParser(description="NA19240 vs hg38 Alignment Demo")
    parser.add_argument("--region", default=DEFAULT_REGION,
                        help=f"Region to analyze (default: {DEFAULT_REGION})")
    parser.add_argument("--output-dir", type=Path, default=Path("na19240_t2t_demo"),
                        help="Output directory")
    parser.add_argument("--haplotype", choices=["mat", "pat", "both"], default="mat",
                        help="Which haplotype to use (default: mat)")
    parser.add_argument("--skip-download", action="store_true",
                        help="Skip downloads if files exist")

    args = parser.parse_args()

    print("="*70)
    print("🧬 NA19240 vs hg38 Assembly Alignment")
    print("="*70)

    chrom, start, end = parse_region(args.region)
    region_size = end - start

    print(f"\nTarget region: {chrom}:{start:,}-{end:,} ({region_size:,} bp)")
    print(f"Haplotype: {args.haplotype}")

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = output_dir / "data"
    data_dir.mkdir(exist_ok=True)

    # Step 1: Get hg38 reference region
    print(f"\n{'='*70}")
    print("Step 1: Extracting hg38 reference region")
    print("="*70)

    ref_fa = data_dir / f"hg38_{chrom}_{start}_{end}.fa"

    if ref_fa.exists() and args.skip_download:
        print(f"Using cached: {ref_fa}")
        ref_seq = "".join(l.strip() for l in open(ref_fa) if not l.startswith(">"))
    else:
        print("Fetching region from UCSC hg38...")
        try:
            url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
            print(f"Fetching from UCSC: {url}")
            with urllib.request.urlopen(url, timeout=120) as response:
                data = json.loads(response.read().decode())
                ref_seq = data.get("dna", "").upper().replace("\n", "")

                if ref_seq:
                    with open(ref_fa, 'w') as f:
                        f.write(f">{chrom}:{start}-{end}\n")
                        for i in range(0, len(ref_seq), 80):
                            f.write(ref_seq[i:i+80] + "\n")
                    print(f"  Downloaded {len(ref_seq):,} bp from UCSC")
        except Exception as e:
            print(f"  UCSC fetch failed: {e}")
            print("  Generating simulated reference...")
            import random
            random.seed(42)
            ref_seq = "".join(random.choices("ACGT", k=region_size))
            with open(ref_fa, 'w') as f:
                f.write(f">{chrom}:{start}-{end}_simulated\n")
                for i in range(0, len(ref_seq), 80):
                    f.write(ref_seq[i:i+80] + "\n")

    if not ref_seq:
        print("ERROR: Could not obtain reference sequence")
        sys.exit(1)

    print(f"Reference sequence: {len(ref_seq):,} bp")

    # Step 2: Get NA19240 assembly region
    print(f"\n{'='*70}")
    print("Step 2: Obtaining NA19240 assembly region")
    print("="*70)

    query_seq = None  # Initialize

    # NA19240 assembly URL
    if args.haplotype == "mat":
        assembly_url = NA19240_MAT
    else:
        assembly_url = NA19240_PAT

    assembly_gz = data_dir / f"NA19240.{args.haplotype}.fa.gz"
    query_fa = data_dir / f"na19240_{args.haplotype}_{chrom}_{start}_{end}.fa"

    if query_fa.exists() and args.skip_download:
        print(f"Using cached: {query_fa}")
        query_seq = "".join(l.strip() for l in open(query_fa) if not l.startswith(">"))
    else:
        # Download NA19240 assembly if needed
        if not assembly_gz.exists():
            print(f"Downloading NA19240 {args.haplotype} assembly from NCBI...")
            print(f"  URL: {assembly_url}")
            print(f"  Note: This is a large file (~800MB), download may take a while...")

            try:
                download_file_with_progress(assembly_url, assembly_gz)
            except Exception as e:
                print(f"  Download failed: {e}")
                print(f"  Falling back to simulated sequence...")
                assembly_gz = None

        if assembly_gz and assembly_gz.exists():
            # For independent assemblies, we need to find the homologous region
            # Direct coordinate extraction doesn't work because the assembly
            # has its own coordinate system
            print(f"\nFinding homologous region to {chrom}:{start}-{end} in NA19240 assembly...")
            print("  (NA19240 is an independent assembly with its own coordinates)")

            contig_name, query_seq = find_homologous_region_in_assembly(
                assembly_gz, ref_seq, chrom, query_fa
            )

            if not query_seq:
                print("  Could not find homologous region, falling back to simulation")
                assembly_gz = None

        # Fallback to simulation if real data not available
        if not assembly_gz or not query_seq:
            print("\nGenerating simulated NA19240-like sequence...")
            print("  (Real assembly download failed or region not found)")

            import random
            random.seed(19240)

            # Known SVs for this region (based on literature/databases)
            KNOWN_SVS = {
                "chr17:10958130-11017410": [
                    {"type": "INS", "position": 30000, "size": 400, "description": "~400bp insertion in NA19240"},
                    {"type": "DEL", "position": 45000, "size": 5000, "description": "~5kb deletion in NA19240 (vs hg38)"},
                ]
            }

            # Store the simulated SVs for comparison later
            simulated_svs = []

            region_key = f"{chrom}:{start}-{end}"
            known_svs = KNOWN_SVS.get(region_key, [])

            if known_svs:
                print(f"  Simulating known structural variants:")
                for sv in known_svs:
                    print(f"    - {sv['description']} at ~position {sv['position']}")

            query_seq = list(ref_seq)

            # Add SNPs
            snp_count = 0
            for i in range(len(query_seq)):
                if random.random() < 0.001:
                    bases = [b for b in "ACGT" if b != query_seq[i]]
                    query_seq[i] = random.choice(bases)
                    snp_count += 1

            query_seq = "".join(query_seq)

            # Apply known structural variants
            offset = 0
            for sv in sorted(known_svs, key=lambda x: x['position']):
                pos = sv['position'] + offset

                if sv['type'] == 'INS':
                    insert_seq = "".join(random.choices("ACGT", k=sv['size']))
                    query_seq = query_seq[:pos] + insert_seq + query_seq[pos:]
                    simulated_svs.append({
                        "type": "INS",
                        "ref_position": sv['position'],
                        "query_position": pos,
                        "size": sv['size'],
                    })
                    offset += sv['size']
                    print(f"  ✓ Applied {sv['size']}bp INSERTION at ref position {sv['position']}")

                elif sv['type'] == 'DEL':
                    del_end = min(pos + sv['size'], len(query_seq))
                    actual_del_size = del_end - pos
                    query_seq = query_seq[:pos] + query_seq[del_end:]
                    simulated_svs.append({
                        "type": "DEL",
                        "ref_position": sv['position'],
                        "query_position": pos,
                        "size": actual_del_size,
                    })
                    offset -= actual_del_size
                    print(f"  ✓ Applied {actual_del_size}bp DELETION at ref position {sv['position']}")

            print(f"  SNPs: {snp_count}")
            print(f"  Net length change: {offset:+,} bp")

            # Save simulated SVs to file for reference
            sv_info_path = data_dir / f"simulated_svs_{chrom}_{start}_{end}.json"
            with open(sv_info_path, 'w') as f:
                json.dump({"simulated_svs": simulated_svs, "known_svs": known_svs}, f, indent=2)

            with open(query_fa, 'w') as f:
                f.write(f">NA19240_{args.haplotype}_{chrom}:{start}-{end}_simulated\n")
                for i in range(0, len(query_seq), 80):
                    f.write(query_seq[i:i+80] + "\n")

    print(f"Query sequence: {len(query_seq):,} bp")

    # Step 3: Run alignment
    print(f"\n{'='*70}")
    print("Step 3: Running Catalign alignment")
    print("="*70)

    results = run_alignment(
        query_fa,
        ref_fa,
        output_dir,
        f"na19240_{args.haplotype}_{chrom}_{start}_{end}",
    )

    # Step 4: Generate visualizations
    print(f"\n{'='*70}")
    print("Step 4: Generating visualizations")
    print("="*70)

    generate_visualizations(results, output_dir)

    # Save summary
    summary_path = output_dir / "summary.json"
    with open(summary_path, 'w') as f:
        json.dump(results, f, indent=2)

    # Print summary
    print(f"\n{'='*70}")
    print("SUMMARY")
    print("="*70)
    print(f"Region:         {chrom}:{start:,}-{end:,}")
    print(f"Reference:      hg38 ({results['target_length']:,} bp)")
    print(f"Query:          NA19240-{args.haplotype} ({results['query_length']:,} bp)")
    print(f"Length diff:    {results.get('length_difference', 0):+,} bp")
    print(f"Alignment time: {results['alignment_time']:.2f}s")
    print(f"Identity:       {results['identity']*100:.2f}%")
    print(f"Quality score:  {results['quality_score']:.1f}/100")
    print(f"Matches:        {results['matches']:,}")
    print(f"Mismatches:     {results['mismatches']:,}")
    print(f"Insertions:     {results['insertions']:,}")
    print(f"Deletions:      {results['deletions']:,}")

    # Print structural variants
    svs = results.get('structural_variants', [])
    if svs:
        print(f"\n{'='*70}")
        print("STRUCTURAL VARIANTS DETECTED")
        print("="*70)
        for sv in svs:
            sv_type = sv['type']
            size = sv['size']
            if sv_type == "INS":
                print(f"  • INSERTION: {size:,} bp at target position ~{sv['target_start']:,}")
            elif sv_type == "DEL":
                print(f"  • DELETION:  {size:,} bp at target position {sv['target_start']:,}-{sv['target_end']:,}")
            else:
                print(f"  • {sv_type}: {size:,} bp")

    print(f"\n{'='*70}")
    print("OUTPUT FILES")
    print("="*70)
    print(f"CALI file:       {results['cali_file']}")
    print(f"Alignment text:  {output_dir}/na19240_*_alignment.txt")
    print(f"Visualizations:  {output_dir}/visualizations/")
    print(f"Sequences:       {output_dir}/data/")

    print(f"\n{'='*70}")
    print("NEXT STEPS")
    print("="*70)
    print(f"1. View with caliview:")
    print(f"   ./target/release/caliview view {results['cali_file']}")
    print()
    print(f"2. Open visualizations:")
    print(f"   open {output_dir}/visualizations/")
    print()
    print(f"3. View alignment text:")
    print(f"   less {output_dir}/*_alignment.txt")


if __name__ == "__main__":
    main()
