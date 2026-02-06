"""CLI entry point for Catalign."""

from __future__ import annotations

import argparse
import sys

from catalign.align import CatalignAligner
from catalign.energy import EnergyModel
from catalign.io import read_fasta
from catalign.quality import evaluate_quality


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="catalign",
        description="Catalign – energy-based DNA sequence alignment",
    )
    sub = parser.add_subparsers(dest="command")

    # align sub-command
    align_p = sub.add_parser("align", help="Align query against target FASTA")
    align_p.add_argument("query", help="Query FASTA file")
    align_p.add_argument("target", help="Target FASTA file")
    align_p.add_argument("--kmer-size", type=int, default=15)
    align_p.add_argument("--window-size", type=int, default=50)
    align_p.add_argument("--bandwidth", type=int, default=100)
    align_p.add_argument("--match-energy", type=float, default=-2.0)
    align_p.add_argument("--mismatch-energy", type=float, default=3.0)
    align_p.add_argument("--gap-open", type=float, default=5.0)
    align_p.add_argument("--gap-extend", type=float, default=1.0)
    align_p.add_argument("--output", choices=["text", "paf", "json"], default="text")

    # quality sub-command
    qual_p = sub.add_parser("quality", help="Evaluate alignment quality")
    qual_p.add_argument("query", help="Query FASTA file")
    qual_p.add_argument("target", help="Target FASTA file")
    qual_p.add_argument("--kmer-size", type=int, default=15)
    qual_p.add_argument("--window-size", type=int, default=50)

    # viz sub-command
    viz_p = sub.add_parser("viz", help="Launch interactive visualization dashboard")
    viz_p.add_argument("--port", type=int, default=8501, help="Port for dashboard server")
    viz_p.add_argument("--host", type=str, default="localhost", help="Host address")

    # benchmark sub-command
    bench_p = sub.add_parser("benchmark", help="Run benchmark suite against test data")
    bench_p.add_argument("--resources-dir", type=str, default="tests/resources",
                         help="Path to test resources directory")
    bench_p.add_argument("--output", type=str, help="Output file for benchmark report")

    # metrics sub-command
    metrics_p = sub.add_parser("metrics", help="Generate metrics report")
    metrics_p.add_argument("query", help="Query FASTA file")
    metrics_p.add_argument("target", help="Target FASTA file")
    metrics_p.add_argument("--output", type=str, help="Output file for metrics report")
    metrics_p.add_argument("--json", action="store_true", help="Output as JSON")

    return parser


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        sys.exit(0)

    if args.command == "align":
        _cmd_align(args)
    elif args.command == "quality":
        _cmd_quality(args)
    elif args.command == "viz":
        _cmd_viz(args)
    elif args.command == "benchmark":
        _cmd_benchmark(args)
    elif args.command == "metrics":
        _cmd_metrics(args)


def _cmd_align(args) -> None:
    em = EnergyModel(
        match_energy=args.match_energy,
        mismatch_energy=args.mismatch_energy,
        gap_open_energy=args.gap_open,
        gap_extend_energy=args.gap_extend,
    )
    aligner = CatalignAligner(energy_model=em, k=args.kmer_size, w=args.window_size)

    queries = list(read_fasta(args.query))
    targets = list(read_fasta(args.target))

    for qname, qseq in queries:
        for tname, tseq in targets:
            aln = aligner.align(qseq, tseq, query_name=qname, target_name=tname)
            if args.output == "paf":
                _print_paf(aln, len(qseq), len(tseq))
            elif args.output == "json":
                _print_json(aln, qseq, tseq)
            else:
                _print_text(aln)


def _cmd_quality(args) -> None:
    aligner = CatalignAligner(k=args.kmer_size, w=args.window_size)
    queries = list(read_fasta(args.query))
    targets = list(read_fasta(args.target))

    for qname, qseq in queries:
        for tname, tseq in targets:
            aln = aligner.align(qseq, tseq, query_name=qname, target_name=tname)
            qual = evaluate_quality(aln, qseq, tseq)
            print(f"Query: {qname}  Target: {tname}")
            print(f"  Identity:  {qual.overall_identity:.4f}")
            print(f"  Matches:   {qual.total_matches}")
            print(f"  Mismatches:{qual.total_mismatches}")
            print(f"  Insertions:{qual.total_insertions}")
            print(f"  Deletions: {qual.total_deletions}")
            print(f"  Energy:    {qual.total_energy:.2f}")
            print(f"  Q-score:   {qual.quality_score:.1f}")
            if qual.region_quality:
                rq = qual.region_quality
                print(f"  Blocks:    {rq.num_blocks}")
                print(f"  Coverage Q:{rq.coverage_query:.4f}")
                print(f"  Coverage T:{rq.coverage_target:.4f}")
            print()


def _print_text(aln) -> None:
    print(
        f"{aln.query_name}\t{aln.query_start}-{aln.query_end}\t"
        f"{aln.target_name}\t{aln.target_start}-{aln.target_end}\t"
        f"strand={aln.strand}\tenergy={aln.energy_score:.2f}\t"
        f"cigar={aln.cigar}"
    )


def _print_paf(aln, qlen: int, tlen: int) -> None:
    n_matches = sum(1 for _, _, op in aln.aligned_pairs if op == "M")
    block_len = len(aln.aligned_pairs)
    print(
        f"{aln.query_name}\t{qlen}\t{aln.query_start}\t{aln.query_end}\t"
        f"{aln.strand}\t"
        f"{aln.target_name}\t{tlen}\t{aln.target_start}\t{aln.target_end}\t"
        f"{n_matches}\t{block_len}\t255\tcg:Z:{aln.cigar}"
    )


def _print_json(aln, qseq: str, tseq: str) -> None:
    import json
    from catalign.metrics import compute_metrics

    qual = evaluate_quality(aln, qseq, tseq)
    metrics = compute_metrics(aln, qseq, tseq, qual)

    data = {
        "query_name": aln.query_name,
        "target_name": aln.target_name,
        "query_start": aln.query_start,
        "query_end": aln.query_end,
        "target_start": aln.target_start,
        "target_end": aln.target_end,
        "strand": aln.strand,
        "cigar": aln.cigar,
        "energy_score": aln.energy_score,
        "metrics": metrics.to_dict(),
    }
    print(json.dumps(data, indent=2))


def _cmd_viz(args) -> None:
    try:
        from catalign.viz.dashboard import launch_dashboard
        print(f"Launching Catalign dashboard at http://{args.host}:{args.port}")
        launch_dashboard(port=args.port, host=args.host)
    except ImportError as e:
        print(f"Error: {e}")
        print("Install visualization dependencies with: pip install catalign[viz]")
        sys.exit(1)


def _cmd_benchmark(args) -> None:
    import json
    from pathlib import Path
    from catalign.metrics import run_benchmark_suite

    resources_dir = Path(args.resources_dir)
    manifest_path = resources_dir / "ground_truth" / "manifest.json"

    if not manifest_path.exists():
        print(f"Error: Ground truth manifest not found at {manifest_path}")
        print("Generate test data first with: python scripts/generate_test_data.py")
        sys.exit(1)

    with open(manifest_path) as f:
        manifest = json.load(f)

    test_cases = manifest.get("test_cases", [])
    print(f"Running {len(test_cases)} benchmark tests...")

    aligner = CatalignAligner()
    suite = run_benchmark_suite(test_cases, aligner, resources_dir)

    report = suite.summary()
    print(report)

    if args.output:
        Path(args.output).write_text(report)
        print(f"\nReport saved to: {args.output}")

    sys.exit(0 if suite.failed == 0 else 1)


def _cmd_metrics(args) -> None:
    from pathlib import Path
    from catalign.metrics import compute_metrics
    import json

    aligner = CatalignAligner()
    queries = list(read_fasta(args.query))
    targets = list(read_fasta(args.target))

    all_metrics = []

    for qname, qseq in queries:
        for tname, tseq in targets:
            aln = aligner.align(qseq, tseq, query_name=qname, target_name=tname)
            qual = evaluate_quality(aln, qseq, tseq)
            metrics = compute_metrics(aln, qseq, tseq, qual)
            all_metrics.append({
                "query_name": qname,
                "target_name": tname,
                "metrics": metrics.to_dict(),
            })

    if args.json:
        output = json.dumps(all_metrics, indent=2)
    else:
        lines = []
        for item in all_metrics:
            m = item["metrics"]
            lines.extend([
                f"Query: {item['query_name']}  Target: {item['target_name']}",
                f"  Identity:      {m['identity']:.4f} ({m['identity']*100:.2f}%)",
                f"  Matches:       {m['matches']}",
                f"  Mismatches:    {m['mismatches']}",
                f"  Insertions:    {m['insertions']}",
                f"  Deletions:     {m['deletions']}",
                f"  Gap rate:      {m['gap_rate']:.4f}",
                f"  Query cov:     {m['query_coverage']:.4f}",
                f"  Target cov:    {m['target_coverage']:.4f}",
                f"  Energy:        {m['total_energy']:.2f}",
                f"  Quality score: {m['quality_score']:.1f}",
                "",
            ])
        output = "\n".join(lines)

    print(output)

    if args.output:
        Path(args.output).write_text(output)
        print(f"\nReport saved to: {args.output}")

