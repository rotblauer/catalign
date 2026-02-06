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
    align_p.add_argument("--output", choices=["text", "paf"], default="text")

    # quality sub-command
    qual_p = sub.add_parser("quality", help="Evaluate alignment quality")
    qual_p.add_argument("query", help="Query FASTA file")
    qual_p.add_argument("target", help="Target FASTA file")
    qual_p.add_argument("--kmer-size", type=int, default=15)
    qual_p.add_argument("--window-size", type=int, default=50)

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
