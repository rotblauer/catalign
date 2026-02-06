#!/usr/bin/env python3
"""Generate synthetic test data for Catalign alignment validation.

All data is generated with fixed random seeds for reproducibility.
Run this script to regenerate the test resources directory.

Usage:
    python scripts/generate_test_data.py [--output-dir tests/resources]
"""

from __future__ import annotations

import argparse
import json
import random
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List, Tuple, Optional


BASES = "ACGT"
DEFAULT_OUTPUT_DIR = Path(__file__).parent.parent / "tests" / "resources"


@dataclass
class TestCase:
    """Describes a test case for ground truth validation."""

    name: str
    query_file: str
    target_file: str
    description: str
    expected_cigar: Optional[str] = None
    expected_identity: Optional[float] = None
    expected_matches: Optional[int] = None
    expected_mismatches: Optional[int] = None
    expected_insertions: Optional[int] = None
    expected_deletions: Optional[int] = None
    mutation_positions: Optional[List[int]] = None
    indel_positions: Optional[List[int]] = None
    sv_type: Optional[str] = None
    sv_size: Optional[int] = None


def random_sequence(length: int, seed: int) -> str:
    """Generate a random DNA sequence with fixed seed."""
    random.seed(seed)
    return "".join(random.choice(BASES) for _ in range(length))


def mutate_base(base: str) -> str:
    """Return a different base (transition/transversion)."""
    others = [b for b in BASES if b != base]
    return random.choice(others)


def write_fasta(filepath: Path, name: str, sequence: str, line_width: int = 80) -> None:
    """Write a single sequence to a FASTA file."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(sequence), line_width):
            f.write(sequence[i:i + line_width] + "\n")


def write_ground_truth(filepath: Path, test_case: TestCase) -> None:
    """Write ground truth JSON for a test case."""
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w") as f:
        data = {k: v for k, v in asdict(test_case).items() if v is not None}
        json.dump(data, f, indent=2)


# =============================================================================
# Test Case Generators
# =============================================================================

def generate_identical_cases(output_dir: Path) -> List[TestCase]:
    """Generate identical sequence pairs of various lengths."""
    cases = []
    lengths = [20, 100, 500, 1000, 5000]

    for length in lengths:
        seq = random_sequence(length, seed=42 + length)
        name = f"identical_{length}bp"

        query_path = output_dir / "synthetic" / "identical" / f"{name}_query.fa"
        target_path = output_dir / "synthetic" / "identical" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", seq)
        write_fasta(target_path, f"{name}_target", seq)

        case = TestCase(
            name=name,
            query_file=f"synthetic/identical/{name}_query.fa",
            target_file=f"synthetic/identical/{name}_target.fa",
            description=f"Identical {length}bp sequences",
            expected_cigar=f"{length}M",
            expected_identity=1.0,
            expected_matches=length,
            expected_mismatches=0,
            expected_insertions=0,
            expected_deletions=0,
        )
        cases.append(case)

    return cases


def generate_snp_cases(output_dir: Path) -> List[TestCase]:
    """Generate sequences with single nucleotide polymorphisms."""
    cases = []

    # Single SNP at various positions
    for pos in [0, 50, 99]:
        seq_len = 100
        seed = 100 + pos
        random.seed(seed)

        target = random_sequence(seq_len, seed)
        query = list(target)
        query[pos] = mutate_base(target[pos])
        query = "".join(query)

        name = f"snp_pos{pos}"
        query_path = output_dir / "synthetic" / "mutations" / f"{name}_query.fa"
        target_path = output_dir / "synthetic" / "mutations" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", query)
        write_fasta(target_path, f"{name}_target", target)

        case = TestCase(
            name=name,
            query_file=f"synthetic/mutations/{name}_query.fa",
            target_file=f"synthetic/mutations/{name}_target.fa",
            description=f"Single SNP at position {pos}",
            expected_identity=0.99,
            expected_matches=seq_len - 1,
            expected_mismatches=1,
            expected_insertions=0,
            expected_deletions=0,
            mutation_positions=[pos],
        )
        cases.append(case)

    # Multiple SNPs
    for num_snps in [5, 10, 20]:
        seq_len = 200
        seed = 200 + num_snps
        random.seed(seed)

        target = random_sequence(seq_len, seed)
        query = list(target)

        # Pick random positions for SNPs
        positions = sorted(random.sample(range(seq_len), num_snps))
        for pos in positions:
            query[pos] = mutate_base(target[pos])
        query = "".join(query)

        name = f"snp_multi_{num_snps}"
        query_path = output_dir / "synthetic" / "mutations" / f"{name}_query.fa"
        target_path = output_dir / "synthetic" / "mutations" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", query)
        write_fasta(target_path, f"{name}_target", target)

        case = TestCase(
            name=name,
            query_file=f"synthetic/mutations/{name}_query.fa",
            target_file=f"synthetic/mutations/{name}_target.fa",
            description=f"{num_snps} SNPs spread across {seq_len}bp",
            expected_identity=(seq_len - num_snps) / seq_len,
            expected_matches=seq_len - num_snps,
            expected_mismatches=num_snps,
            expected_insertions=0,
            expected_deletions=0,
            mutation_positions=positions,
        )
        cases.append(case)

    return cases


def generate_indel_cases(output_dir: Path) -> List[TestCase]:
    """Generate sequences with insertions and deletions."""
    cases = []

    # Single insertion of various sizes
    for indel_size in [1, 3, 10, 50]:
        seq_len = 200
        seed = 300 + indel_size
        random.seed(seed)

        target = random_sequence(seq_len, seed)
        insert_pos = seq_len // 2
        insert_seq = random_sequence(indel_size, seed + 1000)

        # Query has insertion relative to target
        query = target[:insert_pos] + insert_seq + target[insert_pos:]

        name = f"insertion_{indel_size}bp"
        query_path = output_dir / "synthetic" / "indels" / f"{name}_query.fa"
        target_path = output_dir / "synthetic" / "indels" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", query)
        write_fasta(target_path, f"{name}_target", target)

        case = TestCase(
            name=name,
            query_file=f"synthetic/indels/{name}_query.fa",
            target_file=f"synthetic/indels/{name}_target.fa",
            description=f"{indel_size}bp insertion at position {insert_pos}",
            expected_matches=seq_len,
            expected_insertions=indel_size,
            expected_deletions=0,
            indel_positions=[insert_pos],
        )
        cases.append(case)

    # Single deletion of various sizes
    for indel_size in [1, 3, 10, 50]:
        seq_len = 200
        seed = 400 + indel_size
        random.seed(seed)

        target = random_sequence(seq_len, seed)
        del_pos = seq_len // 2

        # Query has deletion relative to target
        query = target[:del_pos] + target[del_pos + indel_size:]

        name = f"deletion_{indel_size}bp"
        query_path = output_dir / "synthetic" / "indels" / f"{name}_query.fa"
        target_path = output_dir / "synthetic" / "indels" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", query)
        write_fasta(target_path, f"{name}_target", target)

        case = TestCase(
            name=name,
            query_file=f"synthetic/indels/{name}_query.fa",
            target_file=f"synthetic/indels/{name}_target.fa",
            description=f"{indel_size}bp deletion at position {del_pos}",
            expected_matches=seq_len - indel_size,
            expected_insertions=0,
            expected_deletions=indel_size,
            indel_positions=[del_pos],
        )
        cases.append(case)

    # Complex: multiple indels
    for num_indels in [3, 5]:
        seq_len = 500
        seed = 500 + num_indels
        random.seed(seed)

        target = random_sequence(seq_len, seed)
        query = target

        # Apply indels at random positions (spacing them out)
        indel_positions = []
        offset = 0
        spacing = seq_len // (num_indels + 1)

        for i in range(num_indels):
            pos = spacing * (i + 1) + offset
            if i % 2 == 0:
                # Insertion
                insert_seq = random_sequence(random.randint(1, 5), seed + i * 100)
                query = query[:pos] + insert_seq + query[pos:]
                offset += len(insert_seq)
            else:
                # Deletion
                del_size = random.randint(1, 5)
                query = query[:pos] + query[pos + del_size:]
                offset -= del_size
            indel_positions.append(pos - offset + (offset if i % 2 == 0 else 0))

        name = f"complex_indels_{num_indels}"
        query_path = output_dir / "synthetic" / "indels" / f"{name}_query.fa"
        target_path = output_dir / "synthetic" / "indels" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", query)
        write_fasta(target_path, f"{name}_target", target)

        case = TestCase(
            name=name,
            query_file=f"synthetic/indels/{name}_query.fa",
            target_file=f"synthetic/indels/{name}_target.fa",
            description=f"Complex case with {num_indels} indels",
            indel_positions=indel_positions,
        )
        cases.append(case)

    return cases


def generate_structural_variant_cases(output_dir: Path) -> List[TestCase]:
    """Generate sequences with structural variants (inversions, duplications, translocations)."""
    cases = []

    # Inversion
    for inv_size in [20, 100, 500]:
        seq_len = max(1000, inv_size * 3)
        seed = 600 + inv_size
        random.seed(seed)

        target = random_sequence(seq_len, seed)
        inv_start = seq_len // 3
        inv_end = inv_start + inv_size

        # Invert the middle region
        inverted = target[inv_start:inv_end][::-1]
        # Complement
        comp = {"A": "T", "T": "A", "C": "G", "G": "C"}
        inverted = "".join(comp.get(b, b) for b in inverted)

        query = target[:inv_start] + inverted + target[inv_end:]

        name = f"inversion_{inv_size}bp"
        query_path = output_dir / "synthetic" / "structural" / f"{name}_query.fa"
        target_path = output_dir / "synthetic" / "structural" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", query)
        write_fasta(target_path, f"{name}_target", target)

        case = TestCase(
            name=name,
            query_file=f"synthetic/structural/{name}_query.fa",
            target_file=f"synthetic/structural/{name}_target.fa",
            description=f"{inv_size}bp inversion at position {inv_start}",
            sv_type="inversion",
            sv_size=inv_size,
        )
        cases.append(case)

    # Tandem duplication
    for dup_size in [10, 50, 200]:
        seq_len = max(500, dup_size * 5)
        seed = 700 + dup_size
        random.seed(seed)

        target = random_sequence(seq_len, seed)
        dup_start = seq_len // 3
        dup_end = dup_start + dup_size

        # Duplicate the region
        duplicated_region = target[dup_start:dup_end]
        query = target[:dup_end] + duplicated_region + target[dup_end:]

        name = f"tandem_dup_{dup_size}bp"
        query_path = output_dir / "synthetic" / "structural" / f"{name}_query.fa"
        target_path = output_dir / "synthetic" / "structural" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", query)
        write_fasta(target_path, f"{name}_target", target)

        case = TestCase(
            name=name,
            query_file=f"synthetic/structural/{name}_query.fa",
            target_file=f"synthetic/structural/{name}_target.fa",
            description=f"{dup_size}bp tandem duplication at position {dup_start}",
            sv_type="duplication",
            sv_size=dup_size,
        )
        cases.append(case)

    # Large deletion (SV-scale)
    for del_size in [100, 500, 1000]:
        seq_len = max(2000, del_size * 3)
        seed = 800 + del_size
        random.seed(seed)

        target = random_sequence(seq_len, seed)
        del_start = seq_len // 3

        query = target[:del_start] + target[del_start + del_size:]

        name = f"sv_deletion_{del_size}bp"
        query_path = output_dir / "synthetic" / "structural" / f"{name}_query.fa"
        target_path = output_dir / "synthetic" / "structural" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", query)
        write_fasta(target_path, f"{name}_target", target)

        case = TestCase(
            name=name,
            query_file=f"synthetic/structural/{name}_query.fa",
            target_file=f"synthetic/structural/{name}_target.fa",
            description=f"{del_size}bp structural deletion at position {del_start}",
            sv_type="deletion",
            sv_size=del_size,
        )
        cases.append(case)

    return cases


def generate_repetitive_cases(output_dir: Path) -> List[TestCase]:
    """Generate sequences with repetitive elements."""
    cases = []

    # Simple tandem repeats
    for unit in ["AT", "ATG", "AATG", "AAATGG"]:
        repeats = 50
        target = unit * repeats
        query = unit * repeats  # Identical

        name = f"str_{unit}x{repeats}"
        query_path = output_dir / "synthetic" / "repetitive" / f"{name}_query.fa"
        target_path = output_dir / "synthetic" / "repetitive" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", query)
        write_fasta(target_path, f"{name}_target", target)

        case = TestCase(
            name=name,
            query_file=f"synthetic/repetitive/{name}_query.fa",
            target_file=f"synthetic/repetitive/{name}_target.fa",
            description=f"Simple tandem repeat: {unit} x {repeats}",
            expected_identity=1.0,
        )
        cases.append(case)

    # Repeat expansion/contraction
    for unit in ["CAG", "CTG"]:
        target_repeats = 30
        query_repeats = 50  # Expansion

        seed = 900 + ord(unit[0])
        random.seed(seed)

        flank_5 = random_sequence(100, seed)
        flank_3 = random_sequence(100, seed + 1)

        target = flank_5 + (unit * target_repeats) + flank_3
        query = flank_5 + (unit * query_repeats) + flank_3

        name = f"repeat_expansion_{unit}_{target_repeats}_to_{query_repeats}"
        query_path = output_dir / "synthetic" / "repetitive" / f"{name}_query.fa"
        target_path = output_dir / "synthetic" / "repetitive" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", query)
        write_fasta(target_path, f"{name}_target", target)

        case = TestCase(
            name=name,
            query_file=f"synthetic/repetitive/{name}_query.fa",
            target_file=f"synthetic/repetitive/{name}_target.fa",
            description=f"Repeat expansion: {unit} from {target_repeats} to {query_repeats}",
            sv_type="repeat_expansion",
            sv_size=(query_repeats - target_repeats) * len(unit),
        )
        cases.append(case)

    # Interspersed repeats (Alu-like)
    seed = 950
    random.seed(seed)
    alu_like = random_sequence(300, seed)  # Simulated Alu element

    seq_len = 2000
    target = random_sequence(seq_len, seed + 1)

    # Insert Alu-like elements at multiple positions
    positions = [300, 800, 1400]
    query = target
    offset = 0
    for pos in positions:
        insert_pos = pos + offset
        query = query[:insert_pos] + alu_like + query[insert_pos:]
        offset += len(alu_like)

    name = "interspersed_repeats"
    query_path = output_dir / "synthetic" / "repetitive" / f"{name}_query.fa"
    target_path = output_dir / "synthetic" / "repetitive" / f"{name}_target.fa"

    write_fasta(query_path, f"{name}_query", query)
    write_fasta(target_path, f"{name}_target", target)

    case = TestCase(
        name=name,
        query_file=f"synthetic/repetitive/{name}_query.fa",
        target_file=f"synthetic/repetitive/{name}_target.fa",
        description=f"Interspersed repeats (Alu-like) at positions {positions}",
    )
    cases.append(case)

    return cases


def generate_benchmark_cases(output_dir: Path) -> List[TestCase]:
    """Generate larger sequences for benchmarking."""
    cases = []

    for length in [10_000, 50_000, 100_000]:
        seed = 1000 + length
        random.seed(seed)

        target = random_sequence(length, seed)

        # Add ~1% mutations
        num_mutations = length // 100
        positions = sorted(random.sample(range(length), num_mutations))
        query = list(target)
        for pos in positions:
            query[pos] = mutate_base(target[pos])
        query = "".join(query)

        name = f"benchmark_{length // 1000}kb"
        query_path = output_dir / "benchmarks" / f"{name}_query.fa"
        target_path = output_dir / "benchmarks" / f"{name}_target.fa"

        write_fasta(query_path, f"{name}_query", query)
        write_fasta(target_path, f"{name}_target", target)

        case = TestCase(
            name=name,
            query_file=f"benchmarks/{name}_query.fa",
            target_file=f"benchmarks/{name}_target.fa",
            description=f"Benchmark {length // 1000}kb with ~1% divergence",
            expected_identity=0.99,
        )
        cases.append(case)

    return cases


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic test data for Catalign")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Output directory for test resources"
    )
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    print(f"Generating test data in: {output_dir}")

    # Collect all test cases
    all_cases: List[TestCase] = []

    print("Generating identical sequence cases...")
    all_cases.extend(generate_identical_cases(output_dir))

    print("Generating SNP cases...")
    all_cases.extend(generate_snp_cases(output_dir))

    print("Generating indel cases...")
    all_cases.extend(generate_indel_cases(output_dir))

    print("Generating structural variant cases...")
    all_cases.extend(generate_structural_variant_cases(output_dir))

    print("Generating repetitive sequence cases...")
    all_cases.extend(generate_repetitive_cases(output_dir))

    print("Generating benchmark cases...")
    all_cases.extend(generate_benchmark_cases(output_dir))

    # Write ground truth manifest
    ground_truth_dir = output_dir / "ground_truth"
    ground_truth_dir.mkdir(parents=True, exist_ok=True)

    manifest = {
        "version": "1.0",
        "generated_with_seed": True,
        "test_cases": [asdict(c) for c in all_cases],
    }

    manifest_path = ground_truth_dir / "manifest.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    # Write individual ground truth files
    for case in all_cases:
        gt_path = ground_truth_dir / f"{case.name}.json"
        write_ground_truth(gt_path, case)

    print(f"\nGenerated {len(all_cases)} test cases")
    print(f"Ground truth manifest: {manifest_path}")


if __name__ == "__main__":
    main()
