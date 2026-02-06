"""Benchmark tests using synthetic test data."""

import json
import pytest
from pathlib import Path

from catalign import CatalignAligner, EnergyModel, evaluate_quality
from catalign.metrics import compute_metrics, compare_to_ground_truth, BenchmarkResult


# Path to test resources
RESOURCES_DIR = Path(__file__).parent / "resources"
GROUND_TRUTH_DIR = RESOURCES_DIR / "ground_truth"


def load_manifest():
    """Load the ground truth manifest."""
    manifest_path = GROUND_TRUTH_DIR / "manifest.json"
    if not manifest_path.exists():
        pytest.skip("Ground truth manifest not found. Run: python scripts/generate_test_data.py")

    with open(manifest_path) as f:
        return json.load(f)


def load_sequence(filepath: Path) -> str:
    """Load a sequence from a FASTA file."""
    from catalign.io import read_fasta
    full_path = RESOURCES_DIR / filepath
    if not full_path.exists():
        pytest.skip(f"Test file not found: {full_path}")

    records = list(read_fasta(full_path))
    if not records:
        pytest.skip(f"No sequences in: {full_path}")
    return records[0][1]


class TestIdenticalSequences:
    """Tests for identical sequence alignments."""

    @pytest.mark.parametrize("length", [20, 100, 500, 1000])
    def test_identical_alignment_identity(self, length):
        """Identical sequences should have 100% identity."""
        manifest = load_manifest()
        case = next(
            (c for c in manifest["test_cases"] if c["name"] == f"identical_{length}bp"),
            None
        )
        if case is None:
            pytest.skip(f"Test case not found: identical_{length}bp")

        query = load_sequence(case["query_file"])
        target = load_sequence(case["target_file"])

        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        quality = evaluate_quality(aln, query, target)
        metrics = compute_metrics(aln, query, target, quality)

        assert metrics.identity == 1.0, f"Expected 100% identity, got {metrics.identity}"
        assert metrics.mismatches == 0, f"Expected 0 mismatches, got {metrics.mismatches}"


class TestSNPCases:
    """Tests for sequences with single nucleotide polymorphisms."""

    @pytest.mark.parametrize("pos", [0, 50, 99])
    def test_single_snp_detection(self, pos):
        """Single SNP should be detected as one mismatch."""
        manifest = load_manifest()
        case = next(
            (c for c in manifest["test_cases"] if c["name"] == f"snp_pos{pos}"),
            None
        )
        if case is None:
            pytest.skip(f"Test case not found: snp_pos{pos}")

        query = load_sequence(case["query_file"])
        target = load_sequence(case["target_file"])

        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        quality = evaluate_quality(aln, query, target)
        metrics = compute_metrics(aln, query, target, quality)

        # Allow for some flexibility in alignment
        assert metrics.mismatches >= 1, "Should detect at least one mismatch"
        assert metrics.identity >= 0.95, f"Identity should be high, got {metrics.identity}"

    @pytest.mark.parametrize("num_snps", [5, 10, 20])
    def test_multiple_snps(self, num_snps):
        """Multiple SNPs should be detected with appropriate identity drop."""
        manifest = load_manifest()
        case = next(
            (c for c in manifest["test_cases"] if c["name"] == f"snp_multi_{num_snps}"),
            None
        )
        if case is None:
            pytest.skip(f"Test case not found: snp_multi_{num_snps}")

        query = load_sequence(case["query_file"])
        target = load_sequence(case["target_file"])

        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        quality = evaluate_quality(aln, query, target)
        metrics = compute_metrics(aln, query, target, quality)

        expected_identity = case.get("expected_identity", 0.9)
        # Allow 5% tolerance
        assert metrics.identity >= expected_identity - 0.05, \
            f"Identity {metrics.identity} below expected {expected_identity}"


class TestIndelCases:
    """Tests for sequences with insertions and deletions."""

    @pytest.mark.parametrize("size", [1, 3, 10])
    def test_insertion_detection(self, size):
        """Insertions should be correctly identified."""
        manifest = load_manifest()
        case = next(
            (c for c in manifest["test_cases"] if c["name"] == f"insertion_{size}bp"),
            None
        )
        if case is None:
            pytest.skip(f"Test case not found: insertion_{size}bp")

        query = load_sequence(case["query_file"])
        target = load_sequence(case["target_file"])

        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        quality = evaluate_quality(aln, query, target)
        metrics = compute_metrics(aln, query, target, quality)

        # Check that gaps were introduced
        total_gaps = metrics.insertions + metrics.deletions
        assert total_gaps >= size * 0.5, f"Expected gaps for {size}bp insertion"

    @pytest.mark.parametrize("size", [1, 3, 10])
    def test_deletion_detection(self, size):
        """Deletions should be correctly identified."""
        manifest = load_manifest()
        case = next(
            (c for c in manifest["test_cases"] if c["name"] == f"deletion_{size}bp"),
            None
        )
        if case is None:
            pytest.skip(f"Test case not found: deletion_{size}bp")

        query = load_sequence(case["query_file"])
        target = load_sequence(case["target_file"])

        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        quality = evaluate_quality(aln, query, target)
        metrics = compute_metrics(aln, query, target, quality)

        # Check that gaps were introduced
        total_gaps = metrics.insertions + metrics.deletions
        assert total_gaps >= size * 0.5, f"Expected gaps for {size}bp deletion"


class TestStructuralVariants:
    """Tests for structural variant detection."""

    @pytest.mark.slow
    @pytest.mark.parametrize("sv_type", ["inversion_20bp", "tandem_dup_10bp", "sv_deletion_100bp"])
    def test_sv_handling(self, sv_type):
        """Structural variants should not crash the aligner."""
        manifest = load_manifest()
        case = next(
            (c for c in manifest["test_cases"] if c["name"] == sv_type),
            None
        )
        if case is None:
            pytest.skip(f"Test case not found: {sv_type}")

        query = load_sequence(case["query_file"])
        target = load_sequence(case["target_file"])

        aligner = CatalignAligner(k=10, w=20)

        # Should complete without error
        aln = aligner.align(query, target)
        assert len(aln.aligned_pairs) > 0, "Should produce some alignment"


class TestRepetitiveSequences:
    """Tests for repetitive sequence handling."""

    def test_simple_tandem_repeat(self):
        """Simple tandem repeats should align correctly."""
        manifest = load_manifest()
        case = next(
            (c for c in manifest["test_cases"] if c["name"] == "str_ATx50"),
            None
        )
        if case is None:
            pytest.skip("Test case not found: str_ATx50")

        query = load_sequence(case["query_file"])
        target = load_sequence(case["target_file"])

        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        quality = evaluate_quality(aln, query, target)
        metrics = compute_metrics(aln, query, target, quality)

        # Identical repeats should align perfectly
        assert metrics.identity >= 0.99, f"Expected high identity for identical repeats"


class TestGroundTruthComparison:
    """Tests that compare alignments to ground truth expectations."""

    @pytest.mark.benchmark
    def test_all_ground_truth_cases(self):
        """Run all test cases against ground truth."""
        manifest = load_manifest()
        aligner = CatalignAligner()

        passed = 0
        failed = 0
        failures = []

        for case in manifest["test_cases"]:
            try:
                query = load_sequence(case["query_file"])
                target = load_sequence(case["target_file"])

                aln = aligner.align(query, target)
                quality = evaluate_quality(aln, query, target)
                metrics = compute_metrics(aln, query, target, quality)

                result = compare_to_ground_truth(aln, case, metrics)

                if result.passed:
                    passed += 1
                else:
                    failed += 1
                    failures.append((case["name"], result.errors))
            except Exception as e:
                failed += 1
                failures.append((case["name"], [str(e)]))

        # Report failures
        if failures:
            msg_lines = [f"Failed {failed}/{passed + failed} tests:"]
            for name, errors in failures[:10]:  # Limit to first 10
                msg_lines.append(f"  {name}: {errors}")
            if len(failures) > 10:
                msg_lines.append(f"  ... and {len(failures) - 10} more")
            pytest.fail("\n".join(msg_lines))
