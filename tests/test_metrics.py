"""Tests for metrics module."""

import pytest
from catalign import CatalignAligner, evaluate_quality
from catalign.metrics import (
    compute_metrics,
    compare_to_ground_truth,
    AlignmentMetrics,
    BenchmarkResult,
    BenchmarkSuite,
)


class TestAlignmentMetrics:
    """Tests for AlignmentMetrics class."""

    def test_metrics_from_identical_alignment(self):
        """Compute metrics for identical sequences."""
        query = "ACGTACGTACGT"
        target = "ACGTACGTACGT"

        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        quality = evaluate_quality(aln, query, target)

        metrics = compute_metrics(aln, query, target, quality)

        assert metrics.identity == 1.0
        assert metrics.mismatches == 0
        assert metrics.gap_rate == 0.0

    def test_metrics_from_mismatched_alignment(self):
        """Compute metrics for sequences with mismatches."""
        query = "ACGTACGTACGT"
        target = "ACGTAAGTACGT"  # One mismatch

        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        quality = evaluate_quality(aln, query, target)

        metrics = compute_metrics(aln, query, target, quality)

        assert metrics.identity < 1.0
        assert metrics.mismatches >= 1

    def test_metrics_to_dict(self):
        """Test metrics serialization to dict."""
        metrics = AlignmentMetrics(
            alignment_length=100,
            matches=90,
            mismatches=10,
            identity=0.9,
        )

        d = metrics.to_dict()

        assert d["alignment_length"] == 100
        assert d["matches"] == 90
        assert d["identity"] == 0.9

    def test_metrics_to_json(self, tmp_path):
        """Test metrics export to JSON."""
        metrics = AlignmentMetrics(
            alignment_length=100,
            matches=90,
            mismatches=10,
        )

        output_path = tmp_path / "metrics.json"
        json_str = metrics.to_json(str(output_path))

        assert output_path.exists()
        assert '"alignment_length": 100' in json_str

    def test_metrics_from_dict(self):
        """Test metrics creation from dict."""
        data = {
            "alignment_length": 100,
            "matches": 90,
            "mismatches": 10,
            "identity": 0.9,
        }

        metrics = AlignmentMetrics.from_dict(data)

        assert metrics.alignment_length == 100
        assert metrics.identity == 0.9


class TestBenchmarkResult:
    """Tests for BenchmarkResult class."""

    def test_passed_result(self):
        """Test a passing benchmark result."""
        result = BenchmarkResult(
            name="test_case",
            passed=True,
            errors=[],
        )

        assert result.passed
        assert len(result.errors) == 0

    def test_failed_result(self):
        """Test a failing benchmark result."""
        result = BenchmarkResult(
            name="test_case",
            passed=False,
            errors=["Identity too low"],
        )

        assert not result.passed
        assert len(result.errors) == 1


class TestBenchmarkSuite:
    """Tests for BenchmarkSuite class."""

    def test_empty_suite(self):
        """Test empty benchmark suite."""
        suite = BenchmarkSuite()

        assert suite.total == 0
        assert suite.passed == 0
        assert suite.failed == 0

    def test_suite_with_results(self):
        """Test suite with mixed results."""
        suite = BenchmarkSuite(results=[
            BenchmarkResult(name="test1", passed=True),
            BenchmarkResult(name="test2", passed=True),
            BenchmarkResult(name="test3", passed=False, errors=["error"]),
        ])

        assert suite.total == 3
        assert suite.passed == 2
        assert suite.failed == 1
        assert suite.pass_rate == 2/3

    def test_suite_summary(self):
        """Test suite summary generation."""
        suite = BenchmarkSuite(results=[
            BenchmarkResult(name="test1", passed=True),
            BenchmarkResult(name="test2", passed=False, errors=["error"]),
        ])

        summary = suite.summary()

        assert "Total tests: 2" in summary
        assert "Passed: 1" in summary
        assert "Failed: 1" in summary
        assert "test2" in summary


class TestGroundTruthComparison:
    """Tests for ground truth comparison."""

    def test_cigar_match(self):
        """Test CIGAR string comparison."""
        query = "ACGT"
        target = "ACGT"

        aligner = CatalignAligner()
        aln = aligner.align(query, target)

        expected = {
            "name": "test",
            "expected_cigar": "4M",
        }

        result = compare_to_ground_truth(aln, expected)

        assert result.cigar_match

    def test_identity_check(self):
        """Test identity comparison."""
        query = "ACGTACGT"
        target = "ACGTACGT"

        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        quality = evaluate_quality(aln, query, target)
        metrics = compute_metrics(aln, query, target, quality)

        expected = {
            "name": "test",
            "expected_identity": 1.0,
        }

        result = compare_to_ground_truth(aln, expected, metrics)

        assert result.passed

    def test_identity_failure(self):
        """Test identity comparison failure."""
        query = "ACGTACGT"
        target = "ACGTACGT"

        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        quality = evaluate_quality(aln, query, target)
        metrics = compute_metrics(aln, query, target, quality)

        expected = {
            "name": "test",
            "expected_identity": 0.5,  # Way off
        }

        result = compare_to_ground_truth(aln, expected, metrics, identity_tolerance=0.01)

        assert not result.passed
        assert len(result.errors) > 0
