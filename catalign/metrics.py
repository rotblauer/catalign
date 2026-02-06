"""Consolidated alignment metrics and benchmarking utilities."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple, TYPE_CHECKING
import json
from pathlib import Path

if TYPE_CHECKING:
    from catalign.align import Alignment
    from catalign.quality import AlignmentQuality


@dataclass
class AlignmentMetrics:
    """Comprehensive metrics for a single alignment."""

    # Basic metrics
    alignment_length: int = 0
    matches: int = 0
    mismatches: int = 0
    insertions: int = 0
    deletions: int = 0

    # Derived metrics
    identity: float = 0.0  # matches / (matches + mismatches)
    similarity: float = 0.0  # (matches + mismatches) / alignment_length
    gap_rate: float = 0.0  # (insertions + deletions) / alignment_length

    # Coverage metrics
    query_coverage: float = 0.0
    target_coverage: float = 0.0

    # Energy metrics
    total_energy: float = 0.0
    energy_per_base: float = 0.0

    # Quality metrics
    quality_score: float = 0.0

    # Alignment coordinates
    query_start: int = 0
    query_end: int = 0
    target_start: int = 0
    target_end: int = 0

    def to_dict(self) -> Dict:
        """Convert metrics to dictionary."""
        return {
            "alignment_length": self.alignment_length,
            "matches": self.matches,
            "mismatches": self.mismatches,
            "insertions": self.insertions,
            "deletions": self.deletions,
            "identity": self.identity,
            "similarity": self.similarity,
            "gap_rate": self.gap_rate,
            "query_coverage": self.query_coverage,
            "target_coverage": self.target_coverage,
            "total_energy": self.total_energy,
            "energy_per_base": self.energy_per_base,
            "quality_score": self.quality_score,
            "query_start": self.query_start,
            "query_end": self.query_end,
            "target_start": self.target_start,
            "target_end": self.target_end,
        }

    def to_json(self, filepath: Optional[str] = None) -> str:
        """Export metrics to JSON."""
        json_str = json.dumps(self.to_dict(), indent=2)
        if filepath:
            Path(filepath).write_text(json_str)
        return json_str

    @classmethod
    def from_dict(cls, data: Dict) -> "AlignmentMetrics":
        """Create from dictionary."""
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass
class BenchmarkResult:
    """Result of comparing alignment against ground truth."""

    name: str
    passed: bool = True
    metrics: Optional[AlignmentMetrics] = None
    expected: Optional[Dict] = None
    errors: List[str] = field(default_factory=list)

    # Precision/recall for alignment positions
    position_precision: float = 0.0
    position_recall: float = 0.0
    position_f1: float = 0.0

    # CIGAR comparison
    cigar_match: bool = False
    cigar_similarity: float = 0.0


@dataclass
class BenchmarkSuite:
    """Collection of benchmark results."""

    results: List[BenchmarkResult] = field(default_factory=list)

    @property
    def passed(self) -> int:
        return sum(1 for r in self.results if r.passed)

    @property
    def failed(self) -> int:
        return sum(1 for r in self.results if not r.passed)

    @property
    def total(self) -> int:
        return len(self.results)

    @property
    def pass_rate(self) -> float:
        return self.passed / max(self.total, 1)

    def summary(self) -> str:
        """Generate summary report."""
        lines = [
            "=" * 60,
            "Benchmark Summary",
            "=" * 60,
            f"Total tests: {self.total}",
            f"Passed: {self.passed}",
            f"Failed: {self.failed}",
            f"Pass rate: {self.pass_rate:.1%}",
            "",
        ]

        if self.failed > 0:
            lines.append("Failed tests:")
            for r in self.results:
                if not r.passed:
                    lines.append(f"  - {r.name}")
                    for err in r.errors:
                        lines.append(f"      {err}")

        return "\n".join(lines)


def compute_metrics(
    alignment: "Alignment",
    query_seq: str,
    target_seq: str,
    quality: Optional["AlignmentQuality"] = None,
) -> AlignmentMetrics:
    """Compute comprehensive metrics for an alignment.

    Parameters
    ----------
    alignment : Alignment
        The alignment to analyze
    query_seq : str
        Query sequence
    target_seq : str
        Target sequence
    quality : AlignmentQuality, optional
        Pre-computed quality assessment

    Returns
    -------
    AlignmentMetrics
        Comprehensive alignment metrics
    """
    pairs = alignment.aligned_pairs

    matches = sum(1 for _, _, op in pairs if op == "M")
    mismatches = sum(1 for _, _, op in pairs if op == "X")
    insertions = sum(1 for _, _, op in pairs if op == "I")
    deletions = sum(1 for _, _, op in pairs if op == "D")

    aln_length = len(pairs)
    aligned_bases = matches + mismatches

    # Compute coverage
    q_positions = {qp for qp, _, op in pairs if qp is not None}
    t_positions = {tp for _, tp, op in pairs if tp is not None}

    q_coverage = len(q_positions) / max(len(query_seq), 1)
    t_coverage = len(t_positions) / max(len(target_seq), 1)

    # Get energy from quality if available
    total_energy = 0.0
    q_score = 0.0
    if quality:
        total_energy = quality.total_energy
        q_score = quality.quality_score

    return AlignmentMetrics(
        alignment_length=aln_length,
        matches=matches,
        mismatches=mismatches,
        insertions=insertions,
        deletions=deletions,
        identity=matches / max(aligned_bases, 1),
        similarity=aligned_bases / max(aln_length, 1),
        gap_rate=(insertions + deletions) / max(aln_length, 1),
        query_coverage=q_coverage,
        target_coverage=t_coverage,
        total_energy=total_energy,
        energy_per_base=total_energy / max(aln_length, 1),
        quality_score=q_score,
        query_start=alignment.query_start,
        query_end=alignment.query_end,
        target_start=alignment.target_start,
        target_end=alignment.target_end,
    )


def compare_to_ground_truth(
    alignment: "Alignment",
    expected: Dict,
    metrics: Optional[AlignmentMetrics] = None,
    identity_tolerance: float = 0.02,
    position_tolerance: int = 5,
) -> BenchmarkResult:
    """Compare alignment against expected ground truth.

    Parameters
    ----------
    alignment : Alignment
        The alignment to evaluate
    expected : Dict
        Expected values (from ground truth JSON)
    metrics : AlignmentMetrics, optional
        Pre-computed metrics
    identity_tolerance : float
        Tolerance for identity comparison
    position_tolerance : int
        Tolerance for position comparison (bp)

    Returns
    -------
    BenchmarkResult
        Benchmark comparison result
    """
    result = BenchmarkResult(
        name=expected.get("name", "unnamed"),
        metrics=metrics,
        expected=expected,
    )

    errors = []

    # Check CIGAR if expected
    if "expected_cigar" in expected and expected["expected_cigar"] is not None:
        expected_cigar = expected["expected_cigar"]
        result.cigar_match = alignment.cigar == expected_cigar
        if not result.cigar_match:
            # Compute CIGAR similarity
            result.cigar_similarity = _cigar_similarity(
                alignment.cigar, expected_cigar
            )
            if result.cigar_similarity < 0.9:
                errors.append(
                    f"CIGAR mismatch: got {alignment.cigar[:30]}..., "
                    f"expected {expected_cigar[:30]}..."
                )

    # Check identity if expected
    if "expected_identity" in expected and expected["expected_identity"] is not None and metrics:
        expected_id = expected["expected_identity"]
        actual_id = metrics.identity
        if abs(actual_id - expected_id) > identity_tolerance:
            errors.append(
                f"Identity mismatch: got {actual_id:.4f}, "
                f"expected {expected_id:.4f} (±{identity_tolerance})"
            )

    # Check matches if expected (with 5% tolerance)
    if "expected_matches" in expected and expected["expected_matches"] is not None and metrics:
        expected_m = expected["expected_matches"]
        actual_m = metrics.matches
        tolerance = max(2, int(expected_m * 0.05))  # 5% or at least 2
        if abs(actual_m - expected_m) > tolerance:
            errors.append(
                f"Match count mismatch: got {actual_m}, expected {expected_m} (±{tolerance})"
            )

    # Check mismatches if expected (with tolerance)
    if "expected_mismatches" in expected and expected["expected_mismatches"] is not None and metrics:
        expected_mm = expected["expected_mismatches"]
        actual_mm = metrics.mismatches
        tolerance = max(2, int(expected_mm * 0.1) if expected_mm > 0 else 2)  # 10% or at least 2
        if abs(actual_mm - expected_mm) > tolerance:
            errors.append(
                f"Mismatch count mismatch: got {actual_mm}, expected {expected_mm} (±{tolerance})"
            )

    result.errors = errors
    result.passed = len(errors) == 0

    return result


def _cigar_similarity(cigar1: str, cigar2: str) -> float:
    """Compute similarity between two CIGAR strings."""
    import re

    # Handle None or empty strings
    if not cigar1 or not cigar2:
        return 0.0 if (cigar1 or cigar2) else 1.0

    def parse(cigar: str) -> List[Tuple[int, str]]:
        return [(int(m.group(1)), m.group(2))
                for m in re.finditer(r"(\d+)([MIDNSHP=X])", cigar)]

    ops1 = parse(cigar1)
    ops2 = parse(cigar2)

    # Expand to base-level operations
    expand1 = "".join(op * length for length, op in ops1)
    expand2 = "".join(op * length for length, op in ops2)

    # Compare character by character
    matches = sum(1 for a, b in zip(expand1, expand2) if a == b)
    max_len = max(len(expand1), len(expand2))

    return matches / max(max_len, 1)


def run_benchmark_suite(
    test_cases: List[Dict],
    aligner,
    resources_dir: Path,
) -> BenchmarkSuite:
    """Run a suite of benchmark tests.

    Parameters
    ----------
    test_cases : List[Dict]
        List of test case definitions
    aligner : CatalignAligner
        Configured aligner instance
    resources_dir : Path
        Directory containing test resources

    Returns
    -------
    BenchmarkSuite
        Results from all benchmarks
    """
    from catalign.io import read_fasta
    from catalign.quality import evaluate_quality

    suite = BenchmarkSuite()

    for case in test_cases:
        try:
            # Load sequences
            query_path = resources_dir / case["query_file"]
            target_path = resources_dir / case["target_file"]

            query_records = list(read_fasta(query_path))
            target_records = list(read_fasta(target_path))

            if not query_records or not target_records:
                suite.results.append(BenchmarkResult(
                    name=case.get("name", "unnamed"),
                    passed=False,
                    errors=["Could not load sequences"],
                ))
                continue

            _, query_seq = query_records[0]
            _, target_seq = target_records[0]

            # Run alignment
            alignment = aligner.align(query_seq, target_seq)
            quality = evaluate_quality(alignment, query_seq, target_seq)

            # Compute metrics
            metrics = compute_metrics(alignment, query_seq, target_seq, quality)

            # Compare to ground truth
            result = compare_to_ground_truth(alignment, case, metrics)
            suite.results.append(result)

        except Exception as e:
            suite.results.append(BenchmarkResult(
                name=case.get("name", "unnamed"),
                passed=False,
                errors=[f"Exception: {str(e)}"],
            ))

    return suite


def generate_metrics_report(
    alignments: List[Tuple["Alignment", str, str]],
    output_path: Optional[Path] = None,
) -> str:
    """Generate a comprehensive metrics report for multiple alignments.

    Parameters
    ----------
    alignments : List[Tuple[Alignment, str, str]]
        List of (alignment, query_seq, target_seq) tuples
    output_path : Path, optional
        Path to save the report

    Returns
    -------
    str
        Formatted report
    """
    from catalign.quality import evaluate_quality

    lines = [
        "=" * 80,
        "Catalign Alignment Metrics Report",
        "=" * 80,
        "",
    ]

    all_metrics = []

    for i, (aln, q_seq, t_seq) in enumerate(alignments, 1):
        quality = evaluate_quality(aln, q_seq, t_seq)
        metrics = compute_metrics(aln, q_seq, t_seq, quality)
        all_metrics.append(metrics)

        lines.extend([
            f"Alignment {i}: {aln.query_name} vs {aln.target_name}",
            "-" * 60,
            f"  Length:        {metrics.alignment_length:,} bp",
            f"  Matches:       {metrics.matches:,}",
            f"  Mismatches:    {metrics.mismatches:,}",
            f"  Insertions:    {metrics.insertions:,}",
            f"  Deletions:     {metrics.deletions:,}",
            f"  Identity:      {metrics.identity:.4f} ({metrics.identity*100:.2f}%)",
            f"  Gap rate:      {metrics.gap_rate:.4f}",
            f"  Query cov:     {metrics.query_coverage:.4f}",
            f"  Target cov:    {metrics.target_coverage:.4f}",
            f"  Total energy:  {metrics.total_energy:.2f}",
            f"  Quality score: {metrics.quality_score:.1f}/100",
            "",
        ])

    # Summary statistics
    if len(all_metrics) > 1:
        lines.extend([
            "=" * 60,
            "Summary Statistics",
            "=" * 60,
            f"  Total alignments: {len(all_metrics)}",
            f"  Mean identity:    {sum(m.identity for m in all_metrics) / len(all_metrics):.4f}",
            f"  Mean quality:     {sum(m.quality_score for m in all_metrics) / len(all_metrics):.1f}",
            f"  Total energy:     {sum(m.total_energy for m in all_metrics):.2f}",
            "",
        ])

    report = "\n".join(lines)

    if output_path:
        output_path.write_text(report)

    return report
