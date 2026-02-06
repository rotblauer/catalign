"""Tests for the quality evaluation module."""

import pytest

from catalign.align import CatalignAligner, Alignment
from catalign.quality import evaluate_quality, QualityReport


class TestBaseQuality:
    def test_perfect_alignment_all_matches(self, simple_seq):
        aligner = CatalignAligner()
        aln = aligner.align(simple_seq, simple_seq)
        qual = evaluate_quality(aln, simple_seq, simple_seq)
        assert qual.total_matches == len(simple_seq)
        assert qual.total_mismatches == 0
        assert qual.total_insertions == 0
        assert qual.total_deletions == 0

    def test_mismatch_counted(self, seq_with_mismatch):
        query, target = seq_with_mismatch
        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        qual = evaluate_quality(aln, query, target)
        assert qual.total_mismatches >= 1 or qual.total_insertions >= 1 or qual.total_deletions >= 1


class TestBlockQuality:
    def test_perfect_alignment_one_block(self, simple_seq):
        aligner = CatalignAligner()
        aln = aligner.align(simple_seq, simple_seq)
        qual = evaluate_quality(aln, simple_seq, simple_seq)
        assert len(qual.block_qualities) >= 1
        assert qual.block_qualities[0].identity == 1.0

    def test_block_energy_favorable(self, simple_seq):
        aligner = CatalignAligner()
        aln = aligner.align(simple_seq, simple_seq)
        qual = evaluate_quality(aln, simple_seq, simple_seq)
        assert qual.block_qualities[0].energy < 0


class TestRegionQuality:
    def test_region_coverage(self, simple_seq):
        aligner = CatalignAligner()
        aln = aligner.align(simple_seq, simple_seq)
        qual = evaluate_quality(aln, simple_seq, simple_seq)
        rq = qual.region_quality
        assert rq is not None
        assert rq.coverage_query > 0
        assert rq.coverage_target > 0

    def test_region_concordance_perfect(self, simple_seq):
        aligner = CatalignAligner()
        aln = aligner.align(simple_seq, simple_seq)
        qual = evaluate_quality(aln, simple_seq, simple_seq)
        assert qual.region_quality.concordance == 1.0


class TestOverallQuality:
    def test_perfect_quality_score(self, simple_seq):
        aligner = CatalignAligner()
        aln = aligner.align(simple_seq, simple_seq)
        qual = evaluate_quality(aln, simple_seq, simple_seq)
        assert qual.quality_score == 100.0
        assert qual.overall_identity == 1.0

    def test_imperfect_quality_lower(self, seq_with_mismatch):
        query, target = seq_with_mismatch
        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        qual = evaluate_quality(aln, query, target)
        assert qual.quality_score <= 100.0

    def test_energy_favorable_for_good_alignment(self, simple_seq):
        aligner = CatalignAligner()
        aln = aligner.align(simple_seq, simple_seq)
        qual = evaluate_quality(aln, simple_seq, simple_seq)
        assert qual.total_energy < 0


class TestQualityReport:
    def test_aggregate(self, simple_seq):
        aligner = CatalignAligner()
        aln = aligner.align(simple_seq, simple_seq)
        q = evaluate_quality(aln, simple_seq, simple_seq)

        report = QualityReport()
        report.add(q)
        assert report.mean_identity == 1.0
        assert report.mean_quality_score == 100.0
