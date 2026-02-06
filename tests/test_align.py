"""Tests for the core alignment engine."""

import pytest

from catalign.align import CatalignAligner, Alignment
from catalign.energy import EnergyModel


class TestAlignment:
    def test_cigar_simple(self):
        aln = Alignment(aligned_pairs=[(0, 0, "M"), (1, 1, "M"), (2, 2, "M")])
        assert aln.cigar == "3M"

    def test_cigar_mixed(self):
        aln = Alignment(aligned_pairs=[
            (0, 0, "M"), (1, 1, "X"), (2, None, "D"), (3, 2, "M"),
        ])
        assert aln.cigar == "1M1X1D1M"

    def test_cigar_empty(self):
        aln = Alignment()
        assert aln.cigar == ""


class TestCatalignAligner:
    def test_identical_sequences(self, simple_seq):
        aligner = CatalignAligner()
        aln = aligner.align(simple_seq, simple_seq)
        matches = sum(1 for _, _, op in aln.aligned_pairs if op == "M")
        assert matches == len(simple_seq)
        assert aln.energy_score < 0

    def test_single_mismatch(self, seq_with_mismatch):
        query, target = seq_with_mismatch
        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        mismatches = sum(1 for _, _, op in aln.aligned_pairs if op == "X")
        assert mismatches >= 1

    def test_insertion(self, seq_with_insertion):
        query, target = seq_with_insertion
        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        ops = [op for _, _, op in aln.aligned_pairs]
        # Should contain deletions in query or insertions
        assert "D" in ops or "I" in ops or "X" in ops

    def test_deletion(self, seq_with_deletion):
        query, target = seq_with_deletion
        aligner = CatalignAligner()
        aln = aligner.align(query, target)
        ops = [op for _, _, op in aln.aligned_pairs]
        assert "D" in ops or "I" in ops or "X" in ops

    def test_cigar_string_generated(self, simple_seq):
        aligner = CatalignAligner()
        aln = aligner.align(simple_seq, simple_seq)
        assert len(aln.cigar) > 0
        assert "M" in aln.cigar

    def test_empty_sequences(self):
        aligner = CatalignAligner()
        aln = aligner.align("", "")
        assert aln.aligned_pairs == []

    def test_end_to_end_synthetic(self, longer_seqs):
        query, target = longer_seqs
        aligner = CatalignAligner(k=10, w=20)
        aln = aligner.align(query, target)
        assert len(aln.aligned_pairs) > 0
        matches = sum(1 for _, _, op in aln.aligned_pairs if op == "M")
        # Should align a significant portion of the common region
        assert matches > 50

    def test_query_only_empty(self):
        aligner = CatalignAligner()
        aln = aligner.align("", "ACGT")
        assert all(op == "I" for _, _, op in aln.aligned_pairs)

    def test_target_only_empty(self):
        aligner = CatalignAligner()
        aln = aligner.align("ACGT", "")
        assert all(op == "D" for _, _, op in aln.aligned_pairs)
