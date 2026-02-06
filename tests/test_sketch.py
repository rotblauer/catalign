"""Tests for the sketch module."""

import pytest

from catalign.sketch import kmer_sketch, minimizer_sketch, find_anchors


class TestKmerSketch:
    def test_kmer_extraction_count(self):
        seq = "ACGTACGTACGTACGT"  # length 16
        sk = kmer_sketch(seq, k=5)
        assert len(sk.entries) == 16 - 5 + 1

    def test_kmer_positions(self):
        seq = "ACGTACGT"
        sk = kmer_sketch(seq, k=4)
        positions = [pos for _, pos, _ in sk.entries]
        assert positions == list(range(len(seq) - 4 + 1))

    def test_kmer_skips_n(self):
        seq = "ACGTNACGT"
        sk = kmer_sketch(seq, k=5)
        for _, _, kmer in sk.entries:
            assert "N" not in kmer

    def test_empty_sequence(self):
        sk = kmer_sketch("", k=5)
        assert len(sk.entries) == 0

    def test_short_sequence(self):
        sk = kmer_sketch("ACG", k=5)
        assert len(sk.entries) == 0


class TestMinimizerSketch:
    def test_minimizer_smaller_than_kmer(self):
        seq = "ACGTACGTACGTACGT" * 10  # 160 bases
        sk_full = kmer_sketch(seq, k=10)
        sk_mini = minimizer_sketch(seq, k=10, w=20)
        assert len(sk_mini.entries) <= len(sk_full.entries)
        assert len(sk_mini.entries) > 0

    def test_minimizer_empty(self):
        sk = minimizer_sketch("", k=15, w=50)
        assert len(sk.entries) == 0

    def test_minimizer_positions_sorted(self):
        seq = "ACGTACGTACGTACGT" * 5
        sk = minimizer_sketch(seq, k=10, w=20)
        positions = [pos for _, pos, _ in sk.entries]
        assert positions == sorted(positions)


class TestFindAnchors:
    def test_identical_sequences(self):
        seq = "ACGTACGTACGTACGT"
        sk_a = kmer_sketch(seq, k=5)
        sk_b = kmer_sketch(seq, k=5)
        anchors = find_anchors(sk_a, sk_b)
        assert len(anchors) > 0
        # Diagonal anchors should exist (pos_a == pos_b)
        diag = [(a, b) for a, b, _ in anchors if a == b]
        assert len(diag) > 0

    def test_no_shared_kmers(self):
        seq_a = "AAAAAAAAAA"
        seq_b = "CCCCCCCCCC"
        sk_a = kmer_sketch(seq_a, k=5)
        sk_b = kmer_sketch(seq_b, k=5)
        anchors = find_anchors(sk_a, sk_b)
        assert len(anchors) == 0

    def test_shared_region(self):
        common = "ACGTACGTACGTACGT"
        seq_a = "AAAAAAA" + common + "TTTTTTT"
        seq_b = "CCCCC" + common + "GGGGG"
        sk_a = kmer_sketch(seq_a, k=5)
        sk_b = kmer_sketch(seq_b, k=5)
        anchors = find_anchors(sk_a, sk_b)
        assert len(anchors) > 0
