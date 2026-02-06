"""Tests for I/O module."""

import gzip
import pytest

from catalign.io import read_fasta, write_fasta, Sequence


@pytest.fixture
def fasta_file(tmp_path):
    p = tmp_path / "test.fa"
    p.write_text(">seq1\nACGTACGT\n>seq2\nGGGGAAAA\n")
    return p


@pytest.fixture
def fasta_gz_file(tmp_path):
    p = tmp_path / "test.fa.gz"
    with gzip.open(p, "wt") as f:
        f.write(">seq1\nACGTACGT\n>seq2\nGGGGAAAA\n")
    return p


class TestReadFasta:
    def test_read_plain(self, fasta_file):
        records = list(read_fasta(fasta_file))
        assert len(records) == 2
        assert records[0] == ("seq1", "ACGTACGT")
        assert records[1] == ("seq2", "GGGGAAAA")

    def test_read_gzipped(self, fasta_gz_file):
        records = list(read_fasta(fasta_gz_file))
        assert len(records) == 2
        assert records[0][1] == "ACGTACGT"

    def test_multiline_sequence(self, tmp_path):
        p = tmp_path / "multi.fa"
        p.write_text(">seq1\nACGT\nACGT\n")
        records = list(read_fasta(p))
        assert records[0] == ("seq1", "ACGTACGT")


class TestWriteFasta:
    def test_write_tuples(self, tmp_path):
        p = tmp_path / "out.fa"
        write_fasta(p, [("s1", "ACGT"), ("s2", "GGGG")])
        records = list(read_fasta(p))
        assert records[0] == ("s1", "ACGT")
        assert records[1] == ("s2", "GGGG")

    def test_write_sequence_objects(self, tmp_path):
        p = tmp_path / "out.fa"
        write_fasta(p, [Sequence("s1", "ACGT")])
        records = list(read_fasta(p))
        assert records[0] == ("s1", "ACGT")

    def test_roundtrip(self, tmp_path):
        p = tmp_path / "rt.fa"
        seqs = [("chr1", "ACGTACGTACGT"), ("chr2", "NNNNAAAA")]
        write_fasta(p, seqs)
        result = list(read_fasta(p))
        assert result == seqs

    def test_roundtrip_gz(self, tmp_path):
        p = tmp_path / "rt.fa.gz"
        seqs = [("s1", "ACGT")]
        write_fasta(p, seqs)
        result = list(read_fasta(p))
        assert result == seqs
