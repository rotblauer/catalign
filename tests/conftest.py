"""Shared test fixtures for Catalign tests."""

import pytest


@pytest.fixture
def simple_seq():
    """Short identical sequences."""
    return "ACGTACGTACGTACGT"


@pytest.fixture
def seq_with_mismatch():
    """Two sequences differing by a single base."""
    return ("ACGTACGTACGT", "ACGTACATACGT")


@pytest.fixture
def seq_with_insertion():
    """Query has a 3-base insertion relative to target."""
    return ("ACGTAAAACGT", "ACGTACGT")


@pytest.fixture
def seq_with_deletion():
    """Query has a deletion relative to target."""
    return ("ACGTACGT", "ACGTAAAACGT")


@pytest.fixture
def repetitive_seq():
    """Repetitive sequence."""
    return "ATATAT" * 20


@pytest.fixture
def longer_seqs():
    """Longer sequences sharing a common core."""
    import random

    random.seed(42)
    bases = "ACGT"
    common = "".join(random.choice(bases) for _ in range(200))
    query = "".join(random.choice(bases) for _ in range(50)) + common + "".join(random.choice(bases) for _ in range(50))
    target = "".join(random.choice(bases) for _ in range(30)) + common + "".join(random.choice(bases) for _ in range(30))
    return query, target
