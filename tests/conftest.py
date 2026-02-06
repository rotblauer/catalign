"""Shared test fixtures for Catalign tests."""

import json
import pytest
from pathlib import Path


# Path to test resources
RESOURCES_DIR = Path(__file__).parent / "resources"


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


@pytest.fixture
def resources_dir():
    """Path to test resources directory."""
    return RESOURCES_DIR


@pytest.fixture
def ground_truth_manifest():
    """Load and return the ground truth manifest if available."""
    manifest_path = RESOURCES_DIR / "ground_truth" / "manifest.json"
    if not manifest_path.exists():
        pytest.skip("Ground truth manifest not found. Run: python scripts/generate_test_data.py")
    with open(manifest_path) as f:
        return json.load(f)


@pytest.fixture
def aligner():
    """Default CatalignAligner instance."""
    from catalign import CatalignAligner
    return CatalignAligner()


@pytest.fixture
def energy_model():
    """Default EnergyModel instance."""
    from catalign import EnergyModel
    return EnergyModel()


@pytest.fixture
def custom_energy_model():
    """Custom energy model for testing parameter sensitivity."""
    from catalign import EnergyModel
    return EnergyModel(
        match_energy=-3.0,
        mismatch_energy=4.0,
        gap_open_energy=8.0,
        gap_extend_energy=2.0,
    )


# Fixtures for synthetic test sequences
@pytest.fixture
def identical_100bp():
    """100bp identical sequences."""
    import random
    random.seed(100)
    seq = "".join(random.choice("ACGT") for _ in range(100))
    return seq, seq


@pytest.fixture
def snp_sequence_pair():
    """Sequence pair with a single SNP at position 50."""
    import random
    random.seed(200)
    seq = "".join(random.choice("ACGT") for _ in range(100))
    mutated = list(seq)
    # Introduce SNP at position 50
    original = mutated[50]
    mutated[50] = {"A": "T", "T": "A", "C": "G", "G": "C"}[original]
    return seq, "".join(mutated)


@pytest.fixture
def small_insertion_pair():
    """Sequence pair with a 5bp insertion."""
    import random
    random.seed(300)
    target = "".join(random.choice("ACGT") for _ in range(100))
    insert = "".join(random.choice("ACGT") for _ in range(5))
    query = target[:50] + insert + target[50:]
    return query, target


@pytest.fixture
def small_deletion_pair():
    """Sequence pair with a 5bp deletion."""
    import random
    random.seed(400)
    target = "".join(random.choice("ACGT") for _ in range(100))
    query = target[:50] + target[55:]  # Delete 5bp at position 50
    return query, target

