"""Energy/scoring module mimicking molecular forces for sequence alignment."""

import numpy as np


# Base encoding: A=0, C=1, G=2, T=3, N=4
BASE_TO_INT = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4,
               "a": 0, "c": 1, "g": 2, "t": 3, "n": 4}

# Transition pairs (purine<->purine or pyrimidine<->pyrimidine)
_TRANSITIONS = {(0, 2), (2, 0), (1, 3), (3, 1)}  # A<->G, C<->T


class EnergyModel:
    """Scoring parameters mimicking molecular forces.

    Lower (more negative) energy = more favorable alignment, like
    hydrogen bonding stabilising a DNA duplex.
    """

    def __init__(
        self,
        match_energy: float = -2.0,
        mismatch_energy: float = 3.0,
        gap_open_energy: float = 5.0,
        gap_extend_energy: float = 1.0,
        transition_energy: float = 2.0,
        n_energy: float = 0.0,
        custom_matrix: np.ndarray | None = None,
    ):
        self.match_energy = match_energy
        self.mismatch_energy = mismatch_energy
        self.gap_open_energy = gap_open_energy
        self.gap_extend_energy = gap_extend_energy
        self.transition_energy = transition_energy
        self.n_energy = n_energy

        if custom_matrix is not None:
            self.matrix = np.array(custom_matrix, dtype=float)
        else:
            # Build default 5x5 matrix (A, C, G, T, N)
            self.matrix = np.full((5, 5), mismatch_energy, dtype=float)
            np.fill_diagonal(self.matrix, match_energy)
            for i, j in _TRANSITIONS:
                self.matrix[i, j] = transition_energy
            # N vs anything = neutral
            self.matrix[4, :] = n_energy
            self.matrix[:, 4] = n_energy

    def compute_base_energy(self, base_a: str, base_b: str) -> float:
        """Return the energy for aligning *base_a* against *base_b*."""
        i = BASE_TO_INT.get(base_a, 4)
        j = BASE_TO_INT.get(base_b, 4)
        return float(self.matrix[i, j])

    def compute_block_energy(self, seq_a: str, seq_b: str) -> float:
        """Total energy for an aligned block (same-length strings)."""
        if len(seq_a) != len(seq_b):
            raise ValueError("Aligned blocks must be the same length")
        total = 0.0
        for a, b in zip(seq_a, seq_b):
            total += self.compute_base_energy(a, b)
        return total

    @staticmethod
    def long_range_attraction(distance: float, scale: float = 1000.0) -> float:
        """Inverse-distance attraction for coarse-grained matching.

        Returns a value in (0, 1] that decays with distance.
        """
        if distance < 0:
            distance = -distance
        return scale / (scale + distance)
