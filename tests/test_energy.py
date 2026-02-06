"""Tests for the energy module."""

import numpy as np
import pytest

from catalign.energy import EnergyModel


class TestEnergyModel:
    def test_match_energy_is_favorable(self):
        em = EnergyModel()
        assert em.compute_base_energy("A", "A") < 0
        assert em.compute_base_energy("C", "C") < 0
        assert em.compute_base_energy("G", "G") < 0
        assert em.compute_base_energy("T", "T") < 0

    def test_mismatch_energy_is_unfavorable(self):
        em = EnergyModel()
        assert em.compute_base_energy("A", "C") > 0
        assert em.compute_base_energy("G", "T") > 0

    def test_transition_vs_transversion(self):
        em = EnergyModel()
        # Transitions (A<->G, C<->T) should be less penalised
        transition = em.compute_base_energy("A", "G")
        transversion = em.compute_base_energy("A", "C")
        assert transition < transversion

    def test_gap_penalties(self):
        em = EnergyModel()
        assert em.gap_open_energy > 0
        assert em.gap_extend_energy > 0
        assert em.gap_open_energy > em.gap_extend_energy

    def test_block_energy_identical(self):
        em = EnergyModel()
        energy = em.compute_block_energy("ACGT", "ACGT")
        assert energy < 0  # all matches -> favorable

    def test_block_energy_mismatches(self):
        em = EnergyModel()
        e_match = em.compute_block_energy("AAAA", "AAAA")
        e_mismatch = em.compute_block_energy("AAAA", "CCCC")
        assert e_match < e_mismatch

    def test_block_energy_length_mismatch_raises(self):
        em = EnergyModel()
        with pytest.raises(ValueError):
            em.compute_block_energy("ACG", "AC")

    def test_custom_matrix(self):
        matrix = np.full((5, 5), 10.0)
        np.fill_diagonal(matrix, -5.0)
        em = EnergyModel(custom_matrix=matrix)
        assert em.compute_base_energy("A", "A") == -5.0
        assert em.compute_base_energy("A", "C") == 10.0

    def test_long_range_attraction(self):
        assert EnergyModel.long_range_attraction(0) == 1.0
        close = EnergyModel.long_range_attraction(100)
        far = EnergyModel.long_range_attraction(10000)
        assert close > far
        assert 0 < far < 1

    def test_n_bases_neutral(self):
        em = EnergyModel()
        assert em.compute_base_energy("N", "A") == 0.0
        assert em.compute_base_energy("A", "N") == 0.0
