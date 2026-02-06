"""Tests for the anchor chaining module."""

import pytest

from catalign.chain import chain_anchors, score_chain
from catalign.energy import EnergyModel


class TestChainAnchors:
    def test_colinear_anchors(self):
        anchors = [
            (0, 0, "AAAAA"),
            (10, 10, "CCCCC"),
            (20, 20, "GGGGG"),
            (30, 30, "TTTTT"),
        ]
        chains = chain_anchors(anchors)
        assert len(chains) >= 1
        best = max(chains, key=len)
        assert len(best) >= 3

    def test_non_colinear_separated(self):
        # Two groups on different diagonals
        anchors = [
            (0, 0, "AAAAA"),
            (10, 10, "CCCCC"),
            (0, 5000, "GGGGG"),
            (10, 5010, "TTTTT"),
        ]
        chains = chain_anchors(anchors, band_width=100)
        assert len(chains) >= 1

    def test_empty_anchors(self):
        chains = chain_anchors([])
        assert chains == []

    def test_single_anchor(self):
        chains = chain_anchors([(0, 0, "AAAAA")])
        # Single anchor won't form a chain (min length 2)
        assert len(chains) == 0


class TestScoreChain:
    def test_perfect_chain_score(self):
        em = EnergyModel()
        chain = [
            (0, 0, "AAAAA"),
            (10, 10, "CCCCC"),
            (20, 20, "GGGGG"),
        ]
        score = score_chain(chain, em)
        # All anchors co-linear with equal gaps -> only match rewards
        assert score < 0  # favorable

    def test_gapped_chain_worse(self):
        em = EnergyModel()
        perfect = [
            (0, 0, "AAAAA"),
            (10, 10, "CCCCC"),
        ]
        gapped = [
            (0, 0, "AAAAA"),
            (10, 15, "CCCCC"),  # offset -> gap penalty
        ]
        assert score_chain(perfect, em) < score_chain(gapped, em)

    def test_empty_chain_score(self):
        em = EnergyModel()
        assert score_chain([], em) == 0.0
