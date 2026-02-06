"""Anchor chaining module – groups co-linear anchors via dynamic programming."""

from __future__ import annotations

from typing import List, Tuple

from catalign.energy import EnergyModel


def chain_anchors(
    anchors: List[Tuple[int, int, str]],
    max_gap: int = 10000,
    band_width: int = 500,
) -> List[List[Tuple[int, int, str]]]:
    """Chain compatible anchors into co-linear groups.

    Each anchor is (pos_a, pos_b, kmer).  Two anchors are compatible when
    the second is downstream of the first on both sequences and within
    *max_gap* / *band_width* constraints.

    Returns a list of chains (each a sorted list of anchors).
    """
    if not anchors:
        return []

    n = len(anchors)
    # Sort by position in query
    anchors = sorted(anchors, key=lambda a: (a[0], a[1]))

    dp = [1] * n  # best chain length ending at i
    parent = [-1] * n

    for i in range(1, n):
        for j in range(i - 1, max(-1, i - 500), -1):  # look-back window
            pa, pb, _ = anchors[j]
            ca, cb, _ = anchors[i]
            da = ca - pa
            db = cb - pb
            if da <= 0 or db <= 0:
                continue
            if da > max_gap or db > max_gap:
                continue
            if abs(da - db) > band_width:
                continue
            if dp[j] + 1 > dp[i]:
                dp[i] = dp[j] + 1
                parent[i] = j

    # Traceback: greedily extract chains starting from the best endpoints
    used = [False] * n
    chains: List[List[Tuple[int, int, str]]] = []

    indices = sorted(range(n), key=lambda i: dp[i], reverse=True)
    for idx in indices:
        if used[idx]:
            continue
        chain: List[Tuple[int, int, str]] = []
        cur = idx
        while cur != -1 and not used[cur]:
            chain.append(anchors[cur])
            used[cur] = True
            cur = parent[cur]
        if chain:
            chain.reverse()
            chains.append(chain)

    # Filter tiny chains
    chains = [c for c in chains if len(c) >= 2]
    chains.sort(key=lambda c: c[0][0])
    return chains


def score_chain(chain: List[Tuple[int, int, str]], energy_model: EnergyModel) -> float:
    """Score a chain based on anchor count and consistency.

    Lower (more negative) scores are better, matching the energy convention.
    """
    if not chain:
        return 0.0

    n = len(chain)
    # Reward per anchor
    score = n * energy_model.match_energy

    # Penalise gaps between consecutive anchors
    for i in range(1, n):
        da = chain[i][0] - chain[i - 1][0]
        db = chain[i][1] - chain[i - 1][1]
        gap = abs(da - db)
        if gap > 0:
            score += energy_model.gap_open_energy + (gap - 1) * energy_model.gap_extend_energy

    return score
