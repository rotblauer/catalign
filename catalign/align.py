"""Core alignment engine – multi-scale energy-based sequence alignment."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple, Optional

import numpy as np

from catalign.energy import EnergyModel, BASE_TO_INT
from catalign.sketch import kmer_sketch, minimizer_sketch, find_anchors
from catalign.chain import chain_anchors


@dataclass
class Alignment:
    """Stores the result of an alignment."""

    query_name: str = ""
    target_name: str = ""
    query_start: int = 0
    query_end: int = 0
    target_start: int = 0
    target_end: int = 0
    aligned_pairs: List[Tuple[int, int, str]] = field(default_factory=list)
    chains: List[List[Tuple[int, int, str]]] = field(default_factory=list)
    energy_score: float = 0.0
    strand: str = "+"

    @property
    def cigar(self) -> str:
        """Generate a CIGAR string from aligned_pairs."""
        if not self.aligned_pairs:
            return ""
        ops: list[tuple[str, int]] = []
        for _, _, op in self.aligned_pairs:
            if ops and ops[-1][0] == op:
                ops[-1] = (op, ops[-1][1] + 1)
            else:
                ops.append((op, 1))
        return "".join(f"{count}{op}" for op, count in ops)


class CatalignAligner:
    """Main alignment class orchestrating multi-scale alignment."""

    def __init__(
        self,
        energy_model: Optional[EnergyModel] = None,
        k: int = 15,
        w: int = 50,
    ):
        self.energy_model = energy_model or EnergyModel()
        self.k = k
        self.w = w

    def align(
        self,
        query: str,
        target: str,
        query_name: str = "query",
        target_name: str = "target",
    ) -> Alignment:
        """Full alignment pipeline.

        1. Sketch both sequences
        2. Find anchors (long-range attraction)
        3. Chain anchors
        4. Banded alignment between anchors (short-range forces)
        5. Return Alignment result
        """
        q = query.upper()
        t = target.upper()

        # For short sequences, skip sketching and do direct alignment
        if len(q) < self.k + self.w or len(t) < self.k + self.w:
            pairs, score = self._banded_align(q, t, bandwidth=max(len(q), len(t)))
            aln = Alignment(
                query_name=query_name,
                target_name=target_name,
                aligned_pairs=pairs,
                energy_score=score,
            )
            if pairs:
                q_positions = [p[0] for p in pairs if p[0] is not None]
                t_positions = [p[1] for p in pairs if p[1] is not None]
                aln.query_start = min(q_positions) if q_positions else 0
                aln.query_end = max(q_positions) + 1 if q_positions else 0
                aln.target_start = min(t_positions) if t_positions else 0
                aln.target_end = max(t_positions) + 1 if t_positions else 0
            return aln

        # Long-range: sketch and anchor
        sketch_q = minimizer_sketch(q, k=self.k, w=self.w)
        sketch_t = minimizer_sketch(t, k=self.k, w=self.w)
        anchors = find_anchors(sketch_q, sketch_t)

        # Chain anchors
        chains = chain_anchors(anchors)

        if not chains:
            # Fallback: direct banded alignment
            pairs, score = self._banded_align(q, t, bandwidth=max(100, len(q) // 10))
            aln = Alignment(
                query_name=query_name,
                target_name=target_name,
                aligned_pairs=pairs,
                energy_score=score,
            )
            if pairs:
                q_positions = [p[0] for p in pairs if p[0] is not None]
                t_positions = [p[1] for p in pairs if p[1] is not None]
                aln.query_start = min(q_positions) if q_positions else 0
                aln.query_end = max(q_positions) + 1 if q_positions else 0
                aln.target_start = min(t_positions) if t_positions else 0
                aln.target_end = max(t_positions) + 1 if t_positions else 0
            return aln

        # Use best chain (longest)
        best_chain = max(chains, key=len)

        # Determine alignment region from chain
        q_start = best_chain[0][0]
        q_end = best_chain[-1][0] + self.k
        t_start = best_chain[0][1]
        t_end = best_chain[-1][1] + self.k

        # Add flanking context
        flank = 100
        q_start_ctx = max(0, q_start - flank)
        q_end_ctx = min(len(q), q_end + flank)
        t_start_ctx = max(0, t_start - flank)
        t_end_ctx = min(len(t), t_end + flank)

        sub_q = q[q_start_ctx:q_end_ctx]
        sub_t = t[t_start_ctx:t_end_ctx]

        pairs, score = self._banded_align(sub_q, sub_t, bandwidth=200)

        # Shift coordinates back to global
        global_pairs = []
        for qp, tp, op in pairs:
            gq = qp + q_start_ctx if qp is not None else None
            gt = tp + t_start_ctx if tp is not None else None
            global_pairs.append((gq, gt, op))

        aln = Alignment(
            query_name=query_name,
            target_name=target_name,
            query_start=q_start_ctx,
            query_end=q_end_ctx,
            target_start=t_start_ctx,
            target_end=t_end_ctx,
            aligned_pairs=global_pairs,
            chains=chains,
            energy_score=score,
        )
        return aln

    def _banded_align(
        self, seq_a: str, seq_b: str, bandwidth: int = 100
    ) -> Tuple[List[Tuple[Optional[int], Optional[int], str]], float]:
        """Banded DP alignment with energy-based scoring.

        Returns (aligned_pairs, total_energy).
        Lower energy = better alignment.
        """
        n = len(seq_a)
        m = len(seq_b)

        if n == 0 and m == 0:
            return [], 0.0
        if n == 0:
            pairs = [(None, j, "I") for j in range(m)]
            score = self.energy_model.gap_open_energy + (m - 1) * self.energy_model.gap_extend_energy if m > 0 else 0.0
            return pairs, score
        if m == 0:
            pairs = [(i, None, "D") for i in range(n)]
            score = self.energy_model.gap_open_energy + (n - 1) * self.energy_model.gap_extend_energy if n > 0 else 0.0
            return pairs, score

        em = self.energy_model
        INF = float("inf")

        # DP matrices: H = match/mismatch, E = gap in seq_b (deletion), F = gap in seq_a (insertion)
        H = np.full((n + 1, m + 1), INF)
        E = np.full((n + 1, m + 1), INF)
        F = np.full((n + 1, m + 1), INF)

        H[0, 0] = 0.0
        for j in range(1, m + 1):
            if abs(j) <= bandwidth:
                H[0, j] = em.gap_open_energy + (j - 1) * em.gap_extend_energy
        for i in range(1, n + 1):
            if abs(i) <= bandwidth:
                H[i, 0] = em.gap_open_energy + (i - 1) * em.gap_extend_energy

        # Encode sequences for fast lookup
        enc_a = [BASE_TO_INT.get(c, 4) for c in seq_a]
        enc_b = [BASE_TO_INT.get(c, 4) for c in seq_b]

        for i in range(1, n + 1):
            j_center = int(round(i * m / n)) if n > 0 else 0
            j_lo = max(1, j_center - bandwidth)
            j_hi = min(m, j_center + bandwidth)
            for j in range(j_lo, j_hi + 1):
                # Substitution / match
                sub_cost = float(em.matrix[enc_a[i - 1], enc_b[j - 1]])
                h_diag = H[i - 1, j - 1] if H[i - 1, j - 1] < INF else INF
                diag = h_diag + sub_cost if h_diag < INF else INF

                # Deletion (gap in seq_b)
                e1 = H[i - 1, j] + em.gap_open_energy if H[i - 1, j] < INF else INF
                e2 = E[i - 1, j] + em.gap_extend_energy if E[i - 1, j] < INF else INF
                E[i, j] = min(e1, e2)

                # Insertion (gap in seq_a)
                f1 = H[i, j - 1] + em.gap_open_energy if H[i, j - 1] < INF else INF
                f2 = F[i, j - 1] + em.gap_extend_energy if F[i, j - 1] < INF else INF
                F[i, j] = min(f1, f2)

                H[i, j] = min(diag, E[i, j], F[i, j])

        # Traceback
        pairs: List[Tuple[Optional[int], Optional[int], str]] = []
        i, j = n, m
        total_energy = H[n, m] if H[n, m] < INF else 0.0

        while i > 0 or j > 0:
            if i > 0 and j > 0 and H[i, j] < INF:
                sub_cost = float(em.matrix[enc_a[i - 1], enc_b[j - 1]])
                if H[i - 1, j - 1] < INF and abs(H[i, j] - (H[i - 1, j - 1] + sub_cost)) < 1e-9:
                    op = "M" if enc_a[i - 1] == enc_b[j - 1] else "X"
                    pairs.append((i - 1, j - 1, op))
                    i -= 1
                    j -= 1
                    continue

            if i > 0 and H[i, j] < INF and E[i, j] < INF and abs(H[i, j] - E[i, j]) < 1e-9:
                pairs.append((i - 1, None, "D"))
                i -= 1
                continue

            if j > 0 and H[i, j] < INF and F[i, j] < INF and abs(H[i, j] - F[i, j]) < 1e-9:
                pairs.append((None, j - 1, "I"))
                j -= 1
                continue

            # Fallback: try any direction
            if i > 0 and j > 0:
                op = "M" if enc_a[i - 1] == enc_b[j - 1] else "X"
                pairs.append((i - 1, j - 1, op))
                i -= 1
                j -= 1
            elif i > 0:
                pairs.append((i - 1, None, "D"))
                i -= 1
            else:
                pairs.append((None, j - 1, "I"))
                j -= 1

        pairs.reverse()
        return pairs, total_energy
