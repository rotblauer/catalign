"""Multi-scale quality evaluation for alignments."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from catalign.energy import EnergyModel


@dataclass
class BaseQuality:
    """Per-base alignment quality."""

    position_query: Optional[int]
    position_target: Optional[int]
    operation: str  # M, X, I, D
    energy: float = 0.0


@dataclass
class BlockQuality:
    """Quality for a contiguous aligned block."""

    start_query: int = 0
    end_query: int = 0
    start_target: int = 0
    end_target: int = 0
    length: int = 0
    matches: int = 0
    mismatches: int = 0
    identity: float = 0.0
    energy: float = 0.0


@dataclass
class RegionQuality:
    """Larger-region quality assessment."""

    start_query: int = 0
    end_query: int = 0
    start_target: int = 0
    end_target: int = 0
    num_blocks: int = 0
    total_aligned_bases: int = 0
    coverage_query: float = 0.0
    coverage_target: float = 0.0
    concordance: float = 0.0
    energy: float = 0.0


@dataclass
class AlignmentQuality:
    """Overall alignment quality summary containing all scales."""

    base_qualities: List[BaseQuality] = field(default_factory=list)
    block_qualities: List[BlockQuality] = field(default_factory=list)
    region_quality: Optional[RegionQuality] = None
    total_matches: int = 0
    total_mismatches: int = 0
    total_insertions: int = 0
    total_deletions: int = 0
    overall_identity: float = 0.0
    total_energy: float = 0.0
    quality_score: float = 0.0  # 0-100 composite score


@dataclass
class QualityReport:
    """Aggregates multiple alignment qualities for a genome-wide view."""

    alignment_qualities: List[AlignmentQuality] = field(default_factory=list)
    genome_coverage: float = 0.0
    mean_identity: float = 0.0
    mean_quality_score: float = 0.0
    total_energy: float = 0.0

    def add(self, aq: AlignmentQuality) -> None:
        self.alignment_qualities.append(aq)
        self._recompute()

    def _recompute(self) -> None:
        if not self.alignment_qualities:
            return
        identities = [a.overall_identity for a in self.alignment_qualities]
        scores = [a.quality_score for a in self.alignment_qualities]
        self.mean_identity = sum(identities) / len(identities)
        self.mean_quality_score = sum(scores) / len(scores)
        self.total_energy = sum(a.total_energy for a in self.alignment_qualities)


def evaluate_quality(
    alignment,
    query_seq: str,
    target_seq: str,
    energy_model: Optional[EnergyModel] = None,
) -> AlignmentQuality:
    """Evaluate alignment quality at multiple scales.

    *alignment* is an ``Alignment`` instance (from catalign.align).
    """
    em = energy_model or EnergyModel()
    pairs = alignment.aligned_pairs

    # --- Base-level ---
    base_quals: List[BaseQuality] = []
    matches = mismatches = insertions = deletions = 0
    total_energy = 0.0

    for qp, tp, op in pairs:
        if op == "M":
            e = em.match_energy
            matches += 1
        elif op == "X":
            q_base = query_seq[qp] if qp is not None and qp < len(query_seq) else "N"
            t_base = target_seq[tp] if tp is not None and tp < len(target_seq) else "N"
            e = em.compute_base_energy(q_base, t_base)
            mismatches += 1
        elif op == "I":
            e = em.gap_open_energy
            insertions += 1
        elif op == "D":
            e = em.gap_open_energy
            deletions += 1
        else:
            e = 0.0
        total_energy += e
        base_quals.append(BaseQuality(
            position_query=qp, position_target=tp, operation=op, energy=e
        ))

    # --- Block-level ---
    block_quals = _compute_blocks(pairs, em, query_seq, target_seq)

    # --- Region-level ---
    aligned_bases = matches + mismatches
    total_ops = matches + mismatches + insertions + deletions
    q_positions = {qp for qp, _, op in pairs if qp is not None}
    t_positions = {tp for _, tp, op in pairs if tp is not None}

    region = RegionQuality(
        start_query=alignment.query_start,
        end_query=alignment.query_end,
        start_target=alignment.target_start,
        end_target=alignment.target_end,
        num_blocks=len(block_quals),
        total_aligned_bases=aligned_bases,
        coverage_query=len(q_positions) / max(len(query_seq), 1),
        coverage_target=len(t_positions) / max(len(target_seq), 1),
        concordance=matches / max(total_ops, 1),
        energy=total_energy,
    )

    identity = matches / max(aligned_bases, 1)
    # Quality score: 0-100
    quality_score = identity * 100.0

    return AlignmentQuality(
        base_qualities=base_quals,
        block_qualities=block_quals,
        region_quality=region,
        total_matches=matches,
        total_mismatches=mismatches,
        total_insertions=insertions,
        total_deletions=deletions,
        overall_identity=identity,
        total_energy=total_energy,
        quality_score=quality_score,
    )


def _compute_blocks(pairs, em, query_seq, target_seq) -> List[BlockQuality]:
    """Segment aligned pairs into contiguous match/mismatch blocks."""
    blocks: List[BlockQuality] = []
    if not pairs:
        return blocks

    current: List[Tuple[Optional[int], Optional[int], str]] = []
    for pair in pairs:
        _, _, op = pair
        if op in ("M", "X"):
            current.append(pair)
        else:
            if current:
                blocks.append(_make_block(current, em, query_seq, target_seq))
                current = []
    if current:
        blocks.append(_make_block(current, em, query_seq, target_seq))

    return blocks


def _make_block(pairs, em, query_seq, target_seq) -> BlockQuality:
    q_positions = [p[0] for p in pairs if p[0] is not None]
    t_positions = [p[1] for p in pairs if p[1] is not None]
    m = sum(1 for _, _, op in pairs if op == "M")
    x = sum(1 for _, _, op in pairs if op == "X")
    energy = 0.0
    for qp, tp, op in pairs:
        if op == "M":
            energy += em.match_energy
        else:
            q_base = query_seq[qp] if qp is not None and qp < len(query_seq) else "N"
            t_base = target_seq[tp] if tp is not None and tp < len(target_seq) else "N"
            energy += em.compute_base_energy(q_base, t_base)

    return BlockQuality(
        start_query=min(q_positions) if q_positions else 0,
        end_query=max(q_positions) + 1 if q_positions else 0,
        start_target=min(t_positions) if t_positions else 0,
        end_target=max(t_positions) + 1 if t_positions else 0,
        length=len(pairs),
        matches=m,
        mismatches=x,
        identity=m / max(m + x, 1),
        energy=energy,
    )
