# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Pose-recovery validation for externally generated docking poses."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from .poses import LigandPose
from .rmsd import AtomMapping, symmetry_aware_rmsd


@dataclass(frozen=True)
class PoseRecovery:
    """RMSD-based pose-recovery result with score provenance retained."""

    rmsd: float
    threshold: float
    passed: bool
    best_mapping: AtomMapping
    mapping_index: int
    n_mappings: int
    pose_score: float | None
    score_unit: str | None
    score_source: str


def evaluate_pose_recovery(
    pose: LigandPose,
    reference: LigandPose,
    atom_mappings: Sequence[Sequence[tuple[int, int]]],
    threshold: float,
    *,
    require_complete: bool = True,
) -> PoseRecovery:
    """Compare a pose with a reference using explicit symmetry mappings.

    ``threshold`` is in nm. The docking score is retained as metadata and is
    never interpreted as a binding free energy.
    """

    limit = float(threshold)
    if not np.isfinite(limit) or limit <= 0.0:
        raise ValueError("threshold must be finite and positive in nm.")

    comparison = symmetry_aware_rmsd(
        pose.positions,
        reference.positions,
        atom_mappings,
        require_complete=require_complete,
    )
    return PoseRecovery(
        rmsd=comparison.rmsd,
        threshold=limit,
        passed=comparison.rmsd <= limit,
        best_mapping=comparison.best_mapping,
        mapping_index=comparison.mapping_index,
        n_mappings=comparison.n_mappings,
        pose_score=pose.score,
        score_unit=pose.score_unit,
        score_source=pose.score_source,
    )
