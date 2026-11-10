# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Compare ligand binding modes by geometry and interaction patterns.

Geometry and interaction patterns are reported separately. The module does not
combine them into a single score, because a shared pose and a shared interaction
pattern answer different questions.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from ..fingerprints import InteractionFeature, InteractionFingerprint
from .poses import LigandPose
from .rmsd import AtomMapping, symmetry_aware_rmsd

PatternKey = InteractionFeature | tuple[str, str]


@dataclass(frozen=True)
class PoseGeometryComparison:
    """RMSD and centroid displacement between two mapped poses."""

    rmsd: float
    centroid_distance: float
    centroid_displacement: np.ndarray
    best_mapping: AtomMapping
    mapping_index: int
    n_mappings: int


@dataclass(frozen=True)
class InteractionPatternComparison:
    """Shared and unique interaction features or residue-level patterns."""

    level: str
    shared: tuple[PatternKey, ...]
    only_in_a: tuple[PatternKey, ...]
    only_in_b: tuple[PatternKey, ...]
    jaccard: float


@dataclass(frozen=True)
class BindingModeComparison:
    """Geometry and interaction comparison, deliberately reported separately."""

    geometry: PoseGeometryComparison
    interactions: InteractionPatternComparison


def compare_pose_geometry(
    pose_a: LigandPose,
    pose_b: LigandPose,
    atom_mappings: Sequence[Sequence[tuple[int, int]]],
    *,
    require_complete: bool = True,
) -> PoseGeometryComparison:
    """Compare two poses with explicit symmetry mappings.

    ``atom_mappings`` maps pose A atom indices to pose B atom indices. The poses
    must already be in a common coordinate frame; receptor alignment is a
    separate operation.
    """

    if not isinstance(pose_a, LigandPose) or not isinstance(pose_b, LigandPose):
        raise TypeError("pose_a and pose_b must be LigandPose objects.")

    comparison = symmetry_aware_rmsd(
        pose_a.positions,
        pose_b.positions,
        atom_mappings,
        require_complete=require_complete,
    )
    pose_indices = [pair[0] for pair in comparison.best_mapping]
    reference_indices = [pair[1] for pair in comparison.best_mapping]
    centroid_a = np.mean(pose_a.positions[pose_indices], axis=0)
    centroid_b = np.mean(pose_b.positions[reference_indices], axis=0)
    displacement = centroid_a - centroid_b
    return PoseGeometryComparison(
        rmsd=comparison.rmsd,
        centroid_distance=float(np.linalg.norm(displacement)),
        centroid_displacement=displacement,
        best_mapping=comparison.best_mapping,
        mapping_index=comparison.mapping_index,
        n_mappings=comparison.n_mappings,
    )


def compare_interaction_patterns(
    fingerprint_a: InteractionFingerprint,
    fingerprint_b: InteractionFingerprint,
    *,
    include_ligand_atom: bool = True,
) -> InteractionPatternComparison:
    """Compare two interaction fingerprints.

    With ``include_ligand_atom=True``, exact feature identities are compared,
    including ligand atom and ligand residue. With ``include_ligand_atom=False``,
    features are projected to ``(interaction_type, protein_residue)`` so that
    different ligands can be compared at the residue-interaction level.
    """

    if not isinstance(fingerprint_a, InteractionFingerprint):
        raise TypeError("fingerprint_a must be an InteractionFingerprint.")
    if not isinstance(fingerprint_b, InteractionFingerprint):
        raise TypeError("fingerprint_b must be an InteractionFingerprint.")

    if include_ligand_atom:
        set_a = set(fingerprint_a.present)
        set_b = set(fingerprint_b.present)
        shared = tuple(sorted(set_a & set_b, key=lambda item: item.key))
        only_a = tuple(sorted(set_a - set_b, key=lambda item: item.key))
        only_b = tuple(sorted(set_b - set_a, key=lambda item: item.key))
        level = "feature"
    else:
        set_a = {
            (feature.interaction_type, feature.protein_residue)
            for feature in fingerprint_a.present
        }
        set_b = {
            (feature.interaction_type, feature.protein_residue)
            for feature in fingerprint_b.present
        }
        shared = tuple(sorted(set_a & set_b))
        only_a = tuple(sorted(set_a - set_b))
        only_b = tuple(sorted(set_b - set_a))
        level = "residue_pattern"

    union = set_a | set_b
    jaccard = 1.0 if not union else len(shared) / len(union)
    return InteractionPatternComparison(
        level=level,
        shared=shared,
        only_in_a=only_a,
        only_in_b=only_b,
        jaccard=float(jaccard),
    )


def compare_binding_modes(
    pose_a: LigandPose,
    pose_b: LigandPose,
    atom_mappings: Sequence[Sequence[tuple[int, int]]],
    fingerprint_a: InteractionFingerprint,
    fingerprint_b: InteractionFingerprint,
    *,
    include_ligand_atom: bool = True,
    require_complete: bool = True,
) -> BindingModeComparison:
    """Compare pose geometry and interaction patterns without a composite score."""

    return BindingModeComparison(
        geometry=compare_pose_geometry(
            pose_a,
            pose_b,
            atom_mappings,
            require_complete=require_complete,
        ),
        interactions=compare_interaction_patterns(
            fingerprint_a,
            fingerprint_b,
            include_ligand_atom=include_ligand_atom,
        ),
    )
