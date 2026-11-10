# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Explicit-mapping RMSD for ligand poses.

RMSD depends on atom mapping. This module never guesses which atoms correspond.
The caller supplies one-to-one index mappings, and symmetry-aware comparison is
performed by supplying several valid mappings.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from .poses import LigandPose

AtomMapping = tuple[tuple[int, int], ...]


@dataclass(frozen=True)
class SymmetryAwareRMSD:
    """Best RMSD over an explicit collection of symmetry-equivalent mappings."""

    rmsd: float
    best_mapping: AtomMapping
    mapping_index: int
    n_mappings: int


def validate_atom_mapping(
    atom_mapping: Sequence[tuple[int, int]],
    n_pose_atoms: int,
    n_reference_atoms: int,
    *,
    require_complete: bool = True,
) -> AtomMapping:
    """Validate and normalize a one-to-one pose/reference atom mapping.

    Each pair is ``(pose_index, reference_index)``. Indices must be unique on
    each side. With ``require_complete=True``, the mapping must cover every pose
    and reference atom; otherwise it represents an explicitly partial mapping.
    """

    pairs = tuple((int(pose_index), int(reference_index)) for pose_index, reference_index in atom_mapping)
    if not pairs:
        raise ValueError("atom_mapping must contain at least one pair.")

    pose_indices = [pose_index for pose_index, _ in pairs]
    reference_indices = [reference_index for _, reference_index in pairs]
    if len(set(pose_indices)) != len(pose_indices):
        raise ValueError("pose indices must be unique in atom_mapping.")
    if len(set(reference_indices)) != len(reference_indices):
        raise ValueError("reference indices must be unique in atom_mapping.")

    if any(index < 0 or index >= n_pose_atoms for index in pose_indices):
        raise ValueError("atom_mapping contains an out-of-range pose index.")
    if any(index < 0 or index >= n_reference_atoms for index in reference_indices):
        raise ValueError("atom_mapping contains an out-of-range reference index.")

    if require_complete and (
        len(pairs) != n_pose_atoms or len(pairs) != n_reference_atoms
    ):
        raise ValueError(
            "atom_mapping must cover every pose and reference atom; "
            "set require_complete=False for an explicit partial mapping."
        )
    return pairs


def _as_positions(name: str, positions: np.ndarray) -> np.ndarray:
    coordinates = np.asarray(positions, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(f"{name} must have shape (n_atoms, 3).")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(f"{name} must contain only finite values.")
    return coordinates


def pose_rmsd(
    pose_positions: np.ndarray,
    reference_positions: np.ndarray,
    atom_mapping: Sequence[tuple[int, int]],
    *,
    require_complete: bool = True,
) -> float:
    """Return the RMSD in nm for an explicit atom mapping.

    The two coordinate sets must already be in the same frame. Receptor
    alignment is a separate operation and is not performed here.
    """

    pose = _as_positions("pose_positions", pose_positions)
    reference = _as_positions("reference_positions", reference_positions)
    mapping = validate_atom_mapping(
        atom_mapping,
        pose.shape[0],
        reference.shape[0],
        require_complete=require_complete,
    )

    pose_indices = np.array([pair[0] for pair in mapping], dtype=int)
    reference_indices = np.array([pair[1] for pair in mapping], dtype=int)
    displacements = pose[pose_indices] - reference[reference_indices]
    return float(np.sqrt(np.mean(np.sum(displacements**2, axis=1))))


def symmetry_aware_rmsd(
    pose_positions: np.ndarray,
    reference_positions: np.ndarray,
    atom_mappings: Sequence[Sequence[tuple[int, int]]],
    *,
    require_complete: bool = True,
) -> SymmetryAwareRMSD:
    """Return the minimum RMSD over explicitly supplied symmetry mappings.

    The caller supplies the symmetry-equivalent mappings. This module does not
    perform graph automorphism perception.
    """

    mappings = tuple(atom_mappings)
    if not mappings:
        raise ValueError("atom_mappings must contain at least one mapping.")

    best_rmsd = np.inf
    best_mapping: AtomMapping | None = None
    best_index = -1
    for index, mapping in enumerate(mappings):
        current_rmsd = pose_rmsd(
            pose_positions,
            reference_positions,
            mapping,
            require_complete=require_complete,
        )
        if current_rmsd < best_rmsd:
            best_rmsd = current_rmsd
            best_mapping = validate_atom_mapping(
                mapping,
                np.asarray(pose_positions).shape[0],
                np.asarray(reference_positions).shape[0],
                require_complete=require_complete,
            )
            best_index = index

    if best_mapping is None:
        raise ValueError("atom_mappings must contain at least one mapping.")
    return SymmetryAwareRMSD(
        rmsd=float(best_rmsd),
        best_mapping=best_mapping,
        mapping_index=best_index,
        n_mappings=len(mappings),
    )


def atom_id_mapping(pose: LigandPose, reference: LigandPose) -> AtomMapping:
    """Return a mapping based on identical stable atom identifiers.

    The pose and reference must contain exactly the same set of atom IDs.
    Missing or extra atoms are an explicit error, not something to skip.
    """

    pose_ids = set(pose.atom_ids)
    reference_ids = set(reference.atom_ids)
    if pose_ids != reference_ids:
        missing = sorted(pose_ids - reference_ids)
        extra = sorted(reference_ids - pose_ids)
        raise ValueError(
            "pose and reference atom IDs do not match; "
            f"missing from reference: {missing}; extra in reference: {extra}."
        )

    reference_index = {
        atom_id: index for index, atom_id in enumerate(reference.atom_ids)
    }
    return tuple(
        (pose_index, reference_index[atom_id])
        for pose_index, atom_id in enumerate(pose.atom_ids)
    )
