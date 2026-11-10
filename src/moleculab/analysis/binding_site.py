"""Parser-independent binding-site geometry.

This module works only with Cartesian coordinates in nanometres. It does not
know about PDB files, mmCIF records, residues, alternate locations, or force
fields. The caller is responsible for choosing a structure representation and
for deciding how to handle alternate locations, hydrogens, waters, and
cofactors.

The first operation is deliberately simple: for every protein atom, find the
distance to the nearest ligand atom. Selection then uses a distance cutoff.
Periodic boundary conditions are **not** applied; crystallographic symmetry and
simulation periodicity are separate modelling decisions.
"""

from __future__ import annotations

from collections.abc import Hashable, Sequence
from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class BindingSiteGeometry:
    """Geometric summary of a selected set of protein atoms.

    ``atom_indices`` are indices into the protein position array. ``center`` is
    the arithmetic mean of the selected atom positions. ``lower`` and ``upper``
    are the corners of the axis-aligned bounding box, and ``dimensions`` is
    their difference. All positions and dimensions are in nm.
    """

    atom_indices: np.ndarray
    center: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    dimensions: np.ndarray
    cutoff: float | None = None

    @property
    def n_atoms(self) -> int:
        """Number of selected protein atoms."""

        return int(self.atom_indices.size)


def _as_positions(name: str, positions: np.ndarray) -> np.ndarray:
    coordinates = np.asarray(positions, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(
            f"{name} must have shape (n_atoms, 3), got {coordinates.shape}."
        )
    if not np.all(np.isfinite(coordinates)):
        raise ValueError(f"{name} must contain only finite values.")
    return coordinates


def _as_cutoff(cutoff: float) -> float:
    value = float(cutoff)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("cutoff must be finite and positive.")
    return value


def _as_indices(name: str, values: Sequence[int]) -> np.ndarray:
    try:
        numeric = np.asarray(values, dtype=float).reshape(-1)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain integers.") from exc
    if not np.all(np.isfinite(numeric)) or not np.all(numeric == np.floor(numeric)):
        raise ValueError(f"{name} must contain integers.")
    return numeric.astype(int)


def minimum_distances(
    protein_positions: np.ndarray,
    ligand_positions: np.ndarray,
) -> np.ndarray:
    """Return each protein atom's distance to the nearest ligand atom.

    ``protein_positions`` has shape ``(N, 3)`` and ``ligand_positions`` has
    shape ``(M, 3)``, both in nm. The result has shape ``(N,)``.

    This is an explicit :math:`O(NM)` calculation. It does not build a
    neighbour list or apply periodic minimum images.
    """

    protein = _as_positions("protein_positions", protein_positions)
    ligand = _as_positions("ligand_positions", ligand_positions)
    if ligand.shape[0] == 0:
        raise ValueError("ligand_positions must contain at least one atom.")
    if protein.shape[0] == 0:
        return np.empty(0, dtype=float)

    distances = np.empty(protein.shape[0], dtype=float)
    for index, position in enumerate(protein):
        displacements = ligand - position
        distances[index] = np.min(np.linalg.norm(displacements, axis=1))
    return distances


def select_pocket_atoms(
    protein_positions: np.ndarray,
    ligand_positions: np.ndarray,
    cutoff: float,
    *,
    inclusive: bool = False,
) -> np.ndarray:
    """Return protein atom indices within *cutoff* of any ligand atom.

    The default selection is strict: an atom is selected when its minimum
    distance is less than ``cutoff``. Set ``inclusive=True`` to also select
    atoms exactly at the cutoff. The result is a sorted integer array.

    This function applies no alternate-location, hydrogen, water, cofactor, or
    periodic-boundary policy. Those decisions belong to the structure adapter.
    """

    distance_limit = _as_cutoff(cutoff)
    distances = minimum_distances(protein_positions, ligand_positions)
    mask = (
        distances <= distance_limit
        if inclusive
        else distances < distance_limit
    )
    return np.flatnonzero(mask)


def pocket_geometry(
    protein_positions: np.ndarray,
    atom_indices: Sequence[int],
    *,
    cutoff: float | None = None,
) -> BindingSiteGeometry:
    """Return the centre and bounding box of selected protein atoms.

    The centre is the arithmetic mean of the selected positions, not a
    ligand centre or a cavity centre. The bounding box is axis-aligned.
    """

    protein = _as_positions("protein_positions", protein_positions)
    indices = _as_indices("atom_indices", atom_indices)
    if indices.size == 0:
        raise ValueError("atom_indices must contain at least one atom.")
    if np.any(indices < 0) or np.any(indices >= protein.shape[0]):
        raise ValueError("atom_indices contains an out-of-range index.")
    if np.unique(indices).size != indices.size:
        raise ValueError("atom_indices must be unique.")

    cutoff_value = None if cutoff is None else _as_cutoff(cutoff)
    selected = protein[indices]
    center = np.mean(selected, axis=0)
    lower = np.min(selected, axis=0)
    upper = np.max(selected, axis=0)
    return BindingSiteGeometry(
        atom_indices=indices.copy(),
        center=center,
        lower=lower,
        upper=upper,
        dimensions=upper - lower,
        cutoff=cutoff_value,
    )


def binding_site_geometry(
    protein_positions: np.ndarray,
    ligand_positions: np.ndarray,
    cutoff: float,
    *,
    inclusive: bool = False,
) -> BindingSiteGeometry:
    """Select protein atoms around a ligand and return their geometry."""

    indices = select_pocket_atoms(
        protein_positions,
        ligand_positions,
        cutoff,
        inclusive=inclusive,
    )
    if indices.size == 0:
        raise ValueError("no protein atoms satisfy the binding-site cutoff.")
    return pocket_geometry(protein_positions, indices, cutoff=cutoff)


def group_indices_by_residue(
    atom_indices: Sequence[int],
    residue_ids: Sequence[Hashable],
) -> dict[Hashable, tuple[int, ...]]:
    """Group selected atom indices by residue identifier.

    ``residue_ids`` must have one label per atom index. The labels can be
    strings, tuples, or any other hashable structure chosen by the adapter.
    Residues appear in first-occurrence order, and each atom index appears once
    per residue group.
    """

    indices = list(atom_indices)
    labels = list(residue_ids)
    if len(indices) != len(labels):
        raise ValueError("atom_indices and residue_ids must have the same length.")
    validated = _as_indices("atom_indices", indices)
    if np.unique(validated).size != validated.size:
        raise ValueError("atom_indices must be unique.")

    groups: dict[Hashable, list[int]] = {}
    for index, label in zip(validated, labels, strict=True):
        groups.setdefault(label, []).append(int(index))
    return {label: tuple(values) for label, values in groups.items()}
