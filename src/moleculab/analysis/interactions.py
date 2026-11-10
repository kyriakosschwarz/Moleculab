# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Geometry-based protein-ligand interaction detection.

The first implementation covers hydrogen bonds. It consumes prepared chemical
features, not RDKit or Gemmi objects, and works only with Cartesian positions in
nanometres.

Chemical perception is deliberately outside this module. The caller must decide
which atoms are donors, acceptors, charged, aromatic, or hydrophobic, and must
supply hydrogen positions for donor atoms when angular geometry is required.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class InteractionAtom:
    """A structure atom with the chemical features needed for interactions.

    ``position`` is in nm. ``hydrogen_positions`` contains the positions of
    hydrogens attached to this atom when it acts as a donor. A donor without
    hydrogen positions is retained in the representation but cannot pass the
    angular hydrogen-bond criterion.

    The flags are supplied by a chemistry adapter. This class does not infer
    donor/acceptor status from the element.
    """

    position: np.ndarray
    element: str = ""
    atom_type: str | None = None
    residue_id: str | None = None
    name: str = ""
    is_donor: bool = False
    is_acceptor: bool = False
    formal_charge: int = 0
    is_aromatic: bool = False
    is_hydrophobic: bool = False
    hydrogen_positions: tuple[np.ndarray, ...] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        position = np.asarray(self.position, dtype=float)
        if position.shape != (3,) or not np.all(np.isfinite(position)):
            raise ValueError("position must be a finite three-vector in nm.")
        object.__setattr__(self, "position", position)

        hydrogens = []
        for hydrogen in self.hydrogen_positions:
            coordinates = np.asarray(hydrogen, dtype=float)
            if coordinates.shape != (3,) or not np.all(np.isfinite(coordinates)):
                raise ValueError(
                    "hydrogen_positions must contain finite three-vectors in nm."
                )
            hydrogens.append(coordinates)
        object.__setattr__(self, "hydrogen_positions", tuple(hydrogens))


@dataclass(frozen=True)
class HydrogenBondCriteria:
    """Operational thresholds for hydrogen-bond detection.

    The defaults are conventions, not physical laws. ``max_distance`` is the
    donor-acceptor distance in nm. ``min_angle_degrees`` is the D-H...A angle.
    Both comparisons are inclusive.
    """

    max_distance: float = 0.35
    min_angle_degrees: float = 120.0

    def __post_init__(self) -> None:
        distance = float(self.max_distance)
        if not np.isfinite(distance) or distance <= 0.0:
            raise ValueError("max_distance must be finite and positive.")
        angle = float(self.min_angle_degrees)
        if not np.isfinite(angle) or not 0.0 < angle <= 180.0:
            raise ValueError("min_angle_degrees must be in (0, 180].")
        object.__setattr__(self, "max_distance", distance)
        object.__setattr__(self, "min_angle_degrees", angle)


@dataclass(frozen=True)
class HydrogenBond:
    """One geometrically plausible hydrogen bond.

    ``distance`` is the donor-acceptor distance in nm. ``angle_degrees`` is the
    D-H...A angle at the hydrogen. ``donor_is_protein`` records the direction of
    the interaction. The result is a geometric observation, not an energy.
    """

    protein_index: int
    ligand_index: int
    donor_index: int
    acceptor_index: int
    donor_is_protein: bool
    protein_residue: str | None
    ligand_residue: str | None
    protein_label: str
    ligand_label: str
    donor_label: str
    acceptor_label: str
    distance: float
    angle_degrees: float
    hydrogen_position: np.ndarray
    hydrogen_index: int
    criteria: HydrogenBondCriteria

    @property
    def interaction_type(self) -> str:
        """Stable interaction label for fingerprints and reports."""

        return "hydrogen_bond"


def _as_atoms(name: str, atoms: Sequence[InteractionAtom]) -> tuple[InteractionAtom, ...]:
    values = tuple(atoms)
    for atom in values:
        if not isinstance(atom, InteractionAtom):
            raise TypeError(f"{name} must contain InteractionAtom objects.")
    return values


def _angle_degrees(first: np.ndarray, vertex: np.ndarray, second: np.ndarray) -> float:
    left = first - vertex
    right = second - vertex
    left_norm = float(np.linalg.norm(left))
    right_norm = float(np.linalg.norm(right))
    if left_norm == 0.0 or right_norm == 0.0:
        return 0.0
    cosine = float(np.dot(left, right) / (left_norm * right_norm))
    return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def _atom_label(atom: InteractionAtom, index: int) -> str:
    return atom.name or str(index)


def detect_hydrogen_bonds(
    protein_atoms: Sequence[InteractionAtom],
    ligand_atoms: Sequence[InteractionAtom],
    criteria: HydrogenBondCriteria | None = None,
) -> tuple[HydrogenBond, ...]:
    """Return plausible hydrogen bonds between protein and ligand atoms.

    Both directions are considered: protein donor to ligand acceptor, and
    ligand donor to protein acceptor. A donor without hydrogen positions is
    skipped because its angular geometry cannot be evaluated.

    For a given protein-ligand atom pair, at most one record is returned. If
    both role assignments are possible, the geometry with the larger D-H...A
    angle wins; ties go to the shorter donor-acceptor distance.
    """

    protein = _as_atoms("protein_atoms", protein_atoms)
    ligand = _as_atoms("ligand_atoms", ligand_atoms)
    if criteria is None:
        criteria = HydrogenBondCriteria()
    if not isinstance(criteria, HydrogenBondCriteria):
        raise TypeError("criteria must be a HydrogenBondCriteria instance.")

    records: dict[tuple[int, int], HydrogenBond] = {}

    def consider(
        donor: InteractionAtom,
        acceptor: InteractionAtom,
        donor_is_protein: bool,
        protein_index: int,
        ligand_index: int,
    ) -> None:
        distance = float(np.linalg.norm(donor.position - acceptor.position))
        if distance > criteria.max_distance:
            return
        if not donor.hydrogen_positions:
            return

        best_angle = -1.0
        best_hydrogen = None
        best_hydrogen_index = -1
        for hydrogen_index, hydrogen in enumerate(donor.hydrogen_positions):
            angle = _angle_degrees(donor.position, hydrogen, acceptor.position)
            if angle > best_angle:
                best_angle = angle
                best_hydrogen = hydrogen
                best_hydrogen_index = hydrogen_index
        if best_hydrogen is None or best_angle < criteria.min_angle_degrees:
            return

        if donor_is_protein:
            protein_atom, ligand_atom = donor, acceptor
            donor_index, acceptor_index = protein_index, ligand_index
        else:
            protein_atom, ligand_atom = acceptor, donor
            donor_index, acceptor_index = ligand_index, protein_index

        record = HydrogenBond(
            protein_index=protein_index,
            ligand_index=ligand_index,
            donor_index=donor_index,
            acceptor_index=acceptor_index,
            donor_is_protein=donor_is_protein,
            protein_residue=protein_atom.residue_id,
            ligand_residue=ligand_atom.residue_id,
            protein_label=_atom_label(protein_atom, protein_index),
            ligand_label=_atom_label(ligand_atom, ligand_index),
            donor_label=_atom_label(donor, donor_index),
            acceptor_label=_atom_label(acceptor, acceptor_index),
            distance=distance,
            angle_degrees=best_angle,
            hydrogen_position=best_hydrogen.copy(),
            hydrogen_index=best_hydrogen_index,
            criteria=criteria,
        )

        key = (protein_index, ligand_index)
        previous = records.get(key)
        if previous is None or (
            best_angle,
            -distance,
        ) > (
            previous.angle_degrees,
            -previous.distance,
        ):
            records[key] = record

    for protein_index, protein_atom in enumerate(protein):
        for ligand_index, ligand_atom in enumerate(ligand):
            if protein_atom.is_donor and ligand_atom.is_acceptor:
                consider(
                    protein_atom,
                    ligand_atom,
                    True,
                    protein_index,
                    ligand_index,
                )
            if ligand_atom.is_donor and protein_atom.is_acceptor:
                consider(
                    ligand_atom,
                    protein_atom,
                    False,
                    protein_index,
                    ligand_index,
                )

    return tuple(
        sorted(
            records.values(),
            key=lambda record: (
                record.protein_residue or "",
                record.ligand_residue or "",
                record.protein_label,
                record.ligand_label,
                record.distance,
                -record.angle_degrees,
            ),
        )
    )
