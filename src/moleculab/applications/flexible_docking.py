# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Flexible-receptor docking with ligand rigid-body sampling.

The receptor can contain one or more explicit side-chain groups.  Each group is
moved as a rigid block by a sampled displacement vector.  The calibrated score
adds a harmonic strain penalty to the receptor-ligand interaction energy:

.. math::

    S = E_{\\mathrm{interaction}}
    + \\sum_{\\mathrm{chains}} \\frac{1}{2} k_s |\\Delta \\mathbf r_s|^2.

``k_s`` is a user-provided stiffness in kJ/(mol nm^2).  This is a simple
physical calibration term: it prevents the search from moving side chains
arbitrarily far just to improve the interaction energy.  It is not a trained
empirical scoring function.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from ..composition.system import System
from ..foundation.validation import as_nonnegative_float, as_positive_float, as_vector3
from .docking import (
    RigidLigand,
    _random_rotations,
    _random_translations,
    score_pose,
)


@dataclass(frozen=True)
class SideChain:
    """A group of receptor atoms moved together during docking.

    ``atom_indices`` must refer to receptor particles.  ``max_displacement`` is
    the radius of the translation search sphere in nm, and ``stiffness`` is the
    harmonic strain constant in kJ/(mol nm^2).
    """

    atom_indices: tuple[int, ...]
    max_displacement: float = 0.2
    stiffness: float = 10.0

    def __post_init__(self) -> None:
        indices = tuple(int(index) for index in self.atom_indices)
        if not indices:
            raise ValueError("A side chain must contain at least one atom.")
        if len(set(indices)) != len(indices):
            raise ValueError("Side-chain atom indices must be unique.")
        if any(index < 0 for index in indices):
            raise ValueError("Side-chain atom indices must be non-negative.")
        object.__setattr__(self, "atom_indices", indices)
        object.__setattr__(
            self,
            "max_displacement",
            as_nonnegative_float(self.max_displacement, "max_displacement"),
        )
        object.__setattr__(
            self,
            "stiffness",
            as_positive_float(self.stiffness, "stiffness"),
        )


@dataclass(frozen=True)
class FlexibleDockingPose:
    """A ligand pose plus a receptor side-chain conformation."""

    translation: np.ndarray
    rotation: np.ndarray
    ligand_positions: np.ndarray
    receptor_displacements: tuple[np.ndarray, ...]
    interaction_score: float
    strain_score: float
    score: float


@dataclass(frozen=True)
class FlexibleDockingResult:
    """Ranked flexible-receptor poses, best first."""

    poses: tuple[FlexibleDockingPose, ...]

    @property
    def best(self) -> FlexibleDockingPose:
        return self.poses[0]

    @property
    def scores(self) -> np.ndarray:
        return np.array([pose.score for pose in self.poses], dtype=float)

    def top(self, n: int = 10) -> tuple[FlexibleDockingPose, ...]:
        return self.poses[: int(n)]


def _validate_side_chains(receptor: System, side_chains: Sequence[SideChain]) -> tuple[SideChain, ...]:
    chains = tuple(side_chains)
    for index, chain in enumerate(chains):
        if not isinstance(chain, SideChain):
            raise TypeError(f"side_chains[{index}] must be a SideChain.")
        if any(atom >= receptor.n_particles for atom in chain.atom_indices):
            raise ValueError(
                f"side_chains[{index}] refers to an atom outside the receptor."
            )
    return chains


def make_receptor_conformer(
    receptor: System,
    side_chains: Sequence[SideChain],
    displacements: Sequence[np.ndarray],
) -> System:
    """Return a receptor copy with side-chain groups translated."""

    chains = _validate_side_chains(receptor, side_chains)
    if len(displacements) != len(chains):
        raise ValueError("displacements must match side_chains in length.")
    positions = receptor.positions.copy()
    vectors = [as_vector3(displacement, "displacement") for displacement in displacements]
    for chain, vector in zip(chains, vectors):
        for atom_index in chain.atom_indices:
            positions[atom_index] += vector

    return receptor.with_state(positions, receptor.velocities)


def receptor_strain_score(
    side_chains: Sequence[SideChain],
    displacements: Sequence[np.ndarray],
) -> float:
    """Return the harmonic side-chain strain energy in kJ/mol."""

    chains = tuple(side_chains)
    if len(displacements) != len(chains):
        raise ValueError("displacements must match side_chains in length.")
    total = 0.0
    for chain, displacement in zip(chains, displacements):
        vector = as_vector3(displacement, "displacement")
        total += 0.5 * chain.stiffness * float(np.dot(vector, vector))
    return float(total)


def score_flexible_pose(
    receptor: System,
    ligand: RigidLigand,
    side_chains: Sequence[SideChain],
    displacements: Sequence[np.ndarray],
    translation: np.ndarray,
    rotation: np.ndarray,
) -> FlexibleDockingPose:
    """Score one ligand pose in one explicit receptor conformation."""

    conformer = make_receptor_conformer(receptor, side_chains, displacements)
    pose = score_pose(conformer, ligand, translation, rotation)
    strain = receptor_strain_score(side_chains, displacements)
    return FlexibleDockingPose(
        translation=pose.translation,
        rotation=pose.rotation,
        ligand_positions=pose.positions,
        receptor_displacements=tuple(
            as_vector3(displacement, "displacement").copy()
            for displacement in displacements
        ),
        interaction_score=pose.score,
        strain_score=strain,
        score=pose.score + strain,
    )


def dock_flexible_receptor(
    receptor: System,
    ligand: RigidLigand,
    side_chains: Sequence[SideChain],
    *,
    center: np.ndarray | None = None,
    translation_radius: float = 0.4,
    n_ligand_samples: int = 200,
    n_receptor_samples: int = 20,
    seed: int | None = None,
    max_results: int | None = None,
) -> FlexibleDockingResult:
    """Sample ligand poses and explicit receptor side-chain displacements."""

    if not isinstance(receptor, System):
        raise TypeError("receptor must be a System.")
    if not isinstance(ligand, RigidLigand):
        raise TypeError("ligand must be a RigidLigand.")
    if receptor.nonbonded_model is None:
        raise ValueError(
            "receptor.nonbonded_model is required for typed docking scores."
        )
    chains = _validate_side_chains(receptor, side_chains)
    for name, value in (
        ("n_ligand_samples", n_ligand_samples),
        ("n_receptor_samples", n_receptor_samples),
    ):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
            raise ValueError(f"{name} must be an integer.")
        if int(value) < 1:
            raise ValueError(f"{name} must be at least 1.")
    n_ligand_samples = int(n_ligand_samples)
    n_receptor_samples = int(n_receptor_samples)
    radius = as_nonnegative_float(translation_radius, "translation_radius")

    if center is None:
        center_value = np.mean(receptor.positions, axis=0)
    else:
        center_value = as_vector3(center, "center")

    if max_results is not None:
        if isinstance(max_results, bool) or not isinstance(
            max_results, (int, np.integer)
        ):
            raise ValueError("max_results must be an integer or None.")
        if int(max_results) < 1:
            raise ValueError("max_results must be at least 1.")
        max_results = int(max_results)

    rng = np.random.default_rng(seed)
    poses: list[FlexibleDockingPose] = []
    n_chains = len(chains)

    for _ in range(n_receptor_samples):
        displacements = []
        for chain in chains:
            if chain.max_displacement == 0.0 or n_chains == 0:
                displacements.append(np.zeros(3))
            else:
                direction = rng.normal(size=3)
                direction /= np.linalg.norm(direction)
                distance = chain.max_displacement * rng.random() ** (1.0 / 3.0)
                displacements.append(direction * distance)
        conformer = make_receptor_conformer(receptor, chains, displacements)
        strain = receptor_strain_score(chains, displacements)
        translations = _random_translations(
            rng, center_value, radius, n_ligand_samples
        )
        rotations = _random_rotations(rng, n_ligand_samples)
        for index in range(n_ligand_samples):
            pose = score_pose(
                conformer,
                ligand,
                translations[index],
                rotations[index],
            )
            poses.append(
                FlexibleDockingPose(
                    translation=pose.translation,
                    rotation=pose.rotation,
                    ligand_positions=pose.positions,
                    receptor_displacements=tuple(
                        displacement.copy() for displacement in displacements
                    ),
                    interaction_score=pose.score,
                    strain_score=strain,
                    score=pose.score + strain,
                )
            )

    poses.sort(key=lambda pose: pose.score)
    if max_results is not None:
        poses = poses[:max_results]
    return FlexibleDockingResult(poses=tuple(poses))
