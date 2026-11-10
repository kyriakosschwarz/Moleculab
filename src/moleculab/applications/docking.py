"""Minimal rigid-ligand docking.

Docking asks which position and orientation of a ligand fits best into a
receptor binding site.  This module keeps the problem deliberately small:

* the ligand is **rigid** — its internal geometry never changes;
* only receptor-ligand nonbonded interactions are scored;
* poses are sampled by random translations and rotations;
* the score is the existing typed Lennard-Jones plus reaction-field Coulomb
  energy.

The result is a ranked list of poses.  This is not a production docking
program: it has no flexible receptor, no explicit desolvation term, no
empirical scoring function, and no local optimization.  It is a transparent
way to see how the search and scoring pieces fit together.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from ..composition.system import System
from ..foundation.particles import Particle
from ..foundation.validation import as_nonnegative_float, as_vector3


def _as_rotation_matrix(rotation: np.ndarray) -> np.ndarray:
    matrix = np.asarray(rotation, dtype=float)
    if matrix.shape != (3, 3):
        raise ValueError(f"rotation must have shape (3, 3), got {matrix.shape}.")
    if not np.all(np.isfinite(matrix)):
        raise ValueError("rotation must contain only finite values.")
    if not np.allclose(matrix.T @ matrix, np.eye(3), atol=1e-8):
        raise ValueError("rotation must be an orthonormal matrix.")
    if not np.isclose(np.linalg.det(matrix), 1.0, atol=1e-8):
        raise ValueError("rotation must have determinant +1.")
    return matrix


@dataclass(frozen=True)
class RigidLigand:
    """A rigid collection of particles.

    The constructor computes the center of mass and stores the particle
    positions relative to it.  A pose is defined by translating the center of
    mass and rotating the local coordinates.

    Parameters
    ----------
    particles:
        At least one particle.  Only mass, position, type, and charge are
        used; velocities are ignored.
    """

    particles: Sequence[Particle]

    center_of_mass: np.ndarray = field(init=False, repr=False)
    local_positions: np.ndarray = field(init=False, repr=False)
    masses: np.ndarray = field(init=False, repr=False)
    atom_types: tuple[str, ...] = field(init=False, repr=False)
    charges: np.ndarray = field(init=False, repr=False)

    def __post_init__(self) -> None:
        atoms = tuple(self.particles)
        if not atoms:
            raise ValueError("RigidLigand must contain at least one particle.")
        for index, particle in enumerate(atoms):
            if not isinstance(particle, Particle):
                raise TypeError(f"particles[{index}] must be a Particle.")
        object.__setattr__(self, "particles", atoms)

        masses = np.array([particle.mass for particle in atoms], dtype=float)
        positions = np.array([particle.position for particle in atoms], dtype=float)
        center = np.sum(masses[:, np.newaxis] * positions, axis=0) / np.sum(masses)
        object.__setattr__(self, "center_of_mass", center)
        object.__setattr__(self, "local_positions", positions - center)
        object.__setattr__(self, "masses", masses)
        object.__setattr__(
            self,
            "atom_types",
            tuple(particle.atom_type for particle in atoms),
        )
        object.__setattr__(
            self,
            "charges",
            np.array([particle.charge for particle in atoms], dtype=float),
        )

    @property
    def n_atoms(self) -> int:
        return len(self.particles)

    def transformed_positions(
        self,
        translation: np.ndarray,
        rotation: np.ndarray,
    ) -> np.ndarray:
        """Return local coordinates rotated and translated into the receptor frame."""

        center = as_vector3(translation, "translation")
        matrix = _as_rotation_matrix(rotation)
        return self.local_positions @ matrix.T + center

    def transformed_particles(
        self,
        translation: np.ndarray,
        rotation: np.ndarray,
    ) -> tuple[Particle, ...]:
        """Return particles at the requested rigid-body pose."""

        positions = self.transformed_positions(translation, rotation)
        return tuple(
            Particle(
                mass=particle.mass,
                position=positions[index],
                atom_type=particle.atom_type,
                charge=particle.charge,
            )
            for index, particle in enumerate(self.particles)
        )


@dataclass(frozen=True)
class DockingPose:
    """One scored ligand pose.

    ``translation`` is the ligand center-of-mass position, ``rotation`` is the
    orthogonal matrix applied to the local coordinates, ``score`` is the
    receptor-ligand interaction energy in kJ/mol, and ``positions`` contains
    the transformed ligand atom positions.
    """

    translation: np.ndarray
    rotation: np.ndarray
    score: float
    positions: np.ndarray


@dataclass(frozen=True)
class DockingResult:
    """Ranked docking poses, best first."""

    poses: tuple[DockingPose, ...]

    @property
    def best(self) -> DockingPose:
        return self.poses[0]

    @property
    def scores(self) -> np.ndarray:
        return np.array([pose.score for pose in self.poses], dtype=float)

    def top(self, n: int = 10) -> tuple[DockingPose, ...]:
        """Return the best *n* poses."""

        return self.poses[: int(n)]


def score_pose(
    receptor: System,
    ligand: RigidLigand,
    translation: np.ndarray,
    rotation: np.ndarray,
) -> DockingPose:
    """Score one rigid ligand pose against a receptor.

    Only receptor-ligand pairs are included.  Intramolecular ligand and
    receptor energies are constant for a rigid ligand and are therefore
    omitted.  Periodic minimum images are used when the receptor has a box.
    """

    if not isinstance(receptor, System):
        raise TypeError("receptor must be a System.")
    if not isinstance(ligand, RigidLigand):
        raise TypeError("ligand must be a RigidLigand.")
    if receptor.nonbonded_model is None:
        raise ValueError(
            "receptor.nonbonded_model is required for typed docking scores."
        )
    receptor.nonbonded_model.validate_atom_types(
        [particle.atom_type for particle in receptor.particles]
        + list(ligand.atom_types)
    )

    positions = ligand.transformed_positions(translation, rotation)
    model = receptor.nonbonded_model
    total = 0.0

    for receptor_particle in receptor.particles:
        for ligand_index, ligand_position in enumerate(positions):
            displacement = ligand_position - receptor_particle.position
            if receptor.box is not None:
                displacement = receptor.box.minimum_image(displacement)
            total += model.energy(
                receptor_particle.atom_type,
                ligand.atom_types[ligand_index],
                receptor_particle.charge,
                float(ligand.charges[ligand_index]),
                displacement,
            )

    return DockingPose(
        translation=as_vector3(translation, "translation").copy(),
        rotation=_as_rotation_matrix(rotation).copy(),
        score=float(total),
        positions=positions.copy(),
    )


def _random_translations(
    rng: np.random.Generator,
    center: np.ndarray,
    radius: float,
    n_samples: int,
) -> np.ndarray:
    if radius == 0.0:
        return np.repeat(center[np.newaxis, :], n_samples, axis=0)

    directions = rng.normal(size=(n_samples, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)
    # Uniform in a ball: radius^(1/3) gives the correct volume weighting.
    radii = radius * rng.random(n_samples) ** (1.0 / 3.0)
    return center + directions * radii[:, np.newaxis]


def _random_rotations(
    rng: np.random.Generator,
    n_samples: int,
) -> np.ndarray:
    rotations = []
    for _ in range(n_samples):
        matrix = rng.normal(size=(3, 3))
        q, r = np.linalg.qr(matrix)
        q = q @ np.diag(np.sign(np.diag(r)))
        if np.linalg.det(q) < 0.0:
            q[:, 0] *= -1.0
        rotations.append(q)
    return np.asarray(rotations)


def dock_rigid_ligand(
    receptor: System,
    ligand: RigidLigand,
    *,
    center: np.ndarray | None = None,
    translation_radius: float = 0.4,
    n_samples: int = 500,
    seed: int | None = None,
    max_results: int | None = None,
) -> DockingResult:
    """Sample rigid translations/rotations and return ranked poses.

    Parameters
    ----------
    receptor:
        A :class:`~moleculab.composition.system.System` with a typed ``nonbonded_model``.
    ligand:
        Rigid ligand whose center of mass is placed at each sampled
        translation.
    center:
        Center of the translation search sphere.  Defaults to the geometric
        center of the receptor.
    translation_radius:
        Radius of the translation search sphere in nm.
    n_samples:
        Number of random poses to score.
    seed:
        Random seed for reproducible sampling.
    max_results:
        If given, only the best ``max_results`` poses are returned.

    Returns
    -------
    DockingResult
        Poses sorted from lowest (best) to highest score.
    """

    if not isinstance(receptor, System):
        raise TypeError("receptor must be a System.")
    if not isinstance(ligand, RigidLigand):
        raise TypeError("ligand must be a RigidLigand.")
    if receptor.nonbonded_model is None:
        raise ValueError(
            "receptor.nonbonded_model is required for typed docking scores."
        )
    if isinstance(n_samples, bool) or not isinstance(
        n_samples, (int, np.integer)
    ):
        raise ValueError("n_samples must be an integer.")
    if int(n_samples) < 1:
        raise ValueError("n_samples must be at least 1.")
    n_samples = int(n_samples)

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
    translations = _random_translations(rng, center_value, radius, n_samples)
    rotations = _random_rotations(rng, n_samples)

    poses = [
        score_pose(receptor, ligand, translations[index], rotations[index])
        for index in range(n_samples)
    ]
    poses.sort(key=lambda pose: pose.score)
    if max_results is not None:
        poses = poses[:max_results]
    return DockingResult(poses=tuple(poses))
