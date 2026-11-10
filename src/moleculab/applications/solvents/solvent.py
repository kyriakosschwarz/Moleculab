# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Explicit solvent boxes for molecular simulations.

This module provides a deliberately simple explicit-solvent path.  Each
solvent molecule is represented by one particle with a Lennard-Jones
interaction (and an optional point charge).  That is a **coarse-grained**
solvent model: it captures excluded volume, dispersion, density, and periodic
boundaries, but not the directional hydrogen-bond network of atomistic water.

The model is still useful for learning:

* how a solvent box is built at a target density;
* how periodic boundaries and a thermostat are applied to many solvent
  particles;
* how structural analysis such as the radial distribution function works.

For atomistic water, a later milestone would need per-atom types, charges,
intramolecular exclusions, and likely Ewald summation for long-range
electrostatics.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...composition.system import System
from ...foundation.box import Box
from ...foundation.particles import Particle
from ...foundation.units import BOLTZMANN_CONSTANT
from ...foundation.validation import as_finite_float, as_positive_float
from ...models.potentials import CoulombPotential, LennardJonesPotential, PairPotential


@dataclass(frozen=True)
class SolventModel:
    """A single-site solvent particle model.

    Parameters
    ----------
    name:
        Human-readable model name.
    mass:
        Particle mass in Da.
    epsilon:
        Lennard-Jones well depth in kJ/mol.
    sigma:
        Lennard-Jones diameter in nm.
    charge:
        Optional point charge in elementary-charge units.  A non-zero charge
        adds a Coulomb pair potential.  Long-range electrostatics are not yet
        handled with Ewald summation.
    """

    name: str
    mass: float
    epsilon: float
    sigma: float
    charge: float = 0.0

    def __post_init__(self) -> None:
        name = str(self.name).strip()
        if not name:
            raise ValueError("name must be a non-empty string.")
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "mass", as_positive_float(self.mass, "mass"))
        object.__setattr__(self, "epsilon", as_positive_float(self.epsilon, "epsilon"))
        object.__setattr__(self, "sigma", as_positive_float(self.sigma, "sigma"))
        object.__setattr__(self, "charge", as_finite_float(self.charge, "charge"))

    @property
    def pair_potentials(self) -> tuple[PairPotential, ...]:
        """Return the nonbonded potentials for this solvent."""

        potentials: list[PairPotential] = [
            LennardJonesPotential(epsilon=self.epsilon, sigma=self.sigma)
        ]
        if self.charge != 0.0:
            potentials.append(
                CoulombPotential(charge_i=self.charge, charge_j=self.charge)
            )
        return tuple(potentials)


#: Argon-like monatomic solvent, useful for generic nonpolar solvation tests.
ARGON = SolventModel(
    name="argon",
    mass=39.948,
    epsilon=0.997,
    sigma=0.34,
)

#: United-atom water-like solvent.  This is not atomistic water: the hydrogens
#: and hydrogen bonding are intentionally omitted.
COARSE_GRAINED_WATER = SolventModel(
    name="coarse-grained water",
    mass=18.015,
    epsilon=0.65,
    sigma=0.31,
)


def _validate_count(n_molecules: object) -> int:
    if isinstance(n_molecules, bool) or not isinstance(
        n_molecules, (int, np.integer)
    ):
        raise ValueError("n_molecules must be an integer.")
    count = int(n_molecules)
    if count < 1:
        raise ValueError("n_molecules must be at least 1.")
    return count


def _box_length_from_density(n_molecules: int, density: float) -> float:
    density_value = as_positive_float(density, "density")
    return float((n_molecules / density_value) ** (1.0 / 3.0))


def _check_positions(
    positions: np.ndarray,
    box: Box,
    min_distance: float,
) -> None:
    n_particles = positions.shape[0]
    if n_particles < 2:
        return

    displacement = positions[:, np.newaxis, :] - positions[np.newaxis, :, :]
    displacement = box.minimum_image_array(displacement)
    distances = np.linalg.norm(displacement, axis=-1)
    upper = np.triu_indices(n_particles, k=1)
    if np.any(distances[upper] < min_distance):
        raise ValueError(
            "Solvent placement produced overlapping particles. "
            "Reduce the density, increase box_length, or lower min_distance."
        )


def _lattice_positions(
    n_molecules: int,
    box_length: float,
    box: Box,
    min_distance: float,
) -> np.ndarray:
    n_side = int(np.ceil(n_molecules ** (1.0 / 3.0)))
    spacing = box_length / n_side
    positions = []
    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                if len(positions) == n_molecules:
                    break
                positions.append(
                    [
                        (i + 0.5) * spacing,
                        (j + 0.5) * spacing,
                        (k + 0.5) * spacing,
                    ]
                )
            if len(positions) == n_molecules:
                break
        if len(positions) == n_molecules:
            break

    coordinates = np.asarray(positions, dtype=float)
    _check_positions(coordinates, box, min_distance)
    return coordinates


def _random_positions(
    n_molecules: int,
    box: Box,
    min_distance: float,
    rng: np.random.Generator,
) -> np.ndarray:
    positions: list[np.ndarray] = []
    max_attempts = max(1000, 1000 * n_molecules)

    for _ in range(max_attempts):
        if len(positions) == n_molecules:
            return np.asarray(positions, dtype=float)

        candidate = rng.uniform(0.0, box.lengths)
        if all(
            box.distance(candidate - existing) >= min_distance
            for existing in positions
        ):
            positions.append(candidate)

    raise ValueError(
        "Could not place all solvent particles without overlaps. "
        "Reduce the density or min_distance, or use lattice placement."
    )


def build_solvent_box(
    model: SolventModel,
    n_molecules: int,
    *,
    box_length: float | None = None,
    density: float | None = None,
    placement: str = "lattice",
    temperature: float | None = None,
    seed: int | None = None,
    min_distance: float | None = None,
) -> System:
    """Build an explicit solvent box.

    Exactly one of ``box_length`` and ``density`` must be supplied.  The
    density is in particles per nm^3.

    Parameters
    ----------
    model:
        Solvent particle model.
    n_molecules:
        Number of solvent particles.
    box_length:
        Cubic box edge length in nm.
    density:
        Particle number density in nm^-3.  The box length is chosen as
        ``(n_molecules / density) ** (1/3)``.
    placement:
        ``"lattice"`` (deterministic cubic lattice) or ``"random"``.
    temperature:
        Optional initial temperature in K.  When given, velocities are drawn
        from a Maxwell-Boltzmann distribution and the center-of-mass velocity
        is removed.
    seed:
        Random seed for placement and/or velocities.
    min_distance:
        Minimum allowed particle separation in nm.  Defaults to
        ``0.9 * model.sigma``.

    Returns
    -------
    System
        A periodic system containing the solvent particles.
    """

    if not isinstance(model, SolventModel):
        raise TypeError("model must be a SolventModel.")

    count = _validate_count(n_molecules)
    if (box_length is None) == (density is None):
        raise ValueError("Provide exactly one of box_length or density.")

    if box_length is not None:
        edge = as_positive_float(box_length, "box_length")
    else:
        edge = _box_length_from_density(count, density)

    placement_value = str(placement).lower()
    if placement_value not in {"lattice", "random"}:
        raise ValueError("placement must be 'lattice' or 'random'.")

    threshold = (
        0.9 * model.sigma
        if min_distance is None
        else as_positive_float(min_distance, "min_distance")
    )

    box = Box([edge, edge, edge])
    rng = np.random.default_rng(seed)

    if placement_value == "lattice":
        positions = _lattice_positions(count, edge, box, threshold)
    else:
        positions = _random_positions(count, box, threshold, rng)

    if temperature is None:
        velocities = np.zeros((count, 3), dtype=float)
    else:
        target = as_finite_float(temperature, "temperature")
        if target < 0.0:
            raise ValueError("temperature must be non-negative.")
        if target == 0.0:
            velocities = np.zeros((count, 3), dtype=float)
        else:
            standard_deviation = np.sqrt(
                BOLTZMANN_CONSTANT * target / model.mass
            )
            velocities = rng.normal(
                scale=standard_deviation,
                size=(count, 3),
            )
            if count > 1:
                velocities -= np.mean(velocities, axis=0)

    particles = [
        Particle(
            mass=model.mass,
            position=positions[index],
            velocity=velocities[index],
        )
        for index in range(count)
    ]

    return System(
        particles=particles,
        pair_potentials=model.pair_potentials,
        box=box,
    )
