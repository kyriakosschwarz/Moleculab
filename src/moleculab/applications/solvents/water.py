# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""A minimal three-site water model.

The model is intentionally close to SPC/E:

* one oxygen and two hydrogen particles per molecule;
* O-H bonds and an H-O-H angle represented by harmonic bonded terms;
* Lennard-Jones parameters on oxygen and hydrogen;
* partial charges on all three sites;
* reaction-field electrostatics with a spherical cutoff;
* automatic 1-2 and 1-3 nonbonded exclusions from the topology.

This is still an educational model, not a production water model.  It omits
long-range Ewald summation, virtual sites, and rigid-body constraints.  Its
purpose is to show how atom types, charges, exclusions, and long-range
electrostatics fit together in a solvated molecular system.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ...composition.constraints import DistanceConstraint
from ...composition.system import System
from ...foundation.box import Box
from ...foundation.particles import Particle
from ...foundation.units import BOLTZMANN_CONSTANT
from ...foundation.validation import (
    as_finite_float,
    as_nonnegative_float,
    as_positive_float,
)
from ...models.bonded import HarmonicAnglePotential, HarmonicBondPotential
from ...models.ewald import EwaldElectrostatics
from ...models.nonbonded import (
    LennardJonesParameters,
    NonbondedModel,
    ReactionFieldElectrostatics,
)
from ...models.topology import Angle, Bond, Topology


@dataclass(frozen=True)
class WaterModel:
    """Parameters for a flexible three-site water model.

    Defaults are close to SPC/E.  ``angle`` is in radians.  Hydrogen has no
    Lennard-Jones well depth by default, but a finite ``sigma`` is kept so that
    the Lorentz-Berthelot combining rule is well-defined.
    """

    name: str = "SPC/E-like"
    oxygen_mass: float = 15.999
    hydrogen_mass: float = 1.008
    bond_length: float = 0.1
    angle: float = np.deg2rad(109.47)
    oxygen_epsilon: float = 0.6502
    oxygen_sigma: float = 0.3166
    hydrogen_epsilon: float = 0.0
    hydrogen_sigma: float = 0.1
    oxygen_charge: float = -0.8476
    hydrogen_charge: float = 0.4238
    bond_stiffness: float = 10_000.0
    angle_stiffness: float = 1_000.0
    dielectric: float = 78.3

    def __post_init__(self) -> None:
        name = str(self.name).strip()
        if not name:
            raise ValueError("name must be a non-empty string.")
        object.__setattr__(self, "name", name)
        object.__setattr__(
            self, "oxygen_mass", as_positive_float(self.oxygen_mass, "oxygen_mass")
        )
        object.__setattr__(
            self,
            "hydrogen_mass",
            as_positive_float(self.hydrogen_mass, "hydrogen_mass"),
        )
        object.__setattr__(
            self, "bond_length", as_positive_float(self.bond_length, "bond_length")
        )
        angle = as_finite_float(self.angle, "angle")
        if not 0.0 < angle < np.pi:
            raise ValueError("angle must be between 0 and pi radians.")
        object.__setattr__(self, "angle", angle)
        object.__setattr__(
            self,
            "oxygen_epsilon",
            as_nonnegative_float(self.oxygen_epsilon, "oxygen_epsilon"),
        )
        object.__setattr__(
            self,
            "oxygen_sigma",
            as_positive_float(self.oxygen_sigma, "oxygen_sigma"),
        )
        object.__setattr__(
            self,
            "hydrogen_epsilon",
            as_nonnegative_float(self.hydrogen_epsilon, "hydrogen_epsilon"),
        )
        object.__setattr__(
            self,
            "hydrogen_sigma",
            as_positive_float(self.hydrogen_sigma, "hydrogen_sigma"),
        )
        object.__setattr__(
            self,
            "oxygen_charge",
            as_finite_float(self.oxygen_charge, "oxygen_charge"),
        )
        object.__setattr__(
            self,
            "hydrogen_charge",
            as_finite_float(self.hydrogen_charge, "hydrogen_charge"),
        )
        object.__setattr__(
            self,
            "bond_stiffness",
            as_positive_float(self.bond_stiffness, "bond_stiffness"),
        )
        object.__setattr__(
            self,
            "angle_stiffness",
            as_positive_float(self.angle_stiffness, "angle_stiffness"),
        )
        object.__setattr__(
            self, "dielectric", as_positive_float(self.dielectric, "dielectric")
        )


SPCE_WATER = WaterModel()


def _validate_molecule_count(n_molecules: object) -> int:
    if isinstance(n_molecules, bool) or not isinstance(
        n_molecules, (int, np.integer)
    ):
        raise ValueError("n_molecules must be an integer.")
    count = int(n_molecules)
    if count < 1:
        raise ValueError("n_molecules must be at least 1.")
    return count


def _rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    skew = np.array(
        [
            [0.0, -z, y],
            [z, 0.0, -x],
            [-y, x, 0.0],
        ]
    )
    return (
        np.eye(3)
        + np.sin(angle) * skew
        + (1.0 - np.cos(angle)) * (skew @ skew)
    )


def build_water_box(
    n_molecules: int,
    *,
    model: WaterModel = SPCE_WATER,
    box_length: float | None = None,
    density: float | None = None,
    temperature: float | None = None,
    seed: int | None = None,
    cutoff: float | None = None,
    rigid: bool = True,
    electrostatics: ReactionFieldElectrostatics | EwaldElectrostatics | None = None,
) -> System:
    """Build a periodic box of three-site water molecules.

    Exactly one of ``box_length`` and ``density`` must be provided.  Density is
    in molecules per nm^3.  Oxygen atoms are placed on a cubic lattice and each
    molecule receives a random orientation.  The resulting system contains:

    * ``3 * n_molecules`` particles with types ``"O"`` and ``"H"``;
    * a topology with two O-H bonds and one H-O-H angle per molecule;
    * a :class:`~moleculab.models.nonbonded.NonbondedModel` with reaction-field
      electrostatics by default, or Ewald when an ``EwaldElectrostatics``
      object is supplied;
    * automatic 1-2 and 1-3 nonbonded exclusions.

    When ``rigid`` is true (the default), each molecule also receives three
    SHAKE/RATTLE distance constraints (two O-H bonds and the H-H distance), so
    the constrained Velocity Verlet integrator keeps it geometrically rigid.
    """

    if not isinstance(model, WaterModel):
        raise TypeError("model must be a WaterModel.")
    if not isinstance(rigid, bool):
        raise TypeError("rigid must be a boolean.")
    if electrostatics is not None and not isinstance(
        electrostatics, (ReactionFieldElectrostatics, EwaldElectrostatics)
    ):
        raise TypeError(
            "electrostatics must be ReactionFieldElectrostatics, "
            "EwaldElectrostatics, or None."
        )
    count = _validate_molecule_count(n_molecules)
    if (box_length is None) == (density is None):
        raise ValueError("Provide exactly one of box_length or density.")

    if box_length is not None:
        edge = as_positive_float(box_length, "box_length")
    else:
        density_value = as_positive_float(density, "density")
        edge = float((count / density_value) ** (1.0 / 3.0))

    box = Box([edge, edge, edge])
    cutoff_value = (
        0.45 * float(np.min(box.lengths))
        if cutoff is None
        else as_positive_float(cutoff, "cutoff")
    )
    if cutoff_value > float(np.min(box.lengths)) / 2.0 + 1e-12:
        raise ValueError(
            "cutoff must not exceed half the shortest box edge."
        )

    n_side = int(np.ceil(count ** (1.0 / 3.0)))
    spacing = edge / n_side
    if spacing / 2.0 <= model.bond_length:
        raise ValueError(
            "Box is too small for the water geometry; reduce the density "
            "or use a larger box."
        )

    rng = np.random.default_rng(seed)
    half_angle = model.angle / 2.0
    local_h1 = model.bond_length * np.array(
        [np.cos(half_angle), np.sin(half_angle), 0.0]
    )
    local_h2 = model.bond_length * np.array(
        [np.cos(half_angle), -np.sin(half_angle), 0.0]
    )

    oxygen_positions = []
    for i in range(n_side):
        for j in range(n_side):
            for k in range(n_side):
                if len(oxygen_positions) == count:
                    break
                oxygen_positions.append(
                    [
                        (i + 0.5) * spacing,
                        (j + 0.5) * spacing,
                        (k + 0.5) * spacing,
                    ]
                )
            if len(oxygen_positions) == count:
                break
        if len(oxygen_positions) == count:
            break

    particles: list[Particle] = []
    bonds: list[Bond] = []
    angles: list[Angle] = []

    for molecule_index, oxygen_position in enumerate(oxygen_positions):
        oxygen_position = np.asarray(oxygen_position, dtype=float)
        axis = rng.normal(size=3)
        rotation = _rotation_matrix(axis, rng.uniform(0.0, 2.0 * np.pi))
        hydrogen1 = oxygen_position + rotation @ local_h1
        hydrogen2 = oxygen_position + rotation @ local_h2

        oxygen_index = len(particles)
        hydrogen1_index = oxygen_index + 1
        hydrogen2_index = oxygen_index + 2

        particles.extend(
            [
                Particle(
                    mass=model.oxygen_mass,
                    position=oxygen_position,
                    atom_type="O",
                    charge=model.oxygen_charge,
                ),
                Particle(
                    mass=model.hydrogen_mass,
                    position=hydrogen1,
                    atom_type="H",
                    charge=model.hydrogen_charge,
                ),
                Particle(
                    mass=model.hydrogen_mass,
                    position=hydrogen2,
                    atom_type="H",
                    charge=model.hydrogen_charge,
                ),
            ]
        )

        bond_potential = HarmonicBondPotential(
            stiffness=model.bond_stiffness,
            length0=model.bond_length,
        )
        bonds.append(Bond(oxygen_index, hydrogen1_index, bond_potential))
        bonds.append(Bond(oxygen_index, hydrogen2_index, bond_potential))
        angles.append(
            Angle(
                hydrogen1_index,
                oxygen_index,
                hydrogen2_index,
                HarmonicAnglePotential(
                    stiffness=model.angle_stiffness,
                    angle0=model.angle,
                ),
            )
        )

    if temperature is not None:
        target = as_finite_float(temperature, "temperature")
        if target < 0.0:
            raise ValueError("temperature must be non-negative.")
        if target == 0.0:
            velocities = np.zeros((len(particles), 3), dtype=float)
        else:
            velocities = np.array(
                [
                    rng.normal(
                        scale=np.sqrt(BOLTZMANN_CONSTANT * target / particle.mass),
                        size=3,
                    )
                    for particle in particles
                ]
            )
            masses = np.array([particle.mass for particle in particles])
            center_of_mass_velocity = np.sum(
                masses[:, np.newaxis] * velocities,
                axis=0,
            ) / np.sum(masses)
            velocities -= center_of_mass_velocity
            for particle, velocity in zip(particles, velocities):
                particle.velocity = velocity.copy()
    # Otherwise particles keep their zero default velocities.

    nonbonded_model = NonbondedModel(
        parameters={
            "O": LennardJonesParameters(
                epsilon=model.oxygen_epsilon,
                sigma=model.oxygen_sigma,
            ),
            "H": LennardJonesParameters(
                epsilon=model.hydrogen_epsilon,
                sigma=model.hydrogen_sigma,
            ),
        },
        cutoff=cutoff_value,
        electrostatics=(
            ReactionFieldElectrostatics(model.dielectric)
            if electrostatics is None
            else electrostatics
        ),
    )

    constraints: list[DistanceConstraint] = []
    if rigid:
        hydrogen_hydrogen = 2.0 * model.bond_length * np.sin(model.angle / 2.0)
        for molecule_index in range(count):
            oxygen_index = 3 * molecule_index
            hydrogen1_index = oxygen_index + 1
            hydrogen2_index = oxygen_index + 2
            constraints.extend(
                [
                    DistanceConstraint(
                        oxygen_index,
                        hydrogen1_index,
                        model.bond_length,
                    ),
                    DistanceConstraint(
                        oxygen_index,
                        hydrogen2_index,
                        model.bond_length,
                    ),
                    DistanceConstraint(
                        hydrogen1_index,
                        hydrogen2_index,
                        hydrogen_hydrogen,
                    ),
                ]
            )

    return System(
        particles=particles,
        topology=Topology(bonds=bonds, angles=angles),
        box=box,
        nonbonded_model=nonbonded_model,
        constraints=constraints,
    )
