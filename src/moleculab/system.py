"""A collection of particles and their interactions.

The single-particle examples are useful for learning one equation at a time.
A molecular system, however, contains many particles whose interactions are
sums of pairwise and bonded terms.  :class:`System` makes that idea explicit:

.. math::

    U_{\\mathrm{total}} = \\sum_{i<j} U_{\\mathrm{nonbonded}}(r_{ij})
    + U_{\\mathrm{bonded}}.

For every nonbonded pair, the forces on the two particles are equal and
opposite.  This is Newton's third law, and it is why applying a pair force only
once is enough.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from .box import Box
from .constraints import DistanceConstraint
from .ewald import EwaldElectrostatics
from .nonbonded import NonbondedModel
from .particles import Particle
from .physics import kinetic_energy
from .potentials import PairPotential
from .topology import Topology


@dataclass
class System:
    """A set of particles, nonbonded potentials, and bonded topology.

    Parameters
    ----------
    particles:
        At least one :class:`~moleculab.particles.Particle`.
    pair_potentials:
        Potentials summed over every unique pair of particles.  An empty list
        is allowed: the particles are then non-interacting.
    topology:
        Optional bonded terms.  Atom indices must refer to entries in
        ``particles``.
    box:
        Optional orthorhombic periodic box.  When present, nonbonded
        interactions use the minimum-image convention.  Bonded terms use the
        stored coordinates directly.
    nonbonded_model:
        Optional typed nonbonded model.  It reads each particle's
        ``atom_type`` and ``charge`` and automatically excludes 1-2 and 1-3
        pairs from the topology.
    constraints:
        Optional SHAKE/RATTLE distance constraints applied by the constrained
        Velocity Verlet system integrator.
    indexed_pair_potentials:
        Optional potentials applied only to specific ``(i, j, potential)``
        atom pairs.  Useful for solute-solvent interactions that should not
        be applied to every pair.
    """

    particles: Sequence[Particle]
    pair_potentials: Sequence[PairPotential] = field(default_factory=tuple)
    topology: Topology | None = None
    box: Box | None = None
    nonbonded_model: NonbondedModel | None = None
    constraints: Sequence[DistanceConstraint] = field(default_factory=tuple)
    indexed_pair_potentials: Sequence[
        tuple[int, int, PairPotential]
    ] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        self.particles = tuple(self.particles)
        self.pair_potentials = tuple(self.pair_potentials)
        self.constraints = tuple(self.constraints)

        if not self.particles:
            raise ValueError("System must contain at least one particle.")
        for index, particle in enumerate(self.particles):
            if not isinstance(particle, Particle):
                raise TypeError(f"particles[{index}] must be a Particle.")
        if self.topology is not None:
            if not isinstance(self.topology, Topology):
                raise TypeError("topology must be a Topology or None.")
            self.topology.validate(self.n_particles)
        if self.box is not None and not isinstance(self.box, Box):
            raise TypeError("box must be a Box or None.")
        for index, constraint in enumerate(self.constraints):
            if not isinstance(constraint, DistanceConstraint):
                raise TypeError(f"constraints[{index}] must be a DistanceConstraint.")
            if constraint.i >= self.n_particles or constraint.j >= self.n_particles:
                raise ValueError(
                    f"constraints[{index}] refers to an atom outside the system."
                )
        normalized_pairs = []
        for index, item in enumerate(self.indexed_pair_potentials):
            if len(item) != 3:
                raise ValueError(
                    f"indexed_pair_potentials[{index}] must be (i, j, potential)."
                )
            i, j, potential = item
            i = int(i)
            j = int(j)
            if i == j:
                raise ValueError("An indexed pair potential needs two different atoms.")
            if i < 0 or j < 0 or i >= self.n_particles or j >= self.n_particles:
                raise ValueError(
                    f"indexed_pair_potentials[{index}] refers to an atom outside "
                    "the system."
                )
            normalized_pairs.append((i, j, potential))
        self.indexed_pair_potentials = tuple(normalized_pairs)
        if self.nonbonded_model is not None:
            if not isinstance(self.nonbonded_model, NonbondedModel):
                raise TypeError(
                    "nonbonded_model must be a NonbondedModel or None."
                )
            self.nonbonded_model.validate_atom_types(
                [particle.atom_type for particle in self.particles]
            )
            if (
                self.box is not None
                and self.nonbonded_model.cutoff
                > float(np.min(self.box.lengths)) / 2.0 + 1e-12
            ):
                raise ValueError(
                    "nonbonded cutoff must not exceed half the shortest box edge."
                )
            excluded = set(self.nonbonded_model.exclusions)
            if self.topology is not None:
                excluded.update(self.topology.nonbonded_exclusions())
            self._excluded_pairs = frozenset(excluded)

            self._ewald: EwaldElectrostatics | None = None
            if isinstance(
                self.nonbonded_model.electrostatics,
                EwaldElectrostatics,
            ):
                if self.box is None:
                    raise ValueError("Ewald electrostatics requires a periodic box.")
                self._ewald = self.nonbonded_model.electrostatics
        else:
            self._excluded_pairs = frozenset()
            self._ewald = None

    @property
    def n_particles(self) -> int:
        """Number of particles in the system."""

        return len(self.particles)

    @property
    def positions(self) -> np.ndarray:
        """Positions as an array with shape ``(n_particles, 3)``."""

        return np.array([particle.position for particle in self.particles])

    @property
    def velocities(self) -> np.ndarray:
        """Velocities as an array with shape ``(n_particles, 3)``."""

        return np.array([particle.velocity for particle in self.particles])

    @property
    def charges(self) -> np.ndarray:
        """Particle charges in elementary-charge units."""

        return np.array([particle.charge for particle in self.particles])

    @property
    def atom_types(self) -> tuple[str, ...]:
        """Particle atom-type labels."""

        return tuple(particle.atom_type for particle in self.particles)

    def kinetic_energy(self) -> float:
        """Return total kinetic energy in kJ/mol."""

        return float(sum(kinetic_energy(particle) for particle in self.particles))

    def nonbonded_energy(self) -> float:
        """Return the sum of all pair-potential energies in kJ/mol."""

        total = 0.0
        for i in range(self.n_particles):
            for j in range(i + 1, self.n_particles):
                displacement = self.particles[i].position - self.particles[j].position
                if self.box is not None:
                    displacement = self.box.minimum_image(displacement)
                for potential in self.pair_potentials:
                    total += potential.energy(displacement)
                if (
                    self.nonbonded_model is not None
                    and (i, j) not in self._excluded_pairs
                ):
                    total += self.nonbonded_model.energy(
                        self.particles[i].atom_type,
                        self.particles[j].atom_type,
                        self.particles[i].charge,
                        self.particles[j].charge,
                        displacement,
                    )
        for i, j, potential in self.indexed_pair_potentials:
            displacement = self.particles[i].position - self.particles[j].position
            if self.box is not None:
                displacement = self.box.minimum_image(displacement)
            total += potential.energy(displacement)
        if self._ewald is not None:
            total += self._ewald.global_energy(
                self.positions,
                self.charges,
                self.box,
            )
            total += self._ewald.exclusion_energy(
                self.positions,
                self.charges,
                self.box,
                self._excluded_pairs,
            )
        return float(total)

    def bonded_energy(self) -> float:
        """Return the total bonded energy in kJ/mol."""

        if self.topology is None:
            return 0.0
        return self.topology.energy(self.positions)

    def potential_energy(self) -> float:
        """Return nonbonded plus bonded potential energy in kJ/mol."""

        return self.nonbonded_energy() + self.bonded_energy()

    def forces(self) -> np.ndarray:
        """Return forces on all particles in kJ/(mol nm).

        The result has shape ``(n_particles, 3)``.  Nonbonded and bonded
        contributions are added independently.
        """

        forces = np.zeros((self.n_particles, 3), dtype=float)
        for i in range(self.n_particles):
            for j in range(i + 1, self.n_particles):
                displacement = self.particles[i].position - self.particles[j].position
                if self.box is not None:
                    displacement = self.box.minimum_image(displacement)
                for potential in self.pair_potentials:
                    pair_force = potential.force(displacement)
                    forces[i] += pair_force
                    forces[j] -= pair_force
                if (
                    self.nonbonded_model is not None
                    and (i, j) not in self._excluded_pairs
                ):
                    pair_force = self.nonbonded_model.force(
                        self.particles[i].atom_type,
                        self.particles[j].atom_type,
                        self.particles[i].charge,
                        self.particles[j].charge,
                        displacement,
                    )
                    forces[i] += pair_force
                    forces[j] -= pair_force

        for i, j, potential in self.indexed_pair_potentials:
            displacement = self.particles[i].position - self.particles[j].position
            if self.box is not None:
                displacement = self.box.minimum_image(displacement)
            pair_force = potential.force(displacement)
            forces[i] += pair_force
            forces[j] -= pair_force

        if self._ewald is not None:
            forces += self._ewald.global_forces(
                self.positions,
                self.charges,
                self.box,
            )
            forces += self._ewald.exclusion_forces(
                self.positions,
                self.charges,
                self.box,
                self._excluded_pairs,
            )

        if self.topology is not None:
            forces += self.topology.forces(self.positions)

        return forces

    def total_energy(self) -> float:
        """Return kinetic plus potential energy in kJ/mol."""

        return self.kinetic_energy() + self.potential_energy()

    def __len__(self) -> int:
        return self.n_particles
