"""Time integration for molecular dynamics.

The first integrator in Moleculab is the velocity Verlet algorithm.  It is a
good educational choice because it is simple, reversible in the absence of
rounding error, and conserves energy much better than naive Euler updates.

For a particle with position :math:`\\mathbf r`, velocity :math:`\\mathbf v`,
and acceleration :math:`\\mathbf a = \\mathbf F/m`, one step of size
:math:`\\Delta t` is

.. math::

    \\mathbf r_{n+1} &= \\mathbf r_n
        + \\mathbf v_n \\Delta t
        + \\frac{1}{2}\\mathbf a_n \\Delta t^2,

    \\mathbf v_{n+1} &= \\mathbf v_n
        + \\frac{1}{2}(\\mathbf a_n + \\mathbf a_{n+1})\\Delta t.

Notice that the force (and therefore acceleration) is evaluated at the *new*
position before the velocity is updated.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constraints import rattle_velocities, shake_positions
from .particles import Particle
from .physics import acceleration, kinetic_energy, potential_energy
from .potentials import Potential
from .system import System


@dataclass(frozen=True)
class Trajectory:
    """Snapshots produced by a single-particle simulation.

    All arrays have shape ``(n_steps + 1,)`` for scalar quantities and
    ``(n_steps + 1, 3)`` for vector quantities.  The first entry is the
    initial state.  Energies are in kJ/mol, positions in nm, velocities in
    nm/ps, and time in ps.
    """

    time: np.ndarray
    position: np.ndarray
    velocity: np.ndarray
    kinetic_energy: np.ndarray
    potential_energy: np.ndarray
    total_energy: np.ndarray

    def __len__(self) -> int:
        return int(self.time.size)


@dataclass(frozen=True)
class SystemTrajectory:
    """Snapshots produced by a multi-particle simulation.

    ``position`` and ``velocity`` have shape
    ``(n_steps + 1, n_particles, 3)``.  Energy arrays have shape
    ``(n_steps + 1,)`` and are in kJ/mol.
    """

    time: np.ndarray
    position: np.ndarray
    velocity: np.ndarray
    kinetic_energy: np.ndarray
    potential_energy: np.ndarray
    total_energy: np.ndarray

    def __len__(self) -> int:
        return int(self.time.size)


def _validate_timestep(dt: float) -> float:
    step = float(dt)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError(f"dt must be finite and positive, got {dt!r}.")
    return step


def _validate_steps(steps: int) -> int:
    if not isinstance(steps, (int, np.integer)):
        raise TypeError(f"steps must be an integer, got {type(steps).__name__}.")
    if steps < 0:
        raise ValueError(f"steps must be non-negative, got {steps}.")
    return int(steps)


def velocity_verlet_step(
    particle: Particle,
    potential: Potential,
    dt: float,
) -> Particle:
    """Advance *particle* by one velocity Verlet step of size *dt* (ps).

    A new :class:`~moleculab.particles.Particle` is returned; the input is not
    modified.  Returning a new object makes the update easy to inspect and
    avoids hidden state changes.
    """

    step = _validate_timestep(dt)

    force_now = potential.force(particle.position)
    acceleration_now = acceleration(particle, force_now)

    new_position = (
        particle.position
        + particle.velocity * step
        + 0.5 * acceleration_now * step**2
    )

    force_next = potential.force(new_position)
    acceleration_next = acceleration(particle, force_next)
    new_velocity = particle.velocity + 0.5 * (acceleration_now + acceleration_next) * step

    return Particle(
        mass=particle.mass,
        position=new_position,
        velocity=new_velocity,
        atom_type=particle.atom_type,
        charge=particle.charge,
    )


def simulate_velocity_verlet(
    particle: Particle,
    potential: Potential,
    dt: float,
    steps: int,
) -> Trajectory:
    """Simulate *particle* in *potential* using velocity Verlet.

    Parameters
    ----------
    particle:
        Initial state.  It is copied before the simulation begins.
    potential:
        Any object implementing ``energy`` and ``force``.
    dt:
        Time step in ps.  It must be positive.
    steps:
        Number of integration steps.  Zero is allowed and returns the initial
        snapshot only.

    Returns
    -------
    Trajectory
        State and energy snapshots at every step.  Use the total-energy array
        to study numerical energy conservation.
    """

    step_size = _validate_timestep(dt)
    n_steps = _validate_steps(steps)

    times = np.arange(n_steps + 1, dtype=float) * step_size
    positions = np.empty((n_steps + 1, 3), dtype=float)
    velocities = np.empty((n_steps + 1, 3), dtype=float)
    kinetic = np.empty(n_steps + 1, dtype=float)
    potential_values = np.empty(n_steps + 1, dtype=float)

    current = particle.copy()
    for index in range(n_steps + 1):
        positions[index] = current.position
        velocities[index] = current.velocity
        kinetic[index] = kinetic_energy(current)
        potential_values[index] = potential_energy(current, potential)

        if index < n_steps:
            current = velocity_verlet_step(current, potential, step_size)

    return Trajectory(
        time=times,
        position=positions,
        velocity=velocities,
        kinetic_energy=kinetic,
        potential_energy=potential_values,
        total_energy=kinetic + potential_values,
    )


def velocity_verlet_step_system(system: System, dt: float) -> System:
    """Advance a multi-particle :class:`~moleculab.system.System` by one step.

    The input system is not modified.  A new system is returned, which keeps
    the same pair potentials and particle masses.
    """

    step = _validate_timestep(dt)
    positions = system.positions
    velocities = system.velocities
    masses = np.array([particle.mass for particle in system.particles])
    inverse_masses = (1.0 / masses)[:, np.newaxis]

    forces_now = system.forces()
    acceleration_now = forces_now * inverse_masses

    new_positions = positions + velocities * step + 0.5 * acceleration_now * step**2
    if system.constraints:
        new_positions = shake_positions(
            new_positions,
            masses,
            system.constraints,
        ).positions

    intermediate = System(
        particles=[
            Particle(
                mass=particle.mass,
                position=new_positions[index],
                velocity=particle.velocity,
                atom_type=particle.atom_type,
                charge=particle.charge,
            )
            for index, particle in enumerate(system.particles)
        ],
        pair_potentials=system.pair_potentials,
        topology=system.topology,
        box=system.box,
        nonbonded_model=system.nonbonded_model,
        constraints=system.constraints,
        indexed_pair_potentials=system.indexed_pair_potentials,
    )

    forces_next = intermediate.forces()
    acceleration_next = forces_next * inverse_masses
    new_velocities = velocities + 0.5 * (acceleration_now + acceleration_next) * step
    if system.constraints:
        new_velocities = rattle_velocities(
            new_positions,
            new_velocities,
            masses,
            system.constraints,
        ).velocities

    return System(
        particles=[
            Particle(
                mass=particle.mass,
                position=new_positions[index],
                velocity=new_velocities[index],
                atom_type=particle.atom_type,
                charge=particle.charge,
            )
            for index, particle in enumerate(system.particles)
        ],
        pair_potentials=system.pair_potentials,
        topology=system.topology,
        box=system.box,
        nonbonded_model=system.nonbonded_model,
        constraints=system.constraints,
        indexed_pair_potentials=system.indexed_pair_potentials,
    )


def simulate_velocity_verlet_system(
    system: System,
    dt: float,
    steps: int,
) -> SystemTrajectory:
    """Simulate an interacting system with velocity Verlet.

    The initial system is copied before the simulation begins.  Energy
    snapshots make it possible to check conservation and numerical stability.
    """

    step_size = _validate_timestep(dt)
    n_steps = _validate_steps(steps)
    n_particles = system.n_particles

    times = np.arange(n_steps + 1, dtype=float) * step_size
    positions = np.empty((n_steps + 1, n_particles, 3), dtype=float)
    velocities = np.empty((n_steps + 1, n_particles, 3), dtype=float)
    kinetic = np.empty(n_steps + 1, dtype=float)
    potential_values = np.empty(n_steps + 1, dtype=float)

    current = System(
        particles=[particle.copy() for particle in system.particles],
        pair_potentials=system.pair_potentials,
        topology=system.topology,
        box=system.box,
        nonbonded_model=system.nonbonded_model,
        constraints=system.constraints,
        indexed_pair_potentials=system.indexed_pair_potentials,
    )
    for index in range(n_steps + 1):
        positions[index] = current.positions
        velocities[index] = current.velocities
        kinetic[index] = current.kinetic_energy()
        potential_values[index] = current.potential_energy()

        if index < n_steps:
            current = velocity_verlet_step_system(current, step_size)

    return SystemTrajectory(
        time=times,
        position=positions,
        velocity=velocities,
        kinetic_energy=kinetic,
        potential_energy=potential_values,
        total_energy=kinetic + potential_values,
    )
