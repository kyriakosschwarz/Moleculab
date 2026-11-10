"""Energy minimization by steepest descent.

Molecular dynamics asks, "Where do the atoms go next?"  Energy minimization
asks a different question: "What is the nearest low-energy configuration?"

At a minimum, the net force on every atom is zero:

.. math:: \\mathbf F_i = -\\nabla_i U = 0.

Steepest descent takes a small step along the force (the negative gradient)::

    r_new = r + step_size * F.

The step size has units of nm^2/(kJ/mol), so ``step_size * force`` has units
of nm.  A simple backtracking rule reduces the step whenever it would raise the
energy, which keeps the method stable for stiff bonded terms.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .constraints import assert_constraints_satisfied, shake_positions
from .particles import Particle
from .system import System


@dataclass(frozen=True)
class MinimizationResult:
    """The outcome of an energy-minimization run.

    Attributes
    ----------
    system:
        Final system, including optimized positions.
    energies:
        Potential energy after each accepted state, starting with the input.
    max_forces:
        Largest per-atom force magnitude after each accepted state.
    step_sizes:
        Accepted step size for each move.  This array is one element shorter
        than ``energies``.
    converged:
        True if ``max_forces[-1] <= force_tolerance``.
    """

    system: System
    energies: np.ndarray
    max_forces: np.ndarray
    step_sizes: np.ndarray
    converged: bool

    @property
    def steps(self) -> int:
        """Number of accepted minimization steps."""

        return int(self.step_sizes.size)


def _validate_step_size(step_size: float) -> float:
    value = float(step_size)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"step_size must be finite and positive, got {step_size!r}.")
    return value


def _with_positions(system: System, positions: np.ndarray) -> System:
    coordinates = np.asarray(positions, dtype=float)
    if system.constraints:
        masses = np.array([particle.mass for particle in system.particles])
        coordinates = shake_positions(
            coordinates,
            masses,
            system.constraints,
        ).positions
    updated = System(
        particles=[
            Particle(
                mass=particle.mass,
                position=coordinates[index],
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
    assert_constraints_satisfied(updated)
    return updated


def max_force_magnitude(system: System) -> float:
    """Return the largest per-atom force magnitude in kJ/(mol nm)."""

    forces = system.forces()
    return float(np.max(np.linalg.norm(forces, axis=1)))


def steepest_descent_step(system: System, step_size: float) -> System:
    """Move every atom a small distance along its force.

    The input system is not modified.  Use :func:`minimize_energy` when you
    want automatic step-size control and convergence testing.
    """

    step = _validate_step_size(step_size)
    new_positions = system.positions + step * system.forces()
    return _with_positions(system, new_positions)


def minimize_energy(
    system: System,
    step_size: float = 1e-5,
    max_steps: int = 1000,
    force_tolerance: float = 1e-3,
) -> MinimizationResult:
    """Minimize potential energy with backtracking steepest descent.

    Parameters
    ----------
    system:
        Initial configuration.  It is copied before minimization.
    step_size:
        Initial step size in nm^2/(kJ/mol).
    max_steps:
        Maximum number of accepted moves.
    force_tolerance:
        Convergence threshold for the largest per-atom force magnitude in
        kJ/(mol nm).

    Returns
    -------
    MinimizationResult
        Final system and the energy/force history.
    """

    step = _validate_step_size(step_size)
    if not isinstance(max_steps, (int, np.integer)) or max_steps < 0:
        raise ValueError("max_steps must be a non-negative integer.")
    tolerance = float(force_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("force_tolerance must be finite and non-negative.")

    current = System(
        particles=[particle.copy() for particle in system.particles],
        pair_potentials=system.pair_potentials,
        topology=system.topology,
        box=system.box,
        nonbonded_model=system.nonbonded_model,
        constraints=system.constraints,
        indexed_pair_potentials=system.indexed_pair_potentials,
    )

    energies = [current.potential_energy()]
    force_history = [max_force_magnitude(current)]
    accepted_steps: list[float] = []
    converged = force_history[0] <= tolerance

    for _ in range(int(max_steps)):
        if converged:
            break

        forces = current.forces()
        current_energy = energies[-1]
        trial_step = step
        accepted = False

        for _ in range(20):
            candidate_positions = current.positions + trial_step * forces
            candidate = _with_positions(current, candidate_positions)
            candidate_energy = candidate.potential_energy()

            if candidate_energy < current_energy:
                current = candidate
                energies.append(candidate_energy)
                force_history.append(max_force_magnitude(current))
                accepted_steps.append(trial_step)
                accepted = True
                break

            trial_step *= 0.5

        if not accepted:
            break

        converged = force_history[-1] <= tolerance

    return MinimizationResult(
        system=current,
        energies=np.asarray(energies),
        max_forces=np.asarray(force_history),
        step_sizes=np.asarray(accepted_steps),
        converged=bool(converged),
    )
