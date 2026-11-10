# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""SHAKE and RATTLE distance constraints.

SHAKE corrects positions after an unconstrained integration step so that
selected interatomic distances match their target values.  RATTLE then removes
the radial velocity component so the constraints remain satisfied at the next
step.

For a constraint between atoms i and j with displacement
:math:`\\mathbf r = \\mathbf r_i - \\mathbf r_j`, SHAKE applies

.. math::

    \\Delta \\mathbf r_i = -\\frac{\\mu}{m_i}\\frac{r-d}{r}\\mathbf r,
    \\qquad
    \\Delta \\mathbf r_j = +\\frac{\\mu}{m_j}\\frac{r-d}{r}\\mathbf r,

where :math:`\\mu = (1/m_i + 1/m_j)^{-1}`.  This preserves the center of mass.

RATTLE removes the radial relative velocity:

.. math::

    \\Delta \\mathbf v_i = -\\frac{\\mu}{m_i}
    \\frac{\\mathbf r\\cdot\\mathbf v}{r^2}\\mathbf r,
    \\qquad
    \\Delta \\mathbf v_j = +\\frac{\\mu}{m_j}
    \\frac{\\mathbf r\\cdot\\mathbf v}{r^2}\\mathbf r.

The implementation iterates over all constraints because coupled constraints
in a molecule affect one another.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from ..foundation.validation import as_positive_float


@dataclass(frozen=True)
class DistanceConstraint:
    """A rigid distance between two atoms."""

    i: int
    j: int
    distance: float

    def __post_init__(self) -> None:
        i = int(self.i)
        j = int(self.j)
        if i < 0 or j < 0:
            raise ValueError("Constraint atom indices must be non-negative.")
        if i == j:
            raise ValueError("A distance constraint requires two different atoms.")
        object.__setattr__(self, "i", i)
        object.__setattr__(self, "j", j)
        object.__setattr__(
            self,
            "distance",
            as_positive_float(self.distance, "distance"),
        )


@dataclass(frozen=True)
class ShakeResult:
    """Result of a SHAKE position correction."""

    positions: np.ndarray
    max_error: float
    iterations: int


@dataclass(frozen=True)
class RattleResult:
    """Result of a RATTLE velocity correction."""

    velocities: np.ndarray
    max_error: float
    iterations: int


def _validate_inputs(
    positions: np.ndarray,
    masses: np.ndarray,
    constraints: Sequence[DistanceConstraint],
) -> tuple[np.ndarray, np.ndarray, tuple[DistanceConstraint, ...]]:
    coordinates = np.asarray(positions, dtype=float)
    masses_array = np.asarray(masses, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(
            "positions must have shape (n_atoms, 3), "
            f"got {coordinates.shape}."
        )
    if masses_array.shape != (coordinates.shape[0],):
        raise ValueError(
            "masses must have shape (n_atoms,), "
            f"got {masses_array.shape}."
        )
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("positions must contain only finite values.")
    if not np.all(np.isfinite(masses_array)) or np.any(masses_array <= 0.0):
        raise ValueError("masses must be finite and positive.")

    constraint_tuple = tuple(constraints)
    for constraint in constraint_tuple:
        if not isinstance(constraint, DistanceConstraint):
            raise TypeError("constraints must contain DistanceConstraint objects.")
        if constraint.i >= coordinates.shape[0] or constraint.j >= coordinates.shape[0]:
            raise ValueError(
                f"Constraint atom index out of range for {coordinates.shape[0]} atoms."
            )
    return coordinates, masses_array, constraint_tuple


def constraint_error(
    positions: np.ndarray,
    constraints: Sequence[DistanceConstraint],
) -> float:
    """Return the largest absolute distance-constraint error in nm."""

    maximum = 0.0
    for constraint in constraints:
        distance = float(
            np.linalg.norm(positions[constraint.i] - positions[constraint.j])
        )
        maximum = max(maximum, abs(distance - constraint.distance))
    return maximum


def assert_constraints_satisfied(
    system,
    tolerance: float = 1e-7,
) -> float:
    """Raise if a system's distance constraints are not satisfied.

    This is a safety net for state-advancing algorithms.  Algorithms that move
    positions must apply SHAKE (and RATTLE after a velocity update).  The guard
    turns a forgotten projection into an explicit error instead of a slowly
    drifting, silently wrong trajectory.
    """

    constraints = tuple(system.constraints)
    if not constraints:
        return 0.0
    error = constraint_error(system.positions, constraints)
    if error > tolerance:
        raise RuntimeError(
            "Distance constraints are not satisfied: "
            f"maximum error {error:.3e} nm exceeds tolerance {tolerance:.3e} nm. "
            "Apply SHAKE/RATTLE after advancing the state."
        )
    return error


def shake_positions(
    positions: np.ndarray,
    masses: np.ndarray,
    constraints: Sequence[DistanceConstraint],
    tolerance: float = 1e-10,
    max_iterations: int = 100,
) -> ShakeResult:
    """Apply SHAKE distance constraints to positions.

    Returns a corrected copy; the input array is not modified.
    """

    coordinates, masses_array, constraint_tuple = _validate_inputs(
        positions, masses, constraints
    )
    tolerance_value = as_positive_float(tolerance, "tolerance")
    iterations_limit = int(max_iterations)
    if iterations_limit < 1:
        raise ValueError("max_iterations must be at least 1.")
    corrected = coordinates.copy()
    if not constraint_tuple:
        return ShakeResult(corrected, 0.0, 0)

    iterations = 0
    for iterations in range(1, iterations_limit + 1):
        maximum_error = 0.0
        for constraint in constraint_tuple:
            displacement = corrected[constraint.i] - corrected[constraint.j]
            distance = float(np.linalg.norm(displacement))
            if distance == 0.0:
                raise ValueError(
                    "SHAKE cannot correct two atoms at the same position."
                )
            error = distance - constraint.distance
            maximum_error = max(maximum_error, abs(error))
            inverse_mass_i = 1.0 / masses_array[constraint.i]
            inverse_mass_j = 1.0 / masses_array[constraint.j]
            reduced_mass = 1.0 / (inverse_mass_i + inverse_mass_j)
            scale = reduced_mass * error / distance
            corrected[constraint.i] -= scale * inverse_mass_i * displacement
            corrected[constraint.j] += scale * inverse_mass_j * displacement

        if maximum_error <= tolerance_value:
            break
    else:
        raise RuntimeError(
            "SHAKE did not converge within "
            f"{iterations_limit} iterations (error {maximum_error:.3e})."
        )

    return ShakeResult(
        positions=corrected,
        max_error=constraint_error(corrected, constraint_tuple),
        iterations=iterations,
    )


def rattle_velocities(
    positions: np.ndarray,
    velocities: np.ndarray,
    masses: np.ndarray,
    constraints: Sequence[DistanceConstraint],
    tolerance: float = 1e-10,
    max_iterations: int = 100,
) -> RattleResult:
    """Apply RATTLE velocity constraints.

    ``positions`` should already satisfy the distance constraints, typically
    after :func:`shake_positions`.  Returns a corrected velocity copy.
    """

    coordinates, masses_array, constraint_tuple = _validate_inputs(
        positions, masses, constraints
    )
    velocity_array = np.asarray(velocities, dtype=float)
    if velocity_array.shape != coordinates.shape:
        raise ValueError(
            "velocities must have the same shape as positions, "
            f"got {velocity_array.shape} and {coordinates.shape}."
        )
    if not np.all(np.isfinite(velocity_array)):
        raise ValueError("velocities must contain only finite values.")
    tolerance_value = as_positive_float(tolerance, "tolerance")
    iterations_limit = int(max_iterations)
    if iterations_limit < 1:
        raise ValueError("max_iterations must be at least 1.")
    corrected = velocity_array.copy()
    if not constraint_tuple:
        return RattleResult(corrected, 0.0, 0)

    iterations = 0
    for iterations in range(1, iterations_limit + 1):
        maximum_error = 0.0
        for constraint in constraint_tuple:
            displacement = coordinates[constraint.i] - coordinates[constraint.j]
            squared_distance = float(np.dot(displacement, displacement))
            if squared_distance == 0.0:
                raise ValueError(
                    "RATTLE cannot constrain two atoms at the same position."
                )
            relative_velocity = (
                corrected[constraint.i] - corrected[constraint.j]
            )
            radial_velocity = float(np.dot(displacement, relative_velocity))
            maximum_error = max(
                maximum_error,
                abs(radial_velocity) / np.sqrt(squared_distance),
            )
            inverse_mass_i = 1.0 / masses_array[constraint.i]
            inverse_mass_j = 1.0 / masses_array[constraint.j]
            reduced_mass = 1.0 / (inverse_mass_i + inverse_mass_j)
            scale = reduced_mass * radial_velocity / squared_distance
            corrected[constraint.i] -= scale * inverse_mass_i * displacement
            corrected[constraint.j] += scale * inverse_mass_j * displacement

        if maximum_error <= tolerance_value:
            break
    else:
        raise RuntimeError(
            "RATTLE did not converge within "
            f"{iterations_limit} iterations (error {maximum_error:.3e})."
        )

    return RattleResult(
        velocities=corrected,
        max_error=maximum_error,
        iterations=iterations,
    )
