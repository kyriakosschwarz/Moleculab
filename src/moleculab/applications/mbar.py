# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Multistate Bennett Acceptance Ratio (MBAR) for alchemical windows.

MBAR combines samples from several thermodynamic states and estimates their
relative free energies simultaneously.  For samples :math:`x_n` drawn from
state :math:`s_n`, with reduced potential
:math:`u_i(x) = U_i(x)/(RT)`, the self-consistent equations are

.. math::

    e^{-f_i} = \\sum_{n=1}^{N}
    \\frac{e^{-u_i(x_n)}}
    {\\sum_{k=1}^{K} N_k e^{f_k - u_k(x_n)}},

where :math:`N_k` is the number of samples from state k.  The free energies
:math:`f_i` are defined up to an additive constant; Moleculab fixes
:math:`f_0 = 0`, so :math:`\\Delta G_i = RT f_i`.

The implementation uses fixed-point iteration in log space and a harmonic
alchemical test system with analytically known free energies.  It does not yet
provide MBAR uncertainty estimates or the full overlap matrix.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from ..foundation.units import BOLTZMANN_CONSTANT
from ..models.potentials import HarmonicPotential


@dataclass(frozen=True)
class MBARResult:
    """Result of an MBAR free-energy calculation.

    Attributes
    ----------
    free_energies:
        Dimensionless reduced free energies :math:`f_i`, with ``f[0] = 0``.
    delta_g:
        Free energies relative to state 0 in kJ/mol.
    temperature:
        Temperature in K.
    state_counts:
        Number of samples contributed by each state.
    iterations:
        Fixed-point iterations used.
    converged:
        Whether the self-consistency tolerance was reached.
    """

    free_energies: np.ndarray
    delta_g: np.ndarray
    temperature: float
    state_counts: np.ndarray
    iterations: int
    converged: bool


def _logsumexp(values: np.ndarray, axis: int | None = None) -> np.ndarray | float:
    maximum = np.max(values, axis=axis, keepdims=True)
    result = maximum + np.log(np.sum(np.exp(values - maximum), axis=axis, keepdims=True))
    if axis is None:
        return float(result)
    return np.squeeze(result, axis=axis)


def mbar(
    reduced_potentials: np.ndarray,
    sample_states: np.ndarray,
    temperature: float,
    tolerance: float = 1e-10,
    max_iterations: int = 10_000,
    damping: float = 0.5,
) -> MBARResult:
    """Estimate free energies from reduced potentials.

    Parameters
    ----------
    reduced_potentials:
        Array with shape ``(n_samples, n_states)`` containing
        ``U_k(x_n) / (RT)``.
    sample_states:
        Integer state index for each sample, shape ``(n_samples,)``.
    temperature:
        Temperature in K.
    tolerance:
        Convergence threshold on the maximum change in ``f``.
    max_iterations:
        Maximum fixed-point iterations.
    damping:
        Fixed-point mixing fraction in ``(0, 1]``.

    Raises
    ------
    RuntimeError
        If the self-consistent equations do not converge.
    """

    potentials = np.asarray(reduced_potentials, dtype=float)
    states = np.asarray(sample_states)
    if potentials.ndim != 2:
        raise ValueError(
            "reduced_potentials must have shape (n_samples, n_states), "
            f"got {potentials.shape}."
        )
    n_samples, n_states = potentials.shape
    if n_states < 2:
        raise ValueError("MBAR requires at least two states.")
    if states.shape != (n_samples,):
        raise ValueError(
            "sample_states must have shape (n_samples,), "
            f"got {states.shape}."
        )
    if not np.issubdtype(states.dtype, np.integer):
        raise ValueError("sample_states must contain integer state indices.")
    if np.any(states < 0) or np.any(states >= n_states):
        raise ValueError("sample_states contains an out-of-range state index.")
    if not np.all(np.isfinite(potentials)):
        raise ValueError("reduced_potentials must contain only finite values.")

    temperature_value = float(temperature)
    if not np.isfinite(temperature_value) or temperature_value <= 0.0:
        raise ValueError("temperature must be finite and positive.")
    tolerance_value = float(tolerance)
    if not np.isfinite(tolerance_value) or tolerance_value <= 0.0:
        raise ValueError("tolerance must be finite and positive.")
    if isinstance(max_iterations, bool) or not isinstance(
        max_iterations, (int, np.integer)
    ):
        raise ValueError("max_iterations must be an integer.")
    if int(max_iterations) < 1:
        raise ValueError("max_iterations must be at least 1.")
    max_iterations = int(max_iterations)
    damping_value = float(damping)
    if not 0.0 < damping_value <= 1.0:
        raise ValueError("damping must be in the interval (0, 1].")

    counts = np.bincount(states, minlength=n_states).astype(float)
    if np.any(counts == 0.0):
        raise ValueError("Every MBAR state must contribute at least one sample.")
    log_counts = np.log(counts)
    free_energies = np.zeros(n_states, dtype=float)
    converged = False

    for iteration in range(1, max_iterations + 1):
        log_denominator = _logsumexp(
            log_counts[np.newaxis, :]
            + free_energies[np.newaxis, :]
            - potentials,
            axis=1,
        )
        updated = -_logsumexp(
            -potentials - np.asarray(log_denominator)[:, np.newaxis],
            axis=0,
        )
        updated = np.asarray(updated, dtype=float)
        updated -= updated[0]
        change = float(np.max(np.abs(updated - free_energies)))
        free_energies = (
            damping_value * updated + (1.0 - damping_value) * free_energies
        )
        free_energies -= free_energies[0]
        if change < tolerance_value:
            converged = True
            break
    else:
        iteration = max_iterations

    if not converged:
        raise RuntimeError(
            "MBAR did not converge within "
            f"{max_iterations} iterations (last change {change:.3e})."
        )

    rt = BOLTZMANN_CONSTANT * temperature_value
    return MBARResult(
        free_energies=free_energies,
        delta_g=free_energies * rt,
        temperature=temperature_value,
        state_counts=counts.astype(int),
        iterations=iteration,
        converged=converged,
    )


def mbar_harmonic(
    potentials: Sequence[HarmonicPotential],
    temperature: float,
    n_samples: int = 10_000,
    seed: int | None = None,
    tolerance: float = 1e-10,
    max_iterations: int = 10_000,
) -> MBARResult:
    """Run MBAR on a sequence of harmonic alchemical states.

    Each state contributes the same number of samples.  Samples are drawn from
    the exact Gaussian Boltzmann distribution of each potential.
    """

    states = tuple(potentials)
    if len(states) < 2:
        raise ValueError("MBAR requires at least two harmonic potentials.")
    if not all(isinstance(potential, HarmonicPotential) for potential in states):
        raise TypeError("potentials must contain HarmonicPotential objects.")
    temperature_value = float(temperature)
    if not np.isfinite(temperature_value) or temperature_value <= 0.0:
        raise ValueError("temperature must be finite and positive.")
    if isinstance(n_samples, bool) or not isinstance(n_samples, (int, np.integer)):
        raise ValueError("n_samples must be an integer.")
    if int(n_samples) < 1:
        raise ValueError("n_samples must be at least 1.")
    n_samples = int(n_samples)

    rt = BOLTZMANN_CONSTANT * temperature_value
    rng = np.random.default_rng(seed)
    samples = []
    sample_states = []
    for state_index, potential in enumerate(states):
        standard_deviation = np.sqrt(rt / potential.stiffness)
        samples.append(
            rng.normal(
                loc=potential.center,
                scale=standard_deviation,
                size=(n_samples, 3),
            )
        )
        sample_states.extend([state_index] * n_samples)

    coordinates = np.vstack(samples)
    reduced_potentials = np.empty((coordinates.shape[0], len(states)), dtype=float)
    for state_index, potential in enumerate(states):
        displacement = coordinates - potential.center
        reduced_potentials[:, state_index] = (
            0.5
            * potential.stiffness
            * np.sum(displacement * displacement, axis=1)
            / rt
        )

    return mbar(
        reduced_potentials,
        np.asarray(sample_states, dtype=int),
        temperature_value,
        tolerance=tolerance,
        max_iterations=max_iterations,
    )
