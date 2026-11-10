# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Thermodynamic integration for an alchemically coupled harmonic system.

Thermodynamic integration (TI) computes a free-energy difference by integrating
the ensemble average of the Hamiltonian derivative:

.. math::

    \\Delta G = \\int_0^1
    \\left\\langle \\frac{\\partial U_\\lambda}{\\partial \\lambda}
    \\right\\rangle_\\lambda d\\lambda.

Moleculab uses the linear alchemical path

.. math:: U_\\lambda = (1-\\lambda)U_0 + \\lambda U_1,

so the derivative is simply :math:`U_1-U_0`.  The test system is a harmonic
oscillator whose stiffness and center are interpolated.  At each :math:`\\lambda`
the potential is quadratic, so samples can be drawn exactly from its Boltzmann
distribution.  The exact result is

.. math::

    \\Delta G_{\\mathrm{exact}} =
    \\frac{3}{2}RT\\ln\\left(\\frac{k_1}{k_0}\\right),

independent of the potential centers.

Numerical integration uses the trapezoidal rule or Simpson's rule.  Standard
errors are propagated from the per-window derivatives using the quadrature
weights, assuming independent windows.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..foundation.units import BOLTZMANN_CONSTANT
from ..models.potentials import HarmonicPotential


@dataclass(frozen=True)
class TIResult:
    """Result of a thermodynamic-integration calculation.

    Attributes
    ----------
    lambda_values:
        Alchemical coupling values.
    derivatives:
        Mean :math:`\\partial U_\\lambda/\\partial\\lambda` at each window.
    standard_errors:
        Standard error of each derivative.
    weights:
        Quadrature weights; ``sum(weights) == 1`` over ``[0, 1]``.
    delta_g:
        Integrated free-energy difference in kJ/mol.
    standard_error:
        Approximate propagated standard error in kJ/mol.
    method:
        ``"simpson"`` or ``"trapezoid"``.
    """

    lambda_values: np.ndarray
    derivatives: np.ndarray
    standard_errors: np.ndarray
    weights: np.ndarray
    delta_g: float
    standard_error: float
    method: str


def _validate_lambda_grid(lambda_values: np.ndarray) -> np.ndarray:
    values = np.asarray(lambda_values, dtype=float)
    if values.ndim != 1 or values.size < 2:
        raise ValueError("lambda_values must be a one-dimensional array of size >= 2.")
    if not np.all(np.isfinite(values)):
        raise ValueError("lambda_values must contain only finite values.")
    if np.any(np.diff(values) <= 0.0):
        raise ValueError("lambda_values must be strictly increasing.")
    if not np.isclose(values[0], 0.0) or not np.isclose(values[-1], 1.0):
        raise ValueError("lambda_values must run from 0 to 1.")
    return values


def _quadrature_weights(
    lambda_values: np.ndarray,
    method: str,
) -> np.ndarray:
    values = _validate_lambda_grid(lambda_values)
    spacing = np.diff(values)
    if not np.allclose(spacing, spacing[0]):
        raise ValueError("Thermodynamic-integration nodes must be uniformly spaced.")
    h = float(spacing[0])

    if method == "trapezoid":
        weights = np.ones(values.size, dtype=float)
        weights[0] = 0.5
        weights[-1] = 0.5
        return weights * h

    if method == "simpson":
        if values.size % 2 == 0:
            raise ValueError("Simpson's rule requires an odd number of nodes.")
        weights = np.ones(values.size, dtype=float)
        weights[1:-1:2] = 4.0
        weights[2:-1:2] = 2.0
        return weights * (h / 3.0)

    raise ValueError("method must be 'simpson' or 'trapezoid'.")


def integrate_thermodynamic_derivative(
    lambda_values: np.ndarray,
    derivatives: np.ndarray,
    standard_errors: np.ndarray | None = None,
    method: str = "simpson",
) -> TIResult:
    """Integrate mean Hamiltonian derivatives over the alchemical path."""

    values = _validate_lambda_grid(lambda_values)
    derivative_values = np.asarray(derivatives, dtype=float)
    if derivative_values.shape != values.shape:
        raise ValueError("derivatives must have the same shape as lambda_values.")
    if not np.all(np.isfinite(derivative_values)):
        raise ValueError("derivatives must contain only finite values.")

    weights = _quadrature_weights(values, method)
    delta_g = float(np.dot(weights, derivative_values))

    if standard_errors is None:
        errors = np.zeros_like(derivative_values)
    else:
        errors = np.asarray(standard_errors, dtype=float)
        if errors.shape != values.shape:
            raise ValueError(
                "standard_errors must have the same shape as lambda_values."
            )
        if not np.all(np.isfinite(errors)) or np.any(errors < 0.0):
            raise ValueError("standard_errors must be finite and non-negative.")
    standard_error = float(np.sqrt(np.sum((weights * errors) ** 2)))

    return TIResult(
        lambda_values=values,
        derivatives=derivative_values,
        standard_errors=errors,
        weights=weights,
        delta_g=delta_g,
        standard_error=standard_error,
        method=method,
    )


def thermodynamic_integration_harmonic(
    potential0: HarmonicPotential,
    potential1: HarmonicPotential,
    temperature: float,
    n_windows: int = 21,
    n_samples: int = 20_000,
    seed: int | None = None,
    method: str = "simpson",
) -> TIResult:
    """Run TI for a harmonic stiffness/center transformation.

    Samples at each :math:`\\lambda` are drawn from the exact Gaussian
    Boltzmann distribution of the linearly interpolated potential.
    """

    if not isinstance(potential0, HarmonicPotential):
        raise TypeError("potential0 must be a HarmonicPotential.")
    if not isinstance(potential1, HarmonicPotential):
        raise TypeError("potential1 must be a HarmonicPotential.")
    temperature_value = float(temperature)
    if not np.isfinite(temperature_value) or temperature_value <= 0.0:
        raise ValueError("temperature must be finite and positive.")
    if isinstance(n_windows, bool) or not isinstance(n_windows, (int, np.integer)):
        raise ValueError("n_windows must be an integer.")
    if int(n_windows) < 2:
        raise ValueError("n_windows must be at least 2.")
    n_windows = int(n_windows)
    if method == "simpson" and n_windows % 2 == 0:
        raise ValueError("Simpson's rule requires an odd n_windows.")
    if isinstance(n_samples, bool) or not isinstance(n_samples, (int, np.integer)):
        raise ValueError("n_samples must be an integer.")
    if int(n_samples) < 2:
        raise ValueError("n_samples must be at least 2.")
    n_samples = int(n_samples)

    rt = BOLTZMANN_CONSTANT * temperature_value
    lambda_values = np.linspace(0.0, 1.0, n_windows)
    derivatives = np.empty(n_windows, dtype=float)
    standard_errors = np.empty(n_windows, dtype=float)
    rng = np.random.default_rng(seed)

    for index, lambda_value in enumerate(lambda_values):
        stiffness = (
            (1.0 - lambda_value) * potential0.stiffness
            + lambda_value * potential1.stiffness
        )
        center = (
            (1.0 - lambda_value)
            * potential0.stiffness
            * potential0.center
            + lambda_value
            * potential1.stiffness
            * potential1.center
        ) / stiffness
        standard_deviation = np.sqrt(rt / stiffness)
        positions = rng.normal(
            loc=center,
            scale=standard_deviation,
            size=(n_samples, 3),
        )

        displacement0 = positions - potential0.center
        displacement1 = positions - potential1.center
        work = (
            0.5
            * potential1.stiffness
            * np.sum(displacement1 * displacement1, axis=1)
            - 0.5
            * potential0.stiffness
            * np.sum(displacement0 * displacement0, axis=1)
        )
        derivatives[index] = float(np.mean(work))
        standard_errors[index] = float(np.std(work, ddof=1) / np.sqrt(n_samples))

    return integrate_thermodynamic_derivative(
        lambda_values,
        derivatives,
        standard_errors,
        method=method,
    )
