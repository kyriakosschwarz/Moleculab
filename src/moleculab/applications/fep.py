# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Free-energy perturbation (FEP) for a simple analytically tractable system.

FEP estimates the free-energy difference between two states, 0 and 1, by
sampling from state 0 and reweighting each sample by the work needed to change
the Hamiltonian:

.. math::

    \\Delta G = -k_B T \\ln
    \\left\\langle \\exp\\left(-\\frac{\\Delta U}{k_B T}\\right)
    \\right\\rangle_0,

where :math:`\\Delta U = U_1 - U_0`.  In Moleculab's molar units,
:math:`k_B T` is replaced by :math:`RT`, with
``R = 0.008314462618 kJ/(mol K)``.

The test system is a three-dimensional harmonic oscillator whose stiffness is
changed from :math:`k_0` to :math:`k_1`.  Its partition function is known
exactly, so the FEP estimate can be compared with

.. math::

    \\Delta G_{\\mathrm{exact}} =
    \\frac{3}{2} RT \\ln\\left(\\frac{k_1}{k_0}\\right).

This is the same exponential-averaging idea used in alchemical free-energy
calculations, but with an exact reference result.  It is a safe place to learn
the statistical behavior of FEP before applying it to a ligand in a protein.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..foundation.units import BOLTZMANN_CONSTANT
from ..models.potentials import HarmonicPotential, Potential


@dataclass(frozen=True)
class FEPResult:
    """Result of a free-energy perturbation estimate.

    Attributes
    ----------
    delta_g:
        Estimated free-energy difference in kJ/mol.
    standard_error:
        Standard error of the estimate in kJ/mol.
    work:
        Work values :math:`\\Delta U = U_1 - U_0` for all samples, in kJ/mol.
    temperature:
        Temperature in K.
    n_samples:
        Number of samples used.
    """

    delta_g: float
    standard_error: float
    work: np.ndarray
    temperature: float
    n_samples: int

    @property
    def confidence_interval(self) -> tuple[float, float]:
        """Approximate 95% confidence interval in kJ/mol."""

        margin = 1.96 * self.standard_error
        return (self.delta_g - margin, self.delta_g + margin)


def _validate_temperature(temperature: float) -> float:
    value = float(temperature)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("temperature must be finite and positive.")
    return value


def _validate_samples(work: np.ndarray) -> np.ndarray:
    values = np.asarray(work, dtype=float)
    if values.ndim != 1:
        raise ValueError(f"work must be one-dimensional, got shape {values.shape}.")
    if values.size < 2:
        raise ValueError("FEP requires at least two work samples for an error estimate.")
    if not np.all(np.isfinite(values)):
        raise ValueError("work must contain only finite values.")
    return values


def _logsumexp(values: np.ndarray) -> float:
    maximum = float(np.max(values))
    return maximum + float(np.log(np.sum(np.exp(values - maximum))))


def free_energy_perturbation(
    work: np.ndarray,
    temperature: float,
) -> FEPResult:
    """Estimate :math:`\\Delta G` from precomputed work samples.

    ``work`` contains :math:`\\Delta U = U_1 - U_0` in kJ/mol for samples drawn
    from state 0.  The exponential average is evaluated in log space so that
    large work values do not overflow.
    """

    values = _validate_samples(work)
    temperature_value = _validate_temperature(temperature)
    rt = BOLTZMANN_CONSTANT * temperature_value
    beta = 1.0 / rt
    n_samples = values.size

    log_weights = -beta * values
    log_mean = _logsumexp(log_weights) - np.log(n_samples)
    delta_g = float(-log_mean / beta)

    # The same exponential weights, shifted by their maximum for stability.
    weights = np.exp(log_weights - np.max(log_weights))
    mean_weight = float(np.mean(weights))
    variance_weight = float(np.var(weights, ddof=1))
    standard_error = float(
        np.sqrt(variance_weight / n_samples) / mean_weight / beta
    )

    return FEPResult(
        delta_g=delta_g,
        standard_error=standard_error,
        work=values,
        temperature=temperature_value,
        n_samples=n_samples,
    )


def fep_between_potentials(
    samples: np.ndarray,
    potential0: Potential,
    potential1: Potential,
    temperature: float,
) -> FEPResult:
    """Estimate :math:`\\Delta G` from samples of state 0.

    ``samples`` has shape ``(n_samples, 3)`` and contains positions drawn from
    the Boltzmann distribution of ``potential0``.  This generic helper works
    for any objects implementing ``energy(position)``.
    """

    coordinates = np.asarray(samples, dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError(
            "samples must have shape (n_samples, 3), "
            f"got {coordinates.shape}."
        )
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("samples must contain only finite values.")

    work = np.array(
        [
            potential1.energy(position) - potential0.energy(position)
            for position in coordinates
        ],
        dtype=float,
    )
    return free_energy_perturbation(work, temperature)


def sample_harmonic_positions(
    potential: HarmonicPotential,
    temperature: float,
    n_samples: int,
    seed: int | None = None,
) -> np.ndarray:
    """Draw positions from a harmonic oscillator's Boltzmann distribution.

    For :math:`U = \\frac{1}{2}k|\\mathbf r - \\mathbf r_0|^2`, each Cartesian
    coordinate is Gaussian with mean :math:`\\mathbf r_0` and variance
    :math:`RT/k`.

    Returns an array with shape ``(n_samples, 3)``.
    """

    if not isinstance(potential, HarmonicPotential):
        raise TypeError("potential must be a HarmonicPotential.")
    temperature_value = _validate_temperature(temperature)
    if isinstance(n_samples, bool) or not isinstance(n_samples, (int, np.integer)):
        raise ValueError("n_samples must be an integer.")
    if int(n_samples) < 2:
        raise ValueError("n_samples must be at least 2.")
    n_samples = int(n_samples)

    standard_deviation = np.sqrt(
        BOLTZMANN_CONSTANT * temperature_value / potential.stiffness
    )
    rng = np.random.default_rng(seed)
    return rng.normal(
        loc=potential.center,
        scale=standard_deviation,
        size=(n_samples, 3),
    )


def harmonic_free_energy_difference(
    potential0: HarmonicPotential,
    potential1: HarmonicPotential,
    temperature: float,
) -> float:
    """Return the exact harmonic free-energy difference in kJ/mol.

    The result is independent of the potential centers because a center shift
    changes neither the width nor the normalization of the Gaussian
    distribution.
    """

    if not isinstance(potential0, HarmonicPotential):
        raise TypeError("potential0 must be a HarmonicPotential.")
    if not isinstance(potential1, HarmonicPotential):
        raise TypeError("potential1 must be a HarmonicPotential.")
    temperature_value = _validate_temperature(temperature)
    rt = BOLTZMANN_CONSTANT * temperature_value
    return float(
        1.5 * rt * np.log(potential1.stiffness / potential0.stiffness)
    )


def fep_harmonic_oscillator(
    potential0: HarmonicPotential,
    potential1: HarmonicPotential,
    temperature: float,
    n_samples: int = 10_000,
    seed: int | None = None,
) -> FEPResult:
    """Run FEP for a harmonic stiffness change.

    Samples are drawn from ``potential0`` and reweighted to ``potential1``.
    The work is computed directly from the harmonic energies, so this function
    is fast enough for convergence experiments.
    """

    if not isinstance(potential0, HarmonicPotential):
        raise TypeError("potential0 must be a HarmonicPotential.")
    if not isinstance(potential1, HarmonicPotential):
        raise TypeError("potential1 must be a HarmonicPotential.")
    temperature_value = _validate_temperature(temperature)
    if isinstance(n_samples, bool) or not isinstance(n_samples, (int, np.integer)):
        raise ValueError("n_samples must be an integer.")
    if int(n_samples) < 2:
        raise ValueError("n_samples must be at least 2.")
    n_samples = int(n_samples)

    positions = sample_harmonic_positions(
        potential0,
        temperature_value,
        n_samples,
        seed=seed,
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
    return free_energy_perturbation(work, temperature_value)
