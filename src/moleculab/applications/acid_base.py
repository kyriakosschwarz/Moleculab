# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Acid-base equilibria and discrete protonation-state populations.

This module does the thermodynamic bookkeeping that connects a pH value to
the populations of protonation and tautomeric states. It does **not** predict
pKa values. A pKa is an input: it may come from experiment, from an empirical
model, or from a calculation for a particular environment.

The central equations are the Henderson-Hasselbalch relation and the Boltzmann
distribution. All energies are in kJ/mol and temperatures in K, matching the
rest of Moleculab.
"""

from __future__ import annotations

import numpy as np

from ..foundation.units import MOLAR_GAS_CONSTANT

_LN10 = float(np.log(10.0))


def _finite_float(name: str, value: float) -> float:
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(f"{name} must be finite.")
    return result


def _logistic(value: float) -> float:
    """Return ``1 / (1 + exp(-value))`` without overflow."""

    if value >= 0.0:
        return float(1.0 / (1.0 + np.exp(-value)))
    exponential = float(np.exp(value))
    return float(exponential / (1.0 + exponential))


def deprotonated_fraction(pka: float, ph: float) -> float:
    """Return the fraction of an acid that is deprotonated at *ph*.

    For an acid with dissociation constant :math:`K_a`,

    .. math::

        f_{A^-} = \\frac{1}{1 + 10^{pK_a - pH}}.

    At ``ph == pka`` the result is exactly 0.5. The function uses a logistic
    form so that extreme pH values do not overflow.
    """

    pka_value = _finite_float("pka", pka)
    ph_value = _finite_float("ph", ph)
    return _logistic((ph_value - pka_value) * _LN10)


def protonated_fraction(pka: float, ph: float) -> float:
    """Return the fraction of an acid that remains protonated at *ph*.

    This is :math:`1 - f_{A^-}`, evaluated directly to avoid cancellation.
    For a base, the same expression gives the fraction of the conjugate acid
    :math:`BH^+`.
    """

    pka_value = _finite_float("pka", pka)
    ph_value = _finite_float("ph", ph)
    return _logistic((pka_value - ph_value) * _LN10)


def protonation_microstates(pka_values, ph: float) -> np.ndarray:
    """Return populations of all protonation states for independent sites.

    ``pka_values`` is a sequence of site pKa values. The result has length
    ``2**n`` and is indexed by a bit mask: bit ``i`` is 1 when site ``i`` is
    protonated and 0 when it is deprotonated. For independent sites the
    probability of a microstate is the product of the site fractions.

    The independent-site model ignores electrostatic interactions between
    sites. Use :func:`boltzmann_populations` with explicit microstate free
    energies when coupling matters.
    """

    pkas = np.asarray(pka_values, dtype=float).reshape(-1)
    if pkas.size == 0:
        raise ValueError("pka_values must contain at least one site.")
    if pkas.size > 16:
        raise ValueError("pka_values is limited to 16 independent sites.")
    if not np.all(np.isfinite(pkas)):
        raise ValueError("pka_values must be finite.")
    ph_value = _finite_float("ph", ph)

    n_sites = pkas.size
    probabilities = np.ones(1 << n_sites, dtype=float)
    for index, pka in enumerate(pkas):
        protonated = protonated_fraction(pka, ph_value)
        for state in range(1 << n_sites):
            if (state >> index) & 1:
                probabilities[state] *= protonated
            else:
                probabilities[state] *= 1.0 - protonated
    return probabilities


def mean_protonation(pka_values, ph: float) -> np.ndarray:
    """Return the mean protonated fraction of each independent site."""

    pkas = np.asarray(pka_values, dtype=float).reshape(-1)
    if pkas.size == 0:
        raise ValueError("pka_values must contain at least one site.")
    if not np.all(np.isfinite(pkas)):
        raise ValueError("pka_values must be finite.")
    ph_value = _finite_float("ph", ph)
    return np.array([protonated_fraction(pka, ph_value) for pka in pkas])


def mean_net_charge(
    pka_values,
    ph: float,
    protonated_charges,
    deprotonated_charges,
) -> float:
    """Return the mean net charge of independent ionizable sites.

    ``protonated_charges`` and ``deprotonated_charges`` give the formal charge
    of each site in its two states. For example, a carboxylic acid has
    protonated charge 0 and deprotonated charge -1, while an amine has
    protonated charge +1 and deprotonated charge 0.
    """

    pkas = np.asarray(pka_values, dtype=float).reshape(-1)
    if pkas.size == 0:
        raise ValueError("pka_values must contain at least one site.")
    if not np.all(np.isfinite(pkas)):
        raise ValueError("pka_values must be finite.")

    protonated = np.asarray(protonated_charges, dtype=float).reshape(-1)
    deprotonated = np.asarray(deprotonated_charges, dtype=float).reshape(-1)
    if protonated.size != pkas.size or deprotonated.size != pkas.size:
        raise ValueError("charge arrays must match pka_values in length.")
    if not np.all(np.isfinite(protonated)) or not np.all(np.isfinite(deprotonated)):
        raise ValueError("charges must be finite.")

    fractions = mean_protonation(pkas, ph)
    charges = fractions * protonated + (1.0 - fractions) * deprotonated
    return float(np.sum(charges))


def boltzmann_populations(
    state_energies,
    temperature: float = 298.15,
) -> np.ndarray:
    """Return Boltzmann populations for state energies in kJ/mol.

    Energies are relative free energies. An arbitrary constant offset cancels,
    so only energy differences matter. The temperature must be positive.

    This is appropriate for states with the same composition. For protonation
    equilibria, include the pH-dependent proton chemical potential in the
    state energies; see :func:`protonation_microstates`.
    """

    energies = np.asarray(state_energies, dtype=float).reshape(-1)
    if energies.size == 0:
        raise ValueError("state_energies must contain at least one state.")
    if not np.all(np.isfinite(energies)):
        raise ValueError("state_energies must be finite.")

    kelvin = _finite_float("temperature", temperature)
    if kelvin <= 0.0:
        raise ValueError("temperature must be positive.")

    rt = MOLAR_GAS_CONSTANT * kelvin
    shifted = energies - np.min(energies)
    weights = np.exp(-shifted / rt)
    return weights / np.sum(weights)


def tautomer_populations(
    delta_g: float,
    temperature: float = 298.15,
) -> np.ndarray:
    """Return populations of two tautomers separated by *delta_g* kJ/mol.

    ``delta_g`` is the free-energy difference ``G_B - G_A``. The result is
    ``[p_A, p_B]``. Both tautomers are assumed to have the same composition.
    For states that differ in proton count, use
    :func:`protonation_microstates` or pH-dependent microstate energies
    instead. This function does not describe the rate of tautomer
    interconversion.
    """

    difference = _finite_float("delta_g", delta_g)
    return boltzmann_populations([0.0, difference], temperature)
