"""Tests for thermodynamic integration on a harmonic alchemical path."""

import numpy as np
import pytest

from moleculab import (
    HarmonicPotential,
    harmonic_free_energy_difference,
    integrate_thermodynamic_derivative,
    thermodynamic_integration_harmonic,
)


def harmonic_pair():
    return (
        HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0]),
        HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0]),
    )


def test_simpson_integrates_a_known_function():
    lambda_values = np.linspace(0.0, 1.0, 5)
    derivatives = lambda_values**2

    result = integrate_thermodynamic_derivative(
        lambda_values,
        derivatives,
        method="simpson",
    )

    assert result.delta_g == pytest.approx(1.0 / 3.0)
    assert np.sum(result.weights) == pytest.approx(1.0)
    assert result.method == "simpson"


def test_trapezoid_weights_sum_to_one():
    lambda_values = np.linspace(0.0, 1.0, 6)
    result = integrate_thermodynamic_derivative(
        lambda_values,
        lambda_values,
        method="trapezoid",
    )

    assert np.sum(result.weights) == pytest.approx(1.0)
    assert result.delta_g == pytest.approx(0.5)


def test_ti_harmonic_matches_analytic_result():
    potential0, potential1 = harmonic_pair()
    temperature = 300.0
    exact = harmonic_free_energy_difference(potential0, potential1, temperature)

    result = thermodynamic_integration_harmonic(
        potential0,
        potential1,
        temperature,
        n_windows=21,
        n_samples=20_000,
        seed=123,
    )

    assert result.delta_g == pytest.approx(exact, abs=0.05)
    assert result.standard_error > 0.0
    assert result.lambda_values.shape == (21,)
    assert result.derivatives.shape == (21,)
    assert result.standard_errors.shape == (21,)
    assert np.sum(result.weights) == pytest.approx(1.0)


def test_ti_of_identical_potentials_is_zero():
    potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])

    result = thermodynamic_integration_harmonic(
        potential,
        potential,
        300.0,
        n_windows=11,
        n_samples=5_000,
        seed=1,
    )

    assert result.delta_g == pytest.approx(0.0, abs=1e-12)
    np.testing.assert_allclose(result.derivatives, np.zeros(11), atol=1e-12)


def test_ti_center_shift_has_zero_free_energy_difference():
    potential0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    potential1 = HarmonicPotential(stiffness=100.0, center=[0.2, 0.0, 0.0])

    result = thermodynamic_integration_harmonic(
        potential0,
        potential1,
        300.0,
        n_windows=21,
        n_samples=50_000,
        seed=7,
    )

    assert result.delta_g == pytest.approx(0.0, abs=0.1)


def test_ti_improves_with_more_windows():
    potential0, potential1 = harmonic_pair()
    exact = harmonic_free_energy_difference(potential0, potential1, 300.0)

    coarse = thermodynamic_integration_harmonic(
        potential0,
        potential1,
        300.0,
        n_windows=5,
        n_samples=20_000,
        seed=5,
    )
    fine = thermodynamic_integration_harmonic(
        potential0,
        potential1,
        300.0,
        n_windows=21,
        n_samples=20_000,
        seed=5,
    )

    assert abs(fine.delta_g - exact) <= abs(coarse.delta_g - exact)


def test_ti_validation():
    potential0, potential1 = harmonic_pair()

    with pytest.raises(ValueError, match="temperature"):
        thermodynamic_integration_harmonic(potential0, potential1, 0.0)
    with pytest.raises(ValueError, match="n_windows"):
        thermodynamic_integration_harmonic(
            potential0, potential1, 300.0, n_windows=1
        )
    with pytest.raises(ValueError, match="odd"):
        thermodynamic_integration_harmonic(
            potential0, potential1, 300.0, n_windows=10
        )
    with pytest.raises(ValueError, match="n_samples"):
        thermodynamic_integration_harmonic(
            potential0, potential1, 300.0, n_samples=1
        )
    with pytest.raises(ValueError, match="method"):
        integrate_thermodynamic_derivative(
            np.linspace(0, 1, 5),
            np.ones(5),
            method="midpoint",
        )
    with pytest.raises(ValueError, match="odd"):
        integrate_thermodynamic_derivative(
            np.linspace(0, 1, 4),
            np.ones(4),
            method="simpson",
        )
    with pytest.raises(ValueError, match="uniform"):
        integrate_thermodynamic_derivative(
            np.array([0.0, 0.2, 0.5, 1.0]),
            np.ones(4),
            method="trapezoid",
        )
    with pytest.raises(ValueError, match="0 to 1"):
        integrate_thermodynamic_derivative(
            np.linspace(0.0, 2.0, 5),
            np.ones(5),
            method="simpson",
        )
