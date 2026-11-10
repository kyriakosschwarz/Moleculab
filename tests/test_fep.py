"""Tests for free-energy perturbation on a harmonic oscillator."""

import numpy as np
import pytest

from moleculab import (
    HarmonicPotential,
    fep_between_potentials,
    fep_harmonic_oscillator,
    free_energy_perturbation,
    harmonic_free_energy_difference,
    sample_harmonic_positions,
    units,
)


def potentials():
    return (
        HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0]),
        HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0]),
    )


def test_constant_work_returns_the_constant():
    result = free_energy_perturbation([5.0, 5.0, 5.0, 5.0], temperature=300.0)

    assert result.delta_g == pytest.approx(5.0)
    assert result.standard_error == pytest.approx(0.0)
    assert result.n_samples == 4


def test_zero_work_has_zero_free_energy_change():
    result = free_energy_perturbation(np.zeros(100), temperature=300.0)

    assert result.delta_g == pytest.approx(0.0)
    assert result.standard_error == pytest.approx(0.0)


def test_harmonic_free_energy_difference_matches_analytic_formula():
    potential0, potential1 = potentials()
    temperature = 300.0

    expected = (
        1.5
        * units.BOLTZMANN_CONSTANT
        * temperature
        * np.log(potential1.stiffness / potential0.stiffness)
    )

    assert harmonic_free_energy_difference(
        potential0, potential1, temperature
    ) == pytest.approx(expected)


def test_identical_harmonic_potentials_have_zero_free_energy_difference():
    potential = HarmonicPotential(stiffness=100.0, center=[0.2, -0.1, 0.0])

    assert harmonic_free_energy_difference(potential, potential, 300.0) == pytest.approx(
        0.0
    )


def test_harmonic_center_shift_does_not_change_the_analytic_free_energy():
    potential0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    potential1 = HarmonicPotential(stiffness=100.0, center=[0.5, 0.0, 0.0])

    assert harmonic_free_energy_difference(potential0, potential1, 300.0) == pytest.approx(
        0.0
    )


def test_fep_harmonic_estimate_matches_analytic_result():
    potential0, potential1 = potentials()
    temperature = 300.0
    exact = harmonic_free_energy_difference(potential0, potential1, temperature)

    result = fep_harmonic_oscillator(
        potential0,
        potential1,
        temperature,
        n_samples=100_000,
        seed=123,
    )

    assert result.delta_g == pytest.approx(exact, abs=0.08)
    assert result.work.shape == (100_000,)
    assert result.n_samples == 100_000
    assert result.standard_error > 0.0


def test_fep_generic_potentials_path_agrees_with_harmonic_path():
    potential0, potential1 = potentials()
    temperature = 300.0
    positions = sample_harmonic_positions(
        potential0,
        temperature,
        n_samples=50_000,
        seed=456,
    )

    generic = fep_between_potentials(
        positions,
        potential0,
        potential1,
        temperature,
    )
    exact = harmonic_free_energy_difference(potential0, potential1, temperature)

    assert generic.delta_g == pytest.approx(exact, abs=0.12)
    assert generic.work.shape == (50_000,)


def test_harmonic_samples_have_expected_mean_and_variance():
    potential = HarmonicPotential(stiffness=100.0, center=[0.1, -0.2, 0.3])
    temperature = 300.0
    positions = sample_harmonic_positions(
        potential,
        temperature,
        n_samples=100_000,
        seed=789,
    )

    expected_variance = units.BOLTZMANN_CONSTANT * temperature / potential.stiffness
    np.testing.assert_allclose(
        np.mean(positions, axis=0),
        potential.center,
        atol=0.01,
    )
    np.testing.assert_allclose(
        np.var(positions, axis=0),
        expected_variance,
        rtol=0.05,
    )


def test_forward_and_reverse_fep_are_approximately_opposite():
    potential0, potential1 = potentials()
    temperature = 300.0

    forward = fep_harmonic_oscillator(
        potential0, potential1, temperature, n_samples=100_000, seed=1
    )
    reverse = fep_harmonic_oscillator(
        potential1, potential0, temperature, n_samples=100_000, seed=2
    )

    assert forward.delta_g == pytest.approx(-reverse.delta_g, abs=0.15)


def test_confidence_interval_contains_the_estimate():
    result = free_energy_perturbation([1.0, 2.0, 3.0, 4.0], temperature=300.0)

    lower, upper = result.confidence_interval

    assert lower < result.delta_g < upper


def test_fep_is_numerically_stable_for_large_work_values():
    result = free_energy_perturbation([1.0e6, 1.0e6], temperature=300.0)

    assert np.isfinite(result.delta_g)
    assert result.delta_g == pytest.approx(1.0e6)


def test_fep_validates_work_and_temperature():
    with pytest.raises(ValueError, match="one-dimensional"):
        free_energy_perturbation([[1.0, 2.0]], temperature=300.0)
    with pytest.raises(ValueError, match="at least two"):
        free_energy_perturbation([1.0], temperature=300.0)
    with pytest.raises(ValueError, match="finite"):
        free_energy_perturbation([1.0, np.nan], temperature=300.0)
    with pytest.raises(ValueError, match="temperature"):
        free_energy_perturbation([1.0, 2.0], temperature=0.0)


def test_sampling_and_fep_validate_inputs():
    potential0, potential1 = potentials()

    with pytest.raises(TypeError, match="HarmonicPotential"):
        sample_harmonic_positions("bad", 300.0, 10)
    with pytest.raises(ValueError, match="n_samples"):
        sample_harmonic_positions(potential0, 300.0, 1)
    with pytest.raises(ValueError, match="shape"):
        fep_between_potentials(
            np.zeros((10, 2)), potential0, potential1, 300.0
        )
    with pytest.raises(TypeError, match="HarmonicPotential"):
        fep_harmonic_oscillator("bad", potential1, 300.0, n_samples=10)
