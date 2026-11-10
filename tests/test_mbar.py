"""Tests for MBAR over harmonic alchemical states."""

import numpy as np
import pytest

from moleculab import (
    HarmonicPotential,
    harmonic_free_energy_difference,
    mbar,
    mbar_harmonic,
    units,
)


def harmonic_potentials():
    return [
        HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0]),
        HarmonicPotential(stiffness=150.0, center=[0.0, 0.0, 0.0]),
        HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0]),
    ]


def exact_delta_g(potentials, temperature):
    rt = units.BOLTZMANN_CONSTANT * temperature
    return np.array(
        [
            1.5 * rt * np.log(potential.stiffness / potentials[0].stiffness)
            for potential in potentials
        ]
    )


def test_mbar_harmonic_matches_analytic_free_energies():
    potentials = harmonic_potentials()
    temperature = 300.0
    exact = exact_delta_g(potentials, temperature)

    result = mbar_harmonic(
        potentials,
        temperature,
        n_samples=10_000,
        seed=1,
    )

    assert result.converged
    assert result.free_energies[0] == pytest.approx(0.0)
    np.testing.assert_allclose(result.delta_g, exact, atol=0.1)
    np.testing.assert_array_equal(result.state_counts, [10_000, 10_000, 10_000])
    assert result.delta_g.shape == (3,)
    assert result.free_energies.shape == (3,)


def test_mbar_endpoint_matches_analytic_and_fep_reference():
    potentials = harmonic_potentials()
    temperature = 300.0

    result = mbar_harmonic(
        potentials,
        temperature,
        n_samples=12_000,
        seed=2,
    )
    exact = harmonic_free_energy_difference(
        potentials[0], potentials[-1], temperature
    )

    assert result.delta_g[-1] == pytest.approx(exact, abs=0.1)


def test_mbar_identical_states_have_zero_free_energy_difference():
    potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])

    result = mbar_harmonic(
        [potential, potential],
        300.0,
        n_samples=2_000,
        seed=3,
    )

    np.testing.assert_allclose(result.delta_g, [0.0, 0.0], atol=1e-12)


def test_mbar_is_reproducible_with_a_seed():
    potentials = harmonic_potentials()

    first = mbar_harmonic(potentials, 300.0, n_samples=3_000, seed=4)
    second = mbar_harmonic(potentials, 300.0, n_samples=3_000, seed=4)

    np.testing.assert_allclose(first.delta_g, second.delta_g)


def test_mbar_validation():
    potentials = harmonic_potentials()
    u = np.zeros((6, 3))
    states = np.array([0, 0, 1, 1, 2, 2])

    with pytest.raises(ValueError, match="at least two"):
        mbar_harmonic([potentials[0]], 300.0)
    with pytest.raises(TypeError, match="HarmonicPotential"):
        mbar_harmonic([potentials[0], "bad"], 300.0)
    with pytest.raises(ValueError, match="temperature"):
        mbar_harmonic(potentials, 0.0)
    with pytest.raises(ValueError, match="n_samples"):
        mbar_harmonic(potentials, 300.0, n_samples=0)
    with pytest.raises(ValueError, match="shape"):
        mbar(u[:, 0], states, 300.0)
    with pytest.raises(ValueError, match="sample_states"):
        mbar(u, np.array([0, 0, 1, 1, 2]), 300.0)
    with pytest.raises(ValueError, match="integer"):
        mbar(u, np.array([0.0, 0.0, 1.0, 1.0, 2.0, 2.0]), 300.0)
    with pytest.raises(ValueError, match="out-of-range"):
        mbar(u, np.array([0, 0, 1, 1, 2, 3]), 300.0)
    with pytest.raises(ValueError, match="at least one sample"):
        mbar(u, np.array([0, 0, 1, 1, 1, 1]), 300.0)
    with pytest.raises(ValueError, match="finite"):
        mbar(np.full((6, 3), np.nan), states, 300.0)
    with pytest.raises(ValueError, match="damping"):
        mbar(u, states, 300.0, damping=0.0)


def test_mbar_raises_on_nonconvergence():
    u = np.random.default_rng(0).normal(size=(20, 2))
    states = np.array([0] * 10 + [1] * 10)

    with pytest.raises(RuntimeError, match="converge"):
        mbar(u, states, 300.0, max_iterations=1, damping=1.0)
