"""Tests for soft-core Lennard-Jones alchemical free-energy estimators."""

import numpy as np
import pytest

from moleculab import SoftCoreLennardJones


def model(**overrides):
    values = {
        "epsilon": 1.0,
        "sigma": 0.3,
        "trap_stiffness": 20.0,
        "temperature": 300.0,
        "alpha": 0.5,
        "r_max": 2.0,
        "n_grid": 100_000,
    }
    values.update(overrides)
    return SoftCoreLennardJones(**values)


def test_soft_core_limits_and_reference_convergence():
    pair = model()
    radius = np.array([0.4])

    assert pair.soft_core_energy(0.0, radius)[0] == pytest.approx(0.0)
    expected_lj = 4.0 * 1.0 * ((0.3 / 0.4) ** 12 - (0.3 / 0.4) ** 6)
    assert pair.soft_core_energy(1.0, radius)[0] == pytest.approx(expected_lj)

    coarse = model(n_grid=50_000).reference_delta_g()
    fine = model(n_grid=200_000).reference_delta_g()
    assert fine == pytest.approx(coarse, abs=1e-5)


def test_soft_core_derivative_matches_finite_difference():
    pair = model()
    radius = np.array([0.32, 0.4, 0.6])
    lambda_value = 0.4
    h = 1e-6

    numerical = (
        pair.soft_core_energy(lambda_value + h, radius)
        - pair.soft_core_energy(lambda_value - h, radius)
    ) / (2.0 * h)
    analytical = pair.soft_core_derivative(lambda_value, radius)

    np.testing.assert_allclose(analytical, numerical, atol=1e-6)


def test_sampling_returns_finite_positive_radii_and_positions():
    pair = model()
    radii = pair.sample_radii(0.5, 2000, seed=1)
    positions = pair.sample_positions(0.5, 2000, seed=1)

    assert radii.shape == (2000,)
    assert positions.shape == (2000, 3)
    assert np.all(radii > 0.0)
    assert np.all(np.isfinite(radii))
    np.testing.assert_allclose(
        np.linalg.norm(positions, axis=1),
        radii,
        rtol=1e-12,
    )


def test_multistage_fep_matches_reference():
    pair = model()
    reference = pair.reference_delta_g()

    result = pair.fep(n_windows=11, n_samples=50_000, seed=101)

    assert result.delta_g == pytest.approx(reference, abs=0.012)
    assert result.standard_error > 0.0
    assert result.stage_delta_g.shape == (10,)
    assert result.lambda_values.shape == (11,)


def test_ti_matches_reference():
    pair = model()
    reference = pair.reference_delta_g()

    result = pair.ti(n_windows=21, n_samples=50_000, seed=202)

    assert result.delta_g == pytest.approx(reference, abs=0.012)
    assert result.standard_error > 0.0
    assert result.weights.shape == (21,)


def test_mbar_matches_reference():
    pair = model()
    reference = pair.reference_delta_g()

    result = pair.mbar(n_windows=11, n_samples=20_000, seed=303)

    assert result.delta_g[-1] == pytest.approx(reference, abs=0.015)
    assert result.converged
    np.testing.assert_array_equal(result.state_counts, [20_000] * 11)


def test_fep_ti_mbar_are_consistent():
    pair = model()

    fep = pair.fep(n_windows=11, n_samples=30_000, seed=404)
    ti = pair.ti(n_windows=21, n_samples=30_000, seed=505)
    mbar = pair.mbar(n_windows=11, n_samples=10_000, seed=606)

    values = np.array([fep.delta_g, ti.delta_g, mbar.delta_g[-1]])
    assert np.ptp(values) < 0.02


def test_soft_core_validation():
    with pytest.raises(ValueError, match="epsilon"):
        model(epsilon=0.0)
    with pytest.raises(ValueError, match="sigma"):
        model(sigma=0.0)
    with pytest.raises(ValueError, match="temperature"):
        model(temperature=0.0)
    with pytest.raises(ValueError, match="alpha"):
        model(alpha=0.0)
    with pytest.raises(ValueError, match="n_grid"):
        model(n_grid=10)

    pair = model()
    with pytest.raises(ValueError, match="lambda"):
        pair.soft_core_energy(-0.1, [0.4])
    with pytest.raises(ValueError, match="radius"):
        pair.soft_core_energy(0.5, [0.0])
    with pytest.raises(ValueError, match="n_windows"):
        pair.fep(n_windows=1)
    with pytest.raises(ValueError, match="n_samples"):
        pair.ti(n_samples=1)
    with pytest.raises(ValueError, match="odd"):
        pair.ti(n_windows=10)
