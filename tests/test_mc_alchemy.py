"""Tests for multi-particle Monte Carlo alchemical free energies."""

import numpy as np
import pytest

from moleculab import BondedLennardJonesMC, SoftCoreLennardJonesPair


def model(**overrides):
    values = {
        "epsilon": 1.0,
        "sigma": 0.3,
        "bond_stiffness": 1000.0,
        "bond_length": 0.4,
        "temperature": 300.0,
        "alpha": 0.5,
        "mass": 1.0,
        "max_displacement": 0.015,
        "equilibration_steps": 1000,
        "sample_every": 5,
        "reference_r_max": 2.0,
        "reference_grid": 100_000,
    }
    values.update(overrides)
    return BondedLennardJonesMC(**values)


def test_soft_core_pair_limits_and_derivative():
    zero = SoftCoreLennardJonesPair(1.0, 0.3, 0.5, 0.0)
    full = SoftCoreLennardJonesPair(1.0, 0.3, 0.5, 1.0)
    displacement = np.array([0.4, 0.0, 0.0])

    assert zero.energy(displacement) == pytest.approx(0.0)
    expected = 4.0 * 1.0 * ((0.3 / 0.4) ** 12 - (0.3 / 0.4) ** 6)
    assert full.energy(displacement) == pytest.approx(expected)

    h = 1e-6
    for index in range(3):
        plus = displacement.copy()
        minus = displacement.copy()
        plus[index] += h
        minus[index] -= h
        numerical_force = -(full.energy(plus) - full.energy(minus)) / (2.0 * h)
        assert full.force(displacement)[index] == pytest.approx(
            numerical_force, abs=1e-6
        )


def test_soft_core_derivative_matches_finite_difference():
    radius = np.array([0.32, 0.4, 0.6])
    pair = SoftCoreLennardJonesPair(1.0, 0.3, 0.5, 0.4)
    h = 1e-6
    up = SoftCoreLennardJonesPair(1.0, 0.3, 0.5, 0.4 + h)
    down = SoftCoreLennardJonesPair(1.0, 0.3, 0.5, 0.4 - h)

    numerical = (
        up.energy_from_radius(radius) - down.energy_from_radius(radius)
    ) / (2.0 * h)
    np.testing.assert_allclose(
        pair.derivative_from_radius(radius),
        numerical,
        atol=1e-6,
    )


def test_mc_sample_positions_are_finite_and_bonded():
    pair = model()
    positions = pair.sample_positions(0.5, n_samples=500, seed=1)

    assert positions.shape == (500, 2, 3)
    assert np.all(np.isfinite(positions))
    radii = np.linalg.norm(positions[:, 0] - positions[:, 1], axis=1)
    assert np.all(radii > 0.0)
    assert np.mean(radii) == pytest.approx(pair.bond_length, abs=0.05)
    assert np.std(radii) > 1e-3


def test_reference_converges_with_grid():
    coarse = model(reference_grid=50_000).reference_delta_g()
    fine = model(reference_grid=200_000).reference_delta_g()

    assert fine == pytest.approx(coarse, abs=1e-5)


def test_fep_ti_mbar_match_reference():
    pair = model()
    reference = pair.reference_delta_g()

    fep = pair.fep(n_windows=11, n_samples=1000, seed=101)
    ti = pair.ti(n_windows=21, n_samples=1000, seed=202)
    mbar = pair.mbar(n_windows=11, n_samples=1000, seed=303)

    assert fep.delta_g == pytest.approx(reference, abs=0.02)
    assert ti.delta_g == pytest.approx(reference, abs=0.02)
    assert mbar.delta_g[-1] == pytest.approx(reference, abs=0.02)
    assert fep.standard_error > 0.0
    assert ti.standard_error > 0.0
    assert mbar.converged


def test_fep_ti_mbar_are_consistent():
    pair = model()

    values = np.array(
        [
            pair.fep(n_windows=11, n_samples=1000, seed=404).delta_g,
            pair.ti(n_windows=21, n_samples=1000, seed=505).delta_g,
            pair.mbar(n_windows=11, n_samples=1000, seed=606).delta_g[-1],
        ]
    )

    assert np.ptp(values) < 0.03


def test_make_system_contains_bond_and_pair_interaction():
    pair = model()
    system = pair.make_system(0.5)

    assert system.n_particles == 2
    assert len(system.pair_potentials) == 1
    assert system.topology is not None
    assert system.topology.n_bonds == 1
    assert np.isfinite(system.potential_energy())
    assert np.all(np.isfinite(system.forces()))


def test_alchemy_validation():
    with pytest.raises(ValueError, match="epsilon"):
        model(epsilon=0.0)
    with pytest.raises(ValueError, match="temperature"):
        model(temperature=0.0)
    with pytest.raises(ValueError, match="max_displacement"):
        model(max_displacement=0.0)
    with pytest.raises(ValueError, match="n_samples"):
        model().sample_positions(0.5, n_samples=0)
    with pytest.raises(ValueError, match="n_windows"):
        model().fep(n_windows=1)
    with pytest.raises(ValueError, match="odd"):
        model().ti(n_windows=10)
    with pytest.raises(ValueError, match="lambda"):
        SoftCoreLennardJonesPair(1.0, 0.3, 0.5, 1.1)
