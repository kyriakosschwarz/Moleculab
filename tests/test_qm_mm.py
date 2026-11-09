"""Tests for the finite-difference QM/MM particle-in-a-box model."""

import numpy as np
import pytest

from moleculab import ELECTRON_MASS_DA, MMPointCharge, ParticleInBox
from moleculab.units import (
    AVOGADRO_CONSTANT,
    COULOMB_CONSTANT,
    DALTON_IN_KILOGRAMS,
    H_BAR_J_S,
)


def test_finite_difference_matches_analytic_particle_in_a_box():
    box = ParticleInBox(length=1.0, mass=ELECTRON_MASS_DA, n_grid=400)

    result = box.solve(n_states=3)
    exact = box.analytic_energies(3)

    np.testing.assert_allclose(result.energies, exact, atol=0.05)
    assert result.wavefunctions.shape == (400, 3)
    assert result.probability_density.shape == (3, 400)


def test_energy_converges_with_grid_resolution():
    coarse = ParticleInBox(1.0, ELECTRON_MASS_DA, n_grid=100)
    fine = ParticleInBox(1.0, ELECTRON_MASS_DA, n_grid=400)
    exact = fine.analytic_energies(3)

    coarse_error = np.max(np.abs(coarse.solve(n_states=3).energies - exact))
    fine_error = np.max(np.abs(fine.solve(n_states=3).energies - exact))

    assert fine_error < coarse_error
    assert fine_error < 0.05


def test_constant_potential_shifts_all_levels():
    box = ParticleInBox(1.0, ELECTRON_MASS_DA, n_grid=400)
    shift = 10.0
    potential = np.full(box.n_grid, shift)

    result = box.solve(potential=potential, n_states=3)
    exact = box.analytic_energies(3) + shift

    np.testing.assert_allclose(result.energies, exact, atol=0.05)


def test_mm_potential_matches_point_charge_formula():
    box = ParticleInBox(1.0, ELECTRON_MASS_DA, n_grid=10)
    charge = MMPointCharge(position=-0.2, charge=0.5)

    potential = box.mm_potential([charge], qm_charge=-1.0)

    expected = COULOMB_CONSTANT * (-1.0) * 0.5 / np.abs(box.grid[0] + 0.2)
    assert potential[0] == pytest.approx(expected)


def test_qmmm_weak_charge_matches_first_order_perturbation():
    box = ParticleInBox(1.0, ELECTRON_MASS_DA, n_grid=800)
    charge = MMPointCharge(position=-0.2, charge=0.01)
    exact = box.solve(
        mm_charges=[charge],
        qm_charge=-1.0,
        n_states=3,
    )
    baseline = box.solve(n_states=3)
    shifts = exact.energies - baseline.energies

    potential = box.mm_potential([charge], qm_charge=-1.0)
    first_order = []
    for n in range(1, 4):
        wavefunction = np.sqrt(2.0 / box.length) * np.sin(
            n * np.pi * box.grid / box.length
        )
        first_order.append(
            np.trapezoid(wavefunction**2 * potential, box.grid)
        )

    np.testing.assert_allclose(shifts, first_order, atol=0.02)


def test_harmonic_well_matches_analytic_oscillator():
    length = 4.0
    stiffness = 100.0
    box = ParticleInBox(length, ELECTRON_MASS_DA, n_grid=800)
    potential = box.harmonic_potential(stiffness)

    result = box.solve(potential=potential, n_states=3)

    stiffness_si = stiffness * 1000.0 / AVOGADRO_CONSTANT * 1.0e18
    omega = np.sqrt(stiffness_si / (ELECTRON_MASS_DA * DALTON_IN_KILOGRAMS))
    exact = (
        (np.arange(3) + 0.5) * H_BAR_J_S * omega * AVOGADRO_CONSTANT / 1000.0
    )

    np.testing.assert_allclose(result.energies, exact, atol=0.05)


def test_mm_charge_on_grid_is_rejected():
    box = ParticleInBox(1.0, ELECTRON_MASS_DA, n_grid=10)
    on_grid = MMPointCharge(position=float(box.grid[0]), charge=1.0)

    with pytest.raises(ValueError, match="coincides"):
        box.mm_potential([on_grid], qm_charge=-1.0)


def test_qmmm_validation():
    with pytest.raises(ValueError, match="length"):
        ParticleInBox(0.0, ELECTRON_MASS_DA)
    with pytest.raises(ValueError, match="mass"):
        ParticleInBox(1.0, 0.0)
    with pytest.raises(ValueError, match="n_grid"):
        ParticleInBox(1.0, ELECTRON_MASS_DA, n_grid=5)

    box = ParticleInBox(1.0, ELECTRON_MASS_DA, n_grid=20)
    with pytest.raises(ValueError, match="n_states"):
        box.solve(n_states=0)
    with pytest.raises(ValueError, match="n_states"):
        box.solve(n_states=21)
    with pytest.raises(ValueError, match="potential"):
        box.solve(potential=np.zeros(10))
    with pytest.raises(ValueError, match="finite"):
        box.solve(potential=np.full(box.n_grid, np.nan))
    with pytest.raises(ValueError, match="position"):
        MMPointCharge(float("nan"), 1.0)
    with pytest.raises(ValueError, match="charge"):
        MMPointCharge(0.0, float("inf"))
