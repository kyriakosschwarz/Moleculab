# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Independent physical benchmarks for constants and pair physics.

These tests are deliberately **not** written as "the code agrees with itself".
They use hardcoded literature values, hand-calculable closed forms, absolute
sign laws, and order-of-magnitude bounds.

The four layers are:

1. golden values for fundamental constants;
2. closed-form smoke tests at trivial geometries;
3. sign and direction invariants;
4. dimensional and order-of-magnitude bounds.

A broken operator in a constant formula should fail here even if every
self-consistency test still passes.
"""

import math

import numpy as np
import pytest

from moleculab import CoulombPotential, LennardJonesPotential, units


def test_golden_si_constants():
    """Fundamental constants against independent literature values."""

    assert units.AVOGADRO_CONSTANT == pytest.approx(6.02214076e23, rel=1e-12)
    assert units.BOLTZMANN_CONSTANT_J_PER_K == pytest.approx(1.380649e-23, rel=1e-12)
    assert units.ELEMENTARY_CHARGE_COULOMBS == pytest.approx(
        1.602176634e-19,
        rel=1e-12,
    )
    assert units.VACUUM_PERMITTIVITY_F_PER_M == pytest.approx(
        8.8541878128e-12,
        rel=1e-12,
    )
    assert units.PLANCK_CONSTANT_J_S == pytest.approx(6.62607015e-34, rel=1e-12)
    assert units.H_BAR_J_S == pytest.approx(1.054571817e-34, rel=1e-12)
    assert units.ELECTRON_MASS_KG == pytest.approx(9.1093837015e-31, rel=1e-12)
    assert units.DALTON_IN_KILOGRAMS == pytest.approx(
        1.66053906660e-27,
        rel=1e-12,
    )


def test_golden_moleculab_constants():
    """Moleculab-unit constants against independent reference values."""

    # R = 8.314462618 J/(mol K) = 0.008314462618 kJ/(mol K)
    assert units.MOLAR_GAS_CONSTANT == pytest.approx(0.008314462618, rel=1e-12)

    # k_e = (1/(4 pi eps0)) e^2 (1/nm) N_A / 1000
    assert units.COULOMB_CONSTANT == pytest.approx(138.935456, rel=1e-7)

    # 1 kJ/(mol nm^3) = 16.6054 bar
    assert units.KILOJOULES_PER_MOLE_PER_NM3_TO_BAR == pytest.approx(
        16.6054,
        rel=1e-5,
    )


def test_physical_constant_bounds():
    """Order-of-magnitude guards against algebraic corruption."""

    assert 100.0 < units.COULOMB_CONSTANT < 200.0
    assert 0.008 < units.MOLAR_GAS_CONSTANT < 0.009
    assert 6.0e23 < units.AVOGADRO_CONSTANT < 6.1e23
    assert 1.6e-27 < units.DALTON_IN_KILOGRAMS < 1.7e-27
    assert 1.6e-19 < units.ELEMENTARY_CHARGE_COULOMBS < 1.7e-19
    assert 8.0e-12 < units.VACUUM_PERMITTIVITY_F_PER_M < 1.0e-11


def test_coulomb_closed_form_at_one_nanometre():
    """Two unit charges at 1 nm must give exactly the Coulomb constant."""

    potential = CoulombPotential(charge_i=1.0, charge_j=1.0)

    assert potential.energy([1.0, 0.0, 0.0]) == pytest.approx(
        138.935456,
        rel=1e-7,
    )
    assert potential.force([1.0, 0.0, 0.0])[0] == pytest.approx(
        138.935456,
        rel=1e-7,
    )


def test_coulomb_closed_form_at_0_3_nanometre():
    """A hand-calculable contact-distance reference."""

    like = CoulombPotential(charge_i=1.0, charge_j=1.0)
    opposite = CoulombPotential(charge_i=1.0, charge_j=-1.0)

    assert like.energy([0.3, 0.0, 0.0]) == pytest.approx(
        463.118187,
        rel=1e-6,
    )
    assert opposite.energy([0.3, 0.0, 0.0]) == pytest.approx(
        -463.118187,
        rel=1e-6,
    )
    assert like.force([0.3, 0.0, 0.0])[0] == pytest.approx(
        1543.7273,
        rel=1e-6,
    )
    assert opposite.force([0.3, 0.0, 0.0])[0] == pytest.approx(
        -1543.7273,
        rel=1e-6,
    )


def test_coulomb_sign_and_direction_invariants():
    """Like charges repel; opposite charges attract; forces are central."""

    like = CoulombPotential(charge_i=1.0, charge_j=1.0)
    opposite = CoulombPotential(charge_i=1.0, charge_j=-1.0)
    displacement = np.array([0.7, 0.0, 0.0])

    assert like.energy(displacement) > 0.0
    assert opposite.energy(displacement) < 0.0
    assert like.force(displacement)[0] > 0.0
    assert opposite.force(displacement)[0] < 0.0
    assert like.force(displacement)[1] == pytest.approx(0.0)
    assert like.force(displacement)[2] == pytest.approx(0.0)


def test_lennard_jones_reference_points():
    """Sigma, the well minimum, and the repulsive hard contact."""

    epsilon = 0.997
    sigma = 0.34
    minimum = 2.0 ** (1.0 / 6.0) * sigma
    potential = LennardJonesPotential(epsilon=epsilon, sigma=sigma)

    assert potential.energy([sigma, 0.0, 0.0]) == pytest.approx(0.0, abs=1e-12)
    assert potential.energy([minimum, 0.0, 0.0]) == pytest.approx(-epsilon)
    assert potential.force([minimum, 0.0, 0.0]) == pytest.approx(
        np.zeros(3),
        abs=1e-12,
    )
    assert potential.force([sigma, 0.0, 0.0])[0] == pytest.approx(
        24.0 * epsilon / sigma,
        rel=1e-12,
    )


def test_thermal_energy_reference():
    """RT at 300 K is a hand-checkable reference magnitude."""

    expected = 0.008314462618 * 300.0

    assert units.MOLAR_GAS_CONSTANT * 300.0 == pytest.approx(expected)
    assert expected == pytest.approx(2.4943387854, rel=1e-12)


def test_coulomb_constant_matches_independent_si_formula():
    """Recompute k_e from SI values hardcoded in the test, not imported."""

    charge = 1.602176634e-19
    permittivity = 8.8541878128e-12
    avogadro = 6.02214076e23
    expected = (
        (1.0 / (4.0 * math.pi * permittivity))
        * charge**2
        * (1.0 / 1.0e-9)
        * avogadro
        / 1000.0
    )

    assert units.COULOMB_CONSTANT == pytest.approx(expected, rel=1e-12)
