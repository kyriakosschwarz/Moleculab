# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for the internal unit conventions and constants."""

import math

import pytest

from moleculab import units


def test_unit_labels_match_the_documented_system():
    assert units.DISTANCE_UNIT == "nm"
    assert units.TIME_UNIT == "ps"
    assert units.ENERGY_UNIT == "kJ/mol"
    assert units.MASS_UNIT == "Da"
    assert units.FORCE_UNIT == "kJ/(mol nm)"
    assert units.TEMPERATURE_UNIT == "K"
    assert units.CHARGE_UNIT == "e"


def test_molar_gas_constant_value():
    # R = 8.314462618 J/(mol K) = 0.008314462618 kJ/(mol K)
    assert units.MOLAR_GAS_CONSTANT == pytest.approx(0.008314462618)
    assert units.MOLAR_GAS_CONSTANT == pytest.approx(
        units.AVOGADRO_CONSTANT * units.BOLTZMANN_CONSTANT_J_PER_K / 1000.0
    )
    assert units.BOLTZMANN_CONSTANT_J_PER_K == pytest.approx(1.380649e-23)


def test_conversion_factors_are_consistent():
    assert units.NANOMETRE_IN_METRES == pytest.approx(1e-9)
    assert units.PICOSECOND_IN_SECONDS == pytest.approx(1e-12)
    assert units.DALTON_IN_KILOGRAMS == pytest.approx(1.66053906660e-27)
    assert units.KILOJOULES_PER_MOLE_TO_JOULES_PER_PARTICLE == pytest.approx(
        1000.0 / 6.02214076e23
    )


def test_newton_second_law_conversion_factor_is_one():
    force_to_newtons = (
        units.KILOJOULES_PER_MOLE_TO_JOULES_PER_PARTICLE
        / units.NANOMETRE_IN_METRES
    )
    mass_to_kilograms = units.DALTON_IN_KILOGRAMS
    acceleration_to_nm_per_ps2 = 1.0 / (
        units.NANOMETRE_IN_METRES
        / units.PICOSECOND_IN_SECONDS**2
    )

    factor = (
        force_to_newtons
        / mass_to_kilograms
        * acceleration_to_nm_per_ps2
    )
    # One dalton is approximately, not exactly, one gram per mole.  The
    # residual deviation is about 3.5e-10, far below simulation precision.
    assert factor == pytest.approx(1.0, rel=1e-9)


def test_one_dalton_nanometre_squared_per_picosecond_squared_is_one_kj_per_mol():
    energy_per_particle = units.DALTON_IN_KILOGRAMS * (
        units.NANOMETRE_IN_METRES / units.PICOSECOND_IN_SECONDS
    ) ** 2
    molar_energy_kj = (
        energy_per_particle * units.AVOGADRO_CONSTANT / 1000.0
    )

    assert molar_energy_kj == pytest.approx(1.0, rel=1e-6)


def test_molar_gas_constant_is_avogadro_times_boltzmann():
    expected = (
        units.AVOGADRO_CONSTANT
        * units.BOLTZMANN_CONSTANT_J_PER_K
        / 1000.0
    )

    assert units.MOLAR_GAS_CONSTANT == pytest.approx(expected, rel=1e-12)


def test_coulomb_constant_matches_its_si_derivation():
    expected = (
        (1.0 / (4.0 * math.pi * units.VACUUM_PERMITTIVITY_F_PER_M))
        * units.ELEMENTARY_CHARGE_COULOMBS**2
        * (1.0 / units.NANOMETRE_IN_METRES)
        * units.AVOGADRO_CONSTANT
        / 1000.0
    )

    assert units.COULOMB_CONSTANT == pytest.approx(expected, rel=1e-12)
