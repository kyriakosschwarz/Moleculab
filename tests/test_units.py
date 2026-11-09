"""Tests for the internal unit conventions and constants."""

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
    assert units.BOLTZMANN_CONSTANT == pytest.approx(0.008314462618)
    assert units.BOLTZMANN_CONSTANT_J_PER_K == pytest.approx(1.380649e-23)


def test_conversion_factors_are_consistent():
    assert units.NANOMETRE_IN_METRES == pytest.approx(1e-9)
    assert units.PICOSECOND_IN_SECONDS == pytest.approx(1e-12)
    assert units.DALTON_IN_KILOGRAMS == pytest.approx(1.66053906660e-27)
    assert units.KILOJOULES_PER_MOLE_TO_JOULES_PER_PARTICLE == pytest.approx(
        1000.0 / 6.02214076e23
    )
