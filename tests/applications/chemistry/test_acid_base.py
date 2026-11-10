# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for acid-base equilibria and discrete protonation populations."""

import numpy as np
import pytest

from moleculab import (
    boltzmann_populations,
    deprotonated_fraction,
    mean_net_charge,
    mean_protonation,
    protonated_fraction,
    protonation_microstates,
    tautomer_populations,
    units,
)


def test_at_the_pka_half_of_the_acid_is_deprotonated():
    assert deprotonated_fraction(4.76, 4.76) == pytest.approx(0.5, abs=1e-15)
    assert protonated_fraction(4.76, 4.76) == pytest.approx(0.5, abs=1e-15)


def test_one_ph_unit_each_side_gives_a_ten_to_one_ratio():
    assert deprotonated_fraction(4.76, 5.76) == pytest.approx(10.0 / 11.0)
    assert deprotonated_fraction(4.76, 3.76) == pytest.approx(1.0 / 11.0)


def test_fractions_sum_to_one_over_a_range_of_ph():
    for ph in np.linspace(-2.0, 18.0, 41):
        total = deprotonated_fraction(4.76, ph) + protonated_fraction(4.76, ph)
        assert total == pytest.approx(1.0, abs=1e-15)


def test_acid_fraction_matches_the_independent_henderson_hasselbalch_formula():
    pka = 4.76
    ph = 7.4
    expected = 1.0 / (1.0 + 10.0 ** (pka - ph))
    assert deprotonated_fraction(pka, ph) == pytest.approx(expected, rel=1e-13)
    assert deprotonated_fraction(pka, ph) == pytest.approx(0.997714, abs=1e-6)


def test_base_protonation_uses_the_conjugate_acid_formula():
    pka = 9.25
    ph = 7.4
    expected = 1.0 / (1.0 + 10.0 ** (ph - pka))
    assert protonated_fraction(pka, ph) == pytest.approx(expected, rel=1e-13)
    assert protonated_fraction(pka, ph) == pytest.approx(0.986071, abs=1e-6)


def test_microstate_probabilities_are_products_for_independent_sites():
    pkas = [4.0, 9.0]
    ph = 7.0
    probabilities = protonation_microstates(pkas, ph)

    first = protonated_fraction(4.0, ph)
    second = protonated_fraction(9.0, ph)
    expected = np.array(
        [
            (1.0 - first) * (1.0 - second),
            first * (1.0 - second),
            (1.0 - first) * second,
            first * second,
        ]
    )

    np.testing.assert_allclose(probabilities, expected, atol=1e-15)
    assert np.sum(probabilities) == pytest.approx(1.0, abs=1e-15)


def test_microstate_marginals_match_mean_site_protonation():
    pkas = [3.0, 6.0, 10.0]
    ph = 6.8
    probabilities = protonation_microstates(pkas, ph)

    marginals = np.array(
        [
            sum(
                probability
                for state, probability in enumerate(probabilities)
                if (state >> site) & 1
            )
            for site in range(len(pkas))
        ]
    )

    np.testing.assert_allclose(
        marginals,
        mean_protonation(pkas, ph),
        atol=1e-14,
    )


def test_mean_net_charge_uses_each_sites_protonation_charges():
    pkas = [4.0, 9.0]
    ph = 7.0
    charge = mean_net_charge(
        pkas,
        ph,
        protonated_charges=[0.0, 1.0],
        deprotonated_charges=[-1.0, 0.0],
    )

    acid_protonated = protonated_fraction(4.0, ph)
    base_protonated = protonated_fraction(9.0, ph)
    expected = (1.0 - acid_protonated) * -1.0 + base_protonated * 1.0

    assert charge == pytest.approx(expected, abs=1e-14)
    assert charge == pytest.approx(-0.008902, abs=1e-6)


def test_boltzmann_populations_are_equal_for_degenerate_states():
    populations = boltzmann_populations([0.0, 0.0, 0.0], temperature=300.0)
    np.testing.assert_allclose(populations, np.full(3, 1.0 / 3.0), atol=1e-15)


def test_boltzmann_population_ratio_matches_the_free_energy_difference():
    temperature = 300.0
    delta_g = units.MOLAR_GAS_CONSTANT * temperature * np.log(9.0)

    populations = boltzmann_populations([0.0, delta_g], temperature)

    assert populations[1] / populations[0] == pytest.approx(1.0 / 9.0, rel=1e-13)
    assert populations[1] == pytest.approx(0.1, abs=1e-13)


def test_boltzmann_populations_are_invariant_to_a_constant_offset():
    first = boltzmann_populations([0.0, 2.0, 5.0], temperature=280.0)
    shifted = boltzmann_populations([10.0, 12.0, 15.0], temperature=280.0)
    np.testing.assert_allclose(first, shifted, atol=1e-15)


def test_tautomer_populations_use_delta_g_as_g_b_minus_g_a():
    equal = tautomer_populations(0.0, temperature=300.0)
    np.testing.assert_allclose(equal, [0.5, 0.5], atol=1e-15)

    temperature = 300.0
    delta_g = units.MOLAR_GAS_CONSTANT * temperature * np.log(9.0)
    populations = tautomer_populations(delta_g, temperature)
    assert populations[0] == pytest.approx(0.9, abs=1e-13)
    assert populations[1] == pytest.approx(0.1, abs=1e-13)


def test_acid_base_validation():
    with pytest.raises(ValueError, match="pka"):
        deprotonated_fraction(float("nan"), 7.0)
    with pytest.raises(ValueError, match="ph"):
        deprotonated_fraction(4.76, float("inf"))

    with pytest.raises(ValueError, match="at least one"):
        protonation_microstates([], 7.0)
    with pytest.raises(ValueError, match="16"):
        protonation_microstates([7.0] * 17, 7.0)
    with pytest.raises(ValueError, match="finite"):
        protonation_microstates([4.0, float("nan")], 7.0)

    with pytest.raises(ValueError, match="length"):
        mean_net_charge([4.0], 7.0, [0.0], [])

    with pytest.raises(ValueError, match="at least one"):
        boltzmann_populations([])
    with pytest.raises(ValueError, match="temperature"):
        boltzmann_populations([0.0, 1.0], temperature=0.0)
    with pytest.raises(ValueError, match="delta_g"):
        tautomer_populations(float("nan"))
