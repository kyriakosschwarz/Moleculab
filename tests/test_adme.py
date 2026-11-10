"""Tests for LogP to LogD ionisation corrections."""

import numpy as np
import pytest

from moleculab import logd_from_logp


def test_at_pka_half_the_species_is_ionised():
    expected = 3.0 - np.log10(2.0)

    assert logd_from_logp(3.0, 5.0, 5.0, ionisation="acid") == pytest.approx(expected)
    assert logd_from_logp(3.0, 9.0, 9.0, ionisation="base") == pytest.approx(expected)


def test_acid_logd_drops_above_its_pka():
    at_low_ph = logd_from_logp(3.0, 4.9, 2.0, ionisation="acid")
    at_high_ph = logd_from_logp(3.0, 4.9, 7.4, ionisation="acid")

    assert at_low_ph == pytest.approx(3.0, abs=1e-3)
    assert at_high_ph == pytest.approx(3.0 - np.log10(1.0 + 10.0**2.5))


def test_base_logd_drops_below_its_pka():
    at_high_ph = logd_from_logp(1.0, 9.6, 12.0, ionisation="base")
    at_low_ph = logd_from_logp(1.0, 9.6, 7.4, ionisation="base")

    assert at_high_ph == pytest.approx(
        1.0 - np.log10(1.0 + 10.0 ** (9.6 - 12.0))
    )
    assert at_low_ph == pytest.approx(1.0 - np.log10(1.0 + 10.0**2.2))


def test_ibuprofen_and_atenolol_examples():
    ibuprofen = logd_from_logp(3.0732, 4.9, 7.4, ionisation="acid")
    atenolol = logd_from_logp(0.4521, 9.6, 7.4, ionisation="base")

    assert ibuprofen == pytest.approx(0.572, abs=1e-3)
    assert atenolol == pytest.approx(-1.7506, abs=1e-3)
    assert ibuprofen < 3.0732
    assert atenolol < 0.4521


def test_logd_limiting_behaviour():
    logp = 2.5
    pka = 6.0

    acid_neutral = logd_from_logp(logp, pka, pka - 10.0, ionisation="acid")
    acid_ionised = logd_from_logp(logp, pka, pka + 10.0, ionisation="acid")
    base_neutral = logd_from_logp(logp, pka, pka + 10.0, ionisation="base")
    base_ionised = logd_from_logp(logp, pka, pka - 10.0, ionisation="base")

    assert acid_neutral == pytest.approx(logp, abs=1e-9)
    assert acid_ionised == pytest.approx(logp - 10.0, abs=1e-9)
    assert base_neutral == pytest.approx(logp, abs=1e-9)
    assert base_ionised == pytest.approx(logp - 10.0, abs=1e-9)


def test_logd_is_numerically_stable_for_extreme_exponents():
    for ionisation in ("acid", "base"):
        for ph in (-1000.0, 1000.0):
            value = logd_from_logp(3.0, 5.0, ph, ionisation=ionisation)
            assert np.isfinite(value)


def test_logd_validation():
    with pytest.raises(ValueError, match="ionisation"):
        logd_from_logp(3.0, 5.0, 7.0, ionisation="neutral")
    with pytest.raises(ValueError, match="logp"):
        logd_from_logp(np.nan, 5.0, 7.0, ionisation="acid")
    with pytest.raises(ValueError, match="pka"):
        logd_from_logp(3.0, np.inf, 7.0, ionisation="acid")
    with pytest.raises(ValueError, match="ph"):
        logd_from_logp(3.0, 5.0, np.nan, ionisation="base")
