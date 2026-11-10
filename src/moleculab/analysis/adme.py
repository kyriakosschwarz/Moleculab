# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Ionisation-aware lipophilicity helpers.

The functions here convert a neutral-form LogP into a pH-dependent LogD for a
monoprotic acid or base. They do not predict pKa, permeability, or oral
bioavailability.
"""

from __future__ import annotations

from typing import Literal

import numpy as np

IonisationType = Literal["acid", "base"]


def logd_from_logp(
    logp: float,
    pka: float,
    ph: float,
    *,
    ionisation: IonisationType,
) -> float:
    """Return an approximate LogD for a monoprotic acid or base.

    Only the neutral form is assumed to partition into octanol. For an acid,

    .. math::

        \\log D = \\log P - \\log_{10}(1 + 10^{\\mathrm{pH}-\\mathrm{p}K_a}).

    For a base,

    .. math::

        \\log D = \\log P - \\log_{10}(1 + 10^{\\mathrm{p}K_a-\\mathrm{pH}}).

    The pKa is an input. This is not a pKa predictor.
    """

    logp_value = float(logp)
    pka_value = float(pka)
    ph_value = float(ph)
    if not np.isfinite(logp_value):
        raise ValueError("logp must be finite.")
    if not np.isfinite(pka_value):
        raise ValueError("pka must be finite.")
    if not np.isfinite(ph_value):
        raise ValueError("ph must be finite.")

    if ionisation == "acid":
        exponent = ph_value - pka_value
    elif ionisation == "base":
        exponent = pka_value - ph_value
    else:
        raise ValueError("ionisation must be 'acid' or 'base'.")

    correction = np.logaddexp(0.0, exponent * np.log(10.0)) / np.log(10.0)
    return float(logp_value - correction)
