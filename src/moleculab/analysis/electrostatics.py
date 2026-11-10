"""Dipole moments from formal or partial charges.

A dipole moment is a charge-weighted displacement. For a neutral system it is
independent of the coordinate origin. For a charged system it is origin-
dependent, so an explicit origin is required.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..foundation.units import ELEMENTARY_CHARGE_COULOMBS, NANOMETRE_IN_METRES

_DEBYE_IN_COULOMB_METRES = 3.33564e-30
DEBYE_PER_E_NM = (
    ELEMENTARY_CHARGE_COULOMBS * NANOMETRE_IN_METRES / _DEBYE_IN_COULOMB_METRES
)


@dataclass(frozen=True)
class DipoleMoment:
    """Dipole vector and magnitude in e nm and Debye."""

    vector_e_nm: np.ndarray
    magnitude_e_nm: float
    origin_nm: np.ndarray

    @property
    def vector_debye(self) -> np.ndarray:
        """Dipole vector in Debye."""

        return self.vector_e_nm * DEBYE_PER_E_NM

    @property
    def magnitude_debye(self) -> float:
        """Dipole magnitude in Debye."""

        return self.magnitude_e_nm * DEBYE_PER_E_NM


def dipole_moment(
    charges: np.ndarray,
    positions: np.ndarray,
    *,
    origin: np.ndarray | None = None,
    neutral_tolerance: float = 1e-6,
) -> DipoleMoment:
    """Return the dipole moment for charges at positions in nm.

    For a neutral charge set the result is origin-independent. If the net
    charge is non-zero, an explicit ``origin`` is required because the dipole
    moment then depends on the coordinate origin.

    The result is in elementary-charge nanometres and can be converted to
    Debye through :class:`DipoleMoment`.
    """

    charge_values = np.asarray(charges, dtype=float).reshape(-1)
    coordinates = np.asarray(positions, dtype=float)
    if charge_values.size == 0:
        raise ValueError("charges must contain at least one atom.")
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("positions must have shape (n_atoms, 3).")
    if coordinates.shape[0] != charge_values.size:
        raise ValueError("charges and positions must have the same length.")
    if not np.all(np.isfinite(charge_values)):
        raise ValueError("charges must contain only finite values.")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("positions must contain only finite values.")

    tolerance = float(neutral_tolerance)
    if not np.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError("neutral_tolerance must be finite and non-negative.")

    net_charge = float(np.sum(charge_values))
    if origin is None:
        if abs(net_charge) > tolerance:
            raise ValueError(
                "dipole moment is origin-dependent for a charged system; "
                "provide an explicit origin."
            )
        origin_value = np.mean(coordinates, axis=0)
    else:
        origin_value = np.asarray(origin, dtype=float)
        if origin_value.shape != (3,) or not np.all(np.isfinite(origin_value)):
            raise ValueError("origin must be a finite three-vector in nm.")

    vector = np.sum(
        charge_values[:, np.newaxis] * (coordinates - origin_value),
        axis=0,
    )
    return DipoleMoment(
        vector_e_nm=vector,
        magnitude_e_nm=float(np.linalg.norm(vector)),
        origin_nm=origin_value.copy(),
    )
