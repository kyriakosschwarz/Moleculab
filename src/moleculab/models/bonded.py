"""Potentials for bonds and angles.

Nonbonded potentials act between every pair of particles.  Bonded potentials
act only between atoms that are connected by the molecular graph, so they are
stored together with explicit atom indices in :mod:`moleculab.models.topology`.

The displacement conventions are:

``BondPotential``
    ``displacement = r_i - r_j`` and ``force`` returns the force on atom i.

``AnglePotential``
    ``r_ji = r_i - r_j`` and ``r_jk = r_k - r_j`` for the angle i-j-k.
    ``forces`` returns ``(force_i, force_j, force_k)``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

import numpy as np

from ..foundation.validation import (
    as_displacement,
    as_finite_float,
    as_positive_float,
    as_vector3,
)


@runtime_checkable
class BondPotential(Protocol):
    """The interface for a two-atom bonded interaction."""

    def energy(self, displacement: np.ndarray) -> float:
        """Return bond energy in kJ/mol."""

    def force(self, displacement: np.ndarray) -> np.ndarray:
        """Return the force on atom i in kJ/(mol nm)."""


@runtime_checkable
class AnglePotential(Protocol):
    """The interface for a three-atom angle interaction."""

    def energy(self, r_ji: np.ndarray, r_jk: np.ndarray) -> float:
        """Return angle energy in kJ/mol."""

    def forces(
        self,
        r_ji: np.ndarray,
        r_jk: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(force_i, force_j, force_k)`` in kJ/(mol nm)."""


@runtime_checkable
class TorsionPotential(Protocol):
    """The interface for a four-atom torsional interaction."""

    def energy(
        self,
        r_ij: np.ndarray,
        r_jk: np.ndarray,
        r_kl: np.ndarray,
    ) -> float:
        """Return torsion energy in kJ/mol.

        The vectors are ``r_ij = r_j - r_i``, ``r_jk = r_k - r_j``, and
        ``r_kl = r_l - r_k``.
        """

    def forces(
        self,
        r_ij: np.ndarray,
        r_jk: np.ndarray,
        r_kl: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(force_i, force_j, force_k, force_l)`` in kJ/(mol nm)."""



@dataclass
class HarmonicBondPotential:
    """A harmonic bond-stretching potential.

    .. math:: U(r) = \\frac{1}{2} k (r - r_0)^2

    ``stiffness`` is :math:`k` in kJ/(mol nm^2) and ``length0`` is the
    equilibrium bond length :math:`r_0` in nm.
    """

    stiffness: float
    length0: float

    def __post_init__(self) -> None:
        self.stiffness = as_positive_float(self.stiffness, "stiffness")
        self.length0 = as_positive_float(self.length0, "length0")

    def energy(self, displacement: np.ndarray) -> float:
        """Return bond energy in kJ/mol."""

        _, distance = as_displacement(displacement)
        stretch = distance - self.length0
        return float(0.5 * self.stiffness * stretch**2)

    def force(self, displacement: np.ndarray) -> np.ndarray:
        """Return the force on atom i in kJ/(mol nm).

        For a stretched bond this points from i toward j; for a compressed
        bond it points away from j.
        """

        vector, distance = as_displacement(displacement)
        return -(self.stiffness * (distance - self.length0) / distance) * vector


@dataclass
class HarmonicAnglePotential:
    """A harmonic angle-bending potential.

    .. math:: U(\\theta) = \\frac{1}{2} k (\\theta - \\theta_0)^2

    ``stiffness`` is :math:`k` in kJ/(mol rad^2) and ``angle0`` is the
    equilibrium angle :math:`\\theta_0` in **radians**.  Use
    ``np.deg2rad`` when converting a value in degrees.
    """

    stiffness: float
    angle0: float

    def __post_init__(self) -> None:
        self.stiffness = as_positive_float(self.stiffness, "stiffness")
        self.angle0 = as_finite_float(self.angle0, "angle0")
        if not 0.0 <= self.angle0 <= np.pi:
            raise ValueError("angle0 must be between 0 and pi radians.")

    def energy(self, r_ji: np.ndarray, r_jk: np.ndarray) -> float:
        """Return angle energy in kJ/mol."""

        _, _, _, angle, _ = self._geometry(r_ji, r_jk)
        return float(0.5 * self.stiffness * (angle - self.angle0) ** 2)

    def forces(
        self,
        r_ji: np.ndarray,
        r_jk: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return forces on atoms i, j, and k in kJ/(mol nm).

        The force on the central atom j follows from translational
        invariance: the three forces must sum to zero.
        """

        a, b, cos_theta, angle, sin_theta = self._geometry(r_ji, r_jk)
        if sin_theta < 1e-12:
            raise ValueError(
                "Angle force is undefined for a linear (0 or pi) geometry."
            )

        length_a = float(np.linalg.norm(a))
        length_b = float(np.linalg.norm(b))
        unit_product = 1.0 / (length_a * length_b)

        # dc/da and dc/db for c = cos(theta) = a . b / (|a||b|).
        dc_da = b * unit_product - cos_theta * a / length_a**2
        dc_db = a * unit_product - cos_theta * b / length_b**2

        # F = -dU/dr; dtheta = -dc/sin(theta).
        prefactor = self.stiffness * (angle - self.angle0) / sin_theta
        force_i = prefactor * dc_da
        force_k = prefactor * dc_db
        force_j = -(force_i + force_k)
        return force_i, force_j, force_k

    def _geometry(
        self,
        r_ji: np.ndarray,
        r_jk: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, float, float, float]:
        a, length_a = as_displacement(r_ji, "r_ji")
        b, length_b = as_displacement(r_jk, "r_jk")
        cos_theta = float(np.dot(a, b) / (length_a * length_b))
        cos_theta = float(np.clip(cos_theta, -1.0, 1.0))
        angle = float(np.arccos(cos_theta))
        sin_theta = float(np.sqrt(max(0.0, 1.0 - cos_theta**2)))
        return a, b, cos_theta, angle, sin_theta


def _torsion_vectors(
    r_i: np.ndarray,
    r_j: np.ndarray,
    r_k: np.ndarray,
    r_l: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return bond vectors and the signed dihedral angle."""

    b1, _ = as_displacement(as_vector3(r_j, "r_j") - as_vector3(r_i, "r_i"), "b1")
    b2, _ = as_displacement(as_vector3(r_k, "r_k") - as_vector3(r_j, "r_j"), "b2")
    b3, _ = as_displacement(as_vector3(r_l, "r_l") - as_vector3(r_k, "r_k"), "b3")

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    if np.dot(n1, n1) < 1e-24 or np.dot(n2, n2) < 1e-24:
        raise ValueError("Dihedral angle is undefined for collinear atoms.")

    x = float(np.dot(n1, n2))
    y = float(np.dot(np.cross(n1, n2), b2 / np.linalg.norm(b2)))
    return b1, b2, b3, float(np.arctan2(y, x))


def dihedral_angle(
    r_i: np.ndarray,
    r_j: np.ndarray,
    r_k: np.ndarray,
    r_l: np.ndarray,
) -> float:
    """Return the signed i-j-k-l dihedral angle in radians.

    The angle lies in ``(-pi, pi]``.  Collinear atoms make the dihedral
    undefined and raise :class:`ValueError`.
    """

    return _torsion_vectors(r_i, r_j, r_k, r_l)[3]


@dataclass
class PeriodicTorsionPotential:
    """A periodic proper-torsion potential.

    .. math::

        U(\\phi) = k \\left[1 + \\cos(n\\phi - \\delta)\\right]

    ``barrier`` is :math:`k` in kJ/mol (the full barrier height is
    :math:`2k`), ``periodicity`` is the integer :math:`n`, and ``phase`` is
    :math:`\\delta` in radians.  The dihedral angle is defined by atoms
    i-j-k-l and is returned by :func:`dihedral_angle`.
    """

    barrier: float
    periodicity: int
    phase: float = 0.0

    def __post_init__(self) -> None:
        self.barrier = as_positive_float(self.barrier, "barrier")
        if isinstance(self.periodicity, bool) or not isinstance(
            self.periodicity, (int, np.integer)
        ):
            raise ValueError("periodicity must be an integer.")
        if int(self.periodicity) < 1:
            raise ValueError("periodicity must be at least 1.")
        self.periodicity = int(self.periodicity)
        self.phase = as_finite_float(self.phase, "phase")

    def energy(
        self,
        r_ij: np.ndarray,
        r_jk: np.ndarray,
        r_kl: np.ndarray,
    ) -> float:
        """Return torsion energy in kJ/mol."""

        r_ij = as_vector3(r_ij, "r_ij")
        r_jk = as_vector3(r_jk, "r_jk")
        r_kl = as_vector3(r_kl, "r_kl")
        angle = dihedral_angle(
            np.zeros(3),
            r_ij,
            r_ij + r_jk,
            r_ij + r_jk + r_kl,
        )
        return self.energy_from_angle(angle)

    def energy_from_angle(self, angle: float) -> float:
        """Return the torsion energy for a dihedral angle in radians.

        This is useful for conformational scans, where the geometry is varied
        directly rather than reconstructed from three bond vectors.
        """

        return float(
            self.barrier
            * (1.0 + np.cos(self.periodicity * float(angle) - self.phase))
        )

    def force_from_angle(self, angle: float) -> float:
        """Return the generalized torsional force in kJ/(mol rad).

        This is :math:`-dU/d\\phi`, not a Cartesian force in kJ/(mol nm).
        """

        return float(
            self.barrier
            * self.periodicity
            * np.sin(self.periodicity * float(angle) - self.phase)
        )

    def forces(
        self,
        r_ij: np.ndarray,
        r_jk: np.ndarray,
        r_kl: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Return forces on atoms i, j, k, and l in kJ/(mol nm)."""

        r_ij = as_vector3(r_ij, "r_ij")
        r_jk = as_vector3(r_jk, "r_jk")
        r_kl = as_vector3(r_kl, "r_kl")
        b1, b2, b3, angle = _torsion_vectors(
            np.zeros(3),
            r_ij,
            r_ij + r_jk,
            r_ij + r_jk + r_kl,
        )

        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        n1_squared = float(np.dot(n1, n1))
        n2_squared = float(np.dot(n2, n2))
        b2_squared = float(np.dot(b2, b2))
        b2_length = float(np.sqrt(b2_squared))

        d1 = -(b2_length / n1_squared) * n1
        d4 = (b2_length / n2_squared) * n2
        projection_1 = float(np.dot(b1, b2)) / b2_squared
        projection_2 = -float(np.dot(b2, b3)) / b2_squared
        d2 = -(1.0 + projection_1) * d1 - projection_2 * d4
        d3 = projection_1 * d1 + (projection_2 - 1.0) * d4

        influence = (
            self.barrier
            * self.periodicity
            * np.sin(self.periodicity * angle - self.phase)
        )
        return influence * d1, influence * d2, influence * d3, influence * d4
