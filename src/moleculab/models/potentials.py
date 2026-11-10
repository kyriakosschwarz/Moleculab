"""Potential-energy models.

A potential is a *physical model*: it says how energy and force depend on
position.  The functions in :mod:`moleculab.models.physics` then use that model
without knowing which particular potential is active.

Every potential must provide two mathematically consistent methods:

``energy(position)``
    Return the potential energy in kJ/mol.

``force(position)``
    Return the force vector in kJ/(mol nm), defined by

    .. math:: \\mathbf F = -\\nabla U.

Keeping the energy and force derived from the same model is one of the most
important correctness requirements in molecular simulation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol, runtime_checkable

import numpy as np

from ..foundation.units import COULOMB_CONSTANT
from ..foundation.validation import (
    as_displacement,
    as_finite_float,
    as_positive_float,
    as_vector3,
)


@runtime_checkable
class Potential(Protocol):
    """The interface expected by the generic physics functions."""

    def energy(self, position: np.ndarray) -> float:
        """Return potential energy in kJ/mol at *position* (nm)."""

    def force(self, position: np.ndarray) -> np.ndarray:
        """Return force in kJ/(mol nm) at *position* (nm)."""


@runtime_checkable
class PairPotential(Protocol):
    """The interface for an interaction between two particles.

    ``displacement`` is the vector :math:`\\mathbf r_i - \\mathbf r_j` in nm.
    ``energy`` returns the pair energy in kJ/mol.  ``force`` returns the force
    *on particle i* in kJ/(mol nm); the force on particle j is its negative.
    """

    def energy(self, displacement: np.ndarray) -> float:
        """Return pair energy in kJ/mol for a displacement in nm."""

    def force(self, displacement: np.ndarray) -> np.ndarray:
        """Return the force on particle i in kJ/(mol nm)."""



@dataclass
class HarmonicPotential:
    """A three-dimensional harmonic potential.

    The potential energy is

    .. math:: U(\\mathbf r) = \\frac{1}{2} k |\\mathbf r - \\mathbf r_0|^2,

    where ``stiffness`` is :math:`k` in kJ/(mol nm^2) and ``center`` is
    :math:`\\mathbf r_0` in nm.  The corresponding force is

    .. math:: \\mathbf F(\\mathbf r) = -k(\\mathbf r - \\mathbf r_0).

    This is an *external* potential for one particle, not a chemical bond
    between two atoms.  Bonded force-field terms will be introduced as a
    separate interaction model later in the roadmap.
    """

    stiffness: float
    center: np.ndarray = field(default_factory=lambda: np.zeros(3))

    def __post_init__(self) -> None:
        self.stiffness = as_positive_float(self.stiffness, "stiffness")
        self.center = as_vector3(self.center, "center").copy()

    def energy(self, position: np.ndarray) -> float:
        """Potential energy in kJ/mol at *position* (nm)."""

        displacement = self._displacement(position)
        return float(0.5 * self.stiffness * np.dot(displacement, displacement))

    def force(self, position: np.ndarray) -> np.ndarray:
        """Force in kJ/(mol nm) at *position* (nm)."""

        return -self.stiffness * self._displacement(position)

    def _displacement(self, position: np.ndarray) -> np.ndarray:
        return as_vector3(position, "position") - self.center


@dataclass
class DoubleWellPotential:
    """A one-dimensional double-well potential along x.

    The energy is

    .. math::

        U(x) = B \\left[\\left(\\frac{x-x_0}{w}\\right)^4
        - 2\\left(\\frac{x-x_0}{w}\\right)^2\\right],

    with minima at :math:`x_0 \\pm w` and an unstable barrier at
    :math:`x_0`.  The well depth is :math:`B` and the barrier height from a
    minimum to the top is also :math:`B`.

    This is an *external* reaction-coordinate model: the potential is flat in
    y and z.  It is intended for teaching transitions between two states, not
    for representing a full molecular force field.
    """

    barrier: float
    width: float
    center: np.ndarray = field(default_factory=lambda: np.zeros(3))

    def __post_init__(self) -> None:
        self.barrier = as_positive_float(self.barrier, "barrier")
        self.width = as_positive_float(self.width, "width")
        self.center = as_vector3(self.center, "center").copy()

    def energy(self, position: np.ndarray) -> float:
        """Potential energy in kJ/mol at *position* (nm)."""

        coordinate = (as_vector3(position, "position")[0] - self.center[0]) / self.width
        return float(self.barrier * (coordinate**4 - 2.0 * coordinate**2))

    def force(self, position: np.ndarray) -> np.ndarray:
        """Force in kJ/(mol nm) at *position* (nm)."""

        coordinate = (as_vector3(position, "position")[0] - self.center[0]) / self.width
        gradient = self.barrier * (4.0 * coordinate**3 - 4.0 * coordinate) / self.width
        return np.array([-gradient, 0.0, 0.0])


@dataclass
class LennardJonesPotential:
    """The 12-6 Lennard-Jones pair potential.

    .. math::

        U(r) = 4\\epsilon
        \\left[\\left(\\frac{\\sigma}{r}\\right)^{12}
        - \\left(\\frac{\\sigma}{r}\\right)^{6}\\right]

    ``epsilon`` is the depth of the attractive well in kJ/mol and ``sigma`` is
    the distance in nm at which the energy is zero.  The minimum is at
    :math:`r = 2^{1/6}\\sigma`, where ``U = -epsilon``.

    ``displacement`` is :math:`\\mathbf r_i - \\mathbf r_j`; the force returned
    by :meth:`force` acts on particle i.
    """

    epsilon: float
    sigma: float

    def __post_init__(self) -> None:
        self.epsilon = as_positive_float(self.epsilon, "epsilon")
        self.sigma = as_positive_float(self.sigma, "sigma")

    def energy(self, displacement: np.ndarray) -> float:
        """Return the pair energy in kJ/mol."""

        _, distance = as_displacement(displacement)
        sigma_over_r_sixth = (self.sigma / distance) ** 6
        return float(
            4.0
            * self.epsilon
            * (sigma_over_r_sixth**2 - sigma_over_r_sixth)
        )

    def force(self, displacement: np.ndarray) -> np.ndarray:
        """Return the force on particle i in kJ/(mol nm)."""

        vector, distance = as_displacement(displacement)
        sigma_over_r_sixth = (self.sigma / distance) ** 6
        return (
            24.0
            * self.epsilon
            * (2.0 * sigma_over_r_sixth**2 - sigma_over_r_sixth)
            / distance**2
            * vector
        )


@dataclass
class CoulombPotential:
    """Electrostatic interaction between two point charges.

    .. math:: U(r) = k_e \\frac{q_i q_j}{r}

    Charges are in units of the elementary charge ``e``, distance is in nm,
    and energy is in kJ/mol.  The constant :data:`moleculab.foundation.units.COULOMB_CONSTANT`
    performs the SI-to-internal-unit conversion, so the equation above can be
    used directly.

    Like charges repel (positive force along
    :math:`\\mathbf r_i - \\mathbf r_j`); opposite charges attract.
    """

    charge_i: float
    charge_j: float

    def __post_init__(self) -> None:
        self.charge_i = as_finite_float(self.charge_i, "charge_i")
        self.charge_j = as_finite_float(self.charge_j, "charge_j")

    def energy(self, displacement: np.ndarray) -> float:
        """Return the pair energy in kJ/mol."""

        _, distance = as_displacement(displacement)
        return float(
            COULOMB_CONSTANT * self.charge_i * self.charge_j / distance
        )

    def force(self, displacement: np.ndarray) -> np.ndarray:
        """Return the force on particle i in kJ/(mol nm)."""

        vector, distance = as_displacement(displacement)
        return (
            COULOMB_CONSTANT
            * self.charge_i
            * self.charge_j
            / distance**3
            * vector
        )
