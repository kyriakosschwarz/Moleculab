# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Minimal QM/MM: a quantum particle in a 1D box coupled to MM point charges.

The time-independent Schrödinger equation is solved by finite differences on
a uniform grid:

.. math::

    -\\frac{\\hbar^2}{2m}\\frac{d^2\\psi}{dx^2}
    + V(x)\\psi = E\\psi,

with Dirichlet boundary conditions at the box edges.  The Hamiltonian is
tridiagonal and is diagonalized with :func:`numpy.linalg.eigh`.

The QM particle carries a charge :math:`q_{\\mathrm{QM}}` and interacts with
fixed MM point charges through

.. math::

    V_{\\mathrm{MM}}(x) = k_e q_{\\mathrm{QM}}
    \\sum_i \\frac{q_i}{|x-x_i|}.

Lengths are in nm, masses in Da, charges in elementary charge, and energies
in kJ/mol.  This is a deliberately one-dimensional teaching model; it is not
a production QM/MM engine.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from ..foundation.units import (
    AVOGADRO_CONSTANT,
    COULOMB_CONSTANT,
    DALTON_IN_KILOGRAMS,
    H_BAR_J_S,
)
from ..foundation.validation import as_finite_float, as_positive_float


@dataclass(frozen=True)
class MMPointCharge:
    """A fixed MM point charge on the 1D axis."""

    position: float
    charge: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "position", as_finite_float(self.position, "position")
        )
        object.__setattr__(
            self, "charge", as_finite_float(self.charge, "charge")
        )


@dataclass(frozen=True)
class QMResult:
    """Eigenvalues and wavefunctions from a finite-difference solve."""

    grid: np.ndarray
    potential: np.ndarray
    energies: np.ndarray
    wavefunctions: np.ndarray

    @property
    def probability_density(self) -> np.ndarray:
        """Probability density for each state, shape ``(n_states, n_grid)``."""

        return np.abs(self.wavefunctions.T) ** 2


@dataclass(frozen=True)
class ParticleInBox:
    """A quantum particle in a 1D box with optional QM/MM coupling.

    Parameters
    ----------
    length:
        Box length in nm.
    mass:
        Particle mass in Da.  Use ``moleculab.foundation.units.ELECTRON_MASS_DA`` for an
        electron.
    n_grid:
        Number of interior finite-difference grid points.
    """

    length: float
    mass: float
    n_grid: int = 400

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "length", as_positive_float(self.length, "length")
        )
        object.__setattr__(self, "mass", as_positive_float(self.mass, "mass"))
        if isinstance(self.n_grid, bool) or not isinstance(
            self.n_grid, (int, np.integer)
        ):
            raise ValueError("n_grid must be an integer.")
        if int(self.n_grid) < 10:
            raise ValueError("n_grid must be at least 10.")
        object.__setattr__(self, "n_grid", int(self.n_grid))

    @property
    def grid(self) -> np.ndarray:
        """Interior grid positions in nm."""

        spacing = self.length / (self.n_grid + 1)
        return np.arange(1, self.n_grid + 1, dtype=float) * spacing

    @property
    def kinetic_coefficient(self) -> float:
        """Coefficient of :math:`-d^2/dx^2` in kJ/mol nm^2."""

        mass_kg = self.mass * DALTON_IN_KILOGRAMS
        return (
            H_BAR_J_S**2
            / (2.0 * mass_kg)
            * 1.0e18
            * AVOGADRO_CONSTANT
            / 1000.0
        )

    def analytic_energies(self, n_states: int = 3) -> np.ndarray:
        """Return exact particle-in-a-box energies in kJ/mol."""

        states = self._validate_states(n_states)
        quantum_numbers = np.arange(1, states + 1, dtype=float)
        return self.kinetic_coefficient * (
            quantum_numbers * np.pi / self.length
        ) ** 2

    def mm_potential(
        self,
        mm_charges: Sequence[MMPointCharge],
        qm_charge: float = 0.0,
    ) -> np.ndarray:
        """Return the MM electrostatic potential on the QM grid in kJ/mol."""

        charge_qm = as_finite_float(qm_charge, "qm_charge")
        potential = np.zeros(self.n_grid, dtype=float)
        for index, mm_charge in enumerate(mm_charges):
            if not isinstance(mm_charge, MMPointCharge):
                raise TypeError(f"mm_charges[{index}] must be an MMPointCharge.")
            distance = np.abs(self.grid - mm_charge.position)
            if np.any(distance < 1e-8):
                raise ValueError(
                    "An MM point charge coincides with a QM grid point."
                )
            potential += (
                COULOMB_CONSTANT
                * charge_qm
                * mm_charge.charge
                / distance
            )
        return potential

    def harmonic_potential(
        self,
        stiffness: float,
        center: float | None = None,
    ) -> np.ndarray:
        """Return a harmonic potential on the QM grid in kJ/mol."""

        k = as_positive_float(stiffness, "stiffness")
        center_value = (
            self.length / 2.0
            if center is None
            else as_finite_float(center, "center")
        )
        return 0.5 * k * (self.grid - center_value) ** 2

    def solve(
        self,
        potential: np.ndarray | None = None,
        mm_charges: Sequence[MMPointCharge] = (),
        qm_charge: float = 0.0,
        n_states: int = 3,
    ) -> QMResult:
        """Solve the finite-difference Schrödinger equation."""

        states = self._validate_states(n_states)
        if potential is None:
            potential_values = np.zeros(self.n_grid, dtype=float)
        else:
            potential_values = np.asarray(potential, dtype=float)
            if potential_values.shape != (self.n_grid,):
                raise ValueError(
                    f"potential must have shape ({self.n_grid},), "
                    f"got {potential_values.shape}."
                )
            if not np.all(np.isfinite(potential_values)):
                raise ValueError("potential must contain only finite values.")
        potential_values = potential_values + self.mm_potential(
            mm_charges, qm_charge=qm_charge
        )

        spacing = self.length / (self.n_grid + 1)
        kinetic = self.kinetic_coefficient
        diagonal = 2.0 * kinetic / spacing**2 + potential_values
        off_diagonal = -kinetic / spacing**2
        hamiltonian = np.diag(diagonal)
        hamiltonian += np.diag(
            np.full(self.n_grid - 1, off_diagonal),
            k=1,
        )
        hamiltonian += np.diag(
            np.full(self.n_grid - 1, off_diagonal),
            k=-1,
        )
        energies, wavefunctions = np.linalg.eigh(hamiltonian)

        return QMResult(
            grid=self.grid,
            potential=potential_values,
            energies=energies[:states],
            wavefunctions=wavefunctions[:, :states],
        )

    def _validate_states(self, n_states: int) -> int:
        if isinstance(n_states, bool) or not isinstance(
            n_states, (int, np.integer)
        ):
            raise ValueError("n_states must be an integer.")
        states = int(n_states)
        if states < 1 or states > self.n_grid:
            raise ValueError(
                f"n_states must be between 1 and {self.n_grid}, got {states}."
            )
        return states
