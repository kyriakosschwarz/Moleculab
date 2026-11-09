"""Ewald summation for periodic point charges.

For a periodic system, the Coulomb interaction cannot be truncated naively:
every charge interacts with every periodic image of every other charge.  Ewald
splits the conditionally convergent lattice sum into two rapidly convergent
parts:

* a short-range real-space sum using :math:`\\operatorname{erfc}(\\alpha r)`;
* a long-range reciprocal-space sum over wave vectors :math:`\\mathbf k`.

The total electrostatic energy is

.. math::

    U = U_{\\mathrm{real}} + U_{\\mathrm{recip}}
        + U_{\\mathrm{self}} + U_{\\mathrm{charge}}.

Moleculab uses the standard tin-foil boundary condition (the surface term is
zero).  The real-space part is evaluated in the normal pair loop of
:class:`~moleculab.nonbonded.NonbondedModel`; the reciprocal, self, and net
charge terms are added once by :class:`~moleculab.system.System`.

This is a minimal educational implementation.  Intramolecular 1-2/1-3
exclusions are supported through a reciprocal-space pair correction: the
real-space pair is skipped by the normal pair loop, and the excluded pair's
reciprocal contribution is subtracted explicitly.  Self-image interactions
remain part of the periodic Ewald energy.
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from ._validation import as_displacement, as_positive_float
from .backend import _load_torch, get_backend, to_numpy, to_torch, torch_dtype
from .box import Box
from .units import COULOMB_CONSTANT


@dataclass(frozen=True)
class EwaldElectrostatics:
    """Ewald electrostatic parameters.

    Parameters
    ----------
    alpha:
        Ewald splitting parameter in nm^-1.  Larger values move more of the
        interaction into reciprocal space.
    reciprocal_cutoff:
        Maximum wave-vector magnitude :math:`k_{\\max}` in nm^-1.  More
        reciprocal vectors improve accuracy at extra cost.
    """

    alpha: float
    reciprocal_cutoff: float

    def __post_init__(self) -> None:
        object.__setattr__(self, "alpha", as_positive_float(self.alpha, "alpha"))
        object.__setattr__(
            self,
            "reciprocal_cutoff",
            as_positive_float(self.reciprocal_cutoff, "reciprocal_cutoff"),
        )

    def energy(
        self,
        charge_i: float,
        charge_j: float,
        distance: float,
        cutoff: float,
    ) -> float:
        """Return the real-space Ewald pair energy in kJ/mol."""

        cutoff_value = as_positive_float(cutoff, "cutoff")
        if distance == 0.0:
            raise ValueError("Ewald real-space energy is undefined at r=0.")
        if distance >= cutoff_value:
            return 0.0
        return float(
            COULOMB_CONSTANT
            * float(charge_i)
            * float(charge_j)
            * math.erfc(self.alpha * distance)
            / distance
        )

    def force(
        self,
        charge_i: float,
        charge_j: float,
        displacement: np.ndarray,
        cutoff: float,
    ) -> np.ndarray:
        """Return the real-space Ewald force on particle i in kJ/(mol nm)."""

        vector, distance = as_displacement(displacement)
        cutoff_value = as_positive_float(cutoff, "cutoff")
        if distance >= cutoff_value:
            return np.zeros(3)

        prefactor = COULOMB_CONSTANT * float(charge_i) * float(charge_j)
        erfc_term = math.erfc(self.alpha * distance)
        exponential_term = np.exp(-(self.alpha**2) * distance**2)
        factor = prefactor * (
            erfc_term / distance**3
            + 2.0
            * self.alpha
            / np.sqrt(np.pi)
            * exponential_term
            / distance**2
        )
        return factor * vector

    def reciprocal_k_vectors(self, box: Box) -> tuple[np.ndarray, np.ndarray]:
        """Return reciprocal wave vectors and squared magnitudes.

        Only half of the symmetric :math:`\\pm\\mathbf k` pairs are needed for
        the energy, but the implementation keeps both signs because the force
        expression is simpler and less error-prone with the full set.
        """

        if not isinstance(box, Box):
            raise TypeError("box must be a Box.")

        lengths = box.lengths
        n_max = np.ceil(
            self.reciprocal_cutoff * lengths / (2.0 * np.pi)
        ).astype(int)
        ranges = [np.arange(-n, n + 1) for n in n_max]
        integers = np.stack(
            np.meshgrid(*ranges, indexing="ij"),
            axis=-1,
        ).reshape(-1, 3)
        wave_vectors = 2.0 * np.pi * integers / lengths
        squared = np.sum(wave_vectors * wave_vectors, axis=1)
        mask = (squared > 0.0) & (squared <= self.reciprocal_cutoff**2 + 1e-12)
        return wave_vectors[mask], squared[mask]

    def reciprocal_energy(
        self,
        positions: np.ndarray,
        charges: np.ndarray,
        box: Box,
    ) -> float:
        """Return the reciprocal-space energy in kJ/mol."""

        coordinates, charge_values = self._validate_system(positions, charges, box)
        wave_vectors, squared = self.reciprocal_k_vectors(box)
        if wave_vectors.size == 0:
            return 0.0

        if get_backend() == "torch":
            return self._reciprocal_energy_torch(
                coordinates,
                charge_values,
                wave_vectors,
                squared,
                box.volume,
            )

        phases = np.exp(1j * (coordinates @ wave_vectors.T))
        structure_factor = charge_values @ phases
        prefactor = (
            4.0
            * np.pi
            / squared
            * np.exp(-squared / (4.0 * self.alpha**2))
        )
        return float(
            COULOMB_CONSTANT
            / (2.0 * box.volume)
            * np.sum(prefactor * np.abs(structure_factor) ** 2)
        )

    def _reciprocal_energy_torch(
        self,
        coordinates: np.ndarray,
        charge_values: np.ndarray,
        wave_vectors: np.ndarray,
        squared: np.ndarray,
        volume: float,
    ) -> float:
        torch = _load_torch()
        dtype = torch_dtype()
        t_positions = to_torch(coordinates, dtype=dtype)
        t_charges = to_torch(charge_values, dtype=dtype)
        t_wave_vectors = to_torch(wave_vectors, dtype=dtype)
        t_squared = to_torch(squared, dtype=dtype)

        angles = t_positions @ t_wave_vectors.T
        phases = torch.complex(torch.cos(angles), torch.sin(angles))
        structure_factor = torch.complex(
            t_charges, torch.zeros_like(t_charges)
        ) @ phases
        prefactor = (
            4.0
            * np.pi
            / t_squared
            * torch.exp(-t_squared / (4.0 * self.alpha**2))
        )
        energy = (
            COULOMB_CONSTANT
            / (2.0 * volume)
            * torch.sum(prefactor * torch.abs(structure_factor) ** 2)
        )
        return float(energy.detach().cpu().item())

    def self_energy(self, charges: np.ndarray) -> float:
        """Return the Ewald self-energy correction in kJ/mol."""

        charge_values = np.asarray(charges, dtype=float)
        if charge_values.ndim != 1 or not np.all(np.isfinite(charge_values)):
            raise ValueError("charges must be a finite one-dimensional array.")
        return float(
            -COULOMB_CONSTANT
            * self.alpha
            / np.sqrt(np.pi)
            * np.sum(charge_values * charge_values)
        )

    def charge_correction(self, charges: np.ndarray, box: Box) -> float:
        """Return the net-charge correction in kJ/mol.

        This term accounts for a uniform neutralizing background and vanishes
        for a neutral system.
        """

        charge_values = np.asarray(charges, dtype=float)
        if charge_values.ndim != 1 or not np.all(np.isfinite(charge_values)):
            raise ValueError("charges must be a finite one-dimensional array.")
        if not isinstance(box, Box):
            raise TypeError("box must be a Box.")
        total_charge = float(np.sum(charge_values))
        return float(
            -COULOMB_CONSTANT
            * np.pi
            / (2.0 * self.alpha**2 * box.volume)
            * total_charge**2
        )

    def global_energy(
        self,
        positions: np.ndarray,
        charges: np.ndarray,
        box: Box,
    ) -> float:
        """Return reciprocal + self + net-charge energy in kJ/mol.

        The real-space part is intentionally excluded; it is evaluated by the
        pair loop in :class:`~moleculab.nonbonded.NonbondedModel`.
        """

        return float(
            self.reciprocal_energy(positions, charges, box)
            + self.self_energy(charges)
            + self.charge_correction(charges, box)
        )

    def global_forces(
        self,
        positions: np.ndarray,
        charges: np.ndarray,
        box: Box,
    ) -> np.ndarray:
        """Return reciprocal-space forces with shape ``(n_particles, 3)``.

        The self and net-charge terms are position-independent and therefore
        contribute no force.
        """

        coordinates, charge_values = self._validate_system(positions, charges, box)
        wave_vectors, squared = self.reciprocal_k_vectors(box)
        forces = np.zeros_like(coordinates)
        if wave_vectors.size == 0:
            return forces

        if get_backend() == "torch":
            return self._global_forces_torch(
                coordinates,
                charge_values,
                wave_vectors,
                squared,
                box.volume,
            )

        phases = np.exp(1j * (coordinates @ wave_vectors.T))
        structure_factor = charge_values @ phases
        prefactor = (
            4.0
            * np.pi
            / squared
            * np.exp(-squared / (4.0 * self.alpha**2))
        )
        imaginary = np.imag(
            structure_factor[np.newaxis, :] * np.conj(phases)
        )
        forces = (
            -COULOMB_CONSTANT
            / box.volume
            * charge_values[:, np.newaxis]
            * np.sum(
                prefactor[np.newaxis, :, np.newaxis]
                * wave_vectors[np.newaxis, :, :]
                * imaginary[:, :, np.newaxis],
                axis=1,
            )
        )
        return forces

    def _global_forces_torch(
        self,
        coordinates: np.ndarray,
        charge_values: np.ndarray,
        wave_vectors: np.ndarray,
        squared: np.ndarray,
        volume: float,
    ) -> np.ndarray:
        torch = _load_torch()
        dtype = torch_dtype()
        t_positions = to_torch(coordinates, dtype=dtype)
        t_charges = to_torch(charge_values, dtype=dtype)
        t_wave_vectors = to_torch(wave_vectors, dtype=dtype)
        t_squared = to_torch(squared, dtype=dtype)

        angles = t_positions @ t_wave_vectors.T
        phases = torch.complex(torch.cos(angles), torch.sin(angles))
        structure_factor = torch.complex(
            t_charges, torch.zeros_like(t_charges)
        ) @ phases
        prefactor = (
            4.0
            * np.pi
            / t_squared
            * torch.exp(-t_squared / (4.0 * self.alpha**2))
        )
        imaginary = torch.imag(
            structure_factor.unsqueeze(0) * torch.conj(phases)
        )
        forces = (
            -COULOMB_CONSTANT
            / volume
            * t_charges.unsqueeze(1)
            * torch.sum(
                prefactor.unsqueeze(0).unsqueeze(-1)
                * t_wave_vectors.unsqueeze(0)
                * imaginary.unsqueeze(-1),
                dim=1,
            )
        )
        return to_numpy(forces)

    def exclusion_energy(
        self,
        positions: np.ndarray,
        charges: np.ndarray,
        box: Box,
        excluded_pairs: Sequence[tuple[int, int]],
    ) -> float:
        """Return the reciprocal-space correction for excluded pairs.

        The real-space pair terms are omitted by the pair loop.  This method
        subtracts the remaining reciprocal contribution of each excluded pair:

        .. math::

            U_{\\mathrm{excl}} = -\\sum_{(i,j)}
            q_i q_j \\frac{k_e}{V}\\sum_{\\mathbf k \\ne 0}
            \\frac{4\\pi}{k^2}e^{-k^2/(4\\alpha^2)}
            \\cos(\\mathbf k\\cdot\\mathbf r_{ij}).

        Self-image interactions are not removed; they are part of the Ewald
        energy even when intramolecular pairs are excluded.
        """

        coordinates, charge_values = self._validate_system(positions, charges, box)
        pairs = tuple(excluded_pairs)
        if not pairs:
            return 0.0

        pair_indices = np.asarray(pairs, dtype=int)
        if pair_indices.ndim != 2 or pair_indices.shape[1] != 2:
            raise ValueError("excluded_pairs must contain (i, j) index pairs.")
        if np.any(pair_indices < 0) or np.any(pair_indices >= coordinates.shape[0]):
            raise ValueError("excluded_pairs contains an out-of-range atom index.")

        displacement = (
            coordinates[pair_indices[:, 0]] - coordinates[pair_indices[:, 1]]
        )
        displacement = box.minimum_image_array(displacement)
        wave_vectors, squared = self.reciprocal_k_vectors(box)
        if wave_vectors.size == 0:
            return 0.0

        prefactor = (
            4.0
            * np.pi
            / squared
            * np.exp(-squared / (4.0 * self.alpha**2))
        )
        cosine = np.cos(displacement @ wave_vectors.T)
        green = (COULOMB_CONSTANT / box.volume) * (cosine @ prefactor)
        pair_charges = (
            charge_values[pair_indices[:, 0]] * charge_values[pair_indices[:, 1]]
        )
        return float(-np.sum(pair_charges * green))

    def exclusion_forces(
        self,
        positions: np.ndarray,
        charges: np.ndarray,
        box: Box,
        excluded_pairs: Sequence[tuple[int, int]],
    ) -> np.ndarray:
        """Return reciprocal-space force corrections for excluded pairs."""

        coordinates, charge_values = self._validate_system(positions, charges, box)
        pairs = tuple(excluded_pairs)
        forces = np.zeros_like(coordinates)
        if not pairs:
            return forces

        pair_indices = np.asarray(pairs, dtype=int)
        if pair_indices.ndim != 2 or pair_indices.shape[1] != 2:
            raise ValueError("excluded_pairs must contain (i, j) index pairs.")
        if np.any(pair_indices < 0) or np.any(pair_indices >= coordinates.shape[0]):
            raise ValueError("excluded_pairs contains an out-of-range atom index.")

        displacement = (
            coordinates[pair_indices[:, 0]] - coordinates[pair_indices[:, 1]]
        )
        displacement = box.minimum_image_array(displacement)
        wave_vectors, squared = self.reciprocal_k_vectors(box)
        if wave_vectors.size == 0:
            return forces

        prefactor = (
            4.0
            * np.pi
            / squared
            * np.exp(-squared / (4.0 * self.alpha**2))
        )
        sine = np.sin(displacement @ wave_vectors.T)
        weighted_sine = sine * prefactor[np.newaxis, :]
        pair_charges = (
            charge_values[pair_indices[:, 0]] * charge_values[pair_indices[:, 1]]
        )
        pair_force = (
            -COULOMB_CONSTANT
            / box.volume
            * pair_charges[:, np.newaxis]
            * np.sum(
                weighted_sine[:, :, np.newaxis]
                * wave_vectors[np.newaxis, :, :],
                axis=1,
            )
        )
        np.add.at(forces, pair_indices[:, 0], pair_force)
        np.add.at(forces, pair_indices[:, 1], -pair_force)
        return forces

    def _validate_system(
        self,
        positions: np.ndarray,
        charges: np.ndarray,
        box: Box,
    ) -> tuple[np.ndarray, np.ndarray]:
        if not isinstance(box, Box):
            raise TypeError("box must be a Box.")
        coordinates = np.asarray(positions, dtype=float)
        charge_values = np.asarray(charges, dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError(
                "positions must have shape (n_particles, 3), "
                f"got {coordinates.shape}."
            )
        if charge_values.shape != (coordinates.shape[0],):
            raise ValueError(
                "charges must have shape (n_particles,), "
                f"got {charge_values.shape}."
            )
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("positions must contain only finite values.")
        if not np.all(np.isfinite(charge_values)):
            raise ValueError("charges must contain only finite values.")
        return coordinates, charge_values
