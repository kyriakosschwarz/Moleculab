"""Molecular connectivity and bonded energy terms.

A force field needs to know which atoms are connected.  That information is
the *topology*: bonds, angles, and later torsions.  Keeping topology separate
from particle coordinates makes it possible to reuse the same molecular graph
with different conformations.

Bonded terms are indexed by atom number.  For example, an angle ``(0, 1, 2)``
means atom 1 is the central atom of the angle formed by atoms 0-1-2.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field

import numpy as np

from .bonded import AnglePotential, BondPotential, TorsionPotential


def _atom_index(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{name} must be an integer atom index.")
    index = int(value)
    if index < 0:
        raise ValueError(f"{name} must be non-negative, got {index}.")
    return index


@dataclass(frozen=True)
class Bond:
    """A bonded interaction between atoms ``i`` and ``j``."""

    i: int
    j: int
    potential: BondPotential

    def __post_init__(self) -> None:
        i = _atom_index(self.i, "i")
        j = _atom_index(self.j, "j")
        if i == j:
            raise ValueError("A bond must connect two different atoms.")
        object.__setattr__(self, "i", i)
        object.__setattr__(self, "j", j)


@dataclass(frozen=True)
class Angle:
    """An angle interaction with central atom ``j``.

    The three atoms must be distinct; ``j`` is the vertex.
    """

    i: int
    j: int
    k: int
    potential: AnglePotential

    def __post_init__(self) -> None:
        i = _atom_index(self.i, "i")
        j = _atom_index(self.j, "j")
        k = _atom_index(self.k, "k")
        if len({i, j, k}) != 3:
            raise ValueError("An angle must involve three different atom indices.")
        object.__setattr__(self, "i", i)
        object.__setattr__(self, "j", j)
        object.__setattr__(self, "k", k)


@dataclass(frozen=True)
class Torsion:
    """A proper torsion interaction along the bonded path i-j-k-l."""

    i: int
    j: int
    k: int
    l: int
    potential: TorsionPotential

    def __post_init__(self) -> None:
        i = _atom_index(self.i, "i")
        j = _atom_index(self.j, "j")
        k = _atom_index(self.k, "k")
        l = _atom_index(self.l, "l")
        if len({i, j, k, l}) != 4:
            raise ValueError("A torsion must involve four different atom indices.")
        object.__setattr__(self, "i", i)
        object.__setattr__(self, "j", j)
        object.__setattr__(self, "k", k)
        object.__setattr__(self, "l", l)


@dataclass(frozen=True)
class Topology:
    """The bonded terms of a molecular system.

    Examples
    --------
    >>> import numpy as np
    >>> from moleculab import (
    ...     Angle, Bond, HarmonicAnglePotential, HarmonicBondPotential,
    ...     PeriodicTorsionPotential, Torsion
    ... )
    >>> topology = Topology(
    ...     bonds=[Bond(0, 1, HarmonicBondPotential(1000.0, 0.15))],
    ...     angles=[Angle(0, 1, 2, HarmonicAnglePotential(100.0, np.deg2rad(109.5)))],
    ...     torsions=[Torsion(0, 1, 2, 3, PeriodicTorsionPotential(2.0, 3))],
    ... )
    """

    bonds: Sequence[Bond] = field(default_factory=tuple)
    angles: Sequence[Angle] = field(default_factory=tuple)
    torsions: Sequence[Torsion] = field(default_factory=tuple)

    def __post_init__(self) -> None:
        object.__setattr__(self, "bonds", tuple(self.bonds))
        object.__setattr__(self, "angles", tuple(self.angles))
        object.__setattr__(self, "torsions", tuple(self.torsions))
        for index, bond in enumerate(self.bonds):
            if not isinstance(bond, Bond):
                raise TypeError(f"bonds[{index}] must be a Bond.")
        for index, angle in enumerate(self.angles):
            if not isinstance(angle, Angle):
                raise TypeError(f"angles[{index}] must be an Angle.")
        for index, torsion in enumerate(self.torsions):
            if not isinstance(torsion, Torsion):
                raise TypeError(f"torsions[{index}] must be a Torsion.")

    @property
    def n_bonds(self) -> int:
        return len(self.bonds)

    @property
    def n_angles(self) -> int:
        return len(self.angles)

    @property
    def n_torsions(self) -> int:
        return len(self.torsions)

    @property
    def max_atom_index(self) -> int:
        """Largest atom index referenced by the topology, or -1 if empty."""

        indices = [bond.i for bond in self.bonds]
        indices.extend(bond.j for bond in self.bonds)
        indices.extend(angle.i for angle in self.angles)
        indices.extend(angle.j for angle in self.angles)
        indices.extend(angle.k for angle in self.angles)
        indices.extend(torsion.i for torsion in self.torsions)
        indices.extend(torsion.j for torsion in self.torsions)
        indices.extend(torsion.k for torsion in self.torsions)
        indices.extend(torsion.l for torsion in self.torsions)
        return max(indices, default=-1)

    def validate(self, n_atoms: int) -> None:
        """Raise if any term refers to an atom outside ``range(n_atoms)``."""

        if self.max_atom_index >= n_atoms:
            raise ValueError(
                f"Topology refers to atom {self.max_atom_index}, but the system "
                f"only has {n_atoms} atoms."
            )

    def nonbonded_exclusions(self) -> frozenset[tuple[int, int]]:
        """Return 1-2 and 1-3 atom pairs that should not interact nonbonded.

        Bonds exclude their two endpoints (1-2).  Angles exclude the two outer
        atoms (1-3) because they are connected through the central atom.
        Torsion 1-4 pairs are not excluded here; force fields often scale them
        rather than removing them completely.
        """

        excluded: set[tuple[int, int]] = set()
        for bond in self.bonds:
            excluded.add((min(bond.i, bond.j), max(bond.i, bond.j)))
        for angle in self.angles:
            excluded.add((min(angle.i, angle.k), max(angle.i, angle.k)))
        return frozenset(excluded)

    def energy(self, positions: np.ndarray) -> float:
        """Return total bonded energy in kJ/mol."""

        coordinates = self._as_positions(positions)
        total = 0.0

        for bond in self.bonds:
            displacement = coordinates[bond.i] - coordinates[bond.j]
            total += bond.potential.energy(displacement)

        for angle in self.angles:
            r_ji = coordinates[angle.i] - coordinates[angle.j]
            r_jk = coordinates[angle.k] - coordinates[angle.j]
            total += angle.potential.energy(r_ji, r_jk)

        for torsion in self.torsions:
            r_ij = coordinates[torsion.j] - coordinates[torsion.i]
            r_jk = coordinates[torsion.k] - coordinates[torsion.j]
            r_kl = coordinates[torsion.l] - coordinates[torsion.k]
            total += torsion.potential.energy(r_ij, r_jk, r_kl)

        return float(total)

    def forces(self, positions: np.ndarray) -> np.ndarray:
        """Return bonded forces in kJ/(mol nm) with shape ``(n_atoms, 3)``."""

        coordinates = self._as_positions(positions)
        forces = np.zeros_like(coordinates, dtype=float)

        for bond in self.bonds:
            displacement = coordinates[bond.i] - coordinates[bond.j]
            pair_force = bond.potential.force(displacement)
            forces[bond.i] += pair_force
            forces[bond.j] -= pair_force

        for angle in self.angles:
            r_ji = coordinates[angle.i] - coordinates[angle.j]
            r_jk = coordinates[angle.k] - coordinates[angle.j]
            force_i, force_j, force_k = angle.potential.forces(r_ji, r_jk)
            forces[angle.i] += force_i
            forces[angle.j] += force_j
            forces[angle.k] += force_k

        for torsion in self.torsions:
            r_ij = coordinates[torsion.j] - coordinates[torsion.i]
            r_jk = coordinates[torsion.k] - coordinates[torsion.j]
            r_kl = coordinates[torsion.l] - coordinates[torsion.k]
            force_i, force_j, force_k, force_l = torsion.potential.forces(
                r_ij, r_jk, r_kl
            )
            forces[torsion.i] += force_i
            forces[torsion.j] += force_j
            forces[torsion.k] += force_k
            forces[torsion.l] += force_l

        return forces

    def _as_positions(self, positions: np.ndarray) -> np.ndarray:
        coordinates = np.asarray(positions, dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError(
                "positions must have shape (n_atoms, 3), "
                f"got {coordinates.shape}."
            )
        if not np.all(np.isfinite(coordinates)):
            raise ValueError("positions must contain only finite values.")
        self.validate(coordinates.shape[0])
        return coordinates
