"""Typed nonbonded interactions with intramolecular exclusions.

A molecular force field cannot use one pair of Lennard-Jones parameters for
every atom.  Oxygen and hydrogen in water, for example, need different
parameters and charges, and bonded atoms must not interact through the
nonbonded model.

This module adds the minimal machinery for that:

* typed Lennard-Jones parameters with Lorentz-Berthelot combining;
* per-particle charges (stored on :class:`~moleculab.particles.Particle`);
* 1-2 and 1-3 exclusions derived from the bonded topology;
* reaction-field electrostatics with a spherical cutoff.

The reaction-field potential inside the cutoff is

.. math::

    U(r) = k_e q_i q_j
    \\left[\\frac{1}{r}
    + \\frac{\\epsilon_{rf}-1}{2\\epsilon_{rf}+1}
      \\frac{r^2}{R_c^3}
    - \\frac{3\\epsilon_{rf}}{2\\epsilon_{rf}+1}
      \\frac{1}{R_c}\\right],

which is shifted so that :math:`U(R_c)=0`.  The force follows from
:math:`\\mathbf F=-\\nabla U`.  This is a common minimal treatment for
explicit-solvent simulations.  It is not Ewald summation: long-range
electrostatics beyond the cutoff are represented only through the reaction
field.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field

import numpy as np

from ._validation import (
    as_displacement,
    as_nonnegative_float,
    as_positive_float,
)
from .ewald import EwaldElectrostatics
from .units import COULOMB_CONSTANT


@dataclass(frozen=True)
class LennardJonesParameters:
    """Lennard-Jones parameters for one atom type.

    ``epsilon`` is in kJ/mol and ``sigma`` is in nm.  A zero ``epsilon`` is
    allowed for united-atom or hydrogen sites that have no dispersion term.
    """

    epsilon: float
    sigma: float

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "epsilon", as_nonnegative_float(self.epsilon, "epsilon")
        )
        object.__setattr__(self, "sigma", as_positive_float(self.sigma, "sigma"))

    def energy(self, distance: float) -> float:
        """Return the Lennard-Jones energy in kJ/mol."""

        if self.epsilon == 0.0:
            return 0.0
        sigma_over_r_sixth = (self.sigma / distance) ** 6
        return float(
            4.0
            * self.epsilon
            * (sigma_over_r_sixth**2 - sigma_over_r_sixth)
        )

    def force(self, displacement: np.ndarray) -> np.ndarray:
        """Return the force on particle i in kJ/(mol nm)."""

        vector, distance = as_displacement(displacement)
        if self.epsilon == 0.0:
            return np.zeros(3)
        sigma_over_r_sixth = (self.sigma / distance) ** 6
        factor = (
            24.0
            * self.epsilon
            * (2.0 * sigma_over_r_sixth**2 - sigma_over_r_sixth)
            / distance**2
        )
        return factor * vector


def combine_lorentz_berthelot(
    first: LennardJonesParameters,
    second: LennardJonesParameters,
) -> LennardJonesParameters:
    """Combine two atom types with the Lorentz-Berthelot rules.

    .. math::

        \\sigma_{ij} = \\frac{\\sigma_i + \\sigma_j}{2},
        \\qquad
        \\epsilon_{ij} = \\sqrt{\\epsilon_i \\epsilon_j}.
    """

    return LennardJonesParameters(
        epsilon=float(np.sqrt(first.epsilon * second.epsilon)),
        sigma=0.5 * (first.sigma + second.sigma),
    )


@dataclass(frozen=True)
class ReactionFieldElectrostatics:
    """Reaction-field electrostatics with a spherical cutoff.

    ``dielectric`` is the dielectric constant outside the cutoff
    :math:`R_c`.  A value of 1 gives ordinary Coulomb interactions shifted to
    zero at the cutoff.
    """

    dielectric: float = 78.3

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "dielectric", as_positive_float(self.dielectric, "dielectric")
        )

    def energy(
        self,
        charge_i: float,
        charge_j: float,
        distance: float,
        cutoff: float,
    ) -> float:
        """Return the reaction-field pair energy in kJ/mol."""

        cutoff_value = as_positive_float(cutoff, "cutoff")
        if distance == 0.0:
            raise ValueError("Electrostatic energy is undefined at zero distance.")
        if distance >= cutoff_value:
            return 0.0

        prefactor = COULOMB_CONSTANT * float(charge_i) * float(charge_j)
        reaction = (
            (self.dielectric - 1.0)
            / (2.0 * self.dielectric + 1.0)
            / cutoff_value**3
        )
        shift = (
            3.0
            * self.dielectric
            / (2.0 * self.dielectric + 1.0)
            / cutoff_value
        )
        return float(prefactor * (1.0 / distance + reaction * distance**2 - shift))

    def force(
        self,
        charge_i: float,
        charge_j: float,
        displacement: np.ndarray,
        cutoff: float,
    ) -> np.ndarray:
        """Return the force on particle i in kJ/(mol nm)."""

        vector, distance = as_displacement(displacement)
        cutoff_value = as_positive_float(cutoff, "cutoff")
        if distance >= cutoff_value:
            return np.zeros(3)

        prefactor = COULOMB_CONSTANT * float(charge_i) * float(charge_j)
        reaction = (
            2.0
            * (self.dielectric - 1.0)
            / (2.0 * self.dielectric + 1.0)
            / cutoff_value**3
        )
        factor = prefactor * (1.0 / distance**3 - reaction)
        return factor * vector


@dataclass(frozen=True)
class NonbondedModel:
    """Typed Lennard-Jones plus optional reaction-field electrostatics.

    Parameters
    ----------
    parameters:
        Mapping from atom type name to :class:`LennardJonesParameters`.
    cutoff:
        Nonbonded cutoff in nm.  Both Lennard-Jones and electrostatic
        interactions are truncated at this distance.  The cutoff should not
        exceed half the shortest periodic box edge.
    electrostatics:
        Optional electrostatic model: :class:`ReactionFieldElectrostatics` or
        :class:`EwaldElectrostatics`.
    exclusions:
        Additional atom-index pairs to exclude, on top of the 1-2 and 1-3
        exclusions derived from the topology.
    """

    parameters: Mapping[str, LennardJonesParameters]
    cutoff: float
    electrostatics: ReactionFieldElectrostatics | EwaldElectrostatics | None = None
    exclusions: frozenset[tuple[int, int]] = field(default_factory=frozenset)

    def __post_init__(self) -> None:
        normalized: dict[str, LennardJonesParameters] = {}
        for name, parameters in self.parameters.items():
            type_name = str(name).strip()
            if not type_name:
                raise ValueError("Atom type names must be non-empty strings.")
            if not isinstance(parameters, LennardJonesParameters):
                raise TypeError(
                    f"parameters[{type_name!r}] must be LennardJonesParameters."
                )
            normalized[type_name] = parameters
        if not normalized:
            raise ValueError("parameters must contain at least one atom type.")
        object.__setattr__(self, "parameters", normalized)
        object.__setattr__(self, "cutoff", as_positive_float(self.cutoff, "cutoff"))

        if self.electrostatics is not None and not isinstance(
            self.electrostatics,
            (ReactionFieldElectrostatics, EwaldElectrostatics),
        ):
            raise TypeError(
                "electrostatics must be ReactionFieldElectrostatics, "
                "EwaldElectrostatics, or None."
            )

        normalized_exclusions: set[tuple[int, int]] = set()
        for pair in self.exclusions:
            if len(pair) != 2:
                raise ValueError("Each exclusion must contain two atom indices.")
            i, j = int(pair[0]), int(pair[1])
            if i == j:
                raise ValueError("An exclusion cannot pair an atom with itself.")
            if i < 0 or j < 0:
                raise ValueError("Exclusion atom indices must be non-negative.")
            normalized_exclusions.add((min(i, j), max(i, j)))
        object.__setattr__(self, "exclusions", frozenset(normalized_exclusions))

    def validate_atom_types(self, atom_types: list[str]) -> None:
        """Raise if a particle type has no Lennard-Jones parameters."""

        missing = sorted(set(atom_types) - set(self.parameters))
        if missing:
            raise ValueError(
                "Missing Lennard-Jones parameters for atom types: "
                + ", ".join(missing)
            )

    def combined_parameters(
        self,
        type_i: str,
        type_j: str,
    ) -> LennardJonesParameters:
        """Return the Lorentz-Berthelot parameters for a pair of types."""

        try:
            first = self.parameters[type_i]
            second = self.parameters[type_j]
        except KeyError as exc:
            raise ValueError(f"No parameters for atom type {exc.args[0]!r}.") from exc
        return combine_lorentz_berthelot(first, second)

    def energy(
        self,
        type_i: str,
        type_j: str,
        charge_i: float,
        charge_j: float,
        displacement: np.ndarray,
    ) -> float:
        """Return nonbonded pair energy in kJ/mol."""

        _, distance = as_displacement(displacement)
        if distance >= self.cutoff:
            return 0.0

        total = self.combined_parameters(type_i, type_j).energy(distance)
        if self.electrostatics is not None:
            total += self.electrostatics.energy(
                charge_i,
                charge_j,
                distance,
                self.cutoff,
            )
        return float(total)

    def force(
        self,
        type_i: str,
        type_j: str,
        charge_i: float,
        charge_j: float,
        displacement: np.ndarray,
    ) -> np.ndarray:
        """Return the force on particle i in kJ/(mol nm)."""

        vector, distance = as_displacement(displacement)
        if distance >= self.cutoff:
            return np.zeros(3)

        force = self.combined_parameters(type_i, type_j).force(vector)
        if self.electrostatics is not None:
            force = force + self.electrostatics.force(
                charge_i,
                charge_j,
                vector,
                self.cutoff,
            )
        return force
