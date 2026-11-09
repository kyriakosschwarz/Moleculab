"""The mechanical state of a single particle.

A particle in Moleculab is deliberately small: it stores the three quantities
needed to describe classical motion at one instant,

    state = (mass, position, velocity).

Acceleration is *not* stored.  Newton's second law lets us calculate it from
the force and mass whenever it is needed.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ._validation import as_finite_float, as_positive_float, as_vector3


@dataclass
class Particle:
    """A point particle with mass, position, and velocity.

    Parameters
    ----------
    mass:
        Mass in Da.  Must be positive.
    position:
        Three-dimensional position in nm.
    velocity:
        Three-dimensional velocity in nm/ps.  Defaults to zero.

    Notes
    -----
    ``position`` and ``velocity`` are copied into NumPy arrays during
    construction, so later mutation of the input list cannot silently change a
    particle's state.
    """

    mass: float
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))
    velocity: np.ndarray = field(default_factory=lambda: np.zeros(3))
    atom_type: str = "default"
    charge: float = 0.0

    def __post_init__(self) -> None:
        self.mass = as_positive_float(self.mass, "mass")
        self.position = as_vector3(self.position, "position").copy()
        self.velocity = as_vector3(self.velocity, "velocity").copy()
        atom_type = str(self.atom_type).strip()
        if not atom_type:
            raise ValueError("atom_type must be a non-empty string.")
        self.atom_type = atom_type
        self.charge = as_finite_float(self.charge, "charge")

    @property
    def speed(self) -> float:
        """Speed in nm/ps."""

        return float(np.linalg.norm(self.velocity))

    def copy(self) -> "Particle":
        """Return an independent copy of the particle."""

        return Particle(
            self.mass,
            self.position.copy(),
            self.velocity.copy(),
            atom_type=self.atom_type,
            charge=self.charge,
        )

    def __repr__(self) -> str:
        position = np.array2string(self.position, precision=3)
        velocity = np.array2string(self.velocity, precision=3)
        return (
            f"Particle(mass={self.mass:g} Da, "
            f"position={position} nm, velocity={velocity} nm/ps, "
            f"atom_type={self.atom_type!r}, charge={self.charge:g} e)"
        )
