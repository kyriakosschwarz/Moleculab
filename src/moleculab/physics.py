"""General physics calculations.

These functions implement textbook equations independently of any particular
potential.  They answer questions such as:

* How much energy is stored in motion?
* How much energy is stored in a configuration?
* What is the total mechanical energy?
* What acceleration does a force produce?

The separation between these functions and the potential models is the main
architectural idea in Moleculab.
"""

from __future__ import annotations

import numpy as np

from ._validation import as_vector3
from .particles import Particle
from .potentials import Potential


def kinetic_energy(particle: Particle) -> float:
    """Return kinetic energy in kJ/mol.

    .. math:: K = \\frac{1}{2} m \\mathbf v \\cdot \\mathbf v
    """

    velocity_squared = float(np.dot(particle.velocity, particle.velocity))
    return 0.5 * particle.mass * velocity_squared


def potential_energy(particle: Particle, potential: Potential) -> float:
    """Return the potential energy of *particle* in kJ/mol."""

    return potential.energy(particle.position)


def total_energy(particle: Particle, potential: Potential) -> float:
    """Return total mechanical energy, :math:`E = K + U`, in kJ/mol."""

    return kinetic_energy(particle) + potential_energy(particle, potential)


def acceleration(particle: Particle, force: np.ndarray) -> np.ndarray:
    """Return acceleration in nm/ps^2 from Newton's second law.

    .. math:: \\mathbf a = \\frac{\\mathbf F}{m}
    """

    force_vector = as_vector3(force, "force")
    return force_vector / particle.mass


def force_from_potential(particle: Particle, potential: Potential) -> np.ndarray:
    """Return the force on *particle* in kJ/(mol nm)."""

    return potential.force(particle.position)
