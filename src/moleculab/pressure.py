"""Pressure and the Berendsen barostat.

The instantaneous pressure has an ideal-gas kinetic contribution and a virial
contribution from interactions:

.. math::

    P = \\frac{2K}{3V} + P_{\\mathrm{virial}}.

Moleculab estimates :math:`P_{\\mathrm{virial}}` by an isotropic finite-
difference volume change.  All particle positions and the periodic box are
scaled together, so bond lengths, angles, nonbonded distances, and reciprocal
vectors all respond consistently:

.. math::

    P_{\\mathrm{virial}} = -\\frac{\\partial U}{\\partial V}.

This is a transparent and general way to include every interaction in the
virial, including Ewald electrostatics.  It is more expensive than an analytic
virial but much less error-prone for a library that is still growing.

The Berendsen barostat applies the volume scale factor

.. math::

    \\mu = \\left[1 - \\frac{\\kappa \\Delta t}{\\tau_P}
    (P_0 - P)\\right]^{1/3},

where :math:`\\kappa` is the isothermal compressibility in bar^-1,
:math:`\\tau_P` is the pressure coupling time in ps, and pressures are in bar.
"""

from __future__ import annotations

import numpy as np

from ._validation import as_nonnegative_float, as_positive_float
from .box import Box
from .constraints import assert_constraints_satisfied, shake_positions
from .particles import Particle
from .system import System
from .units import KILOJOULES_PER_MOLE_PER_NM3_TO_BAR


def _require_box(system: System) -> Box:
    if not isinstance(system, System):
        raise TypeError("system must be a System.")
    if system.box is None:
        raise ValueError("Pressure coupling requires a periodic Box.")
    return system.box


def _scaled_system(system: System, scale: float) -> System:
    """Return a copy with positions and box scaled isotropically."""

    box = _require_box(system)
    scaled_positions = np.array(
        [particle.position * scale for particle in system.particles],
        dtype=float,
    )
    if system.constraints:
        masses = np.array([particle.mass for particle in system.particles])
        scaled_positions = shake_positions(
            scaled_positions,
            masses,
            system.constraints,
        ).positions
    scaled_particles = [
        Particle(
            mass=particle.mass,
            position=scaled_positions[index],
            velocity=particle.velocity,
            atom_type=particle.atom_type,
            charge=particle.charge,
        )
        for index, particle in enumerate(system.particles)
    ]
    updated = System(
        particles=scaled_particles,
        pair_potentials=system.pair_potentials,
        topology=system.topology,
        box=Box(box.lengths * scale),
        nonbonded_model=system.nonbonded_model,
        constraints=system.constraints,
        indexed_pair_potentials=system.indexed_pair_potentials,
    )
    assert_constraints_satisfied(updated)
    return updated


def kinetic_pressure(system: System) -> float:
    """Return the kinetic pressure in kJ/(mol nm^3)."""

    box = _require_box(system)
    return float(2.0 * system.kinetic_energy() / (3.0 * box.volume))


def virial_pressure(system: System, strain: float = 1e-5) -> float:
    """Return the potential virial pressure in kJ/(mol nm^3).

    A central finite difference is used:

    .. math::

        P_{\\mathrm{virial}} \\approx
        -\\frac{U(V+\\Delta V)-U(V-\\Delta V)}{2\\Delta V}.
    """

    box = _require_box(system)
    strain_value = as_positive_float(strain, "strain")
    if strain_value >= 1.0:
        raise ValueError("strain must be smaller than 1.")

    scale_plus = (1.0 + strain_value) ** (1.0 / 3.0)
    scale_minus = (1.0 - strain_value) ** (1.0 / 3.0)
    energy_plus = _scaled_system(system, scale_plus).potential_energy()
    volume_plus = box.volume * (1.0 + strain_value)
    volume_minus = box.volume * (1.0 - strain_value)

    # If shrinking would put the cutoff above half the box edge, fall back to
    # a one-sided expansion difference.  This happens for systems that already
    # sit exactly at cutoff = L/2, such as the default water box.
    scaled_minimum_edge = float(min(box.lengths * scale_minus))
    if (
        system.nonbonded_model is not None
        and system.nonbonded_model.cutoff
        > scaled_minimum_edge / 2.0 + 1e-12
    ):
        energy_zero = system.potential_energy()
        return float(-(energy_plus - energy_zero) / (volume_plus - box.volume))

    energy_minus = _scaled_system(system, scale_minus).potential_energy()
    return float(-(energy_plus - energy_minus) / (volume_plus - volume_minus))


def pressure(system: System, strain: float = 1e-5) -> float:
    """Return total pressure in kJ/(mol nm^3)."""

    return kinetic_pressure(system) + virial_pressure(system, strain=strain)


def pressure_bar(system: System, strain: float = 1e-5) -> float:
    """Return total pressure in bar."""

    return float(
        pressure(system, strain=strain)
        * KILOJOULES_PER_MOLE_PER_NM3_TO_BAR
    )


def berendsen_barostat(
    system: System,
    target_pressure: float,
    coupling_time: float,
    dt: float,
    compressibility: float = 4.5e-5,
    strain: float = 1e-5,
) -> System:
    """Return a copy of *system* after one Berendsen pressure-coupling step.

    Parameters
    ----------
    system:
        Periodic system to scale.
    target_pressure:
        Target pressure in bar.
    coupling_time:
        Pressure coupling time :math:`\\tau_P` in ps.
    dt:
        Integration timestep in ps.
    compressibility:
        Isothermal compressibility :math:`\\kappa` in bar^-1.
    strain:
        Finite-difference volume strain used by :func:`virial_pressure`.

    Returns
    -------
    System
        A new system with positions and box scaled by the barostat factor.
    """

    _require_box(system)
    target = as_nonnegative_float(target_pressure, "target_pressure")
    tau = as_positive_float(coupling_time, "coupling_time")
    step = as_positive_float(dt, "dt")
    kappa = as_positive_float(compressibility, "compressibility")
    current_pressure = pressure_bar(system, strain=strain)

    scale_cubed = 1.0 - (kappa / tau) * step * (target - current_pressure)
    if scale_cubed <= 0.0:
        raise ValueError(
            "Barostat scale factor is non-positive; reduce dt or increase "
            "coupling_time."
        )

    scale = scale_cubed ** (1.0 / 3.0)
    return _scaled_system(system, scale)
