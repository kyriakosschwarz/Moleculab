"""Temperature and simple velocity-rescaling thermostats.

In the internal unit system, kinetic energy is in kJ/mol and the molar gas
constant is 0.008314462618 kJ/(mol K).  The instantaneous temperature is

.. math::

    T = \\frac{2K}{N_{\\mathrm{dof}} R},

where :math:`N_{\\mathrm{dof}}` is the number of degrees of freedom.  For a
system without external constraints, removing the three center-of-mass
translation modes gives :math:`N_{\\mathrm{dof}} = 3N - 3`.

A thermostat is an *algorithm* for controlling temperature.  It is not part of
the physical potential and is deliberately kept in its own module.
"""

from __future__ import annotations

import numpy as np

from ..composition.system import System
from ..foundation.units import BOLTZMANN_CONSTANT


def degrees_of_freedom(system: System, remove_com: bool = True) -> int:
    """Return the number of thermal degrees of freedom.

    The center of mass of an isolated system moves at constant velocity, so
    its three translational modes are normally removed from the temperature
    count.  For a single particle there is no meaningful relative motion, and
    all three components are retained.
    """

    n_particles = system.n_particles
    if remove_com and n_particles > 1:
        return 3 * n_particles - 3
    return 3 * n_particles


def kinetic_temperature(system: System, remove_com: bool = True) -> float:
    """Return the instantaneous kinetic temperature in K."""

    dof = degrees_of_freedom(system, remove_com=remove_com)
    kinetic = system.kinetic_energy()
    return float(2.0 * kinetic / (dof * BOLTZMANN_CONSTANT))


def _with_velocities(system: System, velocities: np.ndarray) -> System:
    return system.with_state(system.positions, velocities)


def berendsen_thermostat(
    system: System,
    temperature: float,
    coupling_time: float,
    dt: float,
    remove_com: bool = True,
) -> System:
    """Return a copy of *system* with velocities weakly coupled to *temperature*.

    The Berendsen scale factor is

    .. math::

        \\lambda = \\sqrt{1 + \\frac{\\Delta t}{\\tau}
        \\left(\\frac{T_0}{T} - 1\\right)}.

    ``temperature`` is the target :math:`T_0` in K, ``coupling_time`` is
    :math:`\\tau` in ps, and ``dt`` is the integration step in ps.  The input
    system is not modified.

    A target temperature of zero is handled by setting all velocities to zero.
    """

    target = float(temperature)
    if not np.isfinite(target) or target < 0.0:
        raise ValueError("temperature must be finite and non-negative.")

    tau = float(coupling_time)
    if not np.isfinite(tau) or tau <= 0.0:
        raise ValueError("coupling_time must be finite and positive.")

    step = float(dt)
    if not np.isfinite(step) or step <= 0.0:
        raise ValueError("dt must be finite and positive.")

    current_temperature = kinetic_temperature(system, remove_com=remove_com)

    if target == 0.0:
        return _with_velocities(system, np.zeros_like(system.velocities))

    if current_temperature == 0.0:
        raise ValueError(
            "Cannot couple a system with zero kinetic energy to a non-zero "
            "target temperature."
        )

    scale_squared = 1.0 + (step / tau) * (target / current_temperature - 1.0)
    scale = float(np.sqrt(max(0.0, scale_squared)))
    return _with_velocities(system, system.velocities * scale)
