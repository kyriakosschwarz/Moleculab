# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Physical units and constants used throughout Moleculab.

Moleculab uses one consistent internal unit system.  Keeping the numbers in
these units avoids conversion factors in the physics code and makes the
equations look the way they are usually written in molecular simulation.

======================  ======================  ==========================
Quantity                Unit                    Meaning
======================  ======================  ==========================
distance                nm                      nanometre
time                    ps                      picosecond
energy                  kJ/mol                  kilojoule per mole
mass                    Da                      dalton (atomic mass unit)
force                   kJ/(mol nm)             energy per distance
temperature             K                       kelvin
charge                  e                       elementary charge
======================  ======================  ==========================

A useful consequence of this choice is that Newton's second law can be
evaluated numerically as ``a = F / m`` without an extra constant.  This is
why the examples and tests can use forces in kJ/(mol nm), masses in Da, and
obtain accelerations in nm/ps^2 directly.
"""

from __future__ import annotations

DISTANCE_UNIT = "nm"
TIME_UNIT = "ps"
ENERGY_UNIT = "kJ/mol"
MASS_UNIT = "Da"
FORCE_UNIT = "kJ/(mol nm)"
TEMPERATURE_UNIT = "K"
CHARGE_UNIT = "e"

#: Molar gas constant, R = N_A * k_B, in kJ/(mol K).
#:
#: The historical name in this project is ``BOLTZMANN_CONSTANT``.  The value
#: is the *molar* constant because Moleculab expresses energy per mole.  It is
#: numerically equal to 0.008314462618 kJ/(mol K).
BOLTZMANN_CONSTANT = 0.008314462618

#: Descriptive alias for :data:`BOLTZMANN_CONSTANT`.
#:
#: ``BOLTZMANN_CONSTANT`` is retained for compatibility, but its value is the
#: *molar* gas constant ``R`` in kJ/(mol K).  New code should prefer this name.
MOLAR_GAS_CONSTANT = BOLTZMANN_CONSTANT

#: Avogadro constant in 1/mol.
AVOGADRO_CONSTANT = 6.02214076e23

#: Single-particle Boltzmann constant in J/K.
BOLTZMANN_CONSTANT_J_PER_K = 1.380649e-23

#: Conversion factors for moving between SI and the internal unit system.
NANOMETRE_IN_METRES = 1.0e-9
PICOSECOND_IN_SECONDS = 1.0e-12
DALTON_IN_KILOGRAMS = 1.66053906660e-27
KILOJOULES_PER_MOLE_TO_JOULES_PER_PARTICLE = 1000.0 / AVOGADRO_CONSTANT

#: Planck constant and reduced Planck constant in J s.
PLANCK_CONSTANT_J_S = 6.62607015e-34
H_BAR_J_S = 1.054571817e-34

#: Electron mass in kg and Da.
ELECTRON_MASS_KG = 9.1093837015e-31
ELECTRON_MASS_DA = ELECTRON_MASS_KG / DALTON_IN_KILOGRAMS

#: Pressure unit used by the virial expression: kJ/(mol nm^3).
PRESSURE_UNIT = "kJ/(mol nm^3)"

#: Convert one internal pressure unit to bar.
#:
#: 1 kJ/mol/nm^3 = (1000 J / N_A) / (1e-27 m^3) / 1e5 Pa/bar.
KILOJOULES_PER_MOLE_PER_NM3_TO_BAR = (
    (1000.0 / AVOGADRO_CONSTANT)
    / NANOMETRE_IN_METRES**3
    / 1.0e5
)

#: Elementary charge in coulomb.
ELEMENTARY_CHARGE_COULOMBS = 1.602176634e-19

#: Vacuum permittivity in F/m.
VACUUM_PERMITTIVITY_F_PER_M = 8.8541878128e-12

#: Coulomb constant in Moleculab units: kJ nm / (mol e^2).
#:
#: The SI expression for the energy of two charges is
#:
#:     U = (1 / (4 pi eps0)) * (q1 e)(q2 e) / r_SI.
#:
#: Converting r from nm to m and the result from J per pair to kJ per mole
#: gives the factor below.  With it, the internal-unit equation is simply
#:
#:     U = COULOMB_CONSTANT * q1 * q2 / r
#:
#: where r is in nm and the charges are in units of the elementary charge.
COULOMB_CONSTANT = (
    (1.0 / (4.0 * 3.141592653589793 * VACUUM_PERMITTIVITY_F_PER_M))
    * ELEMENTARY_CHARGE_COULOMBS**2
    * (1.0 / NANOMETRE_IN_METRES)
    * AVOGADRO_CONSTANT
    / 1000.0
)
