"""Moleculab: an educational molecular physics and simulation library.

The package is intentionally small and explicit.  Physical models (potentials)
are kept separate from general physics functions and from numerical
algorithms (integrators).

Examples
--------
>>> from moleculab import HarmonicPotential, Particle, total_energy
>>> particle = Particle(mass=12.0, position=[0.1, 0.0, 0.0], velocity=[0.1, 0.0, 0.0])
>>> potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
>>> round(total_energy(particle, potential), 6)
0.56
"""

from .analysis import (
    MSDResult,
    RDFResult,
    mean_squared_displacement,
    radial_distribution_function,
)
from .backend import (
    available_backends,
    backend_info,
    get_backend,
    get_precision,
    set_backend,
    set_precision,
    to_numpy,
    to_torch,
    torch_available,
    torch_device,
    torch_dtype,
    use_backend,
)
from .bonded import (
    AnglePotential,
    BondPotential,
    HarmonicAnglePotential,
    HarmonicBondPotential,
    PeriodicTorsionPotential,
    TorsionPotential,
    dihedral_angle,
)
from .box import Box
from .constraints import (
    DistanceConstraint,
    RattleResult,
    ShakeResult,
    assert_constraints_satisfied,
    constraint_error,
    rattle_velocities,
    shake_positions,
)
from .docking import (
    DockingPose,
    DockingResult,
    RigidLigand,
    dock_rigid_ligand,
    score_pose,
)
from .ewald import EwaldElectrostatics
from .fep import (
    FEPResult,
    fep_between_potentials,
    fep_harmonic_oscillator,
    free_energy_perturbation,
    harmonic_free_energy_difference,
    sample_harmonic_positions,
)
from .flexible_docking import (
    FlexibleDockingPose,
    FlexibleDockingResult,
    SideChain,
    dock_flexible_receptor,
    make_receptor_conformer,
    receptor_strain_score,
    score_flexible_pose,
)
from .integrators import (
    SystemTrajectory,
    Trajectory,
    simulate_velocity_verlet,
    simulate_velocity_verlet_system,
    velocity_verlet_step,
    velocity_verlet_step_system,
)
from .lj_alchemy import MultistageFEPResult, SoftCoreLennardJones
from .minimization import (
    MinimizationResult,
    max_force_magnitude,
    minimize_energy,
    steepest_descent_step,
)
from .mbar import MBARResult, mbar, mbar_harmonic
from .mc_alchemy import BondedLennardJonesMC, SoftCoreLennardJonesPair
from .neighbour import NeighbourList, build_neighbour_list
from .nonbonded import (
    LennardJonesParameters,
    NonbondedModel,
    ReactionFieldElectrostatics,
    combine_lorentz_berthelot,
)
from .particles import Particle
from .physics import (
    acceleration,
    force_from_potential,
    kinetic_energy,
    potential_energy,
    total_energy,
)
from .potentials import (
    CoulombPotential,
    DoubleWellPotential,
    HarmonicPotential,
    LennardJonesPotential,
    PairPotential,
    Potential,
)
from .pressure import (
    berendsen_barostat,
    kinetic_pressure,
    pressure,
    pressure_bar,
    virial_pressure,
)
from .qm_mm import MMPointCharge, ParticleInBox, QMResult
from .system import System
from .solvent import (
    ARGON,
    COARSE_GRAINED_WATER,
    SolventModel,
    build_solvent_box,
)
from .solvated_alchemy import SolvatedLigandAlchemy
from .ti import (
    TIResult,
    integrate_thermodynamic_derivative,
    thermodynamic_integration_harmonic,
)
from .thermostats import berendsen_thermostat, degrees_of_freedom, kinetic_temperature
from .topology import Angle, Bond, Torsion, Topology
from .units import ELECTRON_MASS_DA
from .water import SPCE_WATER, WaterModel, build_water_box

__all__ = [
    "ARGON",
    "Angle",
    "AnglePotential",
    "Bond",
    "BondPotential",
    "BondedLennardJonesMC",
    "Box",
    "COARSE_GRAINED_WATER",
    "CoulombPotential",
    "DistanceConstraint",
    "DockingPose",
    "DockingResult",
    "DoubleWellPotential",
    "EwaldElectrostatics",
    "ELECTRON_MASS_DA",
    "FEPResult",
    "FlexibleDockingPose",
    "FlexibleDockingResult",
    "HarmonicAnglePotential",
    "HarmonicBondPotential",
    "HarmonicPotential",
    "LennardJonesPotential",
    "LennardJonesParameters",
    "MBARResult",
    "MMPointCharge",
    "MSDResult",
    "MinimizationResult",
    "MultistageFEPResult",
    "NeighbourList",
    "NonbondedModel",
    "PairPotential",
    "Particle",
    "ParticleInBox",
    "PeriodicTorsionPotential",
    "Potential",
    "QMResult",
    "RDFResult",
    "RattleResult",
    "ReactionFieldElectrostatics",
    "RigidLigand",
    "SPCE_WATER",
    "ShakeResult",
    "SideChain",
    "SoftCoreLennardJones",
    "SoftCoreLennardJonesPair",
    "SolvatedLigandAlchemy",
    "SolventModel",
    "System",
    "SystemTrajectory",
    "TIResult",
    "Topology",
    "Torsion",
    "TorsionPotential",
    "Trajectory",
    "WaterModel",
    "acceleration",
    "assert_constraints_satisfied",
    "available_backends",
    "backend_info",
    "berendsen_barostat",
    "berendsen_thermostat",
    "build_neighbour_list",
    "build_solvent_box",
    "build_water_box",
    "combine_lorentz_berthelot",
    "constraint_error",
    "degrees_of_freedom",
    "dihedral_angle",
    "dock_flexible_receptor",
    "dock_rigid_ligand",
    "fep_between_potentials",
    "fep_harmonic_oscillator",
    "force_from_potential",
    "free_energy_perturbation",
    "get_backend",
    "get_precision",
    "harmonic_free_energy_difference",
    "integrate_thermodynamic_derivative",
    "kinetic_energy",
    "kinetic_pressure",
    "kinetic_temperature",
    "make_receptor_conformer",
    "max_force_magnitude",
    "mbar",
    "mbar_harmonic",
    "mean_squared_displacement",
    "minimize_energy",
    "potential_energy",
    "pressure",
    "pressure_bar",
    "radial_distribution_function",
    "rattle_velocities",
    "receptor_strain_score",
    "sample_harmonic_positions",
    "score_flexible_pose",
    "score_pose",
    "set_backend",
    "set_precision",
    "shake_positions",
    "simulate_velocity_verlet",
    "simulate_velocity_verlet_system",
    "steepest_descent_step",
    "thermodynamic_integration_harmonic",
    "to_numpy",
    "to_torch",
    "torch_available",
    "torch_device",
    "torch_dtype",
    "total_energy",
    "use_backend",
    "velocity_verlet_step",
    "velocity_verlet_step_system",
    "virial_pressure",
]
