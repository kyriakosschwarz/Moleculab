"""Moleculab: an educational molecular physics and simulation library.

The package is intentionally small and explicit.  Physical models (potentials)
are kept separate from general physics functions and from numerical
algorithms (integrators).

The source is organised by architecture layer:

``foundation``
    units, validation, particles, and box geometry;
``models``
    potentials, bonded terms, topology, nonbonded interactions, Ewald, physics;
``composition``
    :class:`System`, neighbour screening, and constraints;
``algorithms``
    Velocity Verlet and minimisation;
``protocols``
    thermostats and pressure coupling;
``analysis``
    trajectory measurements, contact persistence, binding-site geometry,
    interaction records, fingerprints, external pose analysis, and SAR comparison;
``applications``
    docking, solvent, water, free energy, acid-base, and QM/MM.

The public names below are re-exported so tutorials and users can write
``from moleculab import ...`` without depending on the internal layout.

Examples
--------
>>> from moleculab import HarmonicPotential, Particle, total_energy
>>> particle = Particle(mass=12.0, position=[0.1, 0.0, 0.0], velocity=[0.1, 0.0, 0.0])
>>> potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
>>> round(total_energy(particle, potential), 6)
0.56
"""

from .algorithms.integrators import (
    SystemTrajectory,
    Trajectory,
    simulate_velocity_verlet,
    simulate_velocity_verlet_system,
    velocity_verlet_step,
    velocity_verlet_step_system,
)
from .algorithms.minimization import (
    MinimizationResult,
    max_force_magnitude,
    minimize_energy,
    steepest_descent_step,
)
from .analysis.binding_site import (
    BindingSiteGeometry,
    binding_site_geometry,
    group_indices_by_residue,
    minimum_distances,
    pocket_geometry,
    select_pocket_atoms,
)
from .analysis.docking import (
    BenchmarkCase,
    BenchmarkSummary,
    BindingModeComparison,
    DockingBenchmark,
    InteractionPatternComparison,
    LigandPose,
    PoseCluster,
    PoseClusteringResult,
    PoseGeometryComparison,
    PoseRecovery,
    SymmetryAwareRMSD,
    atom_id_mapping,
    cluster_poses,
    compare_binding_modes,
    compare_interaction_patterns,
    compare_pose_geometry,
    evaluate_pose_recovery,
    pose_distance,
    pose_rmsd,
    summarize_benchmark,
    symmetry_aware_rmsd,
    validate_atom_mapping,
)
from .analysis.fingerprints import (
    FingerprintSchema,
    InteractionFeature,
    InteractionFingerprint,
    OccupancyResult,
    build_fingerprint,
    calculate_occupancy,
)
from .analysis.interactions import (
    HydrogenBond,
    HydrogenBondCriteria,
    InteractionAtom,
    detect_hydrogen_bonds,
)
from .analysis.sar import (
    AssayMeasurement,
    CompoundRecord,
    MatchedMolecularPair,
    MatchedPairComparison,
    compare_matched_pair,
    ligand_efficiency,
    pIC50,
)
from .analysis.trajectory import (
    MSDResult,
    RDFResult,
    mean_squared_displacement,
    radial_distribution_function,
)
from .analysis.trajectory_contacts import (
    ContactLifetimeResult,
    ContactOccupancyResult,
    ResidueContactMatrix,
    aggregate_contacts_by_residue,
    calculate_contact_lifetimes,
    calculate_contact_matrix,
    calculate_contact_occupancy,
)
from .applications.acid_base import (
    boltzmann_populations,
    deprotonated_fraction,
    mean_net_charge,
    mean_protonation,
    protonated_fraction,
    protonation_microstates,
    tautomer_populations,
)
from .applications.docking import (
    DockingPose,
    DockingResult,
    RigidLigand,
    dock_rigid_ligand,
    score_pose,
)
from .applications.fep import (
    FEPResult,
    fep_between_potentials,
    fep_harmonic_oscillator,
    free_energy_perturbation,
    harmonic_free_energy_difference,
    sample_harmonic_positions,
)
from .applications.flexible_docking import (
    FlexibleDockingPose,
    FlexibleDockingResult,
    SideChain,
    dock_flexible_receptor,
    make_receptor_conformer,
    receptor_strain_score,
    score_flexible_pose,
)
from .applications.lj_alchemy import MultistageFEPResult, SoftCoreLennardJones
from .applications.mbar import MBARResult, mbar, mbar_harmonic
from .applications.mc_alchemy import BondedLennardJonesMC, SoftCoreLennardJonesPair
from .applications.qm_mm import MMPointCharge, ParticleInBox, QMResult
from .applications.solvated_alchemy import SolvatedLigandAlchemy
from .applications.solvent import (
    ARGON,
    COARSE_GRAINED_WATER,
    SolventModel,
    build_solvent_box,
)
from .applications.ti import (
    TIResult,
    integrate_thermodynamic_derivative,
    thermodynamic_integration_harmonic,
)
from .applications.water import SPCE_WATER, WaterModel, build_water_box
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
from .composition.constraints import (
    DistanceConstraint,
    RattleResult,
    ShakeResult,
    assert_constraints_satisfied,
    constraint_error,
    rattle_velocities,
    shake_positions,
)
from .composition.neighbour import NeighbourList, build_neighbour_list
from .composition.system import System
from .foundation import units as units
from .foundation.box import Box
from .foundation.particles import Particle
from .foundation.units import ELECTRON_MASS_DA
from .models.bonded import (
    AnglePotential,
    BondPotential,
    HarmonicAnglePotential,
    HarmonicBondPotential,
    PeriodicTorsionPotential,
    TorsionPotential,
    dihedral_angle,
)
from .models.ewald import EwaldElectrostatics
from .models.nonbonded import (
    LennardJonesParameters,
    NonbondedModel,
    ReactionFieldElectrostatics,
    combine_lorentz_berthelot,
)
from .models.physics import (
    acceleration,
    force_from_potential,
    kinetic_energy,
    potential_energy,
    total_energy,
)
from .models.potentials import (
    CoulombPotential,
    DoubleWellPotential,
    HarmonicPotential,
    LennardJonesPotential,
    PairPotential,
    Potential,
)
from .models.topology import Angle, Bond, Topology, Torsion
from .protocols.pressure import (
    berendsen_barostat,
    kinetic_pressure,
    pressure,
    pressure_bar,
    virial_pressure,
)
from .protocols.thermostats import (
    berendsen_thermostat,
    degrees_of_freedom,
    kinetic_temperature,
)

__all__ = [
    "ARGON",
    "COARSE_GRAINED_WATER",
    "CompoundRecord",
    "ContactLifetimeResult",
    "ContactOccupancyResult",
    "ELECTRON_MASS_DA",
    "SPCE_WATER",
    "Angle",
    "AnglePotential",
    "AssayMeasurement",
    "BenchmarkCase",
    "BenchmarkSummary",
    "BindingModeComparison",
    "BindingSiteGeometry",
    "Bond",
    "BondPotential",
    "BondedLennardJonesMC",
    "Box",
    "CoulombPotential",
    "DistanceConstraint",
    "DockingBenchmark",
    "DockingPose",
    "DockingResult",
    "DoubleWellPotential",
    "EwaldElectrostatics",
    "FEPResult",
    "FingerprintSchema",
    "FlexibleDockingPose",
    "FlexibleDockingResult",
    "HarmonicAnglePotential",
    "HarmonicBondPotential",
    "HarmonicPotential",
    "HydrogenBond",
    "HydrogenBondCriteria",
    "InteractionAtom",
    "InteractionFeature",
    "InteractionFingerprint",
    "InteractionPatternComparison",
    "LennardJonesParameters",
    "LennardJonesPotential",
    "LigandPose",
    "MatchedMolecularPair",
    "MatchedPairComparison",
    "MBARResult",
    "MMPointCharge",
    "MSDResult",
    "MinimizationResult",
    "MultistageFEPResult",
    "NeighbourList",
    "NonbondedModel",
    "OccupancyResult",
    "PairPotential",
    "Particle",
    "ParticleInBox",
    "PeriodicTorsionPotential",
    "PoseCluster",
    "PoseClusteringResult",
    "PoseGeometryComparison",
    "PoseRecovery",
    "Potential",
    "QMResult",
    "RDFResult",
    "RattleResult",
    "ReactionFieldElectrostatics",
    "ResidueContactMatrix",
    "RigidLigand",
    "ShakeResult",
    "SideChain",
    "SoftCoreLennardJones",
    "SoftCoreLennardJonesPair",
    "SolvatedLigandAlchemy",
    "SolventModel",
    "SymmetryAwareRMSD",
    "System",
    "SystemTrajectory",
    "TIResult",
    "Topology",
    "Torsion",
    "TorsionPotential",
    "Trajectory",
    "WaterModel",
    "acceleration",
    "aggregate_contacts_by_residue",
    "assert_constraints_satisfied",
    "atom_id_mapping",
    "available_backends",
    "backend_info",
    "berendsen_barostat",
    "berendsen_thermostat",
    "binding_site_geometry",
    "boltzmann_populations",
    "build_fingerprint",
    "build_neighbour_list",
    "build_solvent_box",
    "build_water_box",
    "calculate_contact_lifetimes",
    "calculate_contact_matrix",
    "calculate_contact_occupancy",
    "calculate_occupancy",
    "cluster_poses",
    "combine_lorentz_berthelot",
    "compare_binding_modes",
    "compare_interaction_patterns",
    "compare_matched_pair",
    "compare_pose_geometry",
    "constraint_error",
    "degrees_of_freedom",
    "deprotonated_fraction",
    "detect_hydrogen_bonds",
    "dihedral_angle",
    "dock_flexible_receptor",
    "dock_rigid_ligand",
    "evaluate_pose_recovery",
    "fep_between_potentials",
    "fep_harmonic_oscillator",
    "force_from_potential",
    "free_energy_perturbation",
    "get_backend",
    "get_precision",
    "group_indices_by_residue",
    "harmonic_free_energy_difference",
    "integrate_thermodynamic_derivative",
    "kinetic_energy",
    "kinetic_pressure",
    "kinetic_temperature",
    "ligand_efficiency",
    "make_receptor_conformer",
    "max_force_magnitude",
    "mbar",
    "mbar_harmonic",
    "mean_net_charge",
    "mean_protonation",
    "mean_squared_displacement",
    "minimize_energy",
    "minimum_distances",
    "pIC50",
    "pocket_geometry",
    "pose_distance",
    "pose_rmsd",
    "potential_energy",
    "pressure",
    "pressure_bar",
    "protonated_fraction",
    "protonation_microstates",
    "radial_distribution_function",
    "rattle_velocities",
    "receptor_strain_score",
    "sample_harmonic_positions",
    "score_flexible_pose",
    "score_pose",
    "select_pocket_atoms",
    "set_backend",
    "set_precision",
    "shake_positions",
    "simulate_velocity_verlet",
    "simulate_velocity_verlet_system",
    "steepest_descent_step",
    "summarize_benchmark",
    "symmetry_aware_rmsd",
    "tautomer_populations",
    "thermodynamic_integration_harmonic",
    "to_numpy",
    "to_torch",
    "torch_available",
    "torch_device",
    "torch_dtype",
    "total_energy",
    "use_backend",
    "validate_atom_mapping",
    "velocity_verlet_step",
    "velocity_verlet_step_system",
    "virial_pressure",
]
