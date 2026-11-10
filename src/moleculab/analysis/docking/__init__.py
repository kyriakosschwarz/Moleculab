"""Pose records, RMSD, clustering, and binding-mode comparison.

This package analyses poses produced elsewhere. It does not contain a docking
search engine or a docking-program file parser. Scores are metadata with an
explicit unit and source; they are not reinterpreted as binding free energies.
"""

from .benchmark import (
    BenchmarkCase,
    BenchmarkSummary,
    DockingBenchmark,
    summarize_benchmark,
)
from .clustering import (
    PoseCluster,
    PoseClusteringResult,
    cluster_poses,
    pose_distance,
)
from .comparison import (
    BindingModeComparison,
    InteractionPatternComparison,
    PoseGeometryComparison,
    compare_binding_modes,
    compare_interaction_patterns,
    compare_pose_geometry,
)
from .poses import LigandPose
from .rmsd import (
    SymmetryAwareRMSD,
    atom_id_mapping,
    pose_rmsd,
    symmetry_aware_rmsd,
    validate_atom_mapping,
)
from .validation import PoseRecovery, evaluate_pose_recovery

__all__ = [
    "BenchmarkCase",
    "BenchmarkSummary",
    "BindingModeComparison",
    "DockingBenchmark",
    "InteractionPatternComparison",
    "LigandPose",
    "PoseCluster",
    "PoseClusteringResult",
    "PoseGeometryComparison",
    "PoseRecovery",
    "SymmetryAwareRMSD",
    "atom_id_mapping",
    "cluster_poses",
    "compare_binding_modes",
    "compare_interaction_patterns",
    "compare_pose_geometry",
    "evaluate_pose_recovery",
    "pose_distance",
    "pose_rmsd",
    "summarize_benchmark",
    "symmetry_aware_rmsd",
    "validate_atom_mapping",
]
