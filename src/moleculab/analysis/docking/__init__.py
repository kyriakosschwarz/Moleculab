"""Pose records and docking pose analysis.

This package analyses poses produced elsewhere. It does not contain a docking
search engine or a docking-program file parser. Scores are metadata with an
explicit unit and source; they are not reinterpreted as binding free energies.
"""

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
    "LigandPose",
    "PoseRecovery",
    "SymmetryAwareRMSD",
    "atom_id_mapping",
    "evaluate_pose_recovery",
    "pose_rmsd",
    "symmetry_aware_rmsd",
    "validate_atom_mapping",
]
