"""Pose clustering by explicit geometric distance.

The default distance is symmetry-unaware RMSD based on stable atom identifiers.
For symmetric ligands, pass a distance function that uses explicit equivalent
atom mappings. Clustering is single-linkage connected components, which is
deterministic but can chain through intermediate poses.
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass

import numpy as np

from .poses import LigandPose
from .rmsd import AtomMapping, atom_id_mapping, symmetry_aware_rmsd

PoseDistance = Callable[[LigandPose, LigandPose], float]


@dataclass(frozen=True)
class PoseCluster:
    """A connected component of poses within a distance threshold."""

    member_indices: tuple[int, ...]
    representative_index: int

    @property
    def size(self) -> int:
        """Number of poses in the cluster."""

        return len(self.member_indices)


@dataclass(frozen=True)
class PoseClusteringResult:
    """Single-linkage clustering result and its pairwise distance matrix."""

    threshold: float
    clusters: tuple[PoseCluster, ...]
    distance_matrix: np.ndarray
    n_poses: int

    @property
    def n_clusters(self) -> int:
        """Number of clusters."""

        return len(self.clusters)

    @property
    def labels(self) -> np.ndarray:
        """Cluster label for each pose index."""

        labels = np.empty(self.n_poses, dtype=int)
        for label, cluster in enumerate(self.clusters):
            labels[list(cluster.member_indices)] = label
        return labels


def pose_distance(
    pose_a: LigandPose,
    pose_b: LigandPose,
    *,
    atom_mappings: Sequence[Sequence[tuple[int, int]]] | None = None,
    require_complete: bool = True,
) -> float:
    """Return RMSD between two poses.

    If ``atom_mappings`` is omitted, atoms are mapped by identical stable atom
    identifiers. If it is supplied, the minimum RMSD over those explicit
    mappings is returned.
    """

    if not isinstance(pose_a, LigandPose) or not isinstance(pose_b, LigandPose):
        raise TypeError("pose_a and pose_b must be LigandPose objects.")

    if atom_mappings is None:
        mapping: AtomMapping = atom_id_mapping(pose_a, pose_b)
        mappings = (mapping,)
    else:
        mappings = tuple(atom_mappings)
        if not mappings:
            raise ValueError("atom_mappings must contain at least one mapping.")

    return symmetry_aware_rmsd(
        pose_a.positions,
        pose_b.positions,
        mappings,
        require_complete=require_complete,
    ).rmsd


def cluster_poses(
    poses: Sequence[LigandPose],
    threshold: float,
    *,
    distance_function: PoseDistance | None = None,
) -> PoseClusteringResult:
    """Cluster poses using single-linkage connected components.

    ``threshold`` is in nm. The distance matrix is symmetric; each off-diagonal
    pair is evaluated once with the earlier pose as the first argument. The
    default distance maps atoms by stable identifier and is not symmetry-aware.
    Pass a custom ``distance_function`` for explicit symmetry mappings.
    """

    values = tuple(poses)
    if not values:
        raise ValueError("poses must contain at least one LigandPose.")
    for pose in values:
        if not isinstance(pose, LigandPose):
            raise TypeError("poses must contain LigandPose objects.")

    limit = float(threshold)
    if not np.isfinite(limit) or limit <= 0.0:
        raise ValueError("threshold must be finite and positive in nm.")

    if distance_function is None:
        distance_function = pose_distance
    if not callable(distance_function):
        raise TypeError("distance_function must be callable.")

    n_poses = len(values)
    distances = np.zeros((n_poses, n_poses), dtype=float)
    for first in range(n_poses):
        for second in range(first + 1, n_poses):
            value = float(distance_function(values[first], values[second]))
            if not np.isfinite(value) or value < 0.0:
                raise ValueError("distance_function must return finite non-negative values.")
            distances[first, second] = value
            distances[second, first] = value

    unassigned = set(range(n_poses))
    clusters: list[PoseCluster] = []
    while unassigned:
        start = min(unassigned)
        unassigned.remove(start)
        stack = [start]
        members = []
        while stack:
            current = stack.pop()
            members.append(current)
            neighbours = [
                candidate
                for candidate in sorted(unassigned)
                if distances[current, candidate] <= limit
            ]
            for neighbour in neighbours:
                unassigned.remove(neighbour)
                stack.append(neighbour)
        members.sort()
        clusters.append(
            PoseCluster(
                member_indices=tuple(members),
                representative_index=members[0],
            )
        )

    clusters.sort(key=lambda cluster: cluster.representative_index)
    return PoseClusteringResult(
        threshold=limit,
        clusters=tuple(clusters),
        distance_matrix=distances,
        n_poses=n_poses,
    )
