"""Trajectory analysis for molecular simulations.

Two of the most useful first analyses of a trajectory are:

Mean-squared displacement (MSD)
    How far particles move after a lag time :math:`\\tau`.  For diffusive
    motion, :math:`\\mathrm{MSD}(\\tau) \\sim 6D\\tau` in three dimensions.

Radial distribution function (RDF)
    How the local density around a particle compares with the average density.
    For an ideal gas, :math:`g(r) = 1`; peaks reveal shells of neighbours in a
    liquid or crystal.

Both functions accept positions with shape ``(n_frames, n_particles, 3)``.
A single frame with shape ``(n_particles, 3)`` is also accepted.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..foundation.box import Box


@dataclass(frozen=True)
class MSDResult:
    """Result of :func:`mean_squared_displacement`.

    ``lag`` contains the elapsed times (or frame indices when no times are
    supplied), and ``msd`` contains the corresponding mean-squared
    displacement in nm^2.
    """

    lag: np.ndarray
    msd: np.ndarray


@dataclass(frozen=True)
class RDFResult:
    """Result of :func:`radial_distribution_function`.

    ``r`` is the midpoint of each radial bin in nm, ``g`` is the radial
    distribution function, ``counts`` is the number of pair distances in each
    bin across all frames, and ``edges`` are the bin boundaries.
    """

    r: np.ndarray
    g: np.ndarray
    counts: np.ndarray
    edges: np.ndarray


def _as_trajectory_positions(positions: np.ndarray) -> np.ndarray:
    coordinates = np.asarray(positions, dtype=float)
    if coordinates.ndim == 2:
        coordinates = coordinates[np.newaxis, ...]
    if coordinates.ndim != 3 or coordinates.shape[2] != 3:
        raise ValueError(
            "positions must have shape (n_frames, n_particles, 3) or "
            f"(n_particles, 3), got {coordinates.shape}."
        )
    if coordinates.shape[0] == 0 or coordinates.shape[1] == 0:
        raise ValueError("positions must contain at least one frame and one particle.")
    if not np.all(np.isfinite(coordinates)):
        raise ValueError("positions must contain only finite values.")
    return coordinates


def _as_times(times: np.ndarray | None, n_frames: int) -> np.ndarray:
    if times is None:
        return np.arange(n_frames, dtype=float)

    values = np.asarray(times, dtype=float)
    if values.shape != (n_frames,):
        raise ValueError(
            f"times must have shape ({n_frames},), got {values.shape}."
        )
    if not np.all(np.isfinite(values)):
        raise ValueError("times must contain only finite values.")
    if n_frames > 1 and np.any(np.diff(values) <= 0.0):
        raise ValueError("times must be strictly increasing.")
    return values - values[0]


def mean_squared_displacement(
    positions: np.ndarray,
    box: Box | None = None,
    times: np.ndarray | None = None,
) -> MSDResult:
    """Compute the mean-squared displacement of a trajectory.

    The MSD at lag index ``k`` is averaged over every time origin and every
    particle:

    .. math::

        \\mathrm{MSD}(k) = \\left\\langle
        |\\mathbf r_i(t+k) - \\mathbf r_i(t)|^2
        \\right\\rangle_{i,t}.

    When ``box`` is supplied, displacements use the minimum-image convention.
    This is appropriate for wrapped coordinates and short lag times.  For long
    diffusive trajectories, unwrap the coordinates first or MSD will saturate
    at the size of the box.

    Parameters
    ----------
    positions:
        Array with shape ``(n_frames, n_particles, 3)`` in nm.
    box:
        Optional periodic box.
    times:
        Optional frame times.  When omitted, the lag is the frame index.

    Returns
    -------
    MSDResult
        Lag values and MSD values (nm^2).
    """

    coordinates = _as_trajectory_positions(positions)
    if box is not None and not isinstance(box, Box):
        raise TypeError("box must be a Box or None.")
    lag = _as_times(times, coordinates.shape[0])

    n_frames = coordinates.shape[0]
    msd = np.empty(n_frames, dtype=float)
    for offset in range(n_frames):
        displacement = coordinates[offset:] - coordinates[: n_frames - offset]
        if box is not None:
            displacement = box.minimum_image_array(displacement)
        msd[offset] = float(np.mean(np.sum(displacement * displacement, axis=-1)))

    return MSDResult(lag=lag, msd=msd)


def radial_distribution_function(
    positions: np.ndarray,
    box: Box,
    r_max: float | None = None,
    n_bins: int = 100,
) -> RDFResult:
    """Compute the radial distribution function :math:`g(r)`.

    For each frame, every unique particle pair contributes one distance.  The
    histogram is normalized by the number of frames, the number of pairs, and
    the ideal-gas number of pairs expected in each spherical shell:

    .. math::

        g(r_k) = \\frac{N_k}
        {N_{\\mathrm{frames}} N_{\\mathrm{pairs}}
         (4\\pi r_k^2 \\Delta r / V)}.

    Parameters
    ----------
    positions:
        Array with shape ``(n_frames, n_particles, 3)`` in nm.
    box:
        Periodic box used for the minimum-image convention and the density
        normalization.
    r_max:
        Largest distance to include in nm.  Defaults to half the shortest box
        edge.  Larger values are rejected because pair distances beyond that
        range are affected by periodic images.
    n_bins:
        Number of radial bins.

    Returns
    -------
    RDFResult
        Bin midpoints, :math:`g(r)`, raw pair counts, and bin edges.
    """

    coordinates = _as_trajectory_positions(positions)
    if not isinstance(box, Box):
        raise TypeError("box must be a Box; RDF normalization requires a volume.")
    if isinstance(n_bins, bool) or not isinstance(n_bins, (int, np.integer)):
        raise ValueError("n_bins must be an integer.")
    if int(n_bins) < 1:
        raise ValueError("n_bins must be at least 1.")
    n_bins = int(n_bins)

    n_frames, n_particles, _ = coordinates.shape
    if n_particles < 2:
        raise ValueError("RDF requires at least two particles.")

    shortest_half_edge = float(np.min(box.lengths) / 2.0)
    if r_max is None:
        r_max = shortest_half_edge
    r_max = float(r_max)
    if not np.isfinite(r_max) or r_max <= 0.0:
        raise ValueError("r_max must be finite and positive.")
    if r_max > shortest_half_edge + 1e-12:
        raise ValueError(
            "r_max must not exceed half the shortest box edge "
            f"({shortest_half_edge:.6g} nm)."
        )

    edges = np.linspace(0.0, r_max, n_bins + 1)
    counts = np.zeros(n_bins, dtype=float)
    upper_triangle = np.triu_indices(n_particles, k=1)

    for frame in coordinates:
        displacement = frame[:, np.newaxis, :] - frame[np.newaxis, :, :]
        displacement = box.minimum_image_array(displacement)
        distances = np.linalg.norm(displacement, axis=-1)
        pair_distances = distances[upper_triangle]
        histogram, _ = np.histogram(pair_distances, bins=edges)
        counts += histogram

    shell_volumes = (4.0 * np.pi / 3.0) * (edges[1:] ** 3 - edges[:-1] ** 3)
    n_pairs = n_particles * (n_particles - 1) / 2.0
    expected_counts = (
        n_frames * n_pairs * shell_volumes / box.volume
    )
    g = counts / expected_counts

    return RDFResult(
        r=0.5 * (edges[1:] + edges[:-1]),
        g=g,
        counts=counts,
        edges=edges,
    )
