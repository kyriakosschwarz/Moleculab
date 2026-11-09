"""Tests for trajectory analysis: MSD and the radial distribution function."""

import numpy as np
import pytest

from moleculab import (
    Box,
    mean_squared_displacement,
    radial_distribution_function,
)


def test_msd_is_zero_at_zero_lag():
    positions = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[0.1, 0.0, 0.0], [1.0, 0.0, 0.0]],
            [[0.2, 0.0, 0.0], [1.0, 0.0, 0.0]],
        ]
    )

    result = mean_squared_displacement(positions)

    assert result.lag[0] == 0.0
    assert result.msd[0] == 0.0


def test_msd_matches_constant_velocity_motion():
    times = np.array([0.0, 1.0, 2.0, 3.0])
    velocity = np.array([0.1, 0.0, 0.0])
    positions = times[:, np.newaxis, np.newaxis] * velocity

    result = mean_squared_displacement(positions, times=times)

    np.testing.assert_allclose(result.lag, [0.0, 1.0, 2.0, 3.0])
    np.testing.assert_allclose(result.msd, [0.0, 0.01, 0.04, 0.09])


def test_msd_times_are_relative_to_the_first_frame():
    times = np.array([5.0, 6.0, 7.0])
    positions = np.zeros((3, 1, 3))

    result = mean_squared_displacement(positions, times=times)

    np.testing.assert_allclose(result.lag, [0.0, 1.0, 2.0])


def test_msd_uses_minimum_image_when_box_is_given():
    positions = np.array([0.0, 0.6, 0.8])[:, np.newaxis, np.newaxis] * np.array(
        [1.0, 0.0, 0.0]
    )
    box = Box([1.0, 1.0, 1.0])

    without_box = mean_squared_displacement(positions)
    with_box = mean_squared_displacement(positions, box=box)

    np.testing.assert_allclose(without_box.msd, [0.0, 0.20, 0.64])
    np.testing.assert_allclose(with_box.msd, [0.0, 0.10, 0.04])


def test_msd_accepts_a_single_frame():
    positions = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])

    result = mean_squared_displacement(positions)

    assert len(result.lag) == 1
    assert result.msd[0] == 0.0


def test_msd_validates_inputs():
    with pytest.raises(ValueError, match="positions"):
        mean_squared_displacement(np.zeros((2, 2)))
    with pytest.raises(ValueError, match="times"):
        mean_squared_displacement(np.zeros((3, 1, 3)), times=[0.0, 1.0])
    with pytest.raises(ValueError, match="strictly increasing"):
        mean_squared_displacement(np.zeros((3, 1, 3)), times=[0.0, 0.0, 1.0])
    with pytest.raises(TypeError, match="box"):
        mean_squared_displacement(np.zeros((3, 1, 3)), box="not a box")


def test_rdf_single_pair_has_expected_normalization():
    box = Box([1.0, 1.0, 1.0])
    positions = np.array([[[0.0, 0.0, 0.0], [0.2, 0.0, 0.0]]])

    result = radial_distribution_function(
        positions, box=box, r_max=0.5, n_bins=10
    )

    # The pair distance 0.2 falls in bin 4: [0.20, 0.25).
    assert result.counts[4] == 1.0
    assert np.sum(result.counts) == 1.0

    expected_shell_volume = (4.0 * np.pi / 3.0) * (0.25**3 - 0.20**3)
    assert result.g[4] == pytest.approx(1.0 / expected_shell_volume)


def test_rdf_counts_each_unique_pair_once_per_frame():
    box = Box([1.0, 1.0, 1.0])
    frame = np.array(
        [
            [0.0, 0.0, 0.0],
            [0.1, 0.0, 0.0],
            [0.0, 0.2, 0.0],
        ]
    )
    positions = np.stack([frame, frame + np.array([0.01, 0.0, 0.0])])

    result = radial_distribution_function(
        positions, box=box, r_max=0.5, n_bins=20
    )

    assert np.sum(result.counts) == 2 * 3


def test_rdf_uses_minimum_image_for_pairs_across_the_boundary():
    box = Box([1.0, 1.0, 1.0])
    positions = np.array([[[0.05, 0.0, 0.0], [0.95, 0.0, 0.0]]])

    result = radial_distribution_function(
        positions, box=box, r_max=0.5, n_bins=10
    )

    # Minimum-image distance is 0.1, which lies in bin 2: [0.10, 0.15).
    assert result.counts[2] == 1.0
    assert np.sum(result.counts) == 1.0


def test_rdf_of_a_random_ideal_gas_is_close_to_one():
    rng = np.random.default_rng(1234)
    n_frames = 200
    n_particles = 30
    box = Box([10.0, 10.0, 10.0])
    positions = rng.uniform(0.0, 10.0, size=(n_frames, n_particles, 3))

    result = radial_distribution_function(
        positions, box=box, r_max=4.0, n_bins=20
    )

    # Away from the smallest r bins, a uniform ideal gas has g(r) = 1.
    assert np.mean(result.g[2:]) == pytest.approx(1.0, abs=0.12)
    assert np.all(result.g >= 0.0)


def test_rdf_result_shapes_and_edges():
    box = Box([2.0, 2.0, 2.0])
    positions = np.array([[[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]]])

    result = radial_distribution_function(
        positions, box=box, r_max=1.0, n_bins=5
    )

    assert result.r.shape == (5,)
    assert result.g.shape == (5,)
    assert result.counts.shape == (5,)
    assert result.edges.shape == (6,)
    assert result.edges[0] == 0.0
    assert result.edges[-1] == pytest.approx(1.0)


def test_rdf_validates_box_and_geometry():
    box = Box([1.0, 1.0, 1.0])
    positions = np.array([[[0.0, 0.0, 0.0], [0.2, 0.0, 0.0]]])

    with pytest.raises(TypeError, match="box"):
        radial_distribution_function(positions, box=None)
    with pytest.raises(ValueError, match="at least two"):
        radial_distribution_function(
            np.array([[[0.0, 0.0, 0.0]]]), box=box, r_max=0.4
        )
    with pytest.raises(ValueError, match="r_max"):
        radial_distribution_function(positions, box=box, r_max=0.6)
    with pytest.raises(ValueError, match="n_bins"):
        radial_distribution_function(positions, box=box, n_bins=0)
