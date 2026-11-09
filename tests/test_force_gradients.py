"""Force-gradient audit for every core potential.

Each potential is numerically differentiated and compared with its analytical
force.  This is the recurring invariant behind the tutorials:

    F = -grad U.

The numerical derivative calls only the energy function.  It never calls the
analytical force, so a sign error or missing factor cannot be hidden by reusing
the same implementation twice.

The configurations are deliberately asymmetric and away from equilibrium unless
an equilibrium case is explicitly included.  The step h = 1e-6 nm is small
compared with the coordinate scale but large enough to avoid catastrophic
cancellation for double-precision energies.  The tolerance atol = 1e-6
absorbs both truncation error and roundoff.
"""

import numpy as np
import pytest

from moleculab import (
    CoulombPotential,
    HarmonicAnglePotential,
    HarmonicBondPotential,
    HarmonicPotential,
    LennardJonesParameters,
    LennardJonesPotential,
    NonbondedModel,
    PeriodicTorsionPotential,
    ReactionFieldElectrostatics,
)

H = 1e-6
ATOL = 1e-6


def numerical_force(energy, coordinates, h=H):
    """Return -grad U at *coordinates* by central differences."""

    coordinates = np.asarray(coordinates, dtype=float)
    gradient = np.zeros_like(coordinates)
    for index in np.ndindex(coordinates.shape):
        plus = coordinates.copy()
        minus = coordinates.copy()
        plus[index] += h
        minus[index] -= h
        gradient[index] = (energy(plus) - energy(minus)) / (2.0 * h)
    return -gradient


def assert_force_matches_gradient(analytical, energy, coordinates):
    np.testing.assert_allclose(
        analytical,
        numerical_force(energy, coordinates),
        atol=ATOL,
    )


@pytest.mark.parametrize(
    "position",
    [
        [0.12, -0.03, 0.04],
        [-0.08, 0.11, -0.06],
        [0.16, 0.02, 0.09],
    ],
)
def test_harmonic_potential_force_matches_gradient(position):
    potential = HarmonicPotential(stiffness=100.0, center=[0.05, -0.02, 0.01])
    position = np.asarray(position, dtype=float)

    assert_force_matches_gradient(
        potential.force(position),
        potential.energy,
        position,
    )


@pytest.mark.parametrize(
    "displacement",
    [
        [0.45, -0.12, 0.08],
        [0.36, 0.01, -0.02],
        [0.30, 0.00, 0.00],
        [0.80, -0.30, 0.20],
    ],
)
def test_lennard_jones_force_matches_gradient(displacement):
    potential = LennardJonesPotential(epsilon=0.997, sigma=0.34)
    displacement = np.asarray(displacement, dtype=float)

    assert_force_matches_gradient(
        potential.force(displacement),
        potential.energy,
        displacement,
    )


@pytest.mark.parametrize(
    "displacement",
    [
        [0.45, -0.12, 0.08],
        [0.25, 0.10, -0.03],
        [0.70, 0.20, 0.10],
    ],
)
def test_typed_lennard_jones_parameters_force_matches_gradient(displacement):
    parameters = LennardJonesParameters(epsilon=0.997, sigma=0.34)
    displacement = np.asarray(displacement, dtype=float)

    def energy(vector):
        return parameters.energy(float(np.linalg.norm(vector)))

    assert_force_matches_gradient(
        parameters.force(displacement),
        energy,
        displacement,
    )


@pytest.mark.parametrize(
    "displacement",
    [
        [0.7, 0.2, -0.1],
        [0.35, -0.2, 0.1],
        [1.3, 0.4, -0.2],
    ],
)
def test_coulomb_force_matches_gradient(displacement):
    potential = CoulombPotential(charge_i=1.5, charge_j=-0.75)
    displacement = np.asarray(displacement, dtype=float)

    assert_force_matches_gradient(
        potential.force(displacement),
        potential.energy,
        displacement,
    )


@pytest.mark.parametrize(
    "displacement",
    [
        [0.18, 0.02, -0.01],
        [0.12, -0.03, 0.02],
        [0.20, 0.10, -0.05],
    ],
)
def test_harmonic_bond_force_matches_gradient(displacement):
    potential = HarmonicBondPotential(stiffness=1000.0, length0=0.15)
    displacement = np.asarray(displacement, dtype=float)

    assert_force_matches_gradient(
        potential.force(displacement),
        potential.energy,
        displacement,
    )


@pytest.mark.parametrize(
    "coordinates",
    [
        [
            [0.10, 0.00, 0.00],
            [0.00, 0.00, 0.00],
            [0.02, 0.11, 0.03],
        ],
        [
            [0.08, -0.03, 0.02],
            [-0.01, 0.02, 0.00],
            [-0.04, 0.10, -0.03],
        ],
        [
            [0.14, 0.02, -0.01],
            [0.01, -0.01, 0.02],
            [0.03, 0.15, 0.04],
        ],
    ],
)
def test_harmonic_angle_forces_match_gradient(coordinates):
    potential = HarmonicAnglePotential(
        stiffness=100.0,
        angle0=np.deg2rad(104.5),
    )
    coordinates = np.asarray(coordinates, dtype=float)

    def energy(points):
        points = points.reshape(3, 3)
        r_i, r_j, r_k = points
        return potential.energy(r_i - r_j, r_k - r_j)

    r_i, r_j, r_k = coordinates
    analytical = np.concatenate(potential.forces(r_i - r_j, r_k - r_j))

    assert_force_matches_gradient(
        analytical,
        energy,
        coordinates.reshape(-1),
    )


@pytest.mark.parametrize(
    "coordinates",
    [
        [
            [0.00, 0.00, 0.00],
            [0.15, 0.00, 0.00],
            [0.15, 0.16, 0.00],
            [0.14, 0.17, 0.15],
        ],
        [
            [0.00, 0.00, 0.00],
            [0.15, 0.02, 0.01],
            [0.13, 0.16, 0.03],
            [0.18, 0.18, 0.15],
        ],
        [
            [0.01, -0.02, 0.03],
            [0.16, 0.00, 0.01],
            [0.14, 0.17, -0.02],
            [0.12, 0.19, 0.16],
        ],
    ],
)
def test_periodic_torsion_forces_match_gradient(coordinates):
    potential = PeriodicTorsionPotential(
        barrier=2.0,
        periodicity=3,
        phase=0.4,
    )
    coordinates = np.asarray(coordinates, dtype=float)

    def energy(points):
        points = points.reshape(4, 3)
        r_i, r_j, r_k, r_l = points
        return potential.energy(r_j - r_i, r_k - r_j, r_l - r_k)

    points = coordinates.reshape(4, 3)
    r_i, r_j, r_k, r_l = points
    analytical = np.concatenate(
        potential.forces(r_j - r_i, r_k - r_j, r_l - r_k)
    )

    assert_force_matches_gradient(
        analytical,
        energy,
        coordinates.reshape(-1),
    )


@pytest.mark.parametrize(
    "displacement",
    [
        [0.4, 0.1, -0.05],
        [0.25, -0.1, 0.15],
        [0.7, 0.2, 0.1],
    ],
)
def test_reaction_field_force_matches_gradient(displacement):
    electrostatics = ReactionFieldElectrostatics(dielectric=78.3)
    displacement = np.asarray(displacement, dtype=float)
    cutoff = 1.0

    def energy(vector):
        distance = float(np.linalg.norm(vector))
        return electrostatics.energy(1.0, -0.5, distance, cutoff)

    assert_force_matches_gradient(
        electrostatics.force(1.0, -0.5, displacement, cutoff),
        energy,
        displacement,
    )


@pytest.mark.parametrize(
    "type_i,type_j,charge_i,charge_j,displacement",
    [
        ("C", "H", -0.2, 0.2, [0.4, 0.1, -0.05]),
        ("O", "H", -0.8, 0.4, [0.3, -0.2, 0.1]),
        ("O", "O", -0.8, -0.8, [0.5, 0.05, -0.1]),
    ],
)
def test_nonbonded_model_force_matches_gradient(
    type_i,
    type_j,
    charge_i,
    charge_j,
    displacement,
):
    model = NonbondedModel(
        parameters={
            "C": LennardJonesParameters(epsilon=0.5, sigma=0.34),
            "H": LennardJonesParameters(epsilon=0.1, sigma=0.25),
            "O": LennardJonesParameters(epsilon=0.65, sigma=0.30),
        },
        cutoff=1.0,
        electrostatics=ReactionFieldElectrostatics(dielectric=78.3),
    )
    displacement = np.asarray(displacement, dtype=float)

    def energy(vector):
        return model.energy(
            type_i,
            type_j,
            charge_i,
            charge_j,
            vector,
        )

    assert_force_matches_gradient(
        model.force(type_i, type_j, charge_i, charge_j, displacement),
        energy,
        displacement,
    )
