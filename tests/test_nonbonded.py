# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for typed nonbonded interactions and reaction-field electrostatics."""

import numpy as np
import pytest

from moleculab import (
    Bond,
    Box,
    HarmonicBondPotential,
    LennardJonesParameters,
    LennardJonesPotential,
    NonbondedModel,
    Particle,
    ReactionFieldElectrostatics,
    System,
    Topology,
    combine_lorentz_berthelot,
    units,
    velocity_verlet_step_system,
)


def numerical_force(energy_function, displacement, h=1e-6):
    force = np.empty(3)
    for axis in range(3):
        plus = np.array(displacement, dtype=float)
        minus = np.array(displacement, dtype=float)
        plus[axis] += h
        minus[axis] -= h
        force[axis] = -(energy_function(plus) - energy_function(minus)) / (2.0 * h)
    return force


def test_lennard_jones_parameters_allow_zero_epsilon():
    parameters = LennardJonesParameters(epsilon=0.0, sigma=0.1)

    assert parameters.energy(0.2) == 0.0
    np.testing.assert_allclose(parameters.force([0.2, 0.0, 0.0]), np.zeros(3))


def test_lorentz_berthelot_combining_rules():
    oxygen = LennardJonesParameters(epsilon=0.6502, sigma=0.3166)
    hydrogen = LennardJonesParameters(epsilon=0.0, sigma=0.1)

    oo = combine_lorentz_berthelot(oxygen, oxygen)
    oh = combine_lorentz_berthelot(oxygen, hydrogen)

    assert oo.epsilon == pytest.approx(0.6502)
    assert oo.sigma == pytest.approx(0.3166)
    assert oh.epsilon == pytest.approx(0.0)
    assert oh.sigma == pytest.approx((0.3166 + 0.1) / 2.0)


def test_reaction_field_with_unit_dielectric_is_shifted_coulomb():
    electrostatics = ReactionFieldElectrostatics(dielectric=1.0)
    cutoff = 1.0

    # U = k q_i q_j (1/r - 1/Rc), so at r=0.5 and Rc=1 the bracket is 1.
    assert electrostatics.energy(1.0, 1.0, 0.5, cutoff) == pytest.approx(
        units.COULOMB_CONSTANT
    )
    np.testing.assert_allclose(
        electrostatics.force(1.0, 1.0, [0.5, 0.0, 0.0], cutoff),
        [4.0 * units.COULOMB_CONSTANT, 0.0, 0.0],
    )


def test_reaction_field_energy_is_zero_at_cutoff():
    electrostatics = ReactionFieldElectrostatics(dielectric=78.3)

    assert electrostatics.energy(1.0, -1.0, 1.0, 1.0) == pytest.approx(0.0)
    np.testing.assert_allclose(
        electrostatics.force(1.0, -1.0, [1.0, 0.0, 0.0], 1.0),
        np.zeros(3),
    )


def test_reaction_field_force_matches_energy_gradient():
    electrostatics = ReactionFieldElectrostatics(dielectric=78.3)
    displacement = np.array([0.4, 0.1, -0.05])
    cutoff = 1.0

    def energy_at(vector):
        distance = float(np.linalg.norm(vector))
        return electrostatics.energy(1.0, -0.5, distance, cutoff)

    np.testing.assert_allclose(
        electrostatics.force(1.0, -0.5, displacement, cutoff),
        numerical_force(energy_at, displacement),
        atol=1e-6,
    )


def test_nonbonded_model_reproduces_typed_lennard_jones_energy():
    parameters = LennardJonesParameters(epsilon=0.6502, sigma=0.3166)
    model = NonbondedModel(parameters={"O": parameters}, cutoff=1.0)
    displacement = [0.4, 0.0, 0.0]

    expected = LennardJonesPotential(epsilon=0.6502, sigma=0.3166).energy(
        displacement
    )

    assert model.energy("O", "O", 0.0, 0.0, displacement) == pytest.approx(expected)


def test_nonbonded_parameters_are_read_only():
    model = NonbondedModel(
        parameters={"C": LennardJonesParameters(epsilon=0.5, sigma=0.34)},
        cutoff=1.0,
    )

    with pytest.raises(TypeError):
        model.parameters["C"] = LennardJonesParameters(epsilon=1.0, sigma=0.4)


def test_nonbonded_model_combines_lennard_jones_and_electrostatics():
    model = NonbondedModel(
        parameters={
            "O": LennardJonesParameters(epsilon=0.65, sigma=0.3166),
            "H": LennardJonesParameters(epsilon=0.0, sigma=0.1),
        },
        cutoff=1.0,
        electrostatics=ReactionFieldElectrostatics(dielectric=78.3),
    )
    displacement = np.array([0.5, 0.0, 0.0])

    expected_electrostatic = model.electrostatics.energy(
        -0.8476, 0.4238, 0.5, 1.0
    )
    expected_lj = model.combined_parameters("O", "H").energy(0.5)

    assert model.energy(
        "O", "H", -0.8476, 0.4238, displacement
    ) == pytest.approx(expected_lj + expected_electrostatic)


def test_nonbonded_model_cuts_off_beyond_cutoff():
    model = NonbondedModel(
        parameters={"O": LennardJonesParameters(epsilon=1.0, sigma=0.3)},
        cutoff=0.8,
        electrostatics=ReactionFieldElectrostatics(dielectric=78.3),
    )

    assert model.energy("O", "O", 1.0, -1.0, [0.9, 0.0, 0.0]) == 0.0
    np.testing.assert_allclose(
        model.force("O", "O", 1.0, -1.0, [0.9, 0.0, 0.0]),
        np.zeros(3),
    )


def test_topology_derives_12_and_13_exclusions():
    topology = Topology(
        bonds=[Bond(0, 1, HarmonicBondPotential(1000.0, 0.1))],
        angles=[],
    )
    assert topology.nonbonded_exclusions() == frozenset({(0, 1)})


def test_bonded_pair_is_excluded_from_typed_nonbonded_model():
    model = NonbondedModel(
        parameters={
            "O": LennardJonesParameters(epsilon=0.65, sigma=0.3166),
            "H": LennardJonesParameters(epsilon=0.0, sigma=0.1),
        },
        cutoff=1.0,
        electrostatics=ReactionFieldElectrostatics(dielectric=78.3),
    )
    particles = [
        Particle(15.999, [0.0, 0.0, 0.0], atom_type="O", charge=-0.8476),
        Particle(1.008, [0.1, 0.0, 0.0], atom_type="H", charge=0.4238),
    ]
    topology = Topology(
        bonds=[Bond(0, 1, HarmonicBondPotential(1000.0, 0.1))]
    )

    connected = System(
        particles=particles,
        topology=topology,
        nonbonded_model=model,
    )
    unconnected = System(
        particles=particles,
        nonbonded_model=model,
    )

    assert connected.nonbonded_energy() == pytest.approx(0.0)
    assert abs(unconnected.nonbonded_energy()) > 0.0


def test_system_validates_atom_types_and_cutoff():
    model = NonbondedModel(
        parameters={"O": LennardJonesParameters(epsilon=0.65, sigma=0.3166)},
        cutoff=0.6,
    )

    with pytest.raises(ValueError, match="Missing"):
        System(
            particles=[Particle(1.0, [0.0, 0.0, 0.0], atom_type="X")],
            nonbonded_model=model,
        )
    with pytest.raises(ValueError, match="cutoff"):
        System(
            particles=[Particle(1.0, [0.0, 0.0, 0.0], atom_type="O")],
            box=Box([1.0, 1.0, 1.0]),
            nonbonded_model=model,
        )


def test_integrator_preserves_atom_types_and_charges():
    model = NonbondedModel(
        parameters={"A": LennardJonesParameters(epsilon=1.0, sigma=0.3)},
        cutoff=1.0,
    )
    system = System(
        particles=[
            Particle(
                1.0,
                [0.0, 0.0, 0.0],
                [0.1, 0.0, 0.0],
                atom_type="A",
                charge=0.5,
            )
        ],
        nonbonded_model=model,
    )

    stepped = velocity_verlet_step_system(system, dt=0.001)

    assert stepped.particles[0].atom_type == "A"
    assert stepped.particles[0].charge == pytest.approx(0.5)
    assert stepped.nonbonded_model is model


def test_nonbonded_parameter_validation():
    with pytest.raises(ValueError, match="epsilon"):
        LennardJonesParameters(epsilon=-1.0, sigma=0.3)
    with pytest.raises(ValueError, match="sigma"):
        LennardJonesParameters(epsilon=1.0, sigma=0.0)
    with pytest.raises(ValueError, match="dielectric"):
        ReactionFieldElectrostatics(dielectric=0.0)
    with pytest.raises(ValueError, match="cutoff"):
        NonbondedModel(
            parameters={"A": LennardJonesParameters(1.0, 0.3)},
            cutoff=0.0,
        )
    with pytest.raises(ValueError, match="Missing"):
        NonbondedModel(
            parameters={"A": LennardJonesParameters(1.0, 0.3)},
            cutoff=1.0,
        ).validate_atom_types(["A", "B"])
