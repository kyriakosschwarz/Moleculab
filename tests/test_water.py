"""Tests for the three-site water model and typed exclusions."""

import numpy as np
import pytest

from moleculab import (
    SPCE_WATER,
    EwaldElectrostatics,
    WaterModel,
    berendsen_thermostat,
    build_water_box,
    kinetic_temperature,
    velocity_verlet_step_system,
)


def test_default_water_model_parameters():
    assert SPCE_WATER.name
    assert SPCE_WATER.oxygen_charge + 2.0 * SPCE_WATER.hydrogen_charge == pytest.approx(
        0.0, abs=1e-12
    )
    assert SPCE_WATER.hydrogen_epsilon == 0.0


def test_water_box_has_three_sites_per_molecule_and_bonded_topology():
    system = build_water_box(2, box_length=3.0)

    assert system.n_particles == 6
    assert [particle.atom_type for particle in system.particles] == [
        "O",
        "H",
        "H",
        "O",
        "H",
        "H",
    ]
    assert system.topology.n_bonds == 4
    assert system.topology.n_angles == 2
    assert system.topology.n_torsions == 0
    assert len(system.constraints) == 6
    assert system.nonbonded_model is not None
    assert system.nonbonded_model.cutoff == pytest.approx(1.35)
    assert system.box.volume == pytest.approx(27.0)


def test_water_charges_are_stored_on_particles():
    system = build_water_box(1, box_length=2.0)
    charges = [particle.charge for particle in system.particles]

    assert charges == pytest.approx(
        [SPCE_WATER.oxygen_charge, SPCE_WATER.hydrogen_charge, SPCE_WATER.hydrogen_charge]
    )


def test_single_water_has_no_intramolecular_nonbonded_energy():
    system = build_water_box(1, box_length=2.0)

    assert len(system.constraints) == 3
    assert system.nonbonded_energy() == pytest.approx(0.0)
    assert system.bonded_energy() == pytest.approx(0.0, abs=1e-9)
    np.testing.assert_allclose(system.forces(), np.zeros((3, 3)), atol=1e-8)


def test_water_exclusions_remove_all_intramolecular_pairs():
    system = build_water_box(1, box_length=2.0)

    assert system.topology.nonbonded_exclusions() == frozenset(
        {(0, 1), (0, 2), (1, 2)}
    )


def test_water_box_from_density():
    n_molecules = 8
    density = 10.0

    system = build_water_box(n_molecules, density=density)

    expected_edge = (n_molecules / density) ** (1.0 / 3.0)
    np.testing.assert_allclose(system.box.lengths, [expected_edge] * 3)


def test_water_temperature_initialization():
    system = build_water_box(
        8,
        box_length=4.0,
        temperature=300.0,
        seed=17,
    )

    assert kinetic_temperature(system) == pytest.approx(300.0, rel=0.2)


def test_water_system_runs_with_thermostat_and_integrator():
    system = build_water_box(
        2,
        box_length=4.0,
        temperature=300.0,
        seed=3,
    )

    heated = berendsen_thermostat(
        system,
        temperature=350.0,
        coupling_time=0.1,
        dt=0.0005,
    )
    stepped = velocity_verlet_step_system(heated, dt=0.0005)

    assert stepped.box is system.box
    assert stepped.nonbonded_model is system.nonbonded_model
    assert stepped.constraints is system.constraints
    assert [particle.atom_type for particle in stepped.particles] == [
        "O",
        "H",
        "H",
        "O",
        "H",
        "H",
    ]
    assert np.all(np.isfinite(stepped.positions))
    assert np.all(np.isfinite(stepped.velocities))
    assert kinetic_temperature(heated) > kinetic_temperature(system)


def test_water_builder_validates_inputs():
    with pytest.raises(ValueError, match="n_molecules"):
        build_water_box(0, box_length=2.0)
    with pytest.raises(ValueError, match="exactly one"):
        build_water_box(1)
    with pytest.raises(ValueError, match="exactly one"):
        build_water_box(1, box_length=2.0, density=10.0)
    with pytest.raises(ValueError, match="cutoff"):
        build_water_box(1, box_length=2.0, cutoff=1.5)
    with pytest.raises(ValueError, match="too small"):
        build_water_box(1, box_length=0.15)


def test_water_model_validates_parameters():
    with pytest.raises(ValueError, match="oxygen_mass"):
        WaterModel(oxygen_mass=0.0)
    with pytest.raises(ValueError, match="angle"):
        WaterModel(angle=0.0)
    with pytest.raises(ValueError, match="oxygen_epsilon"):
        WaterModel(oxygen_epsilon=-1.0)
    with pytest.raises(ValueError, match="dielectric"):
        WaterModel(dielectric=0.0)


def test_water_can_use_ewald_electrostatics():
    ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)
    system = build_water_box(
        1,
        box_length=2.0,
        electrostatics=ewald,
    )

    assert system.nonbonded_model.electrostatics is ewald

    reference = 0.0
    for particle in system.particles:
        position = np.array([particle.position])
        charge = np.array([particle.charge])
        reference += ewald.global_energy(
            position, charge, system.box
        ) - ewald.charge_correction(charge, system.box)

    assert system.potential_energy() == pytest.approx(reference, abs=1e-8)


def test_ewald_water_runs_with_rigid_constraints():
    ewald = EwaldElectrostatics(alpha=6.0, reciprocal_cutoff=40.0)
    system = build_water_box(
        2,
        box_length=3.0,
        temperature=300.0,
        seed=9,
        electrostatics=ewald,
    )

    stepped = velocity_verlet_step_system(system, dt=5e-5)

    assert stepped.nonbonded_model is system.nonbonded_model
    assert stepped.constraints is system.constraints
    assert np.all(np.isfinite(stepped.positions))
    assert np.isfinite(stepped.potential_energy())
    for constraint in system.constraints:
        distance = np.linalg.norm(
            stepped.positions[constraint.i] - stepped.positions[constraint.j]
        )
        assert distance == pytest.approx(constraint.distance, abs=1e-8)


def test_water_builder_rejects_invalid_electrostatics_type():
    with pytest.raises(TypeError, match="electrostatics"):
        build_water_box(1, box_length=2.0, electrostatics="ewald")
