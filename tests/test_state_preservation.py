"""Software-correctness tests for state preservation.

Every rebuild path must keep each particle's chemical identity.  Losing
``atom_type`` or ``charge`` would silently change a nonbonded calculation while
leaving the trajectory numerically plausible.
"""

import pytest

from moleculab import (
    LennardJonesPotential,
    Particle,
    System,
    berendsen_thermostat,
    minimize_energy,
    steepest_descent_step,
    velocity_verlet_step_system,
)


def typed_pair_system():
    return System(
        particles=[
            Particle(
                mass=12.0,
                position=[0.0, 0.0, 0.0],
                velocity=[0.10, 0.0, 0.0],
                atom_type="C",
                charge=-0.20,
            ),
            Particle(
                mass=1.0,
                position=[0.40, 0.0, 0.0],
                velocity=[0.0, 0.05, 0.0],
                atom_type="H",
                charge=0.20,
            ),
        ],
        pair_potentials=[LennardJonesPotential(epsilon=1.0, sigma=0.3)],
    )


def assert_identity_preserved(result, reference):
    assert [particle.atom_type for particle in result.particles] == [
        particle.atom_type for particle in reference.particles
    ]
    assert [particle.charge for particle in result.particles] == pytest.approx(
        [particle.charge for particle in reference.particles]
    )


def test_velocity_verlet_step_system_preserves_identity():
    system = typed_pair_system()

    updated = velocity_verlet_step_system(system, dt=1e-4)

    assert_identity_preserved(updated, system)


def test_berendsen_thermostat_preserves_identity():
    system = typed_pair_system()

    updated = berendsen_thermostat(
        system,
        temperature=300.0,
        coupling_time=0.1,
        dt=0.001,
    )

    assert_identity_preserved(updated, system)


def test_steepest_descent_step_preserves_identity():
    system = typed_pair_system()

    updated = steepest_descent_step(system, step_size=1e-8)

    assert_identity_preserved(updated, system)


def test_minimize_energy_preserves_identity():
    system = typed_pair_system()

    result = minimize_energy(system, step_size=1e-8, max_steps=2)

    assert_identity_preserved(result.system, system)


def test_with_state_is_the_single_preserving_rebuild_path():
    system = typed_pair_system()

    updated = system.with_state(system.positions, system.velocities)

    assert_identity_preserved(updated, system)
    assert updated.pair_potentials == system.pair_potentials
    assert updated.topology is system.topology
    assert updated.box is system.box
    assert updated.nonbonded_model is system.nonbonded_model
    assert updated.constraints == system.constraints
    assert updated.indexed_pair_potentials == system.indexed_pair_potentials


def test_with_state_rejects_wrong_shapes():
    system = typed_pair_system()

    with pytest.raises(ValueError, match="positions"):
        system.with_state([[0.0, 0.0, 0.0]], system.velocities)
    with pytest.raises(ValueError, match="velocities"):
        system.with_state(system.positions, [[0.0, 0.0, 0.0]])
