"""Tests for pressure and the Berendsen barostat."""

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
    System,
    Topology,
    berendsen_barostat,
    build_water_box,
    kinetic_pressure,
    pressure,
    pressure_bar,
    virial_pressure,
)
from moleculab.units import KILOJOULES_PER_MOLE_PER_NM3_TO_BAR


def ideal_gas_system():
    return System(
        particles=[
            Particle(1.0, [0.5, 0.5, 0.5], [1.0, 0.0, 0.0])
        ],
        box=Box([1.0, 1.0, 1.0]),
    )


def lj_pair_system():
    return System(
        particles=[
            Particle(1.0, [0.0, 0.0, 0.0]),
            Particle(1.0, [0.4, 0.0, 0.0]),
        ],
        pair_potentials=[LennardJonesPotential(epsilon=1.0, sigma=0.3)],
        box=Box([2.0, 2.0, 2.0]),
    )


def test_internal_pressure_to_bar_conversion():
    assert KILOJOULES_PER_MOLE_PER_NM3_TO_BAR == pytest.approx(16.6054, rel=1e-5)


def test_ideal_gas_pressure_is_kinetic_only():
    system = ideal_gas_system()

    expected_kinetic = 2.0 * 0.5 / (3.0 * 1.0)
    assert kinetic_pressure(system) == pytest.approx(expected_kinetic)
    assert virial_pressure(system) == pytest.approx(0.0, abs=1e-12)
    assert pressure(system) == pytest.approx(expected_kinetic)
    assert pressure_bar(system) == pytest.approx(
        expected_kinetic * KILOJOULES_PER_MOLE_PER_NM3_TO_BAR
    )


def test_lj_pair_virial_matches_pair_force_dot_product():
    system = lj_pair_system()
    displacement = system.box.minimum_image(
        system.particles[0].position - system.particles[1].position
    )
    force = system.forces()[0]
    expected = float(np.dot(displacement, force)) / (3.0 * system.box.volume)

    assert virial_pressure(system) == pytest.approx(expected, rel=1e-6)


def test_bonded_pair_virial_matches_pair_force_dot_product():
    system = System(
        particles=[
            Particle(1.0, [-0.5, 0.0, 0.0]),
            Particle(1.0, [0.5, 0.0, 0.0]),
        ],
        topology=Topology(
            bonds=[Bond(0, 1, HarmonicBondPotential(1000.0, 0.8))]
        ),
        box=Box([2.0, 2.0, 2.0]),
    )
    displacement = system.box.minimum_image(
        system.particles[0].position - system.particles[1].position
    )
    force = system.forces()[0]
    expected = float(np.dot(displacement, force)) / (3.0 * system.box.volume)

    assert virial_pressure(system) == pytest.approx(expected, rel=1e-6)


def test_pressure_is_kinetic_plus_virial():
    system = lj_pair_system()

    assert pressure(system) == pytest.approx(
        kinetic_pressure(system) + virial_pressure(system)
    )


def test_barostat_compresses_for_higher_target_pressure():
    system = ideal_gas_system()
    current_bar = pressure_bar(system)

    compressed = berendsen_barostat(
        system,
        target_pressure=2.0 * current_bar,
        coupling_time=0.1,
        dt=0.1,
        compressibility=1e-3,
    )

    assert compressed.box.volume < system.box.volume
    np.testing.assert_allclose(
        compressed.positions / compressed.box.lengths,
        system.positions / system.box.lengths,
        atol=1e-12,
    )


def test_barostat_expands_for_lower_target_pressure():
    system = ideal_gas_system()
    current_bar = pressure_bar(system)

    expanded = berendsen_barostat(
        system,
        target_pressure=0.5 * current_bar,
        coupling_time=0.1,
        dt=0.1,
        compressibility=1e-3,
    )

    assert expanded.box.volume > system.box.volume


def test_barostat_preserves_types_charges_topology_and_model():
    system = build_water_box(
        1,
        box_length=2.0,
        temperature=300.0,
        seed=4,
    )
    current_bar = pressure_bar(system)

    scaled = berendsen_barostat(
        system,
        target_pressure=current_bar + 10.0,
        coupling_time=0.5,
        dt=0.001,
    )

    assert scaled.topology is system.topology
    assert scaled.nonbonded_model is system.nonbonded_model
    assert [particle.atom_type for particle in scaled.particles] == [
        particle.atom_type for particle in system.particles
    ]
    assert [particle.charge for particle in scaled.particles] == [
        particle.charge for particle in system.particles
    ]
    assert scaled.box.volume != system.box.volume


def test_barostat_requires_a_box():
    system = System(particles=[Particle(1.0, [0.0, 0.0, 0.0])])

    with pytest.raises(ValueError, match="Box"):
        berendsen_barostat(system, 1.0, 0.5, 0.01)


def test_barostat_validates_inputs():
    system = ideal_gas_system()

    with pytest.raises(ValueError, match="target_pressure"):
        berendsen_barostat(system, -1.0, 0.5, 0.01)
    with pytest.raises(ValueError, match="coupling_time"):
        berendsen_barostat(system, 1.0, 0.0, 0.01)
    with pytest.raises(ValueError, match="dt"):
        berendsen_barostat(system, 1.0, 0.5, 0.0)
    with pytest.raises(ValueError, match="compressibility"):
        berendsen_barostat(system, 1.0, 0.5, 0.01, compressibility=0.0)
    with pytest.raises(ValueError, match="strain"):
        virial_pressure(system, strain=1.0)


def test_barostat_rejects_non_positive_scale_factor():
    system = ideal_gas_system()

    with pytest.raises(ValueError, match="scale factor"):
        berendsen_barostat(
            system,
            target_pressure=1.0e12,
            coupling_time=0.001,
            dt=1.0,
            compressibility=1.0,
        )
