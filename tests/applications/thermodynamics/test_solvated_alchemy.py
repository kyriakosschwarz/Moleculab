# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for solvated-ligand soft-core decoupling in an Ewald water box."""

import numpy as np
import pytest

from moleculab import (
    EwaldElectrostatics,
    Particle,
    SolvatedLigandAlchemy,
    System,
    build_water_box,
    velocity_verlet_step_system,
)


def make_alchemy(n_water=2):
    ewald = EwaldElectrostatics(alpha=4.0, reciprocal_cutoff=8.0)
    water = build_water_box(
        n_water,
        box_length=2.0,
        temperature=300.0,
        seed=3,
        electrostatics=ewald,
    )
    return SolvatedLigandAlchemy(
        water,
        ligand_epsilon=1.0,
        ligand_sigma=0.3,
        temperature=300.0,
        dt=1e-4,
        equilibration_steps=50,
        sample_every=5,
        coupling_time=0.05,
    )


def test_solvated_system_contains_ligand_and_indexed_pairs():
    alchemy = make_alchemy()
    system = alchemy.make_system(0.5)

    assert system.n_particles == alchemy.n_water + 1
    assert system.particles[alchemy.ligand_index].atom_type == "L"
    assert system.particles[alchemy.ligand_index].charge == 0.0
    assert len(system.indexed_pair_potentials) == alchemy.n_water
    assert system.constraints is alchemy.water_system.constraints
    assert system.nonbonded_model.electrostatics is not None
    assert np.isfinite(system.potential_energy())
    assert np.all(np.isfinite(system.forces()))


def test_indexed_pair_potentials_survive_integration():
    alchemy = make_alchemy()
    system = alchemy.make_system(0.5)

    stepped = velocity_verlet_step_system(system, dt=1e-4)

    assert stepped.indexed_pair_potentials == system.indexed_pair_potentials
    assert stepped.constraints is system.constraints


def test_water_sampling_returns_finite_water_configurations():
    alchemy = make_alchemy()
    positions = alchemy.sample_water_positions(0.5, n_samples=5, seed=7)

    assert positions.shape == (5, alchemy.n_water, 3)
    assert np.all(np.isfinite(positions))


def test_fep_ti_mbar_cross_validate_solvated_ligand():
    alchemy = make_alchemy()
    n_windows = 5
    n_samples = 10

    fep = alchemy.fep(n_windows=n_windows, n_samples=n_samples, seed=11)
    ti = alchemy.ti(n_windows=n_windows, n_samples=n_samples, seed=12)
    mbar = alchemy.mbar(n_windows=n_windows, n_samples=n_samples, seed=13)

    values = np.array([fep.delta_g, ti.delta_g, mbar.delta_g[-1]])
    assert np.ptp(values) < 0.005
    assert fep.standard_error > 0.0
    assert ti.standard_error > 0.0
    assert mbar.converged


def test_solvated_alchemy_validation():
    with pytest.raises(ValueError, match="Box"):
        SolvatedLigandAlchemy(
            System(particles=[Particle(1.0, [0.0, 0.0, 0.0])])
        )

    alchemy = make_alchemy()
    with pytest.raises(ValueError, match="ligand_epsilon"):
        SolvatedLigandAlchemy(alchemy.water_system, ligand_epsilon=0.0)
    with pytest.raises(ValueError, match="temperature"):
        SolvatedLigandAlchemy(alchemy.water_system, temperature=0.0)
    with pytest.raises(ValueError, match="ligand_position"):
        SolvatedLigandAlchemy(alchemy.water_system, ligand_position=[0.0, 0.0])
    with pytest.raises(ValueError, match="n_windows"):
        alchemy.fep(n_windows=1)
    with pytest.raises(ValueError, match="n_samples"):
        alchemy.sample_water_positions(0.5, n_samples=0)
    with pytest.raises(ValueError, match="odd"):
        alchemy.ti(n_windows=4)
