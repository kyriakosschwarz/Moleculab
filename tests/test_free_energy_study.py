# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for the harmonic free-energy convergence study."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

from moleculab import (
    FREE_ENERGY_STUDY_SCHEMA_VERSION,
    HarmonicPotential,
    build_harmonic_fep_study_record,
    fep_harmonic_oscillator,
    fep_overlap_diagnostics,
    harmonic_configurational_free_energy,
    harmonic_configurational_partition_function,
    harmonic_free_energy_difference,
    load_harmonic_fep_study,
    run_harmonic_fep_case,
    two_state_free_energy_difference,
    two_state_partition_function,
    two_state_populations,
    two_state_reference,
    validate_harmonic_fep_study,
    write_harmonic_fep_study,
)


def reference_study():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "free_energy"
        / "harmonic-fep-study.json"
    )
    return load_harmonic_fep_study(path)


def test_exact_partition_function_and_free_energy_are_consistent():
    potential = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    temperature = 300.0
    partition_function = harmonic_configurational_partition_function(
        potential, temperature
    )
    free_energy = harmonic_configurational_free_energy(potential, temperature)

    assert partition_function > 0
    assert free_energy == pytest.approx(
        -0.008314462618 * temperature * np.log(partition_function)
    )


def test_two_state_degeneracy_enters_partition_function_and_free_energy():
    partition_function = two_state_partition_function(
        0.0, 5.0, 1, 100, 300.0
    )
    populations = two_state_populations(0.0, 5.0, 1, 100, 300.0)
    delta_g = two_state_free_energy_difference(0.0, 5.0, 1, 100, 300.0)

    rt = 0.008314462618 * 300.0
    assert partition_function == pytest.approx(1.0 + 100.0 * np.exp(-5.0 / rt))
    assert populations[0] + populations[1] == pytest.approx(1.0)
    assert populations[1] > 0.9
    assert delta_g == pytest.approx(5.0 - rt * np.log(100.0))


def test_two_state_reference_reports_populations_and_delta_g():
    reference = two_state_reference(
        energy0_kj_per_mol=0.0,
        energy1_kj_per_mol=5.0,
        degeneracy0=1,
        degeneracy1=100,
        temperature=300.0,
    )

    assert reference.partition_function > 0
    assert reference.population0 + reference.population1 == pytest.approx(1.0)
    assert reference.delta_g_kj_per_mol < 0


def test_exact_partition_function_reproduces_the_harmonic_difference():
    potential0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    potential1 = HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0])
    temperature = 300.0

    z0 = harmonic_configurational_partition_function(potential0, temperature)
    z1 = harmonic_configurational_partition_function(potential1, temperature)
    from_partition_functions = (
        -0.008314462618 * temperature * (np.log(z1) - np.log(z0))
    )

    assert from_partition_functions == pytest.approx(
        harmonic_free_energy_difference(
            potential0, potential1, temperature
        ),
        abs=1e-12,
    )


def test_overlap_diagnostics_distinguish_good_and_poor_overlap():
    temperature = 300.0
    mild0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    mild1 = HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0])
    poor0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    poor1 = HarmonicPotential(stiffness=100.0, center=[0.5, 0.0, 0.0])

    mild = fep_harmonic_oscillator(
        mild0, mild1, temperature, n_samples=20_000, seed=1
    )
    poor = fep_harmonic_oscillator(
        poor0, poor1, temperature, n_samples=20_000, seed=2
    )
    mild_overlap = fep_overlap_diagnostics(mild.work, temperature)
    poor_overlap = fep_overlap_diagnostics(poor.work, temperature)

    assert mild_overlap.effective_sample_fraction > 0.3
    assert poor_overlap.effective_sample_fraction < 0.05
    assert poor_overlap.maximum_weight_fraction > mild_overlap.maximum_weight_fraction
    assert 0 < poor_overlap.effective_sample_fraction <= 1
    assert 0 < poor_overlap.maximum_weight_fraction <= 1


def test_case_study_reports_exact_and_estimated_values():
    potential0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    potential1 = HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0])
    case = run_harmonic_fep_case(
        case_id="test",
        label="test case",
        potential0=potential0,
        potential1=potential1,
        temperature=300.0,
        sample_sizes=(1_000, 10_000),
        n_replicates=5,
        replicate_sample_size=2_000,
        seed=7,
    )

    assert case.exact_delta_g_kj_per_mol == pytest.approx(
        harmonic_free_energy_difference(potential0, potential1, 300.0)
    )
    assert len(case.convergence) == 2
    assert case.replicates.n_replicates == 5
    assert all(
        point.overlap.effective_sample_fraction > 0 for point in case.convergence
    )


def test_study_record_round_trip_and_validation(tmp_path):
    mild0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    mild1 = HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0])
    cases = (
        run_harmonic_fep_case(
            case_id="mild",
            label="mild",
            potential0=mild0,
            potential1=mild1,
            temperature=300.0,
            sample_sizes=(1_000, 5_000),
            n_replicates=4,
            replicate_sample_size=1_000,
            seed=11,
        ),
    )
    record = build_harmonic_fep_study_record(
        temperature_kelvin=300.0,
        two_state_example=two_state_reference(
            energy0_kj_per_mol=0.0,
            energy1_kj_per_mol=5.0,
            degeneracy0=1,
            degeneracy1=100,
            temperature=300.0,
        ),
        sample_sizes=(1_000, 5_000),
        replicate_count=4,
        replicate_sample_size=1_000,
        cases=cases,
    )
    destination = write_harmonic_fep_study(tmp_path / "study.json", record)
    loaded = load_harmonic_fep_study(destination)

    assert validate_harmonic_fep_study(record) == []
    assert loaded == record
    assert loaded.schema_version == FREE_ENERGY_STUDY_SCHEMA_VERSION


def test_reference_study_contains_good_and_poor_overlap_cases():
    record = reference_study()

    assert validate_harmonic_fep_study(record) == []
    assert record.temperature_kelvin == pytest.approx(300.0)
    assert record.two_state_example.population1 > 0.9
    assert record.two_state_example.delta_g_kj_per_mol < 0
    assert tuple(case.case_id for case in record.cases) == (
        "mild-stiffness-change",
        "poor-overlap-center-shift",
    )

    mild, poor = record.cases
    assert mild.exact_delta_g_kj_per_mol > 0
    assert poor.exact_delta_g_kj_per_mol == pytest.approx(0.0)
    assert (
        mild.exact_partition_function0
        != pytest.approx(mild.exact_partition_function1)
    )
    assert poor.exact_partition_function0 == pytest.approx(
        poor.exact_partition_function1
    )

    mild_largest = mild.convergence[-1].overlap.effective_sample_fraction
    poor_largest = poor.convergence[-1].overlap.effective_sample_fraction
    assert mild_largest > 0.3
    assert poor_largest < 0.01
    assert poor.convergence[-1].overlap.maximum_weight_fraction > 0.01


def test_invalid_sample_grid_and_overlap_values_are_rejected():
    record = reference_study()
    case = record.cases[0]
    bad_point = replace(
        case.convergence[0],
        overlap=replace(
            case.convergence[0].overlap,
            effective_sample_fraction=0.0,
        ),
    )
    bad_case = replace(
        case,
        convergence=(bad_point,) + case.convergence[1:],
    )
    bad = replace(record, cases=(bad_case,) + record.cases[1:])
    errors = validate_harmonic_fep_study(bad)

    assert any("effective sample" in error for error in errors)


def test_study_record_has_no_binding_affinity_payload():
    record = reference_study()
    serialized = record.to_dict()

    assert "binding_free_energy_kj_per_mol" not in serialized
    assert "affinity" not in serialized
    assert "rank" not in serialized
    assert "score" not in serialized
    assert any("not a protein-ligand" in limitation for limitation in record.limitations)
