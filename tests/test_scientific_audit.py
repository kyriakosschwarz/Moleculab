"""Tests for the final capstone audit and end-to-end demonstration."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from moleculab import (
    SCIENTIFIC_AUDIT_SCHEMA_VERSION,
    load_scientific_audit,
    run_end_to_end_demonstration,
    validate_scientific_audit,
    write_scientific_audit,
)


def reference_audit():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "free_energy"
        / "capstone-scientific-audit.json"
    )
    return load_scientific_audit(path)


def test_reference_capstone_audit_passes_all_checks():
    record = reference_audit()

    assert validate_scientific_audit(record) == []
    assert record.schema_version == SCIENTIFIC_AUDIT_SCHEMA_VERSION
    assert record.status == "pass"
    assert record.fail_count == 0
    assert record.not_run_count == 0
    assert record.pass_count == len(record.checks)
    assert record.pass_count == 23

    categories = {check.category for check in record.checks}
    assert categories == {
        "documentation",
        "physics",
        "chemistry",
        "drug_discovery",
        "thermodynamics",
        "reproducibility",
    }
    assert len({check.check_id for check in record.checks}) == len(record.checks)


def test_end_to_end_demonstration_connects_ensemble_and_free_energy():
    demonstration = run_end_to_end_demonstration(
        temperature_kelvin=300.0,
        n_samples=20_000,
        seed=123,
    )

    assert demonstration.mean_potential_energy_kj_per_mol == pytest.approx(
        demonstration.analytic_mean_potential_energy_kj_per_mol,
        rel=0.05,
    )
    assert demonstration.fep_delta_g_kj_per_mol == pytest.approx(
        demonstration.exact_delta_g_kj_per_mol,
        abs=0.15,
    )
    assert 0.0 < demonstration.overlap.effective_sample_fraction <= 1.0
    assert 0.0 < demonstration.proposal_effective_sample_fraction <= 1.0
    assert 0.0 < demonstration.proposal_maximum_weight_fraction <= 1.0
    assert (
        demonstration.two_state_reference.population0
        + demonstration.two_state_reference.population1
    ) == pytest.approx(1.0)
    assert demonstration.two_state_reference.delta_g_kj_per_mol < 0.0
    assert demonstration.stop_line
    assert any(
        "not converted into a binding free energy" in item
        for item in demonstration.stop_line
    )


def test_audit_record_round_trip(tmp_path):
    record = reference_audit()
    destination = write_scientific_audit(tmp_path / "audit.json", record)
    loaded = load_scientific_audit(destination)

    assert loaded == record
    assert loaded.to_dict()["status"] == "pass"


def test_audit_status_and_counts_must_match_checks():
    record = reference_audit()
    failed_check = replace(record.checks[0], status="fail")
    bad = replace(
        record,
        checks=(failed_check,) + record.checks[1:],
        pass_count=len(record.checks) - 1,
        fail_count=1,
        status="pass",
    )
    errors = validate_scientific_audit(bad)

    assert any("status must be 'fail'" in error for error in errors)


def test_audit_record_does_not_claim_binding_affinity():
    record = reference_audit()
    serialized = record.to_dict()

    assert "affinity" not in serialized
    assert "binding_free_energy_kj_per_mol" not in serialized
    assert "rank" not in serialized
    assert "score" not in serialized
    assert record.demonstration.exact_delta_g_kj_per_mol == pytest.approx(
        2.5934158446919926
    )
    assert any(
        "not converted into a binding free energy" in item
        for item in record.demonstration.stop_line
    )


def test_invalid_demonstration_overlap_is_rejected():
    record = reference_audit()
    bad_demo = replace(
        record.demonstration,
        overlap=replace(
            record.demonstration.overlap,
            effective_sample_fraction=0.0,
        ),
    )
    errors = validate_scientific_audit(
        replace(record, demonstration=bad_demo)
    )

    assert any("effective sample fraction" in error for error in errors)
