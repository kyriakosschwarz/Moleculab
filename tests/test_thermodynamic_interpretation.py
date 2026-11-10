"""Tests for structural-to-thermodynamic interpretation records."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from moleculab import (
    THERMODYNAMIC_INTERPRETATION_SCHEMA_VERSION,
    RobustnessStatusCount,
    StructuralEvidenceSummary,
    default_free_energy_requirements,
    default_free_energy_terms,
    load_thermodynamic_interpretation,
    structural_thermodynamic_boundaries,
    validate_thermodynamic_interpretation,
    write_thermodynamic_interpretation,
)


def reference_record():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.thermodynamic-interpretation.json"
    )
    return load_thermodynamic_interpretation(path)


def test_reference_2hyy_thermodynamic_record_is_non_quantitative():
    result = reference_record()

    assert validate_thermodynamic_interpretation(result) == []
    assert result.schema_version == THERMODYNAMIC_INTERPRETATION_SCHEMA_VERSION
    assert result.accession == "2HYY"
    assert result.free_energy_status == "not_computed"
    assert result.binding_free_energy_kj_per_mol is None
    assert result.ensemble_status == "not_available"
    assert len(result.boundaries) == 8
    assert len(result.terms) == 8
    assert len(result.requirements) == 10

    summary = result.structural_summary
    assert summary.ligand_state_count == 8
    assert summary.protein_state_count == 2
    assert summary.combination_count == 16
    assert summary.hydrophobic_pair_count == 18
    assert summary.aromatic_pair_count == 1
    assert summary.supported_hydrogen_bond_count == 2
    assert summary.plausible_water_bridge_count == 1
    assert {
        item.status: item.count for item in summary.robustness_conclusion_counts
    } == {
        "invariant": 4,
        "assumption_sensitive": 8,
        "unsupported": 3,
        "unknown": 2,
    }

    assert all(term.status == "not_estimated" for term in result.terms)
    assert {item.status for item in result.requirements} == {
        "partial",
        "not_run",
        "not_available",
    }
    assert all(
        boundary.required_evidence
        for boundary in result.boundaries
        if boundary.status in ("not_established", "not_estimated")
    )


def test_structural_boundaries_keep_counts_away_from_energies():
    summary = StructuralEvidenceSummary(
        ligand_state_count=8,
        protein_state_count=2,
        combination_count=16,
        ligand_contact_atom_count=34,
        protein_contact_atom_count=53,
        hydrophobic_pair_count=18,
        aromatic_pair_count=1,
        supported_hydrogen_bond_count=2,
        plausible_water_bridge_count=1,
        robustness_conclusion_counts=(
            RobustnessStatusCount("invariant", 4),
            RobustnessStatusCount("assumption_sensitive", 8),
            RobustnessStatusCount("unsupported", 3),
            RobustnessStatusCount("unknown", 2),
        ),
    )

    boundaries = structural_thermodynamic_boundaries(summary)
    identifiers = {item.boundary_id for item in boundaries}

    assert identifiers == {
        "hydrophobic_pairs_not_energy",
        "aromatic_pair_not_energy",
        "hydrogen_bonds_not_energy",
        "water_bridge_not_energy",
        "counts_not_ranking",
        "snapshot_not_ensemble",
        "robustness_not_certainty",
        "states_not_populations",
    }
    assert all(
        item.status in ("not_established", "not_estimated") for item in boundaries
    )
    assert all(item.required_evidence for item in boundaries)
    assert all("kJ" not in item.thermodynamic_claim for item in boundaries)


def test_default_terms_and_requirements_are_explicitly_incomplete():
    terms = default_free_energy_terms()
    requirements = default_free_energy_requirements()

    assert {item.term_id for item in terms} == {
        "interaction_enthalpy",
        "desolvation",
        "ligand_conformational_entropy",
        "protein_conformational_entropy",
        "water_reorganisation",
        "translational_rotational_entropy",
        "protonation_tautomer_populations",
        "conformational_strain",
    }
    assert all(item.status == "not_estimated" for item in terms)
    assert all(item.required_evidence for item in terms)

    assert {item.requirement_id for item in requirements} == {
        "molecular_model",
        "protein_preparation",
        "ligand_preparation",
        "water_ion_preparation",
        "force_field_parameters",
        "conformational_sampling",
        "thermodynamic_states",
        "standard_state",
        "convergence_uncertainty",
        "experimental_affinity_validation",
    }
    assert any(item.status != "met" for item in requirements)


def test_record_round_trip(tmp_path):
    result = reference_record()
    destination = write_thermodynamic_interpretation(
        tmp_path / "thermodynamics.json",
        result,
    )
    loaded = load_thermodynamic_interpretation(destination)

    assert loaded == result
    assert loaded.to_dict()["binding_free_energy_kj_per_mol"] is None


def test_non_null_free_energy_is_rejected():
    result = reference_record()
    bad = replace(result, binding_free_energy_kj_per_mol=-42.0)
    errors = validate_thermodynamic_interpretation(bad)

    assert any("binding_free_energy_kj_per_mol" in error for error in errors)


def test_invalid_term_and_requirement_statuses_are_rejected():
    result = reference_record()
    bad_term = replace(result.terms[0], status="estimated")
    bad_requirement = replace(result.requirements[0], status="assumed")
    errors = validate_thermodynamic_interpretation(
        replace(
            result,
            terms=(bad_term,) + result.terms[1:],
            requirements=(bad_requirement,) + result.requirements[1:],
        )
    )

    assert any("terms[0].status" in error for error in errors)
    assert any("requirements[0].status" in error for error in errors)


def test_record_contains_no_affinity_or_ranking_payload():
    result = reference_record()
    serialized = result.to_dict()

    assert "affinity" not in serialized
    assert "score" not in serialized
    assert "rank" not in serialized
    assert serialized["binding_free_energy_kj_per_mol"] is None
