# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for chemistry-aware interaction interpretation."""

from __future__ import annotations

from pathlib import Path

import pytest

from moleculab import (
    INTERACTION_SCHEMA_VERSION,
    AromaticRing,
    AtomFeature,
    InteractionCriteria,
    LigandSelection,
    ProteinSelection,
    interpret_interactions,
    load_interactions,
    validate_interactions,
    write_interactions,
)

CHECKSUM = "e" * 64


def protein_selection() -> ProteinSelection:
    return ProteinSelection(label_asym_id="A", auth_asym_id="A")


def ligand_selection() -> LigandSelection:
    return LigandSelection(
        component_id="LIG",
        label_asym_id="E",
        auth_asym_id="A",
        instance="E",
        auth_seq_id="600",
        available_instances=("E",),
    )


def feature(
    atom_id: str,
    side: str,
    position: tuple[float, float, float],
    *,
    element: str = "C",
    donor: bool | None = None,
    acceptor: bool | None = None,
    charge: int | None = None,
    hydrophobic: bool | None = None,
    ring_id: str | None = None,
    metal: bool | None = None,
    hydrogens: tuple[tuple[float, float, float], ...] = (),
) -> AtomFeature:
    return AtomFeature(
        atom_id=atom_id,
        side=side,
        element=element,
        position_nm=position,
        residue_id=f"{side.upper()} 1",
        is_donor=donor,
        is_acceptor=acceptor,
        formal_charge=charge,
        is_hydrophobic=hydrophobic,
        aromatic_ring_id=ring_id,
        is_metal=metal,
        hydrogen_positions_nm=hydrogens,
        feature_source="test fixture",
    )


def interpret(
    protein,
    ligand,
    *,
    water=(),
    metal=(),
    rings=(),
    criteria=None,
    contact_count: int = 1,
):
    return interpret_interactions(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_identity_record_path="identity.json",
        ligand_identity_record_sha256=CHECKSUM,
        binding_site_record_path="site.json",
        binding_site_record_sha256=CHECKSUM,
        accession="TEST",
        model_number=1,
        protein_selection=protein_selection(),
        ligand_selection=ligand_selection(),
        input_contact_count=contact_count,
        protein_features=tuple(protein),
        ligand_features=tuple(ligand),
        water_features=tuple(water),
        metal_features=tuple(metal),
        aromatic_rings=tuple(rings),
        criteria=criteria,
    )


def status(result, interaction_type: str) -> list[str]:
    return [
        item.status
        for item in result.interactions
        if item.interaction_type == interaction_type
    ]


def test_hydrogen_bond_is_classified_when_features_and_geometry_are_present():
    protein = [
        feature(
            "N1",
            "protein",
            (0.0, 0.0, 0.0),
            element="N",
            donor=True,
            hydrogens=((0.1, 0.0, 0.0),),
        )
    ]
    ligand = [feature("O1", "ligand", (0.3, 0.0, 0.0), element="O", acceptor=True)]
    result = interpret(protein, ligand)

    assert "classified" in status(result, "hydrogen_bond")
    hbond = next(item for item in result.interactions if item.interaction_type == "hydrogen_bond")
    assert hbond.distance_nm == pytest.approx(0.3)
    assert hbond.angle_degrees == pytest.approx(180.0)
    assert validate_interactions(result) == []


def test_hydrogen_bond_remains_unknown_without_hydrogen_positions():
    protein = [
        feature("N1", "protein", (0.0, 0.0, 0.0), element="N", donor=True)
    ]
    ligand = [feature("O1", "ligand", (0.3, 0.0, 0.0), element="O", acceptor=True)]
    result = interpret(protein, ligand)

    assert status(result, "hydrogen_bond") == ["unknown"]
    assert "hydrogen" in result.interactions[0].limitations.lower()


def test_missing_donor_acceptor_assignment_remains_unknown():
    protein = [feature("C1", "protein", (0.0, 0.0, 0.0))]
    ligand = [feature("C2", "ligand", (0.3, 0.0, 0.0))]
    result = interpret(protein, ligand)

    assert status(result, "hydrogen_bond") == ["unknown"]
    assert status(result, "ionic") == ["unknown"]
    assert status(result, "hydrophobic") == ["unknown"]
    assert status(result, "aromatic") == ["unknown"]
    assert status(result, "water_mediated") == ["not_applicable"]
    assert status(result, "metal_coordination") == ["not_applicable"]
    assert all(item.status != "classified" for item in result.interactions)


def test_ionic_classification_requires_explicit_charges():
    protein = [feature("N1", "protein", (0.0, 0.0, 0.0), element="N", charge=1)]
    ligand = [feature("O1", "ligand", (0.3, 0.0, 0.0), element="O", charge=-1)]
    result = interpret(protein, ligand)

    assert "classified" in status(result, "ionic")
    ionic = next(item for item in result.interactions if item.interaction_type == "ionic")
    assert ionic.distance_nm == pytest.approx(0.3)
    assert "+1" in ionic.evidence and "-1" in ionic.evidence


def test_hydrophobic_classification_requires_explicit_flags():
    protein = [feature("C1", "protein", (0.0, 0.0, 0.0), hydrophobic=True)]
    ligand = [feature("C2", "ligand", (0.3, 0.0, 0.0), hydrophobic=True)]
    result = interpret(protein, ligand)

    assert "classified" in status(result, "hydrophobic")
    assert any(
        item.interaction_type == "hydrophobic" and item.status == "classified"
        for item in result.interactions
    )


def test_aromatic_classification_requires_ring_geometry():
    protein = [feature("R1", "protein", (0.0, 0.0, 0.0))]
    ligand = [feature("R2", "ligand", (0.4, 0.0, 0.0))]
    rings = (
        AromaticRing(
            "PHE1",
            "protein",
            ("R1",),
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            "test fixture",
        ),
        AromaticRing(
            "LIG1",
            "ligand",
            ("R2",),
            (0.4, 0.0, 0.0),
            (0.0, 0.0, 1.0),
            "test fixture",
        ),
    )
    result = interpret(protein, ligand, rings=rings)

    assert "classified" in status(result, "aromatic")
    aromatic = next(item for item in result.interactions if item.interaction_type == "aromatic")
    assert aromatic.distance_nm == pytest.approx(0.4)
    assert aromatic.angle_degrees == pytest.approx(0.0)


def test_water_mediated_classification_retains_the_water_participant():
    protein = [feature("O1", "protein", (0.0, 0.0, 0.0), element="O", acceptor=True)]
    ligand = [feature("O2", "ligand", (0.4, 0.0, 0.0), element="O", acceptor=True)]
    water = [
        feature(
            "W1",
            "water",
            (0.2, 0.0, 0.0),
            element="O",
            donor=True,
            hydrogens=((0.1, 0.0, 0.0),),
        )
    ]
    result = interpret(protein, ligand, water=water)

    assert "classified" in status(result, "water_mediated")
    record = next(item for item in result.interactions if item.interaction_type == "water_mediated")
    assert record.water_id == "W1"
    assert "W1" in record.participants


def test_metal_coordination_does_not_claim_geometry_it_lacks():
    protein = [feature("O1", "protein", (0.0, 0.0, 0.0), element="O")]
    ligand = [feature("O2", "ligand", (0.4, 0.0, 0.0), element="O")]
    metal = [feature("ZN", "metal", (0.2, 0.0, 0.0), element="Zn", metal=True)]
    result = interpret(protein, ligand, metal=metal)

    assert status(result, "metal_coordination") == ["unknown", "unknown"]
    records = [
        item for item in result.interactions if item.interaction_type == "metal_coordination"
    ]
    assert all(record.metal_id == "ZN" for record in records)
    assert all("coordination geometry" in record.limitations for record in records)


def test_record_round_trip_and_validation(tmp_path):
    protein = [feature("N1", "protein", (0.0, 0.0, 0.0), element="N", donor=True)]
    ligand = [feature("O1", "ligand", (0.3, 0.0, 0.0), element="O", acceptor=True)]
    result = interpret(protein, ligand)
    path = write_interactions(tmp_path / "interactions.json", result)
    loaded = load_interactions(path)

    assert loaded == result
    assert validate_interactions(loaded) == []
    assert loaded.schema_version == INTERACTION_SCHEMA_VERSION


def test_validation_rejects_bad_criteria():
    protein = [feature("C1", "protein", (0.0, 0.0, 0.0))]
    ligand = [feature("C2", "ligand", (0.3, 0.0, 0.0))]
    result = interpret(
        protein,
        ligand,
        criteria=InteractionCriteria(
            h_bond_max_distance_nm=-1.0,
            h_bond_min_angle_degrees=200.0,
        ),
    )

    errors = validate_interactions(result)

    assert any("h_bond_max_distance_nm" in error for error in errors)
    assert any("h_bond_min_angle_degrees" in error for error in errors)


def test_reference_2hyy_interaction_record_is_conservative():
    path = Path(__file__).resolve().parents[1] / "data" / "structures" / "2HYY.interactions.json"
    result = load_interactions(path)

    assert validate_interactions(result) == []
    assert result.accession == "2HYY"
    assert result.input_contact_count == 112
    assert all(item.status != "classified" for item in result.interactions)
    assert status(result, "hydrogen_bond") == ["unknown"]
    assert status(result, "ionic") == ["unknown"]
    assert status(result, "hydrophobic") == ["unknown"]
    assert status(result, "aromatic") == ["unknown"]
    assert status(result, "water_mediated") == ["unknown"]
    assert status(result, "metal_coordination") == ["not_applicable"]
