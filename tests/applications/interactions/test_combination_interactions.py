# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for evidence-complete ligand/protein state combinations."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from moleculab import (
    COMBINATION_SCHEMA_VERSION,
    AtomFeature,
    ChemicalState,
    CoordinateAtom,
    InteractionCriteria,
    LigandSelection,
    PreparedAtom,
    PreparedProteinAtom,
    PreparedProteinState,
    PreparedWaterRole,
    ProteinResidueState,
    ProteinSelection,
    assess_combination,
    build_combination_record,
    load_combination_record,
    validate_combination_record,
    write_combination_record,
)

CHECKSUM = "4" * 64


def ligand_state(*, charge: int = 0) -> ChemicalState:
    atoms = tuple(
        PreparedAtom(
            atom_id=f"atom_{index}",
            element="C",
            formal_charge=charge,
            is_donor=True,
            is_acceptor=True,
            is_hydrophobic=True,
            aromatic_ring_ids=(),
            hydrogen_count=0,
            hydrogen_positions_nm=(),
            coordinate_source="generated",
            hydrogen_source=None,
            feature_source="test",
        )
        for index in range(2)
    )
    return ChemicalState(
        state_id="ligand-state",
        kind="reference",
        label="ligand state",
        parent_state_id=None,
        molecular_graph="CC",
        formula="C2H6",
        proton_count=31,
        formal_charge=charge,
        tautomer_of=None,
        ph_assumption=None,
        pka_sites=(),
        enumeration_method="test",
        enumeration_parameters=("test",),
        assumptions=("test",),
        atoms=atoms,
        reference_atom_mapping=(("atom_0", 0), ("atom_1", 1)),
        rings=(),
        water_roles=(PreparedWaterRole("W1", "unassigned", "test"),),
    )


def protein_state(*, charge: int = 0) -> PreparedProteinState:
    atoms = tuple(
        PreparedProteinAtom(
            atom_id=f"A:1:{name}",
            atom_name=name,
            element="C",
            residue_name="ALA",
            label_asym_id="A",
            auth_asym_id="A",
            label_seq_id="1",
            auth_seq_id="1",
            insertion_code=None,
            alternate_location=None,
            position_angstrom=position,
            formal_charge=charge,
            is_donor=True,
            is_acceptor=True,
            is_hydrophobic=True,
            aromatic_ring_ids=(),
            hydrogen_count=0,
            hydrogen_positions_nm=(),
            coordinate_source="deposited",
            hydrogen_source=None,
            feature_source="test template",
        )
        for name, position in (("N", (0.3, 0.0, 0.0)), ("CA", (1.0, 0.0, 0.0)))
    )
    return PreparedProteinState(
        state_id="protein-state",
        label="protein state",
        protonation_model="test",
        hydrogen_method=None,
        assumptions=("test",),
        atoms=atoms,
        residue_states=(
            ProteinResidueState(
                residue_id="A:1",
                residue_name="ALA",
                label_asym_id="A",
                auth_asym_id="A",
                label_seq_id="1",
                auth_seq_id="1",
                insertion_code=None,
                protonation_state="standard",
                formal_charge=charge,
                source="test",
            ),
        ),
        rings=(),
        water_roles=(),
        hetero_roles=(),
    )


def ligand_atoms() -> tuple[CoordinateAtom, ...]:
    return (
        CoordinateAtom(
            atom_name="C1",
            element="C",
            reference_index=0,
            position_angstrom=(0.0, 0.0, 0.0),
            occupancy=1.0,
            b_iso_or_equiv=20.0,
            alternate_location=None,
            model_number=1,
        ),
        CoordinateAtom(
            atom_name="C2",
            element="C",
            reference_index=1,
            position_angstrom=(1.4, 0.0, 0.0),
            occupancy=1.0,
            b_iso_or_equiv=20.0,
            alternate_location=None,
            model_number=1,
        ),
    )


def water_features():
    return (
        AtomFeature(
            atom_id="water:1:O",
            side="water",
            element="O",
            position_nm=(0.2, 0.0, 0.0),
            residue_id="HOH A 1",
            is_donor=None,
            is_acceptor=None,
            formal_charge=None,
            is_hydrophobic=None,
            aromatic_ring_id=None,
            is_metal=None,
            hydrogen_positions_nm=(),
            feature_source=None,
        ),
    )


def assess(ligand_charge: int = 0, protein_charge: int = 0):
    return assess_combination(
        ligand_state=ligand_state(charge=ligand_charge),
        protein_state=protein_state(charge=protein_charge),
        reference_to_deposited={0: "C1", 1: "C2"},
        ligand_atoms=ligand_atoms(),
        protein_atoms=protein_state(charge=protein_charge).atoms,
        contact_ligand_atoms={"C1", "C2"},
        contact_protein_atoms={("A", "1", "N"), ("A", "1", "CA")},
        water_features=water_features(),
        criteria=InteractionCriteria(),
    )


def test_combination_classifies_supported_categories_and_leaves_hydrogen_bonds_unknown():
    result = assess()

    statuses = dict(result.category_statuses)
    assert statuses["hydrophobic"] == "classified"
    assert statuses["hydrogen_bond"] == "unknown"
    assert statuses["water_mediated"] == "unknown"
    assert statuses["metal_coordination"] == "not_applicable"
    assert not any(
        item.interaction_type == "hydrogen_bond"
        for item in result.classified_interactions
    )


def test_opposite_charges_produce_ionic_classification():
    result = assess(ligand_charge=1, protein_charge=-1)
    statuses = dict(result.category_statuses)

    assert statuses["ionic"] == "classified"
    assert any(
        item.interaction_type == "ionic" for item in result.classified_interactions
    )


def test_classification_does_not_imply_ranking_or_affinity():
    result = assess()
    serialized = result.to_dict()

    assert "affinity" not in serialized
    assert "score" not in serialized
    assert "rank" not in serialized


def build_record():
    result = assess()
    return build_combination_record(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_identity_record_path="identity.json",
        ligand_identity_record_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        binding_site_record_path="site.json",
        binding_site_record_sha256=CHECKSUM,
        chemical_state_record_path="states.json",
        chemical_state_record_sha256=CHECKSUM,
        protein_preparation_record_path="protein.json",
        protein_preparation_record_sha256=CHECKSUM,
        accession="TEST",
        protein_selection=ProteinSelection(label_asym_id="A", auth_asym_id="A"),
        ligand_selection=LigandSelection(
            component_id="LIG",
            label_asym_id="E",
            auth_asym_id="A",
            instance="E",
            auth_seq_id="600",
            available_instances=("E",),
        ),
        criteria=InteractionCriteria(),
        ligand_state_ids=("ligand-state",),
        protein_state_ids=("protein-state",),
        combinations=(result,),
    )


def test_record_round_trip_and_validation(tmp_path):
    result = build_record()
    path = write_combination_record(tmp_path / "combinations.json", result)
    loaded = load_combination_record(path)

    assert loaded == result
    assert validate_combination_record(loaded) == []
    assert loaded.schema_version == COMBINATION_SCHEMA_VERSION


def test_matrix_must_be_complete():
    result = build_record()
    bad = replace(result, protein_state_ids=("missing",))
    with pytest.raises(ValueError):
        build_combination_record(
            source_manifest_path=bad.source_manifest_path,
            source_manifest_sha256=bad.source_manifest_sha256,
            ligand_identity_record_path=bad.ligand_identity_record_path,
            ligand_identity_record_sha256=bad.ligand_identity_record_sha256,
            ligand_coordinate_record_path=bad.ligand_coordinate_record_path,
            ligand_coordinate_record_sha256=bad.ligand_coordinate_record_sha256,
            binding_site_record_path=bad.binding_site_record_path,
            binding_site_record_sha256=bad.binding_site_record_sha256,
            chemical_state_record_path=bad.chemical_state_record_path,
            chemical_state_record_sha256=bad.chemical_state_record_sha256,
            protein_preparation_record_path=bad.protein_preparation_record_path,
            protein_preparation_record_sha256=bad.protein_preparation_record_sha256,
            accession=bad.accession,
            protein_selection=bad.protein_selection,
            ligand_selection=bad.ligand_selection,
            criteria=bad.criteria,
            ligand_state_ids=bad.ligand_state_ids,
            protein_state_ids=bad.protein_state_ids,
            combinations=bad.combinations,
        )


def test_reference_2hyy_combination_record_is_conservative():
    path = (
        Path(__file__).resolve().parents[3]
        / "data"
        / "structures"
        / "2HYY.combination-interactions.json"
    )
    result = load_combination_record(path)

    assert validate_combination_record(result) == []
    assert result.accession == "2HYY"
    assert len(result.ligand_state_ids) == 8
    assert len(result.protein_state_ids) == 2
    assert len(result.combinations) == 16
    assert all(
        not any(
            item.interaction_type == "hydrogen_bond"
            for item in combination.classified_interactions
        )
        for combination in result.combinations
    )
    assert all(
        dict(combination.category_statuses)["water_mediated"] == "unknown"
        for combination in result.combinations
    )
