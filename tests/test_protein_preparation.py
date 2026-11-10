"""Tests for protein-pocket chemical preparation records."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from moleculab import (
    PROTEIN_PREPARATION_SCHEMA_VERSION,
    PreparedProteinAtom,
    PreparedProteinRing,
    PreparedProteinState,
    PreparedWaterRole,
    ProteinHeteroRole,
    ProteinPreparationRecord,
    ProteinResidueState,
    ProteinSelection,
    build_protein_preparation_record,
    load_protein_preparation,
    validate_protein_preparation,
    write_protein_preparation,
)

CHECKSUM = "3" * 64


def selection() -> ProteinSelection:
    return ProteinSelection(label_asym_id="A", auth_asym_id="A")


def atom(
    atom_id: str = "A:1:N",
    *,
    donor: bool | None = True,
    acceptor: bool | None = False,
    hydrophobic: bool | None = False,
    charge: int | None = 0,
    hydrogen_count: int = 0,
    hydrogen_source: str | None = None,
) -> PreparedProteinAtom:
    return PreparedProteinAtom(
        atom_id=atom_id,
        atom_name=atom_id.rsplit(":", 1)[-1],
        element="N",
        residue_name="ALA",
        label_asym_id="A",
        auth_asym_id="A",
        label_seq_id="1",
        auth_seq_id="1",
        insertion_code=None,
        alternate_location=None,
        formal_charge=charge,
        is_donor=donor,
        is_acceptor=acceptor,
        is_hydrophobic=hydrophobic,
        aromatic_ring_ids=(),
        hydrogen_count=hydrogen_count,
        hydrogen_positions_nm=tuple(
            (0.0, 0.0, 0.0) for _ in range(hydrogen_count)
        ),
        coordinate_source="deposited",
        hydrogen_source=hydrogen_source,
        feature_source="test template",
    )


def state(state_id: str, *, donor: bool = True, charge: int = 0) -> PreparedProteinState:
    return PreparedProteinState(
        state_id=state_id,
        label=state_id,
        protonation_model="test",
        hydrogen_method=None,
        assumptions=("test assumption",),
        atoms=(atom(donor=donor, charge=charge),),
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
                source="test template",
            ),
        ),
        rings=(),
        water_roles=(
            PreparedWaterRole("A:101:O", "unassigned", "experimental water"),
        ),
        hetero_roles=(
            ProteinHeteroRole("A:700:SO4", "SO4", "unknown", "preserved"),
        ),
    )


def record() -> ProteinPreparationRecord:
    return build_protein_preparation_record(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        binding_site_record_path="site.json",
        binding_site_record_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        scope="test pocket residues",
        toolkit="Moleculab residue templates",
        toolkit_version="test",
        parameters=("test parameters",),
        states=(
            state("standard"),
            state("histidine-protonated", donor=True, charge=1),
        ),
        selected_state_id=None,
    )


def test_complete_protein_preparation_record_is_valid():
    result = record()

    assert validate_protein_preparation(result) == []
    assert result.schema_version == PROTEIN_PREPARATION_SCHEMA_VERSION
    assert len(result.states) == 2
    assert result.selected_state_id is None
    assert any(item.category == "state_selection" for item in result.discrepancies)
    assert any(item.category == "hydrogen_coordinates" for item in result.discrepancies)


def test_json_round_trip(tmp_path):
    result = record()
    path = write_protein_preparation(tmp_path / "protein.json", result)
    loaded = load_protein_preparation(path)

    assert loaded == result
    assert loaded.to_dict()["states"][0]["state_id"] == "standard"


def test_selected_state_must_exist():
    result = record()
    bad = replace(result, selected_state_id="missing")
    errors = validate_protein_preparation(bad)

    assert any("selected_state_id" in error for error in errors)


def test_duplicate_state_ids_are_rejected():
    result = record()
    bad = replace(result, states=(state("same"), state("same")))
    errors = validate_protein_preparation(bad)

    assert any("state_id values must be unique" in error for error in errors)


def test_hydrogen_count_requires_source():
    result = record()
    bad_atom = atom(hydrogen_count=1, hydrogen_source=None)
    bad_state = replace(result.states[0], atoms=(bad_atom,))
    bad = replace(result, states=(bad_state, result.states[1]))
    errors = validate_protein_preparation(bad)

    assert any("generated hydrogens require a source" in error for error in errors)


def test_ring_must_reference_known_atom():
    result = record()
    bad_ring = PreparedProteinRing(
        "ring-1",
        "A:1",
        ("MISSING",),
        (0.0, 0.0, 0.0),
        (0.0, 0.0, 1.0),
        "test",
    )
    bad_state = replace(result.states[0], rings=(bad_ring,))
    bad = replace(result, states=(bad_state, result.states[1]))
    errors = validate_protein_preparation(bad)

    assert any("references an unknown atom" in error for error in errors)


def test_water_and_hetero_roles_must_be_valid():
    result = record()
    bad_state = replace(
        result.states[0],
        water_roles=(PreparedWaterRole("W", "assumed", "test"),),
        hetero_roles=(ProteinHeteroRole("H", "SO4", "assumed", "test"),),
    )
    bad = replace(result, states=(bad_state, result.states[1]))
    errors = validate_protein_preparation(bad)

    assert any("water_roles" in error and "invalid" in error for error in errors)
    assert any("hetero_roles" in error and "invalid" in error for error in errors)


def test_reference_2hyy_protein_preparation_is_reproducible():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.protein-preparation.json"
    )
    result = load_protein_preparation(path)

    assert validate_protein_preparation(result) == []
    assert result.accession == "2HYY"
    assert len(result.states) == 2
    assert result.selected_state_id is None
    standard = next(state for state in result.states if state.state_id == "standard-templates")
    histidine = next(state for state in result.states if state.state_id == "histidine-protonated")
    assert len(standard.atoms) == 186
    assert len(standard.residue_states) == 23
    assert len(standard.rings) == 4
    assert sum(1 for atom in standard.atoms if atom.is_donor) == 25
    assert sum(1 for atom in standard.atoms if atom.is_acceptor) == 25
    assert sum(1 for atom in standard.atoms if atom.is_hydrophobic) == 61
    assert sum(1 for atom in standard.atoms if atom.formal_charge) == 4
    assert sum(1 for atom in histidine.atoms if atom.is_donor) == 26
    assert sum(1 for atom in histidine.atoms if atom.formal_charge) == 5
    assert all(water.role == "unassigned" for water in standard.water_roles)
    assert result.toolkit == "Moleculab residue templates"
