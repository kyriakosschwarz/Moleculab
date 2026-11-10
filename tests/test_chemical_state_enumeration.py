"""Tests for protonation and tautomer state enumeration records."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from moleculab import (
    CHEMICAL_STATE_SCHEMA_VERSION,
    ChemicalState,
    IonizableSite,
    LigandSelection,
    PreparedAtom,
    build_chemical_state_enumeration_record,
    load_chemical_state_enumeration,
    validate_chemical_state_enumeration,
    write_chemical_state_enumeration,
)

CHECKSUM = "1" * 64


def selection() -> LigandSelection:
    return LigandSelection(
        component_id="LIG",
        label_asym_id="E",
        auth_asym_id="A",
        instance="E",
        auth_seq_id="600",
        available_instances=("E",),
    )


def atom() -> PreparedAtom:
    return PreparedAtom(
        atom_id="C1",
        element="C",
        formal_charge=0,
        is_donor=None,
        is_acceptor=None,
        is_hydrophobic=None,
        aromatic_ring_ids=(),
        hydrogen_count=0,
        hydrogen_positions_nm=(),
        coordinate_source="deposited",
        hydrogen_source=None,
        feature_source=None,
    )


def state(
    state_id: str,
    *,
    kind: str = "reference",
    parent: str | None = None,
    protons: int = 31,
    charge: int = 0,
    tautomer_of: str | None = None,
    pka_sites: tuple[IonizableSite, ...] = (),
) -> ChemicalState:
    return ChemicalState(
        state_id=state_id,
        kind=kind,
        label=state_id,
        parent_state_id=parent,
        molecular_graph="CC",
        formula="C2H6",
        proton_count=protons,
        formal_charge=charge,
        tautomer_of=tautomer_of,
        ph_assumption=None,
        pka_sites=pka_sites,
        enumeration_method="test",
        enumeration_parameters=("test",),
        assumptions=("test assumption",),
        atoms=(atom(),),
        rings=(),
        water_roles=(),
    )


def record():
    return build_chemical_state_enumeration_record(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        ligand_identity_record_path="identity.json",
        ligand_identity_record_sha256=CHECKSUM,
        ligand_coordinate_record_path="coordinates.json",
        ligand_coordinate_record_sha256=CHECKSUM,
        binding_site_record_path="site.json",
        binding_site_record_sha256=CHECKSUM,
        chemical_preparation_record_path="preparation.json",
        chemical_preparation_record_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        reference_state_id="reference",
        states=(
            state("reference"),
            state(
                "tautomer-1",
                kind="tautomer",
                parent="reference",
                tautomer_of="reference",
            ),
            state(
                "protonation-1",
                kind="protonation",
                parent="reference",
                protons=32,
                charge=1,
                pka_sites=(
                    IonizableSite(
                        "site-1",
                        "atom-1",
                        "base",
                        pka=None,
                        pka_source="not supplied",
                        pka_uncertainty=None,
                    ),
                ),
            ),
        ),
        toolkit="RDKit",
        toolkit_version="test",
        parameters=("test parameters",),
        selected_state_id=None,
        population_model=None,
    )


def test_complete_state_record_is_valid():
    result = record()

    assert validate_chemical_state_enumeration(result) == []
    assert result.schema_version == CHEMICAL_STATE_SCHEMA_VERSION
    assert result.selected_state_id is None
    assert result.population_model is None
    assert {state.kind for state in result.states} == {
        "reference",
        "tautomer",
        "protonation",
    }


def test_json_round_trip(tmp_path):
    result = record()
    path = write_chemical_state_enumeration(tmp_path / "states.json", result)
    loaded = load_chemical_state_enumeration(path)

    assert loaded == result
    assert loaded.to_dict()["states"][0]["kind"] == "reference"


def test_population_model_is_rejected():
    with pytest.raises(ValueError, match="population"):
        build_chemical_state_enumeration_record(
            source_manifest_path="manifest.json",
            source_manifest_sha256=CHECKSUM,
            ligand_identity_record_path="identity.json",
            ligand_identity_record_sha256=CHECKSUM,
            ligand_coordinate_record_path="coordinates.json",
            ligand_coordinate_record_sha256=CHECKSUM,
            binding_site_record_path="site.json",
            binding_site_record_sha256=CHECKSUM,
            chemical_preparation_record_path="preparation.json",
            chemical_preparation_record_sha256=CHECKSUM,
            accession="TEST",
            selection=selection(),
            reference_state_id="reference",
            states=(state("reference"),),
            toolkit="RDKit",
            toolkit_version="test",
            parameters=("test",),
            population_model="should not be here",
        )


def test_tautomer_must_match_parent_composition():
    result = record()
    bad = replace(
        result,
        states=(
            result.states[0],
            replace(result.states[1], proton_count=30),
            result.states[2],
        ),
    )
    errors = validate_chemical_state_enumeration(bad)

    assert any("tautomer proton count must match" in error for error in errors)


def test_protonation_must_change_proton_count():
    result = record()
    bad = replace(
        result,
        states=(
            result.states[0],
            result.states[1],
            replace(result.states[2], proton_count=31),
        ),
    )
    errors = validate_chemical_state_enumeration(bad)

    assert any("must change proton count" in error for error in errors)


def test_pka_value_requires_source():
    result = record()
    bad_site = IonizableSite(
        "site-1",
        "atom-1",
        "base",
        pka=8.0,
        pka_source=None,
        pka_uncertainty=0.5,
    )
    bad = replace(
        result,
        states=(
            result.states[0],
            result.states[1],
            replace(result.states[2], pka_sites=(bad_site,)),
        ),
    )
    errors = validate_chemical_state_enumeration(bad)

    assert any("requires a pKa source" in error for error in errors)


def test_missing_state_reference_is_rejected():
    result = record()
    bad = replace(result, selected_state_id="missing")
    errors = validate_chemical_state_enumeration(bad)

    assert any("selected_state_id" in error for error in errors)


def test_reference_2hyy_state_record_is_reproducible():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.chemical-states.json"
    )
    result = load_chemical_state_enumeration(path)

    assert validate_chemical_state_enumeration(result) == []
    assert result.accession == "2HYY"
    assert len(result.states) == 8
    assert {state.kind for state in result.states} == {
        "reference",
        "tautomer",
        "protonation",
    }
    assert result.selected_state_id is None
    assert result.population_model is None
    reference = next(state for state in result.states if state.kind == "reference")
    assert reference.proton_count == 31
    assert reference.formal_charge == 0
    for state in result.states:
        if state.kind == "tautomer":
            assert state.proton_count == 31
            assert state.tautomer_of == "reference-neutral"
        if state.kind == "protonation":
            assert state.proton_count == 32
            assert state.formal_charge == 1
            assert state.pka_sites
            assert state.pka_sites[0].pka is None
            assert "not supplied" in state.pka_sites[0].pka_source
