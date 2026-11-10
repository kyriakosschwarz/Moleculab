"""Tests for ligand identity reconciliation."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from moleculab import (
    LIGAND_IDENTITY_SCHEMA_VERSION,
    AtomMapping,
    DepositedAtom,
    DepositedBond,
    DepositedLigand,
    IdentityDiscrepancy,
    LigandSelection,
    ReferenceAtom,
    ReferenceBond,
    ReferenceLigand,
    ValidationCheck,
    ValidationRecord,
    load_ligand_identity,
    reconcile_ligand_identity,
    validate_ligand_identity,
    write_ligand_identity,
)

CHECKSUM = "b" * 64


def selection() -> LigandSelection:
    return LigandSelection(
        component_id="LIG",
        label_asym_id="A",
        auth_asym_id="A",
        instance="A",
        auth_seq_id="1",
        available_instances=("A",),
    )


def deposited_ligand(
    *,
    formula: str = "C2 H6 O",
    atoms: tuple[DepositedAtom, ...] | None = None,
    bonds: tuple[DepositedBond, ...] | None = None,
) -> DepositedLigand:
    if atoms is None:
        atoms = (
            DepositedAtom("A1", "C"),
            DepositedAtom("A2", "C"),
            DepositedAtom("A3", "O"),
        )
    if bonds is None:
        bonds = (
            DepositedBond("A1", "A2", "SINGLE"),
            DepositedBond("A2", "A3", "SINGLE"),
        )
    return DepositedLigand(
        component_id="LIG",
        label_asym_id="A",
        auth_asym_id="A",
        auth_seq_id="1",
        label_seq_id=None,
        insertion_code=None,
        formula=formula,
        atoms=atoms,
        bonds=bonds,
        atoms_are_heavy=True,
        formal_charges_present=False,
        hydrogen_atoms_present=False,
        omitted_bond_records=0,
        omitted_bond_atom_elements=(),
        stereochemistry=(),
    )


def reference_ligand(
    *,
    formula: str = "C2H6O",
    atoms: tuple[ReferenceAtom, ...] | None = None,
    bonds: tuple[ReferenceBond, ...] | None = None,
    stereocenters: tuple[str, ...] = (),
) -> ReferenceLigand:
    if atoms is None:
        atoms = (
            ReferenceAtom(0, "C", 0),
            ReferenceAtom(1, "C", 0),
            ReferenceAtom(2, "O", 0),
        )
    if bonds is None:
        bonds = (
            ReferenceBond(0, 1, "SINGLE"),
            ReferenceBond(1, 2, "SINGLE"),
        )
    return ReferenceLigand(
        name="reference",
        source="test reference",
        smiles="CCO",
        formula=formula,
        formal_charge=0,
        hydrogen_count=6,
        atoms=atoms,
        bonds=bonds,
        atoms_are_heavy=True,
        stereocenters=stereocenters,
    )


def record():
    return reconcile_ligand_identity(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        deposited=deposited_ligand(),
        reference=reference_ligand(),
        candidate_mappings=({"A1": 0, "A2": 1, "A3": 2},),
    )


def check_status(record_value, check_id: str) -> str:
    return next(
        check.status
        for check in record_value.validation.checks
        if check.check_id == check_id
    )


def test_complete_identity_record_is_structurally_valid():
    result = record()

    assert validate_ligand_identity(result) == []
    assert result.schema_version == LIGAND_IDENTITY_SCHEMA_VERSION
    assert result.mapping.candidate_count == 1
    assert result.mapping.ambiguous is False
    assert len(result.mapping.pairs) == 3
    assert check_status(result, "chemical.formula_matches_reference") == "pass"
    assert check_status(result, "chemical.connectivity_matches") == "pass"
    assert check_status(result, "chemical.formal_charge_matches") == "not_run"
    assert check_status(result, "chemical.hydrogens_match") == "not_run"
    assert any(item.category == "formal_charge" and item.status == "open" for item in result.discrepancies)
    assert any(item.category == "hydrogen_atoms" and item.status == "not_applicable" for item in result.discrepancies)


def test_reference_2hyy_identity_record_is_reproducible():
    path = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "structures"
        / "2HYY.ligand-identity.json"
    )
    result = load_ligand_identity(path)

    assert validate_ligand_identity(result) == []
    assert result.accession == "2HYY"
    assert result.selection.component_id == "STI"
    assert result.selection.label_asym_id == "E"
    assert len(result.deposited.atoms) == 37
    assert len(result.deposited.bonds) == 41
    assert result.deposited.omitted_bond_records == 31
    assert result.deposited.omitted_bond_atom_elements == ("H",)
    assert result.mapping.ordering_basis == "deposited.atoms order"
    assert result.normalization.toolkit == "RDKit"
    assert result.normalization.toolkit_version is not None
    assert result.reference.formula == "C29H31N7O"
    assert result.mapping.candidate_count == 4
    assert result.mapping.ambiguous is True
    assert check_status(result, "chemical.formula_matches_reference") == "pass"
    assert check_status(result, "chemical.bond_orders_match") == "pass"
    assert check_status(result, "chemical.formal_charge_matches") == "not_run"
    assert any(
        item.category == "bond_order_representation" and item.status == "resolved"
        for item in result.discrepancies
    )


def test_identity_record_json_round_trip(tmp_path):
    result = record()
    path = write_ligand_identity(tmp_path / "identity.json", result)
    loaded = load_ligand_identity(path)

    assert loaded == result
    assert loaded.to_dict()["mapping"]["candidate_count"] == 1


def test_formula_mismatch_is_reported_not_repaired():
    result = reconcile_ligand_identity(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        deposited=deposited_ligand(formula="C2 H4 O"),
        reference=reference_ligand(),
        candidate_mappings=({"A1": 0, "A2": 1, "A3": 2},),
    )

    assert check_status(result, "chemical.formula_matches_reference") == "fail"
    discrepancy = next(item for item in result.discrepancies if item.category == "formula")
    assert discrepancy.status == "open"
    assert discrepancy.observed == "C2 H4 O"
    assert discrepancy.expected == "C2H6O"
    assert result.deposited.formula == "C2 H4 O"


def test_element_mismatch_is_reported_and_mapping_elements_fail():
    reference = reference_ligand(
        atoms=(
            ReferenceAtom(0, "C", 0),
            ReferenceAtom(1, "C", 0),
            ReferenceAtom(2, "N", 0),
        )
    )
    result = reconcile_ligand_identity(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        deposited=deposited_ligand(),
        reference=reference,
        candidate_mappings=({"A1": 0, "A2": 1, "A3": 2},),
    )

    assert check_status(result, "chemical.element_counts_match") == "fail"
    assert check_status(result, "chemical.atom_mapping_elements_agree") == "fail"
    assert any(item.category == "atom_composition" and item.status == "open" for item in result.discrepancies)


def test_ambiguous_mapping_preserves_a_representative_and_the_rule():
    atoms = (DepositedAtom("A1", "C"), DepositedAtom("A2", "C"))
    bonds = (DepositedBond("A1", "A2", "SINGLE"),)
    reference = ReferenceLigand(
        name="C2",
        source="test reference",
        smiles="CC",
        formula="C2",
        formal_charge=0,
        hydrogen_count=0,
        atoms=(ReferenceAtom(0, "C", 0), ReferenceAtom(1, "C", 0)),
        bonds=(ReferenceBond(0, 1, "SINGLE"),),
        atoms_are_heavy=True,
    )
    result = reconcile_ligand_identity(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        deposited=deposited_ligand(formula="C2", atoms=atoms, bonds=bonds),
        reference=reference,
        candidate_mappings=(
            {"A1": 0, "A2": 1},
            {"A1": 1, "A2": 0},
        ),
    )

    assert result.mapping.candidate_count == 2
    assert result.mapping.ambiguous is True
    assert result.mapping.pairs == (("A1", 0), ("A2", 1))
    assert "lexicographically smallest" in result.mapping.selection_rule
    assert check_status(result, "chemical.atom_mapping_ambiguity_documented") == "pass"


def test_aromatic_normalization_records_the_raw_bond_order_difference():
    atoms = (DepositedAtom("A1", "C"), DepositedAtom("A2", "C"))
    deposited_bonds = (DepositedBond("A1", "A2", "SINGLE", aromatic=True),)
    reference = ReferenceLigand(
        name="C2",
        source="test reference",
        smiles="CC",
        formula="C2",
        formal_charge=0,
        hydrogen_count=0,
        atoms=(ReferenceAtom(0, "C", 0), ReferenceAtom(1, "C", 0)),
        bonds=(ReferenceBond(0, 1, "DOUBLE", aromatic=True),),
        atoms_are_heavy=True,
    )
    result = reconcile_ligand_identity(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        deposited=deposited_ligand(formula="C2", atoms=atoms, bonds=deposited_bonds),
        reference=reference,
        candidate_mappings=({"A1": 0, "A2": 1},),
    )

    assert check_status(result, "chemical.bond_orders_match") == "pass"
    discrepancy = next(
        item for item in result.discrepancies if item.category == "bond_order_representation"
    )
    assert discrepancy.status == "resolved"
    assert discrepancy.observed == "1 raw differences"
    assert discrepancy.resolution is not None
    assert result.deposited.bonds[0].order == "SINGLE"


def test_unresolved_bond_order_mismatch_fails():
    atoms = (DepositedAtom("A1", "C"), DepositedAtom("A2", "C"))
    reference = ReferenceLigand(
        name="C2",
        source="test reference",
        smiles="CC",
        formula="C2",
        formal_charge=0,
        hydrogen_count=0,
        atoms=(ReferenceAtom(0, "C", 0), ReferenceAtom(1, "C", 0)),
        bonds=(ReferenceBond(0, 1, "DOUBLE"),),
        atoms_are_heavy=True,
    )
    result = reconcile_ligand_identity(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        deposited=deposited_ligand(
            formula="C2",
            atoms=atoms,
            bonds=(DepositedBond("A1", "A2", "SINGLE"),),
        ),
        reference=reference,
        candidate_mappings=({"A1": 0, "A2": 1},),
    )

    assert check_status(result, "chemical.bond_orders_match") == "fail"
    discrepancy = next(
        item for item in result.discrepancies if item.category == "bond_order_representation"
    )
    assert discrepancy.status == "open"


def test_present_formal_charge_mismatch_is_open():
    deposited = replace(
        deposited_ligand(),
        atoms=(
            DepositedAtom("A1", "C", 0),
            DepositedAtom("A2", "C", 0),
            DepositedAtom("A3", "O", -1),
        ),
        formal_charges_present=True,
    )
    result = reconcile_ligand_identity(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        deposited=deposited,
        reference=reference_ligand(),
        candidate_mappings=({"A1": 0, "A2": 1, "A3": 2},),
    )

    assert check_status(result, "chemical.formal_charge_matches") == "fail"
    assert any(
        item.category == "formal_charge" and item.status == "open"
        for item in result.discrepancies
    )


def test_stereochemistry_mismatch_fails():
    deposited = replace(deposited_ligand(), stereochemistry=("1:S",))
    result = reconcile_ligand_identity(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        deposited=deposited,
        reference=reference_ligand(stereocenters=("1:R",)),
        candidate_mappings=({"A1": 0, "A2": 1, "A3": 2},),
    )

    assert check_status(result, "chemical.stereochemistry_matches") == "fail"
    assert any(
        item.category == "stereochemistry" and item.status == "open"
        for item in result.discrepancies
    )


def test_missing_mapping_and_invalid_candidate_are_reported():
    missing = reconcile_ligand_identity(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        deposited=deposited_ligand(),
        reference=reference_ligand(),
        candidate_mappings=(),
    )
    invalid = reconcile_ligand_identity(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        deposited=deposited_ligand(),
        reference=reference_ligand(),
        candidate_mappings=({"A1": 0},),
    )

    for result in (missing, invalid):
        assert result.mapping.candidate_count == 0
        assert result.mapping.pairs == ()
        assert check_status(result, "chemical.atom_mapping_complete") == "fail"
        assert any(item.category == "atom_mapping" and item.status == "open" for item in result.discrepancies)
        assert validate_ligand_identity(result) == []


def test_formal_charge_and_stereochemistry_are_not_inferred():
    deposited = replace(deposited_ligand(), stereochemistry=())
    reference = reference_ligand(stereocenters=("1:R",))
    result = reconcile_ligand_identity(
        source_manifest_path="manifest.json",
        source_manifest_sha256=CHECKSUM,
        accession="TEST",
        selection=selection(),
        deposited=deposited,
        reference=reference,
        candidate_mappings=({"A1": 0, "A2": 1, "A3": 2},),
    )

    assert check_status(result, "chemical.formal_charge_matches") == "not_run"
    assert check_status(result, "chemical.stereochemistry_matches") == "not_run"
    assert any(item.category == "formal_charge" and item.status == "open" for item in result.discrepancies)
    assert any(item.category == "stereochemistry" and item.status == "open" for item in result.discrepancies)


def test_validation_rejects_bad_schema_and_incomplete_mapping_records():
    result = record()
    bad_mapping = AtomMapping(
        pairs=(("A1", 0),),
        candidate_count=2,
        ambiguous=False,
        selection_rule="bad",
    )
    bad_discrepancy = IdentityDiscrepancy("x", "resolved", "no resolution")
    duplicate_check = result.validation.checks[0]
    bad = replace(
        result,
        source_manifest_sha256="not-a-checksum",
        mapping=bad_mapping,
        discrepancies=(bad_discrepancy,),
        validation=ValidationRecord(checks=(duplicate_check, duplicate_check)),
    )

    errors = validate_ligand_identity(bad)

    assert any("source_manifest_sha256" in error for error in errors)
    assert any("mapping.ambiguous" in error for error in errors)
    assert any("mapping must cover" in error for error in errors)
    assert any("resolved but has no resolution" in error for error in errors)
    assert any("check_id is duplicated" in error for error in errors)


def test_validation_rejects_unknown_bond_atoms():
    result = record()
    bad = replace(
        result,
        deposited=replace(
            result.deposited,
            bonds=(DepositedBond("A1", "MISSING", "SINGLE"),),
        ),
    )

    errors = validate_ligand_identity(bad)

    assert any("unknown atom name" in error for error in errors)


def test_validation_rejects_invalid_check_status():
    result = record()
    bad_check = ValidationCheck(
        check_id="x",
        category="x",
        status="assumed",
        subject="x",
        criterion="x",
        observed=None,
        expected=None,
        detail="x",
    )
    bad = replace(
        result,
        validation=ValidationRecord(checks=(bad_check,)),
    )

    errors = validate_ligand_identity(bad)

    assert any("invalid status" in error for error in errors)
