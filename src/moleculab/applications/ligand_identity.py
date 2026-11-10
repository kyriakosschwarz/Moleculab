"""Ligand identity reconciliation for experimental structures.

This module compares a deposited ligand graph with a reference chemical graph.
It does not parse coordinate files, build reference molecules, or infer atom
mappings.  Callers supply those explicitly, usually through an mmCIF adapter
and an RDKit adapter.  The module then:

- selects one representative mapping from a set of candidate mappings;
- records mapping ambiguity instead of hiding it;
- compares formula, heavy-atom counts, elements, connectivity, and bond order;
- preserves the raw bond-order comparison before any aromatic normalization;
- records formal-charge and stereochemistry checks as ``not_run`` when the
  deposited record does not provide the information;
- reports discrepancies with an explicit status and resolution.

No coordinates, hydrogens, protonation states, formal charges, or bond orders
are modified by this module.
"""

from __future__ import annotations

import json
import re
from collections import Counter
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    LigandSelection,
    ProvenanceEntry,
    ProvenanceRecord,
    ValidationCheck,
    ValidationRecord,
)

LIGAND_IDENTITY_SCHEMA_VERSION = "1.0"
DISCREPANCY_STATUSES = ("open", "resolved", "not_applicable")

_SHA256_PATTERN = re.compile(r"[0-9a-f]{64}")


def _as_mapping(value: Any, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be a mapping.")
    return value


def _as_list(value: Any, name: str) -> list[Any]:
    if value is None:
        return []
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{name} must be a list.")
    return list(value)


def _required_string(data: Mapping[str, Any], key: str) -> str:
    value = data.get(key)
    if value is None or str(value) == "":
        raise ValueError(f"{key} is required.")
    return str(value)


def _optional_string(data: Mapping[str, Any], key: str) -> str | None:
    value = data.get(key)
    return None if value is None else str(value)


def _required_int(data: Mapping[str, Any], key: str) -> int:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return int(value)


def _optional_int(data: Mapping[str, Any], key: str) -> int | None:
    value = data.get(key)
    return None if value is None else int(value)


def _required_bool(data: Mapping[str, Any], key: str) -> bool:
    value = data.get(key)
    if not isinstance(value, bool):
        raise ValueError(f"{key} must be a boolean.")
    return value


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


@dataclass(frozen=True)
class DepositedAtom:
    atom_name: str
    element: str
    formal_charge: int | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_name": self.atom_name,
            "element": self.element,
            "formal_charge": self.formal_charge,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> DepositedAtom:
        data = _as_mapping(data, "deposited atom")
        return cls(
            atom_name=_required_string(data, "atom_name"),
            element=_required_string(data, "element"),
            formal_charge=_optional_int(data, "formal_charge"),
        )


@dataclass(frozen=True)
class DepositedBond:
    atom_name_1: str
    atom_name_2: str
    order: str
    aromatic: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_name_1": self.atom_name_1,
            "atom_name_2": self.atom_name_2,
            "order": self.order,
            "aromatic": self.aromatic,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> DepositedBond:
        data = _as_mapping(data, "deposited bond")
        return cls(
            atom_name_1=_required_string(data, "atom_name_1"),
            atom_name_2=_required_string(data, "atom_name_2"),
            order=_required_string(data, "order"),
            aromatic=_required_bool(data, "aromatic"),
        )


@dataclass(frozen=True)
class DepositedLigand:
    component_id: str
    label_asym_id: str
    auth_asym_id: str
    auth_seq_id: str
    label_seq_id: str | None
    insertion_code: str | None
    formula: str | None
    atoms: tuple[DepositedAtom, ...]
    bonds: tuple[DepositedBond, ...]
    atoms_are_heavy: bool
    formal_charges_present: bool
    hydrogen_atoms_present: bool
    omitted_hydrogen_bonds: int = 0
    stereochemistry: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "component_id": self.component_id,
            "label_asym_id": self.label_asym_id,
            "auth_asym_id": self.auth_asym_id,
            "auth_seq_id": self.auth_seq_id,
            "label_seq_id": self.label_seq_id,
            "insertion_code": self.insertion_code,
            "formula": self.formula,
            "atoms": [atom.to_dict() for atom in self.atoms],
            "bonds": [bond.to_dict() for bond in self.bonds],
            "atoms_are_heavy": self.atoms_are_heavy,
            "formal_charges_present": self.formal_charges_present,
            "hydrogen_atoms_present": self.hydrogen_atoms_present,
            "omitted_hydrogen_bonds": self.omitted_hydrogen_bonds,
            "stereochemistry": list(self.stereochemistry),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> DepositedLigand:
        data = _as_mapping(data, "deposited ligand")
        return cls(
            component_id=_required_string(data, "component_id"),
            label_asym_id=_required_string(data, "label_asym_id"),
            auth_asym_id=_required_string(data, "auth_asym_id"),
            auth_seq_id=_required_string(data, "auth_seq_id"),
            label_seq_id=_optional_string(data, "label_seq_id"),
            insertion_code=_optional_string(data, "insertion_code"),
            formula=_optional_string(data, "formula"),
            atoms=tuple(
                DepositedAtom.from_dict(_as_mapping(item, "deposited atom"))
                for item in _as_list(data.get("atoms"), "deposited.atoms")
            ),
            bonds=tuple(
                DepositedBond.from_dict(_as_mapping(item, "deposited bond"))
                for item in _as_list(data.get("bonds"), "deposited.bonds")
            ),
            atoms_are_heavy=_required_bool(data, "atoms_are_heavy"),
            formal_charges_present=_required_bool(data, "formal_charges_present"),
            hydrogen_atoms_present=_required_bool(data, "hydrogen_atoms_present"),
            omitted_hydrogen_bonds=_required_int(data, "omitted_hydrogen_bonds"),
            stereochemistry=_string_tuple(data, "stereochemistry"),
        )


@dataclass(frozen=True)
class ReferenceAtom:
    index: int
    element: str
    formal_charge: int

    def to_dict(self) -> dict[str, Any]:
        return {
            "index": self.index,
            "element": self.element,
            "formal_charge": self.formal_charge,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ReferenceAtom:
        data = _as_mapping(data, "reference atom")
        return cls(
            index=_required_int(data, "index"),
            element=_required_string(data, "element"),
            formal_charge=_required_int(data, "formal_charge"),
        )


@dataclass(frozen=True)
class ReferenceBond:
    atom_index_1: int
    atom_index_2: int
    order: str
    aromatic: bool = False

    def to_dict(self) -> dict[str, Any]:
        return {
            "atom_index_1": self.atom_index_1,
            "atom_index_2": self.atom_index_2,
            "order": self.order,
            "aromatic": self.aromatic,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ReferenceBond:
        data = _as_mapping(data, "reference bond")
        return cls(
            atom_index_1=_required_int(data, "atom_index_1"),
            atom_index_2=_required_int(data, "atom_index_2"),
            order=_required_string(data, "order"),
            aromatic=_required_bool(data, "aromatic"),
        )


@dataclass(frozen=True)
class ReferenceLigand:
    name: str
    source: str
    smiles: str
    formula: str
    formal_charge: int
    hydrogen_count: int
    atoms: tuple[ReferenceAtom, ...]
    bonds: tuple[ReferenceBond, ...]
    atoms_are_heavy: bool
    stereocenters: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "source": self.source,
            "smiles": self.smiles,
            "formula": self.formula,
            "formal_charge": self.formal_charge,
            "hydrogen_count": self.hydrogen_count,
            "atoms": [atom.to_dict() for atom in self.atoms],
            "bonds": [bond.to_dict() for bond in self.bonds],
            "atoms_are_heavy": self.atoms_are_heavy,
            "stereocenters": list(self.stereocenters),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ReferenceLigand:
        data = _as_mapping(data, "reference ligand")
        return cls(
            name=_required_string(data, "name"),
            source=_required_string(data, "source"),
            smiles=_required_string(data, "smiles"),
            formula=_required_string(data, "formula"),
            formal_charge=_required_int(data, "formal_charge"),
            hydrogen_count=_required_int(data, "hydrogen_count"),
            atoms=tuple(
                ReferenceAtom.from_dict(_as_mapping(item, "reference atom"))
                for item in _as_list(data.get("atoms"), "reference.atoms")
            ),
            bonds=tuple(
                ReferenceBond.from_dict(_as_mapping(item, "reference bond"))
                for item in _as_list(data.get("bonds"), "reference.bonds")
            ),
            atoms_are_heavy=_required_bool(data, "atoms_are_heavy"),
            stereocenters=_string_tuple(data, "stereocenters"),
        )


@dataclass(frozen=True)
class AtomMapping:
    pairs: tuple[tuple[str, int], ...]
    candidate_count: int
    ambiguous: bool
    selection_rule: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "pairs": [
                {"deposited_atom": name, "reference_index": index}
                for name, index in self.pairs
            ],
            "candidate_count": self.candidate_count,
            "ambiguous": self.ambiguous,
            "selection_rule": self.selection_rule,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> AtomMapping:
        data = _as_mapping(data, "atom mapping")
        pairs = []
        for item in _as_list(data.get("pairs"), "mapping.pairs"):
            item = _as_mapping(item, "mapping pair")
            pairs.append(
                (
                    _required_string(item, "deposited_atom"),
                    _required_int(item, "reference_index"),
                )
            )
        return cls(
            pairs=tuple(pairs),
            candidate_count=_required_int(data, "candidate_count"),
            ambiguous=_required_bool(data, "ambiguous"),
            selection_rule=_required_string(data, "selection_rule"),
        )


@dataclass(frozen=True)
class IdentityDiscrepancy:
    category: str
    status: str
    detail: str
    observed: str | None = None
    expected: str | None = None
    resolution: str | None = None

    def to_dict(self) -> dict[str, Any]:
        return {
            "category": self.category,
            "status": self.status,
            "detail": self.detail,
            "observed": self.observed,
            "expected": self.expected,
            "resolution": self.resolution,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> IdentityDiscrepancy:
        data = _as_mapping(data, "identity discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            resolution=_optional_string(data, "resolution"),
        )


@dataclass(frozen=True)
class LigandIdentityRecord:
    source_manifest_path: str
    source_manifest_sha256: str
    accession: str
    selection: LigandSelection
    deposited: DepositedLigand
    reference: ReferenceLigand
    mapping: AtomMapping
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[IdentityDiscrepancy, ...]
    schema_version: str = LIGAND_IDENTITY_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source_manifest_path": self.source_manifest_path,
            "source_manifest_sha256": self.source_manifest_sha256,
            "accession": self.accession,
            "selection": self.selection.to_dict(),
            "deposited": self.deposited.to_dict(),
            "reference": self.reference.to_dict(),
            "mapping": self.mapping.to_dict(),
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": [item.to_dict() for item in self.discrepancies],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> LigandIdentityRecord:
        data = _as_mapping(data, "ligand identity record")
        return cls(
            schema_version=_required_string(data, "schema_version"),
            source_manifest_path=_required_string(data, "source_manifest_path"),
            source_manifest_sha256=_required_string(data, "source_manifest_sha256"),
            accession=_required_string(data, "accession"),
            selection=LigandSelection.from_dict(_as_mapping(data.get("selection"), "selection")),
            deposited=DepositedLigand.from_dict(_as_mapping(data.get("deposited"), "deposited")),
            reference=ReferenceLigand.from_dict(_as_mapping(data.get("reference"), "reference")),
            mapping=AtomMapping.from_dict(_as_mapping(data.get("mapping"), "mapping")),
            decisions=tuple(
                DecisionRecord.from_dict(_as_mapping(item, "decision record"))
                for item in _as_list(data.get("decisions"), "decisions")
            ),
            provenance=ProvenanceRecord.from_dict(_as_mapping(data.get("provenance"), "provenance")),
            validation=ValidationRecord.from_dict(_as_mapping(data.get("validation"), "validation")),
            discrepancies=tuple(
                IdentityDiscrepancy.from_dict(_as_mapping(item, "identity discrepancy"))
                for item in _as_list(data.get("discrepancies"), "discrepancies")
            ),
        )


def _check(
    category: str,
    check_id: str,
    ok: bool,
    *,
    subject: str,
    criterion: str,
    observed: str,
    expected: str,
    detail: str,
    not_run: bool = False,
) -> ValidationCheck:
    status = "not_run" if not_run else ("pass" if ok else "fail")
    return ValidationCheck(
        check_id=check_id,
        category=category,
        status=status,
        subject=subject,
        criterion=criterion,
        observed=observed,
        expected=expected,
        detail=detail,
    )


def _normalise_formula(formula: str) -> str:
    return "".join(formula.replace(" ", "").split())


def _bond_key(first: Any, second: Any) -> frozenset[Any]:
    return frozenset((first, second))


def _valid_candidate_mappings(
    deposited: DepositedLigand,
    reference: ReferenceLigand,
    candidate_mappings: Sequence[Mapping[str, int]],
) -> tuple[dict[tuple[int, ...], dict[str, int]], int]:
    names = tuple(atom.atom_name for atom in deposited.atoms)
    expected_names = set(names)
    valid: dict[tuple[int, ...], dict[str, int]] = {}
    invalid = 0
    for candidate in candidate_mappings:
        mapping = {str(name): int(index) for name, index in candidate.items()}
        if set(mapping) != expected_names:
            invalid += 1
            continue
        values = tuple(mapping[name] for name in names)
        if len(set(values)) != len(values):
            invalid += 1
            continue
        if any(index < 0 or index >= len(reference.atoms) for index in values):
            invalid += 1
            continue
        valid.setdefault(values, mapping)
    return valid, invalid


def reconcile_ligand_identity(
    *,
    source_manifest_path: str,
    source_manifest_sha256: str,
    accession: str,
    selection: LigandSelection,
    deposited: DepositedLigand,
    reference: ReferenceLigand,
    candidate_mappings: Sequence[Mapping[str, int]],
    mapping_selection_rule: str = (
        "lexicographically smallest reference-index tuple in deposited-atom order"
    ),
) -> LigandIdentityRecord:
    """Compare a deposited ligand with a reference graph.

    ``candidate_mappings`` are explicit deposited-atom-name to reference-index
    mappings.  They are not generated here.  Invalid candidates are skipped and
    reported; a complete, one-to-one candidate is required for a useful record.
    """

    valid_mappings, invalid_mappings = _valid_candidate_mappings(
        deposited, reference, candidate_mappings
    )
    names = tuple(atom.atom_name for atom in deposited.atoms)
    if valid_mappings:
        selected_values, selected_mapping = min(
            valid_mappings.items(), key=lambda item: item[0]
        )
        pairs = tuple((name, selected_mapping[name]) for name in names)
        mapping_pairs = dict(pairs)
    else:
        pairs = ()
        mapping_pairs = {}
    candidate_count = len(valid_mappings)
    ambiguous = candidate_count > 1

    deposited_atoms = {atom.atom_name: atom for atom in deposited.atoms}
    reference_atoms = {atom.index: atom for atom in reference.atoms}
    deposited_bonds = {
        _bond_key(bond.atom_name_1, bond.atom_name_2): bond for bond in deposited.bonds
    }
    reference_bonds = {
        _bond_key(bond.atom_index_1, bond.atom_index_2): bond
        for bond in reference.bonds
    }

    mapped_bonds: dict[frozenset[int], DepositedBond] = {}
    mapping_consistent = bool(mapping_pairs)
    for key, bond in deposited_bonds.items():
        first, second = tuple(key)
        if first not in mapping_pairs or second not in mapping_pairs:
            mapping_consistent = False
            continue
        mapped_bonds[_bond_key(mapping_pairs[first], mapping_pairs[second])] = bond

    connectivity_matches = mapping_consistent and set(mapped_bonds) == set(reference_bonds)
    raw_mismatches = 0
    normalized_mismatches = 0
    if connectivity_matches:
        for key, deposited_bond in mapped_bonds.items():
            reference_bond = reference_bonds[key]
            if deposited_bond.order != reference_bond.order:
                raw_mismatches += 1
            if deposited_bond.aromatic and reference_bond.aromatic:
                continue
            if deposited_bond.order != reference_bond.order:
                normalized_mismatches += 1

    formula_matches = (
        deposited.formula is not None
        and _normalise_formula(deposited.formula) == _normalise_formula(reference.formula)
    )
    heavy_atom_counts_match = len(deposited.atoms) == len(reference.atoms)
    elements_match = Counter(atom.element for atom in deposited.atoms) == Counter(
        atom.element for atom in reference.atoms
    )
    mapping_complete = (
        len(pairs) == len(deposited.atoms)
        and len({index for _, index in pairs}) == len(reference.atoms)
        and {index for _, index in pairs} == set(reference_atoms)
    )
    mapping_elements_agree = bool(pairs) and all(
        name in deposited_atoms
        and index in reference_atoms
        and deposited_atoms[name].element == reference_atoms[index].element
        for name, index in pairs
    )

    if deposited.formal_charges_present:
        formal_charge_matches = mapping_complete and all(
            deposited_atoms[name].formal_charge == reference_atoms[index].formal_charge
            for name, index in pairs
        )
        formal_charge_status = "pass" if formal_charge_matches else "fail"
        formal_charge_observed = "deposited formal charges compared with reference"
        formal_charge_expected = "matching formal charges"
    else:
        formal_charge_matches = False
        formal_charge_status = "not_run"
        formal_charge_observed = "deposited formal charges are not populated"
        formal_charge_expected = f"reference formal charge {reference.formal_charge}"

    if reference.stereocenters:
        if not deposited.stereochemistry:
            stereochemistry_matches = False
            stereochemistry_status = "not_run"
            stereochemistry_observed = "deposited stereochemistry is not specified"
            stereochemistry_expected = ", ".join(reference.stereocenters)
        else:
            stereochemistry_matches = tuple(sorted(deposited.stereochemistry)) == tuple(
                sorted(reference.stereocenters)
            )
            stereochemistry_status = "pass" if stereochemistry_matches else "fail"
            stereochemistry_observed = ", ".join(deposited.stereochemistry)
            stereochemistry_expected = ", ".join(reference.stereocenters)
    else:
        stereochemistry_matches = True
        stereochemistry_status = "pass"
        stereochemistry_observed = "no stereocentres in reference"
        stereochemistry_expected = "no stereocentres"

    if deposited.hydrogen_atoms_present:
        deposited_hydrogens = sum(
            1 for atom in deposited.atoms if atom.element.upper() == "H"
        )
        hydrogens_match = deposited_hydrogens == reference.hydrogen_count
        hydrogens_status = "pass" if hydrogens_match else "fail"
        hydrogens_observed = f"{deposited_hydrogens} deposited hydrogen atoms"
        hydrogens_expected = f"{reference.hydrogen_count} reference hydrogens"
    else:
        hydrogens_match = False
        hydrogens_status = "not_run"
        hydrogens_observed = "no hydrogen atoms in deposited coordinates"
        hydrogens_expected = f"{reference.hydrogen_count} hydrogens in reference formula"

    selection_matches = (
        deposited.component_id == selection.component_id
        and deposited.label_asym_id == (selection.label_asym_id or deposited.label_asym_id)
        and deposited.auth_asym_id == selection.auth_asym_id
        and deposited.auth_seq_id == selection.auth_seq_id
    )

    checks = (
        _check(
            "selection",
            "selection.ligand_matches_record",
            selection_matches,
            subject=f"{deposited.component_id}:{deposited.label_asym_id}",
            criterion="deposited ligand identifiers match the Batch 1 selection",
            observed=(
                f"{deposited.component_id} label {deposited.label_asym_id} auth "
                f"{deposited.auth_asym_id} seq {deposited.auth_seq_id}"
            ),
            expected=(
                f"{selection.component_id} label {selection.label_asym_id} auth "
                f"{selection.auth_asym_id} seq {selection.auth_seq_id}"
            ),
            detail="the identity record was built from the manifest selection",
        ),
        _check(
            "chemical_identity",
            "chemical.formula_matches_reference",
            formula_matches,
            subject=deposited.component_id,
            criterion="deposited formula equals the reference formula",
            observed=deposited.formula or "missing",
            expected=reference.formula,
            detail="formula strings are compared after removing spaces",
        ),
        _check(
            "chemical_identity",
            "chemical.heavy_atom_counts_match",
            heavy_atom_counts_match,
            subject=deposited.component_id,
            criterion="deposited heavy atoms equal reference heavy atoms",
            observed=str(len(deposited.atoms)),
            expected=str(len(reference.atoms)),
            detail="hydrogen atoms are handled by a separate check",
        ),
        _check(
            "chemical_identity",
            "chemical.element_counts_match",
            elements_match,
            subject=deposited.component_id,
            criterion="deposited element counts equal reference element counts",
            observed=str(dict(sorted(Counter(a.element for a in deposited.atoms).items()))),
            expected=str(dict(sorted(Counter(a.element for a in reference.atoms).items()))),
            detail="element symbols are compared without normalization",
        ),
        _check(
            "chemical_identity",
            "chemical.connectivity_matches",
            connectivity_matches,
            subject=deposited.component_id,
            criterion="mapped deposited bonds equal reference bonds",
            observed=f"{len(mapped_bonds)} mapped bonds",
            expected=f"{len(reference_bonds)} reference bonds",
            detail="connectivity is compared as an unordered edge set",
        ),
        _check(
            "chemical_identity",
            "chemical.bond_orders_match",
            connectivity_matches and normalized_mismatches == 0,
            subject=deposited.component_id,
            criterion="bond orders agree after explicit aromatic normalization",
            observed=(
                f"{raw_mismatches} raw differences; "
                f"{normalized_mismatches} after aromatic normalization"
            ),
            expected="0 unresolved differences",
            detail=(
                "raw Kekule orders are recorded before aromatic perception is "
                "applied for the normalized comparison"
            ),
        ),
        _check(
            "chemical_identity",
            "chemical.formal_charge_matches",
            formal_charge_matches,
            subject=deposited.component_id,
            criterion="deposited formal charges equal reference formal charges",
            observed=formal_charge_observed,
            expected=formal_charge_expected,
            detail="missing deposited charges are not inferred",
            not_run=formal_charge_status == "not_run",
        ),
        _check(
            "chemical_identity",
            "chemical.stereochemistry_matches",
            stereochemistry_matches,
            subject=deposited.component_id,
            criterion="deposited stereochemistry agrees with the reference",
            observed=stereochemistry_observed,
            expected=stereochemistry_expected,
            detail="unspecified deposited stereochemistry is not assigned",
            not_run=stereochemistry_status == "not_run",
        ),
        _check(
            "chemical_identity",
            "chemical.hydrogens_match",
            hydrogens_match,
            subject=deposited.component_id,
            criterion="deposited hydrogen count agrees with the reference",
            observed=hydrogens_observed,
            expected=hydrogens_expected,
            detail="no hydrogens are added or removed by this module",
            not_run=hydrogens_status == "not_run",
        ),
        _check(
            "chemical_identity",
            "chemical.atom_mapping_complete",
            mapping_complete,
            subject=deposited.component_id,
            criterion="every deposited atom maps once to every reference heavy atom",
            observed=f"{len(pairs)} mapped pairs",
            expected=f"{len(deposited.atoms)} deposited and {len(reference.atoms)} reference atoms",
            detail="mapping completeness is checked before chemical comparisons",
        ),
        _check(
            "chemical_identity",
            "chemical.atom_mapping_elements_agree",
            mapping_elements_agree,
            subject=deposited.component_id,
            criterion="each mapped atom pair has the same element",
            observed="element agreement checked for selected mapping",
            expected="all mapped pairs agree",
            detail="atom names alone are not treated as chemical identifiers",
        ),
        _check(
            "chemical_identity",
            "chemical.atom_mapping_ambiguity_documented",
            candidate_count > 0,
            subject=deposited.component_id,
            criterion="mapping ambiguity is absent or explicitly recorded",
            observed=f"{candidate_count} valid candidate mappings",
            expected="one or more; ambiguity recorded when more than one",
            detail=(
                "the selected mapping is representative, not a claim that it is unique"
            ),
        ),
        _check(
            "reference",
            "reference.provenance_recorded",
            bool(reference.name and reference.source and reference.smiles and reference.formula),
            subject=reference.name,
            criterion="reference name, source, SMILES, and formula are recorded",
            observed=reference.source,
            expected="documented reference provenance",
            detail="the reference graph is not treated as experimental data",
        ),
    )

    discrepancies: list[IdentityDiscrepancy] = []
    if not formula_matches:
        discrepancies.append(
            IdentityDiscrepancy(
                "formula",
                "open",
                "Deposited and reference formulas differ.",
                observed=deposited.formula,
                expected=reference.formula,
            )
        )
    if not heavy_atom_counts_match or not elements_match:
        discrepancies.append(
            IdentityDiscrepancy(
                "atom_composition",
                "open",
                "Deposited and reference heavy-atom compositions differ.",
                observed=str(dict(sorted(Counter(a.element for a in deposited.atoms).items()))),
                expected=str(dict(sorted(Counter(a.element for a in reference.atoms).items()))),
            )
        )
    if not connectivity_matches:
        discrepancies.append(
            IdentityDiscrepancy(
                "connectivity",
                "open",
                "Mapped deposited connectivity does not equal reference connectivity.",
                observed=f"{len(mapped_bonds)} mapped bonds",
                expected=f"{len(reference_bonds)} reference bonds",
            )
        )
    if raw_mismatches:
        discrepancies.append(
            IdentityDiscrepancy(
                "bond_order_representation",
                "resolved" if normalized_mismatches == 0 else "open",
                (
                    f"{raw_mismatches} bonds differ in raw Kekule assignment; "
                    f"{normalized_mismatches} remain after aromatic normalization."
                ),
                observed=f"{raw_mismatches} raw differences",
                expected="0 raw differences",
                resolution=(
                    "Aromatic perception was applied as a separate comparison; raw "
                    "deposited bond orders remain in the record and were not rewritten."
                    if normalized_mismatches == 0
                    else "A bond-order difference remains after aromatic normalization."
                ),
            )
        )
    if not deposited.formal_charges_present:
        discrepancies.append(
            IdentityDiscrepancy(
                "formal_charge",
                "open",
                "Deposited formal charges are unavailable, so charge identity is unverified.",
                observed="formal charges not populated",
                expected=f"reference formal charge {reference.formal_charge}",
                resolution="No charge assignment was applied in batch 2.",
            )
        )
    elif not formal_charge_matches:
        discrepancies.append(
            IdentityDiscrepancy(
                "formal_charge",
                "open",
                "Deposited formal charges differ from the reference.",
                observed="deposited formal charges do not match",
                expected=f"reference formal charge {reference.formal_charge}",
            )
        )
    if not deposited.hydrogen_atoms_present and reference.hydrogen_count:
        discrepancies.append(
            IdentityDiscrepancy(
                "hydrogen_atoms",
                "not_applicable",
                "The deposited coordinates contain no hydrogen atoms.",
                observed="0 hydrogen atoms",
                expected=f"{reference.hydrogen_count} hydrogens in reference formula",
                resolution="Identity mapping uses heavy atoms only; no hydrogens were added.",
            )
        )
    if reference.stereocenters and not deposited.stereochemistry:
        discrepancies.append(
            IdentityDiscrepancy(
                "stereochemistry",
                "open",
                "Deposited stereochemistry is not specified.",
                observed="unspecified",
                expected=", ".join(reference.stereocenters),
                resolution="No stereochemical assignment was applied in batch 2.",
            )
        )
    elif reference.stereocenters and not stereochemistry_matches:
        discrepancies.append(
            IdentityDiscrepancy(
                "stereochemistry",
                "open",
                "Deposited stereochemistry differs from the reference.",
                observed=", ".join(deposited.stereochemistry),
                expected=", ".join(reference.stereocenters),
            )
        )
    if deposited.hydrogen_atoms_present and not hydrogens_match:
        discrepancies.append(
            IdentityDiscrepancy(
                "hydrogen_atoms",
                "open",
                "Deposited and reference hydrogen counts differ.",
                observed=hydrogens_observed,
                expected=hydrogens_expected,
            )
        )
    if candidate_count == 0:
        discrepancies.append(
            IdentityDiscrepancy(
                "atom_mapping",
                "open",
                "No complete one-to-one atom mapping was supplied.",
                observed=f"{invalid_mappings} invalid candidates",
                expected="at least one complete one-to-one mapping",
                resolution="Invalid candidates were not repaired or merged.",
            )
        )
    elif invalid_mappings:
        discrepancies.append(
            IdentityDiscrepancy(
                "atom_mapping",
                "open",
                f"{invalid_mappings} candidate mappings were rejected as incomplete or invalid.",
                observed=f"{invalid_mappings} invalid candidates",
                expected="complete one-to-one mappings",
                resolution="Invalid candidates were not repaired or merged.",
            )
        )

    decisions = (
        DecisionRecord(
            decision="selected a representative atom mapping",
            actor="ligand identity reconciliation",
            rationale=(
                f"{candidate_count} valid candidate mappings; {mapping_selection_rule}. "
                + (
                    "The record preserves the ambiguity and does not claim uniqueness."
                    if ambiguous
                    else "The mapping is unique among the supplied candidates."
                )
            ),
        ),
        DecisionRecord(
            decision="compared raw bond orders before aromatic normalization",
            actor="ligand identity reconciliation",
            rationale=(
                "raw deposited orders were preserved; aromatic perception was used "
                "only for a separate normalized comparison."
            ),
        ),
        DecisionRecord(
            decision="did not add hydrogens, change protonation, or assign formal charges",
            actor="ligand identity reconciliation",
            rationale=(
                "missing deposited information is reported as not_run or open, not "
                "inferred from the reference."
            ),
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "deposited ligand graph",
                (
                    f"Component {deposited.component_id} from the selected structure "
                    "record; atoms and bonds are taken from the deposited mmCIF."
                ),
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "reference molecular graph",
                (
                    f"{reference.name} from {reference.source}; formula "
                    f"{reference.formula}, formal charge {reference.formal_charge}."
                ),
                actor="reference adapter",
            ),
            ProvenanceEntry(
                "calculated",
                "identity comparison",
                (
                    f"{candidate_count} valid mapping candidates; raw bond-order "
                    f"differences {raw_mismatches}; unresolved normalized differences "
                    f"{normalized_mismatches}."
                ),
                actor="ligand identity reconciliation",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No docking pose, affinity, or free-energy estimate is produced.",
                actor="ligand identity reconciliation",
            ),
            ProvenanceEntry(
                "hypothesised",
                "none",
                "No mechanistic interpretation is made from the comparison.",
                actor="ligand identity reconciliation",
            ),
        )
    )

    return LigandIdentityRecord(
        source_manifest_path=source_manifest_path,
        source_manifest_sha256=source_manifest_sha256,
        accession=accession,
        selection=selection,
        deposited=deposited,
        reference=reference,
        mapping=AtomMapping(
            pairs=pairs,
            candidate_count=candidate_count,
            ambiguous=ambiguous,
            selection_rule=mapping_selection_rule,
        ),
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=tuple(discrepancies),
    )


def write_ligand_identity(path: str | Path, record: LigandIdentityRecord) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_ligand_identity(path: str | Path) -> LigandIdentityRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return LigandIdentityRecord.from_dict(data)


def validate_ligand_identity(record: LigandIdentityRecord) -> list[str]:
    """Return schema and internal-consistency errors for *record*.

    Scientific failures, ambiguous mappings, and open discrepancies belong in
    the record's checks and discrepancy list.  They do not make the JSON record
    structurally invalid.
    """

    errors: list[str] = []
    if record.schema_version != LIGAND_IDENTITY_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {LIGAND_IDENTITY_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    if not record.source_manifest_path:
        errors.append("source_manifest_path is required.")
    if not _SHA256_PATTERN.fullmatch(record.source_manifest_sha256):
        errors.append("source_manifest_sha256 must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")

    deposited = record.deposited
    if not deposited.component_id:
        errors.append("deposited.component_id is required.")
    if not deposited.label_asym_id or not deposited.auth_asym_id:
        errors.append("deposited label and author asymmetry identifiers are required.")
    if not deposited.auth_seq_id:
        errors.append("deposited.auth_seq_id is required.")
    if not deposited.atoms:
        errors.append("deposited.atoms must not be empty.")
    if deposited.omitted_hydrogen_bonds < 0:
        errors.append("deposited.omitted_hydrogen_bonds cannot be negative.")

    atom_names = [atom.atom_name for atom in deposited.atoms]
    if len(set(atom_names)) != len(atom_names):
        errors.append("deposited atom names must be unique.")
    has_deposited_hydrogen = any(
        atom.element.upper() == "H" for atom in deposited.atoms
    )
    if deposited.atoms_are_heavy and has_deposited_hydrogen:
        errors.append("deposited.atoms_are_heavy is true but hydrogen atoms are present.")
    if deposited.hydrogen_atoms_present != has_deposited_hydrogen:
        errors.append(
            "deposited.hydrogen_atoms_present must match the deposited atom elements."
        )
    atom_name_set = set(atom_names)
    for index, bond in enumerate(deposited.bonds):
        if bond.atom_name_1 not in atom_name_set or bond.atom_name_2 not in atom_name_set:
            errors.append(f"deposited.bonds[{index}] references an unknown atom name.")
        if bond.atom_name_1 == bond.atom_name_2:
            errors.append(f"deposited.bonds[{index}] is a self-bond.")

    reference = record.reference
    if not reference.name or not reference.source or not reference.smiles:
        errors.append("reference name, source, and SMILES are required.")
    if not reference.formula:
        errors.append("reference.formula is required.")
    if not reference.atoms:
        errors.append("reference.atoms must not be empty.")
    if reference.atoms_are_heavy and any(
        atom.element.upper() == "H" for atom in reference.atoms
    ):
        errors.append("reference.atoms_are_heavy is true but hydrogen atoms are present.")
    reference_indices = [atom.index for atom in reference.atoms]
    if len(set(reference_indices)) != len(reference_indices):
        errors.append("reference atom indices must be unique.")
    reference_index_set = set(reference_indices)
    for index, bond in enumerate(reference.bonds):
        if (
            bond.atom_index_1 not in reference_index_set
            or bond.atom_index_2 not in reference_index_set
        ):
            errors.append(f"reference.bonds[{index}] references an unknown atom index.")
        if bond.atom_index_1 == bond.atom_index_2:
            errors.append(f"reference.bonds[{index}] is a self-bond.")

    mapping = record.mapping
    pair_names = [name for name, _ in mapping.pairs]
    pair_indices = [index for _, index in mapping.pairs]
    if len(set(pair_names)) != len(pair_names):
        errors.append("mapping deposited atom names must be unique.")
    if len(set(pair_indices)) != len(pair_indices):
        errors.append("mapping reference indices must be unique.")
    if mapping.candidate_count < 0:
        errors.append("mapping.candidate_count cannot be negative.")
    if mapping.ambiguous != (mapping.candidate_count > 1):
        errors.append("mapping.ambiguous must match candidate_count > 1.")
    if mapping.candidate_count > 0 and not mapping.pairs:
        errors.append("mapping.pairs must not be empty when candidates are available.")
    if mapping.candidate_count == 0 and mapping.pairs:
        errors.append("mapping.pairs must be empty when candidate_count is zero.")
    if mapping.pairs:
        if set(pair_names) != atom_name_set:
            errors.append("mapping must cover every deposited atom.")
        if set(pair_indices) != reference_index_set:
            errors.append("mapping must cover every reference heavy atom.")

    if not record.decisions:
        errors.append("decisions must contain at least one record.")
    for index, decision in enumerate(record.decisions):
        if not decision.decision or not decision.actor or not decision.rationale:
            errors.append(f"decisions[{index}] requires decision, actor, and rationale.")

    if not record.provenance.records:
        errors.append("provenance.records must not be empty.")
    for index, entry in enumerate(record.provenance.records):
        if entry.label not in PROVENANCE_LABELS:
            errors.append(f"provenance.records[{index}] has an invalid label.")

    check_ids: set[str] = set()
    for index, check in enumerate(record.validation.checks):
        if check.status not in ("pass", "fail", "not_run"):
            errors.append(f"validation.checks[{index}] has an invalid status.")
        if not check.check_id or not check.category or not check.detail:
            errors.append(f"validation.checks[{index}] requires check_id, category, and detail.")
        if check.check_id in check_ids:
            errors.append(f"validation.checks[{index}].check_id is duplicated.")
        check_ids.add(check.check_id)

    for index, item in enumerate(record.discrepancies):
        if item.status not in DISCREPANCY_STATUSES:
            errors.append(f"discrepancies[{index}] has an invalid status.")
        if not item.category or not item.detail:
            errors.append(f"discrepancies[{index}] requires category and detail.")
        if item.status == "resolved" and not item.resolution:
            errors.append(f"discrepancies[{index}] is resolved but has no resolution.")

    return errors
