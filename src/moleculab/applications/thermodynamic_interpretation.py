# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Thermodynamic interpretation boundaries for structural evidence.

This layer connects the capstone's structural records to the questions a
binding free-energy calculation would have to answer.  It does not estimate an
affinity.  It records which thermodynamic terms and preparation/sampling
requirements are missing, and why interaction counts are not additive energies.

The central invariant is that ``binding_free_energy_kj_per_mol`` is ``None`` and
``free_energy_status`` is ``not_computed``.
"""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from .structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    LigandSelection,
    ProteinSelection,
    ProvenanceEntry,
    ProvenanceRecord,
    ValidationCheck,
    ValidationRecord,
)

THERMODYNAMIC_INTERPRETATION_SCHEMA_VERSION = "1.0"
BOUNDARY_STATUSES = ("not_established", "not_estimated", "not_computed")
TERM_STATUSES = ("not_estimated", "not_available", "not_applicable")
REQUIREMENT_STATUSES = ("met", "partial", "not_run", "not_available")
FREE_ENERGY_STATUSES = ("not_computed",)
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


def _required_int(data: Mapping[str, Any], key: str) -> int:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return int(value)


def _optional_float(data: Mapping[str, Any], key: str) -> float | None:
    value = data.get(key)
    return None if value is None else float(value)


def _optional_string(data: Mapping[str, Any], key: str) -> str | None:
    value = data.get(key)
    return None if value is None else str(value)


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


@dataclass(frozen=True)
class RobustnessStatusCount:
    status: str
    count: int

    def to_dict(self) -> dict[str, Any]:
        return {"status": self.status, "count": self.count}

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> RobustnessStatusCount:
        data = _as_mapping(data, "robustness status count")
        return cls(
            status=_required_string(data, "status"),
            count=_required_int(data, "count"),
        )


@dataclass(frozen=True)
class StructuralEvidenceSummary:
    ligand_state_count: int
    protein_state_count: int
    combination_count: int
    ligand_contact_atom_count: int
    protein_contact_atom_count: int
    hydrophobic_pair_count: int
    aromatic_pair_count: int
    supported_hydrogen_bond_count: int
    plausible_water_bridge_count: int
    robustness_conclusion_counts: tuple[RobustnessStatusCount, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "ligand_state_count": self.ligand_state_count,
            "protein_state_count": self.protein_state_count,
            "combination_count": self.combination_count,
            "ligand_contact_atom_count": self.ligand_contact_atom_count,
            "protein_contact_atom_count": self.protein_contact_atom_count,
            "hydrophobic_pair_count": self.hydrophobic_pair_count,
            "aromatic_pair_count": self.aromatic_pair_count,
            "supported_hydrogen_bond_count": self.supported_hydrogen_bond_count,
            "plausible_water_bridge_count": self.plausible_water_bridge_count,
            "robustness_conclusion_counts": [
                item.to_dict() for item in self.robustness_conclusion_counts
            ],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> StructuralEvidenceSummary:
        data = _as_mapping(data, "structural evidence summary")
        return cls(
            ligand_state_count=_required_int(data, "ligand_state_count"),
            protein_state_count=_required_int(data, "protein_state_count"),
            combination_count=_required_int(data, "combination_count"),
            ligand_contact_atom_count=_required_int(
                data, "ligand_contact_atom_count"
            ),
            protein_contact_atom_count=_required_int(
                data, "protein_contact_atom_count"
            ),
            hydrophobic_pair_count=_required_int(data, "hydrophobic_pair_count"),
            aromatic_pair_count=_required_int(data, "aromatic_pair_count"),
            supported_hydrogen_bond_count=_required_int(
                data, "supported_hydrogen_bond_count"
            ),
            plausible_water_bridge_count=_required_int(
                data, "plausible_water_bridge_count"
            ),
            robustness_conclusion_counts=tuple(
                RobustnessStatusCount.from_dict(
                    _as_mapping(item, "robustness status count")
                )
                for item in _as_list(
                    data.get("robustness_conclusion_counts"),
                    "robustness_conclusion_counts",
                )
            ),
        )


@dataclass(frozen=True)
class StructuralToThermodynamicBoundary:
    boundary_id: str
    structural_observation: str
    thermodynamic_claim: str
    status: str
    reason: str
    required_evidence: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "boundary_id": self.boundary_id,
            "structural_observation": self.structural_observation,
            "thermodynamic_claim": self.thermodynamic_claim,
            "status": self.status,
            "reason": self.reason,
            "required_evidence": list(self.required_evidence),
        }

    @classmethod
    def from_dict(
        cls,
        data: Mapping[str, Any],
    ) -> StructuralToThermodynamicBoundary:
        data = _as_mapping(data, "structural-to-thermodynamic boundary")
        return cls(
            boundary_id=_required_string(data, "boundary_id"),
            structural_observation=_required_string(
                data, "structural_observation"
            ),
            thermodynamic_claim=_required_string(data, "thermodynamic_claim"),
            status=_required_string(data, "status"),
            reason=_required_string(data, "reason"),
            required_evidence=_string_tuple(data, "required_evidence"),
        )


@dataclass(frozen=True)
class FreeEnergyTerm:
    term_id: str
    name: str
    status: str
    detail: str
    required_evidence: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "term_id": self.term_id,
            "name": self.name,
            "status": self.status,
            "detail": self.detail,
            "required_evidence": list(self.required_evidence),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> FreeEnergyTerm:
        data = _as_mapping(data, "free-energy term")
        return cls(
            term_id=_required_string(data, "term_id"),
            name=_required_string(data, "name"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            required_evidence=_string_tuple(data, "required_evidence"),
        )


@dataclass(frozen=True)
class FreeEnergyRequirement:
    requirement_id: str
    name: str
    status: str
    detail: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "requirement_id": self.requirement_id,
            "name": self.name,
            "status": self.status,
            "detail": self.detail,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> FreeEnergyRequirement:
        data = _as_mapping(data, "free-energy requirement")
        return cls(
            requirement_id=_required_string(data, "requirement_id"),
            name=_required_string(data, "name"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
        )


@dataclass(frozen=True)
class ThermodynamicDiscrepancy:
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
    def from_dict(cls, data: Mapping[str, Any]) -> ThermodynamicDiscrepancy:
        data = _as_mapping(data, "thermodynamic discrepancy")
        return cls(
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            detail=_required_string(data, "detail"),
            observed=_optional_string(data, "observed"),
            expected=_optional_string(data, "expected"),
            resolution=_optional_string(data, "resolution"),
        )


@dataclass(frozen=True)
class ThermodynamicInterpretationRecord:
    source_manifest_path: str
    source_manifest_sha256: str
    ligand_identity_record_path: str
    ligand_identity_record_sha256: str
    ligand_coordinate_record_path: str
    ligand_coordinate_record_sha256: str
    binding_site_record_path: str
    binding_site_record_sha256: str
    chemical_state_record_path: str
    chemical_state_record_sha256: str
    protein_preparation_record_path: str
    protein_preparation_record_sha256: str
    hydrogen_evidence_record_path: str
    hydrogen_evidence_record_sha256: str
    water_evidence_record_path: str
    water_evidence_record_sha256: str
    combination_interaction_record_path: str
    combination_interaction_record_sha256: str
    interaction_robustness_record_path: str
    interaction_robustness_record_sha256: str
    accession: str
    protein_selection: ProteinSelection
    ligand_selection: LigandSelection
    reference_ligand_state_id: str
    reference_protein_state_id: str
    structural_summary: StructuralEvidenceSummary
    ensemble_status: str
    ensemble_detail: str
    free_energy_status: str
    binding_free_energy_kj_per_mol: float | None
    boundaries: tuple[StructuralToThermodynamicBoundary, ...]
    terms: tuple[FreeEnergyTerm, ...]
    requirements: tuple[FreeEnergyRequirement, ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[ThermodynamicDiscrepancy, ...]
    schema_version: str = THERMODYNAMIC_INTERPRETATION_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "source_manifest_path": self.source_manifest_path,
            "source_manifest_sha256": self.source_manifest_sha256,
            "ligand_identity_record_path": self.ligand_identity_record_path,
            "ligand_identity_record_sha256": self.ligand_identity_record_sha256,
            "ligand_coordinate_record_path": self.ligand_coordinate_record_path,
            "ligand_coordinate_record_sha256": self.ligand_coordinate_record_sha256,
            "binding_site_record_path": self.binding_site_record_path,
            "binding_site_record_sha256": self.binding_site_record_sha256,
            "chemical_state_record_path": self.chemical_state_record_path,
            "chemical_state_record_sha256": self.chemical_state_record_sha256,
            "protein_preparation_record_path": self.protein_preparation_record_path,
            "protein_preparation_record_sha256": (
                self.protein_preparation_record_sha256
            ),
            "hydrogen_evidence_record_path": self.hydrogen_evidence_record_path,
            "hydrogen_evidence_record_sha256": self.hydrogen_evidence_record_sha256,
            "water_evidence_record_path": self.water_evidence_record_path,
            "water_evidence_record_sha256": self.water_evidence_record_sha256,
            "combination_interaction_record_path": (
                self.combination_interaction_record_path
            ),
            "combination_interaction_record_sha256": (
                self.combination_interaction_record_sha256
            ),
            "interaction_robustness_record_path": (
                self.interaction_robustness_record_path
            ),
            "interaction_robustness_record_sha256": (
                self.interaction_robustness_record_sha256
            ),
            "accession": self.accession,
            "protein_selection": self.protein_selection.to_dict(),
            "ligand_selection": self.ligand_selection.to_dict(),
            "reference_ligand_state_id": self.reference_ligand_state_id,
            "reference_protein_state_id": self.reference_protein_state_id,
            "structural_summary": self.structural_summary.to_dict(),
            "ensemble_status": self.ensemble_status,
            "ensemble_detail": self.ensemble_detail,
            "free_energy_status": self.free_energy_status,
            "binding_free_energy_kj_per_mol": self.binding_free_energy_kj_per_mol,
            "boundaries": [item.to_dict() for item in self.boundaries],
            "terms": [item.to_dict() for item in self.terms],
            "requirements": [item.to_dict() for item in self.requirements],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": [item.to_dict() for item in self.discrepancies],
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ThermodynamicInterpretationRecord:
        data = _as_mapping(data, "thermodynamic interpretation record")
        return cls(
            schema_version=_required_string(data, "schema_version"),
            source_manifest_path=_required_string(data, "source_manifest_path"),
            source_manifest_sha256=_required_string(
                data, "source_manifest_sha256"
            ),
            ligand_identity_record_path=_required_string(
                data, "ligand_identity_record_path"
            ),
            ligand_identity_record_sha256=_required_string(
                data, "ligand_identity_record_sha256"
            ),
            ligand_coordinate_record_path=_required_string(
                data, "ligand_coordinate_record_path"
            ),
            ligand_coordinate_record_sha256=_required_string(
                data, "ligand_coordinate_record_sha256"
            ),
            binding_site_record_path=_required_string(
                data, "binding_site_record_path"
            ),
            binding_site_record_sha256=_required_string(
                data, "binding_site_record_sha256"
            ),
            chemical_state_record_path=_required_string(
                data, "chemical_state_record_path"
            ),
            chemical_state_record_sha256=_required_string(
                data, "chemical_state_record_sha256"
            ),
            protein_preparation_record_path=_required_string(
                data, "protein_preparation_record_path"
            ),
            protein_preparation_record_sha256=_required_string(
                data, "protein_preparation_record_sha256"
            ),
            hydrogen_evidence_record_path=_required_string(
                data, "hydrogen_evidence_record_path"
            ),
            hydrogen_evidence_record_sha256=_required_string(
                data, "hydrogen_evidence_record_sha256"
            ),
            water_evidence_record_path=_required_string(
                data, "water_evidence_record_path"
            ),
            water_evidence_record_sha256=_required_string(
                data, "water_evidence_record_sha256"
            ),
            combination_interaction_record_path=_required_string(
                data, "combination_interaction_record_path"
            ),
            combination_interaction_record_sha256=_required_string(
                data, "combination_interaction_record_sha256"
            ),
            interaction_robustness_record_path=_required_string(
                data, "interaction_robustness_record_path"
            ),
            interaction_robustness_record_sha256=_required_string(
                data, "interaction_robustness_record_sha256"
            ),
            accession=_required_string(data, "accession"),
            protein_selection=ProteinSelection.from_dict(
                _as_mapping(data.get("protein_selection"), "protein_selection")
            ),
            ligand_selection=LigandSelection.from_dict(
                _as_mapping(data.get("ligand_selection"), "ligand_selection")
            ),
            reference_ligand_state_id=_required_string(
                data, "reference_ligand_state_id"
            ),
            reference_protein_state_id=_required_string(
                data, "reference_protein_state_id"
            ),
            structural_summary=StructuralEvidenceSummary.from_dict(
                _as_mapping(data.get("structural_summary"), "structural_summary")
            ),
            ensemble_status=_required_string(data, "ensemble_status"),
            ensemble_detail=_required_string(data, "ensemble_detail"),
            free_energy_status=_required_string(data, "free_energy_status"),
            binding_free_energy_kj_per_mol=_optional_float(
                data, "binding_free_energy_kj_per_mol"
            ),
            boundaries=tuple(
                StructuralToThermodynamicBoundary.from_dict(
                    _as_mapping(item, "structural-to-thermodynamic boundary")
                )
                for item in _as_list(data.get("boundaries"), "boundaries")
            ),
            terms=tuple(
                FreeEnergyTerm.from_dict(_as_mapping(item, "free-energy term"))
                for item in _as_list(data.get("terms"), "terms")
            ),
            requirements=tuple(
                FreeEnergyRequirement.from_dict(
                    _as_mapping(item, "free-energy requirement")
                )
                for item in _as_list(data.get("requirements"), "requirements")
            ),
            decisions=tuple(
                DecisionRecord.from_dict(_as_mapping(item, "decision record"))
                for item in _as_list(data.get("decisions"), "decisions")
            ),
            provenance=ProvenanceRecord.from_dict(
                _as_mapping(data.get("provenance"), "provenance")
            ),
            validation=ValidationRecord.from_dict(
                _as_mapping(data.get("validation"), "validation")
            ),
            discrepancies=tuple(
                ThermodynamicDiscrepancy.from_dict(
                    _as_mapping(item, "thermodynamic discrepancy")
                )
                for item in _as_list(data.get("discrepancies"), "discrepancies")
            ),
        )


def structural_thermodynamic_boundaries(
    summary: StructuralEvidenceSummary,
) -> tuple[StructuralToThermodynamicBoundary, ...]:
    """Return standard boundary statements for the capstone summary."""

    robustness_counts = {
        item.status: item.count for item in summary.robustness_conclusion_counts
    }
    return (
        StructuralToThermodynamicBoundary(
            boundary_id="hydrophobic_pairs_not_energy",
            structural_observation=(
                f"{summary.hydrophobic_pair_count} hydrophobic atom pairs are "
                "classified in the reference state combination."
            ),
            thermodynamic_claim=(
                "Each pair contributes a known favourable binding enthalpy."
            ),
            status="not_established",
            reason=(
                "Interaction classification is geometric and template-conditioned; "
                "no pair energy, desolvation, entropy, or ensemble average is "
                "available."
            ),
            required_evidence=(
                "force-field interaction energy",
                "desolvation free energy",
                "conformational entropy",
                "equilibrium ensemble",
            ),
        ),
        StructuralToThermodynamicBoundary(
            boundary_id="aromatic_pair_not_energy",
            structural_observation=(
                f"{summary.aromatic_pair_count} aromatic ring pair is classified."
            ),
            thermodynamic_claim=(
                "The ring pair contributes a known ring-stacking free energy."
            ),
            status="not_established",
            reason=(
                "Robustness across the tested cutoffs does not supply an energy "
                "or a population."
            ),
            required_evidence=(
                "ring interaction energy",
                "solvation free energy",
                "configurational population",
            ),
        ),
        StructuralToThermodynamicBoundary(
            boundary_id="hydrogen_bonds_not_energy",
            structural_observation=(
                f"{summary.supported_hydrogen_bond_count} ligand-donor "
                "candidates satisfy the tested D-H...A criteria."
            ),
            thermodynamic_claim=(
                "These candidates contribute a known hydrogen-bond free energy."
            ),
            status="not_established",
            reason=(
                "Hydrogen positions are generated, protein-donor hydrogens are "
                "unavailable, and geometric criteria are not energies."
            ),
            required_evidence=(
                "experimental or converged hydrogen positions",
                "hydrogen-bond energy model",
                "solvent and entropy contributions",
            ),
        ),
        StructuralToThermodynamicBoundary(
            boundary_id="water_bridge_not_energy",
            structural_observation=(
                f"{summary.plausible_water_bridge_count} water bridge is "
                "geometrically plausible."
            ),
            thermodynamic_claim=(
                "The water-mediated interaction has a known free-energy "
                "contribution."
            ),
            status="not_established",
            reason=(
                "Water orientation and roles are unassigned, crystallographic "
                "occupancy is not an equilibrium probability, and no displacement "
                "free energy is available."
            ),
            required_evidence=(
                "water orientation and chemical roles",
                "equilibrium water occupancy",
                "water displacement free energy",
            ),
        ),
        StructuralToThermodynamicBoundary(
            boundary_id="counts_not_ranking",
            structural_observation=(
                f"{summary.combination_count} state combinations each classify "
                f"{summary.hydrophobic_pair_count + summary.aromatic_pair_count} "
                "interactions."
            ),
            thermodynamic_claim=(
                "A combination with more classified interactions binds more "
                "strongly."
            ),
            status="not_established",
            reason=(
                "Counts are not additive energies, and chemical-state hypotheses "
                "are not equilibrium populations."
            ),
            required_evidence=(
                "free-energy calculation per state",
                "state populations",
                "experimental affinity validation",
            ),
        ),
        StructuralToThermodynamicBoundary(
            boundary_id="snapshot_not_ensemble",
            structural_observation=(
                "The deposited crystal structure supplies one structural "
                "observation."
            ),
            thermodynamic_claim=(
                "The crystal structure defines the equilibrium Boltzmann "
                "population."
            ),
            status="not_established",
            reason=(
                "A crystal structure is a structural observation; no ensemble or "
                "partition function is available."
            ),
            required_evidence=(
                "equilibrium sampling",
                "configurational partition function",
                "population estimates with uncertainty",
            ),
        ),
        StructuralToThermodynamicBoundary(
            boundary_id="robustness_not_certainty",
            structural_observation=(
                f"Robustness conclusions: "
                f"{robustness_counts.get('invariant', 0)} invariant, "
                f"{robustness_counts.get('assumption_sensitive', 0)} "
                "assumption-sensitive, "
                f"{robustness_counts.get('unsupported', 0)} unsupported, "
                f"{robustness_counts.get('unknown', 0)} unknown."
            ),
            thermodynamic_claim=(
                "Invariant classifications are experimentally validated or "
                "physically certain."
            ),
            status="not_established",
            reason=(
                "Invariance is only with respect to the tested assumption range; "
                "it is not an energy or an experimental validation."
            ),
            required_evidence=(
                "experimental validation",
                "free-energy calculation",
                "tested-range justification",
            ),
        ),
        StructuralToThermodynamicBoundary(
            boundary_id="states_not_populations",
            structural_observation=(
                f"{summary.ligand_state_count} ligand states and "
                f"{summary.protein_state_count} protein states are enumerated."
            ),
            thermodynamic_claim=(
                "Their relative populations or free energies are known."
            ),
            status="not_estimated",
            reason=(
                "No pKa/population model or state free energies are supplied."
            ),
            required_evidence=(
                "state free energies",
                "pKa and tautomer populations",
                "experimental state validation",
            ),
        ),
    )


def default_free_energy_terms() -> tuple[FreeEnergyTerm, ...]:
    """Return the thermodynamic terms that would be needed for a binding ΔG."""

    return (
        FreeEnergyTerm(
            term_id="interaction_enthalpy",
            name="Interaction enthalpy",
            status="not_estimated",
            detail=(
                "Requires a force field and an ensemble average of ligand-protein "
                "interaction energies in the bound state."
            ),
            required_evidence=(
                "force-field parameters",
                "bound-state ensemble",
                "interaction energy decomposition",
            ),
        ),
        FreeEnergyTerm(
            term_id="desolvation",
            name="Desolvation",
            status="not_estimated",
            detail=(
                "Requires solvation free energies for the ligand, protein, and "
                "complex, including polar and nonpolar contributions."
            ),
            required_evidence=(
                "solvent model",
                "solvation free-energy calculation",
                "water and ion preparation",
            ),
        ),
        FreeEnergyTerm(
            term_id="ligand_conformational_entropy",
            name="Ligand conformational entropy",
            status="not_estimated",
            detail=(
                "Requires conformational sampling of the free and bound ligand "
                "and an entropy estimator."
            ),
            required_evidence=(
                "ligand force field",
                "conformational sampling",
                "entropy estimation",
            ),
        ),
        FreeEnergyTerm(
            term_id="protein_conformational_entropy",
            name="Protein conformational entropy",
            status="not_estimated",
            detail=(
                "Requires sampling of protein side chains, loops, and possible "
                "conformational changes."
            ),
            required_evidence=(
                "protein force field",
                "protein conformational sampling",
                "entropy estimation",
            ),
        ),
        FreeEnergyTerm(
            term_id="water_reorganisation",
            name="Water reorganisation",
            status="not_estimated",
            detail=(
                "Requires water roles, residence behaviour, and the free-energy "
                "cost of displacing or reorganising waters."
            ),
            required_evidence=(
                "water model and roles",
                "equilibrium water occupancy",
                "water displacement free energy",
            ),
        ),
        FreeEnergyTerm(
            term_id="translational_rotational_entropy",
            name="Translational and rotational entropy",
            status="not_estimated",
            detail=(
                "Requires a defined standard state and a free-energy method that "
                "treats the loss of translational and rotational freedom."
            ),
            required_evidence=(
                "standard-state definition",
                "free-energy estimator",
                "symmetry and degeneracy treatment",
            ),
        ),
        FreeEnergyTerm(
            term_id="protonation_tautomer_populations",
            name="Protonation and tautomer populations",
            status="not_estimated",
            detail=(
                "The state hypotheses are enumerated, but no state free energies "
                "or populations are available."
            ),
            required_evidence=(
                "state free energies",
                "pKa and tautomer populations",
                "experimental state validation",
            ),
        ),
        FreeEnergyTerm(
            term_id="conformational_strain",
            name="Conformational strain",
            status="not_estimated",
            detail=(
                "Requires a force field and a bound-versus-free strain energy "
                "calculation."
            ),
            required_evidence=(
                "ligand and protein force fields",
                "free and bound conformer ensembles",
                "strain energy calculation",
            ),
        ),
    )


def default_free_energy_requirements() -> tuple[FreeEnergyRequirement, ...]:
    """Return the preparation, sampling, and validation requirements."""

    return (
        FreeEnergyRequirement(
            requirement_id="molecular_model",
            name="Defensible molecular model",
            status="partial",
            detail=(
                "Deposited heavy atoms and explicit chemical-state hypotheses are "
                "available; no force-field-ready model is defined."
            ),
        ),
        FreeEnergyRequirement(
            requirement_id="protein_preparation",
            name="Protein preparation",
            status="partial",
            detail=(
                "Pocket atoms and residue-template chemistry are available; "
                "protein hydrogens, protonation populations, and force-field "
                "parameters are not."
            ),
        ),
        FreeEnergyRequirement(
            requirement_id="ligand_preparation",
            name="Ligand preparation",
            status="partial",
            detail=(
                "The ligand graph and generated hydrogens are available; "
                "force-field charges, parameters, and a validated protonation "
                "population are not."
            ),
        ),
        FreeEnergyRequirement(
            requirement_id="water_ion_preparation",
            name="Water and ion preparation",
            status="not_run",
            detail=(
                "Deposited waters are preserved with unassigned roles; no "
                "equilibrated water or ion environment has been prepared."
            ),
        ),
        FreeEnergyRequirement(
            requirement_id="force_field_parameters",
            name="Force-field parameters",
            status="not_run",
            detail=(
                "No ligand, protein, water, or ion force-field parameter set has "
                "been assigned or validated."
            ),
        ),
        FreeEnergyRequirement(
            requirement_id="conformational_sampling",
            name="Conformational sampling",
            status="not_run",
            detail=(
                "No MD, Monte Carlo, or alchemical sampling has been performed "
                "for this complex."
            ),
        ),
        FreeEnergyRequirement(
            requirement_id="thermodynamic_states",
            name="Defined thermodynamic states",
            status="partial",
            detail=(
                "Ligand and protein chemical-state hypotheses exist, but the "
                "bound/free endpoints and any alchemical path are not defined."
            ),
        ),
        FreeEnergyRequirement(
            requirement_id="standard_state",
            name="Standard state and temperature",
            status="not_run",
            detail=(
                "No standard-state concentration or temperature convention has "
                "been fixed for a binding free-energy calculation."
            ),
        ),
        FreeEnergyRequirement(
            requirement_id="convergence_uncertainty",
            name="Convergence and uncertainty",
            status="not_run",
            detail=(
                "No free-energy estimate, replicate analysis, or uncertainty "
                "propagation is available."
            ),
        ),
        FreeEnergyRequirement(
            requirement_id="experimental_affinity_validation",
            name="Experimental affinity validation",
            status="not_available",
            detail=(
                "No matched assay and affinity record is part of this capstone "
                "record."
            ),
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
    if not_run:
        status = "not_run"
    else:
        status = "pass" if ok else "fail"
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


def build_thermodynamic_interpretation_record(
    *,
    source_manifest_path: str,
    source_manifest_sha256: str,
    ligand_identity_record_path: str,
    ligand_identity_record_sha256: str,
    ligand_coordinate_record_path: str,
    ligand_coordinate_record_sha256: str,
    binding_site_record_path: str,
    binding_site_record_sha256: str,
    chemical_state_record_path: str,
    chemical_state_record_sha256: str,
    protein_preparation_record_path: str,
    protein_preparation_record_sha256: str,
    hydrogen_evidence_record_path: str,
    hydrogen_evidence_record_sha256: str,
    water_evidence_record_path: str,
    water_evidence_record_sha256: str,
    combination_interaction_record_path: str,
    combination_interaction_record_sha256: str,
    interaction_robustness_record_path: str,
    interaction_robustness_record_sha256: str,
    accession: str,
    protein_selection: ProteinSelection,
    ligand_selection: LigandSelection,
    reference_ligand_state_id: str,
    reference_protein_state_id: str,
    structural_summary: StructuralEvidenceSummary,
    boundaries: tuple[StructuralToThermodynamicBoundary, ...] | None = None,
    terms: tuple[FreeEnergyTerm, ...] | None = None,
    requirements: tuple[FreeEnergyRequirement, ...] | None = None,
) -> ThermodynamicInterpretationRecord:
    boundaries = (
        structural_thermodynamic_boundaries(structural_summary)
        if boundaries is None
        else boundaries
    )
    terms = default_free_energy_terms() if terms is None else terms
    requirements = (
        default_free_energy_requirements()
        if requirements is None
        else requirements
    )

    checks = (
        _check(
            "provenance",
            "thermodynamics.source_records",
            bool(
                source_manifest_path
                and ligand_identity_record_path
                and ligand_coordinate_record_path
                and binding_site_record_path
                and chemical_state_record_path
                and protein_preparation_record_path
                and hydrogen_evidence_record_path
                and water_evidence_record_path
                and combination_interaction_record_path
                and interaction_robustness_record_path
            ),
            subject=accession,
            criterion="all structural, interaction, and robustness records are referenced",
            observed="ten input records present",
            expected=(
                "manifest, identity, coordinates, binding site, ligand states, "
                "protein states, hydrogen evidence, water evidence, state matrix, "
                "and robustness record"
            ),
            detail=(
                "the thermodynamic interpretation references prior records "
                "without modifying them"
            ),
        ),
        _check(
            "thermodynamics",
            "thermodynamics.free_energy_not_computed",
            True,
            subject=accession,
            criterion="no binding free energy is estimated from structural evidence",
            observed="free_energy_status=not_computed; value is null",
            expected="no numeric affinity or free-energy estimate",
            detail=(
                "interaction counts, geometric candidates, and robustness "
                "statuses are not thermodynamic quantities"
            ),
        ),
        _check(
            "thermodynamics",
            "thermodynamics.ensemble_not_available",
            True,
            subject=accession,
            criterion="the crystal structure is not treated as an equilibrium ensemble",
            observed="ensemble_status=not_available",
            expected="explicit sampling and a partition function for populations",
            detail="a structural snapshot is not a Boltzmann distribution",
        ),
        _check(
            "thermodynamics",
            "thermodynamics.terms_not_estimated",
            all(item.status == "not_estimated" for item in terms),
            subject=accession,
            criterion="missing thermodynamic terms remain unestimated",
            observed=f"{len(terms)} terms",
            expected="interaction, desolvation, entropy, water, strain, and standard-state terms",
            detail="no term is filled with an illustrative or fitted number",
        ),
        _check(
            "thermodynamics",
            "thermodynamics.requirements_incomplete",
            any(item.status != "met" for item in requirements),
            subject=accession,
            criterion="unmet free-energy requirements are reported",
            observed=f"{len(requirements)} requirements",
            expected="at least one preparation, sampling, or validation requirement",
            detail="missing requirements are not silently treated as satisfied",
        ),
        _check(
            "interpretation",
            "thermodynamics.no_affinity_ranking",
            True,
            subject=accession,
            criterion="states and interactions are not ranked by affinity",
            observed="no rank, score, or affinity field is produced",
            expected="structural and thermodynamic boundaries only",
            detail="interaction counts are not binding free energies",
        ),
    )

    decisions = (
        DecisionRecord(
            decision="did not convert interaction counts into a binding free energy",
            actor="thermodynamic interpretation",
            rationale=(
                "counts lack energies, desolvation, entropy, populations, and a "
                "defined thermodynamic path"
            ),
        ),
        DecisionRecord(
            decision="kept structural evidence and thermodynamic terms separate",
            actor="thermodynamic interpretation",
            rationale=(
                "a geometric classification is not an energy or a free-energy "
                "contribution"
            ),
        ),
        DecisionRecord(
            decision=(
                "treated the crystal structure as one observation, not an "
                "equilibrium ensemble"
            ),
            actor="thermodynamic interpretation",
            rationale="populations require a partition function and adequate sampling",
        ),
        DecisionRecord(
            decision="reported unmet free-energy requirements explicitly",
            actor="thermodynamic interpretation",
            rationale=(
                "a defensible affinity calculation needs a model, parameters, "
                "sampling, and uncertainty"
            ),
        ),
        DecisionRecord(
            decision="did not rank ligand or protein states",
            actor="thermodynamic interpretation",
            rationale="state hypotheses are not populations or free energies",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "deposited structure and waters",
                "Structural observations and water records retain deposited provenance.",
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "structural and robustness summary",
                (
                    "Interaction counts, geometric candidates, and robustness "
                    "statuses are summarized from prior records."
                ),
                actor="thermodynamic interpretation",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No affinity, population, or free-energy estimate is produced.",
                actor="thermodynamic interpretation",
            ),
            ProvenanceEntry(
                "hypothesised",
                "chemical states and thermodynamic requirements",
                (
                    "Chemical-state hypotheses and the requirements checklist "
                    "are explicit assumptions, not measured thermodynamics."
                ),
                actor="thermodynamic interpretation",
            ),
        )
    )

    discrepancies = (
        ThermodynamicDiscrepancy(
            category="ensemble",
            status="not_applicable",
            detail="No equilibrium ensemble or partition function is available.",
            observed="one deposited structural observation",
            expected="sampled equilibrium distribution",
            resolution="Populations remain unestimated.",
        ),
        ThermodynamicDiscrepancy(
            category="energetics",
            status="not_applicable",
            detail="No force field or energy model has been assigned.",
            observed="geometric classifications only",
            expected="interaction and desolvation energies",
            resolution="All free-energy terms remain not_estimated.",
        ),
        ThermodynamicDiscrepancy(
            category="entropy",
            status="not_applicable",
            detail="No conformational-entropy estimate is available.",
            observed="single structural model",
            expected="free and bound conformational ensembles",
            resolution="Entropic terms remain unestimated.",
        ),
        ThermodynamicDiscrepancy(
            category="water_populations",
            status="not_applicable",
            detail=(
                "Crystallographic occupancy is not an equilibrium probability "
                "and water roles remain unassigned."
            ),
            observed="deposited water occupancy and B-factor",
            expected="equilibrium water populations and role assignments",
            resolution="Water reorganisation remains unestimated.",
        ),
        ThermodynamicDiscrepancy(
            category="standard_state",
            status="not_applicable",
            detail="No standard state, temperature, or concentration convention is fixed.",
            observed="not defined",
            expected="explicit thermodynamic standard state",
            resolution="No free-energy comparison is made.",
        ),
    )

    return ThermodynamicInterpretationRecord(
        source_manifest_path=source_manifest_path,
        source_manifest_sha256=source_manifest_sha256,
        ligand_identity_record_path=ligand_identity_record_path,
        ligand_identity_record_sha256=ligand_identity_record_sha256,
        ligand_coordinate_record_path=ligand_coordinate_record_path,
        ligand_coordinate_record_sha256=ligand_coordinate_record_sha256,
        binding_site_record_path=binding_site_record_path,
        binding_site_record_sha256=binding_site_record_sha256,
        chemical_state_record_path=chemical_state_record_path,
        chemical_state_record_sha256=chemical_state_record_sha256,
        protein_preparation_record_path=protein_preparation_record_path,
        protein_preparation_record_sha256=protein_preparation_record_sha256,
        hydrogen_evidence_record_path=hydrogen_evidence_record_path,
        hydrogen_evidence_record_sha256=hydrogen_evidence_record_sha256,
        water_evidence_record_path=water_evidence_record_path,
        water_evidence_record_sha256=water_evidence_record_sha256,
        combination_interaction_record_path=combination_interaction_record_path,
        combination_interaction_record_sha256=(
            combination_interaction_record_sha256
        ),
        interaction_robustness_record_path=interaction_robustness_record_path,
        interaction_robustness_record_sha256=(
            interaction_robustness_record_sha256
        ),
        accession=accession,
        protein_selection=protein_selection,
        ligand_selection=ligand_selection,
        reference_ligand_state_id=reference_ligand_state_id,
        reference_protein_state_id=reference_protein_state_id,
        structural_summary=structural_summary,
        ensemble_status="not_available",
        ensemble_detail=(
            "The deposited crystal structure is one structural observation, not "
            "an equilibrium Boltzmann ensemble."
        ),
        free_energy_status="not_computed",
        binding_free_energy_kj_per_mol=None,
        boundaries=boundaries,
        terms=terms,
        requirements=requirements,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=discrepancies,
    )


def write_thermodynamic_interpretation(
    path: str | Path,
    record: ThermodynamicInterpretationRecord,
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_thermodynamic_interpretation(
    path: str | Path,
) -> ThermodynamicInterpretationRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return ThermodynamicInterpretationRecord.from_dict(data)


def validate_thermodynamic_interpretation(
    record: ThermodynamicInterpretationRecord,
) -> list[str]:
    errors: list[str] = []
    if record.schema_version != THERMODYNAMIC_INTERPRETATION_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be "
            f"{THERMODYNAMIC_INTERPRETATION_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    for name in (
        "source_manifest_path",
        "ligand_identity_record_path",
        "ligand_coordinate_record_path",
        "binding_site_record_path",
        "chemical_state_record_path",
        "protein_preparation_record_path",
        "hydrogen_evidence_record_path",
        "water_evidence_record_path",
        "combination_interaction_record_path",
        "interaction_robustness_record_path",
    ):
        if not getattr(record, name):
            errors.append(f"{name} is required.")
    for name in (
        "source_manifest_sha256",
        "ligand_identity_record_sha256",
        "ligand_coordinate_record_sha256",
        "binding_site_record_sha256",
        "chemical_state_record_sha256",
        "protein_preparation_record_sha256",
        "hydrogen_evidence_record_sha256",
        "water_evidence_record_sha256",
        "combination_interaction_record_sha256",
        "interaction_robustness_record_sha256",
    ):
        if not _SHA256_PATTERN.fullmatch(getattr(record, name)):
            errors.append(f"{name} must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")
    if not record.reference_ligand_state_id or not record.reference_protein_state_id:
        errors.append("reference ligand and protein state IDs are required.")
    if record.free_energy_status not in FREE_ENERGY_STATUSES:
        errors.append("free_energy_status must be 'not_computed' in this batch.")
    if record.binding_free_energy_kj_per_mol is not None:
        errors.append("binding_free_energy_kj_per_mol must be null in this batch.")
    if record.ensemble_status != "not_available":
        errors.append("ensemble_status must be 'not_available' in this batch.")

    boundary_ids = set()
    for index, boundary in enumerate(record.boundaries):
        if boundary.status not in BOUNDARY_STATUSES:
            errors.append(f"boundaries[{index}].status is invalid.")
        if not boundary.boundary_id or not boundary.structural_observation:
            errors.append(
                f"boundaries[{index}] requires boundary_id and structural_observation."
            )
        if boundary.boundary_id in boundary_ids:
            errors.append(f"boundaries[{index}].boundary_id is duplicated.")
        boundary_ids.add(boundary.boundary_id)
        if boundary.status in ("not_established", "not_estimated") and not (
            boundary.required_evidence
        ):
            errors.append(
                f"boundaries[{index}] requires evidence for an unestablished claim."
            )

    term_ids = set()
    for index, term in enumerate(record.terms):
        if term.status not in TERM_STATUSES:
            errors.append(f"terms[{index}].status is invalid.")
        if not term.term_id or not term.name:
            errors.append(f"terms[{index}] requires term_id and name.")
        if term.term_id in term_ids:
            errors.append(f"terms[{index}].term_id is duplicated.")
        term_ids.add(term.term_id)

    requirement_ids = set()
    for index, requirement in enumerate(record.requirements):
        if requirement.status not in REQUIREMENT_STATUSES:
            errors.append(f"requirements[{index}].status is invalid.")
        if not requirement.requirement_id or not requirement.name:
            errors.append(
                f"requirements[{index}] requires requirement_id and name."
            )
        if requirement.requirement_id in requirement_ids:
            errors.append(f"requirements[{index}].requirement_id is duplicated.")
        requirement_ids.add(requirement.requirement_id)

    if not record.boundaries:
        errors.append("boundaries must not be empty.")
    if not record.terms:
        errors.append("terms must not be empty.")
    if not record.requirements:
        errors.append("requirements must not be empty.")

    if not record.decisions:
        errors.append("decisions must contain at least one record.")
    for index, decision in enumerate(record.decisions):
        if not decision.decision or not decision.actor or not decision.rationale:
            errors.append(
                f"decisions[{index}] requires decision, actor, and rationale."
            )
    if not record.provenance.records:
        errors.append("provenance.records must not be empty.")
    for index, entry in enumerate(record.provenance.records):
        if entry.label not in PROVENANCE_LABELS:
            errors.append(f"provenance.records[{index}] has an invalid label.")

    for index, item in enumerate(record.discrepancies):
        if item.status not in ("open", "resolved", "not_applicable"):
            errors.append(f"discrepancies[{index}] has an invalid status.")

    check_ids = set()
    for index, check in enumerate(record.validation.checks):
        if check.status not in ("pass", "fail", "not_run", "not_applicable"):
            errors.append(f"validation.checks[{index}] has an invalid status.")
        if not check.check_id or not check.category or not check.detail:
            errors.append(
                f"validation.checks[{index}] requires check_id, category, and detail."
            )
        if check.check_id in check_ids:
            errors.append(f"validation.checks[{index}].check_id is duplicated.")
        check_ids.add(check.check_id)

    return errors
