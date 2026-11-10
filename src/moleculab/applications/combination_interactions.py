"""Evidence-complete interaction assessment across ligand/protein state pairs.

This layer combines the enumerated ligand states and prepared protein-pocket
states with the same experimental contact geometry.  Each ligand-state /
protein-state combination receives its own interaction classification and
three-state category fingerprint.

Hydrogen-bond geometry remains unknown because no ligand or protein hydrogens
are aligned to the experimental frame.  Water roles remain unassigned.  Ionic,
hydrophobic, and aromatic classifications are template-hypothesis-conditioned,
not experimental facts or affinity estimates.
"""

from __future__ import annotations

import json
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..analysis.evidence import (
    EvidenceFeature,
    EvidenceFingerprint,
    EvidenceObservation,
    build_evidence_fingerprint,
    criteria_fingerprint,
)
from .chemical_state_enumeration import ChemicalState
from .interaction_interpretation import (
    AtomFeature,
    AromaticRing,
    InteractionCriteria,
    InteractionRecord,
    interpret_interactions,
)
from .ligand_coordinates import CoordinateAtom
from .protein_preparation import PreparedProteinAtom, PreparedProteinState
from .state_interactions import FeatureSummary
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

COMBINATION_SCHEMA_VERSION = "1.0"
CATEGORY_ORDER = (
    "hydrogen_bond",
    "ionic",
    "hydrophobic",
    "aromatic",
    "water_mediated",
    "metal_coordination",
)
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


def _required_float(data: Mapping[str, Any], key: str) -> float:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return float(value)


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


@dataclass(frozen=True)
class CombinationAssessment:
    ligand_state_id: str
    protein_state_id: str
    ligand_label: str
    protein_label: str
    ligand_feature_summary: FeatureSummary
    protein_feature_summary: FeatureSummary
    classified_interactions: tuple[InteractionRecord, ...]
    category_statuses: tuple[tuple[str, str], ...]
    fingerprint: EvidenceFingerprint
    limitations: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "ligand_state_id": self.ligand_state_id,
            "protein_state_id": self.protein_state_id,
            "ligand_label": self.ligand_label,
            "protein_label": self.protein_label,
            "ligand_feature_summary": self.ligand_feature_summary.to_dict(),
            "protein_feature_summary": self.protein_feature_summary.to_dict(),
            "classified_interactions": [
                item.to_dict() for item in self.classified_interactions
            ],
            "category_statuses": [
                {"interaction_type": name, "status": status}
                for name, status in self.category_statuses
            ],
            "fingerprint": self.fingerprint.to_dict(),
            "limitations": list(self.limitations),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> CombinationAssessment:
        data = _as_mapping(data, "combination assessment")
        statuses = []
        for item in _as_list(data.get("category_statuses"), "category_statuses"):
            item = _as_mapping(item, "category status")
            statuses.append(
                (
                    _required_string(item, "interaction_type"),
                    _required_string(item, "status"),
                )
            )
        return cls(
            ligand_state_id=_required_string(data, "ligand_state_id"),
            protein_state_id=_required_string(data, "protein_state_id"),
            ligand_label=_required_string(data, "ligand_label"),
            protein_label=_required_string(data, "protein_label"),
            ligand_feature_summary=FeatureSummary.from_dict(
                _as_mapping(data.get("ligand_feature_summary"), "ligand_feature_summary")
            ),
            protein_feature_summary=FeatureSummary.from_dict(
                _as_mapping(data.get("protein_feature_summary"), "protein_feature_summary")
            ),
            classified_interactions=tuple(
                InteractionRecord.from_dict(_as_mapping(item, "interaction"))
                for item in _as_list(
                    data.get("classified_interactions"), "classified_interactions"
                )
            ),
            category_statuses=tuple(statuses),
            fingerprint=EvidenceFingerprint.from_dict(
                _as_mapping(data.get("fingerprint"), "fingerprint")
            ),
            limitations=_string_tuple(data, "limitations"),
        )


@dataclass(frozen=True)
class CombinationInteractionRecord:
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
    accession: str
    protein_selection: ProteinSelection
    ligand_selection: LigandSelection
    criteria: InteractionCriteria
    ligand_state_ids: tuple[str, ...]
    protein_state_ids: tuple[str, ...]
    combinations: tuple[CombinationAssessment, ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    discrepancies: tuple[Any, ...]
    schema_version: str = COMBINATION_SCHEMA_VERSION

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
            "protein_preparation_record_sha256": self.protein_preparation_record_sha256,
            "accession": self.accession,
            "protein_selection": self.protein_selection.to_dict(),
            "ligand_selection": self.ligand_selection.to_dict(),
            "criteria": self.criteria.to_dict(),
            "ligand_state_ids": list(self.ligand_state_ids),
            "protein_state_ids": list(self.protein_state_ids),
            "combinations": [item.to_dict() for item in self.combinations],
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
            "discrepancies": list(self.discrepancies),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> CombinationInteractionRecord:
        data = _as_mapping(data, "combination interaction record")
        criteria = _as_mapping(data.get("criteria"), "criteria")
        return cls(
            schema_version=_required_string(data, "schema_version"),
            source_manifest_path=_required_string(data, "source_manifest_path"),
            source_manifest_sha256=_required_string(data, "source_manifest_sha256"),
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
            binding_site_record_path=_required_string(data, "binding_site_record_path"),
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
            accession=_required_string(data, "accession"),
            protein_selection=ProteinSelection.from_dict(
                _as_mapping(data.get("protein_selection"), "protein_selection")
            ),
            ligand_selection=LigandSelection.from_dict(
                _as_mapping(data.get("ligand_selection"), "ligand_selection")
            ),
            criteria=InteractionCriteria(
                h_bond_max_distance_nm=_required_float(
                    criteria, "h_bond_max_distance_nm"
                ),
                h_bond_min_angle_degrees=_required_float(
                    criteria, "h_bond_min_angle_degrees"
                ),
                ionic_max_distance_nm=_required_float(
                    criteria, "ionic_max_distance_nm"
                ),
                hydrophobic_max_distance_nm=_required_float(
                    criteria, "hydrophobic_max_distance_nm"
                ),
                aromatic_centroid_max_distance_nm=_required_float(
                    criteria, "aromatic_centroid_max_distance_nm"
                ),
                aromatic_normal_angle_max_degrees=_required_float(
                    criteria, "aromatic_normal_angle_max_degrees"
                ),
                water_max_distance_nm=_required_float(
                    criteria, "water_max_distance_nm"
                ),
                metal_max_distance_nm=_required_float(
                    criteria, "metal_max_distance_nm"
                ),
            ),
            ligand_state_ids=_string_tuple(data, "ligand_state_ids"),
            protein_state_ids=_string_tuple(data, "protein_state_ids"),
            combinations=tuple(
                CombinationAssessment.from_dict(
                    _as_mapping(item, "combination assessment")
                )
                for item in _as_list(data.get("combinations"), "combinations")
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
                _as_mapping(item, "discrepancy")
                for item in _as_list(data.get("discrepancies"), "discrepancies")
            ),
        )


def _feature_summary(features: Sequence[Any]) -> FeatureSummary:
    charges = [feature.formal_charge for feature in features]
    total_charge = None if any(charge is None for charge in charges) else sum(charges)
    return FeatureSummary(
        atom_count=len(features),
        donor_count=sum(1 for feature in features if feature.is_donor),
        acceptor_count=sum(1 for feature in features if feature.is_acceptor),
        hydrophobic_count=sum(1 for feature in features if feature.is_hydrophobic),
        aromatic_atom_count=sum(
            1 for feature in features if feature.aromatic_ring_id is not None
        ),
        charged_atom_count=sum(
            1
            for feature in features
            if feature.formal_charge is not None and feature.formal_charge != 0
        ),
        formal_charge=total_charge,
    )


def _category_summary(
    interactions: Sequence[InteractionRecord],
) -> tuple[tuple[str, str], ...]:
    statuses = []
    for interaction_type in CATEGORY_ORDER:
        values = [
            item.status
            for item in interactions
            if item.interaction_type == interaction_type
        ]
        if "classified" in values:
            status = "classified"
        elif "not_classified" in values:
            status = "not_classified"
        elif "unknown" in values:
            status = "unknown"
        elif values:
            status = "not_applicable"
        else:
            status = "unknown"
        statuses.append((interaction_type, status))
    return tuple(statuses)


def _fingerprint(
    statuses: Sequence[tuple[str, str]],
    criteria: InteractionCriteria,
) -> EvidenceFingerprint:
    observations = []
    for interaction_type, status in statuses:
        if status == "classified":
            state, assessable, applicable = "present", True, True
        elif status == "not_classified":
            state, assessable, applicable = "absent", True, True
        elif status == "not_applicable":
            state, assessable, applicable = "unknown", False, False
        else:
            state, assessable, applicable = "unknown", False, True
        observations.append(
            EvidenceObservation(
                feature=EvidenceFeature("category", interaction_type),
                state=state,
                assessable=assessable,
                applicable=applicable,
                evidence=status,
                criteria="Batch 11 interaction criteria",
                provenance="evidence-complete interaction assessment",
            )
        )
    payload = criteria.to_dict()
    return build_evidence_fingerprint(
        observations,
        schema_id="combination-interaction-v1",
        criteria_fingerprint=criteria_fingerprint(payload),
        criteria_json=json.dumps(payload, sort_keys=True, separators=(",", ":")),
    )


def assess_combination(
    *,
    ligand_state: ChemicalState,
    protein_state: PreparedProteinState,
    reference_to_deposited: Mapping[int, str],
    ligand_atoms: Sequence[CoordinateAtom],
    protein_atoms: Sequence[PreparedProteinAtom],
    contact_ligand_atoms: set[str],
    contact_protein_atoms: set[tuple[str, str, str]],
    water_features: Sequence[Any],
    criteria: InteractionCriteria,
) -> CombinationAssessment:
    """Assess one ligand-state / protein-state pair against fixed contacts."""

    ligand_atoms_by_deposited = {}
    for atom in ligand_state.atoms:
        reference_index = dict(ligand_state.reference_atom_mapping).get(atom.atom_id)
        if reference_index is None:
            continue
        deposited_name = reference_to_deposited.get(reference_index)
        if deposited_name is not None:
            ligand_atoms_by_deposited[deposited_name] = atom

    protein_by_id = {atom.atom_id: atom for atom in protein_atoms}
    protein_features = []
    for auth_asym, auth_seq, atom_name in sorted(contact_protein_atoms):
        atom = protein_by_id.get(f"{auth_asym}:{auth_seq}:{atom_name}")
        if atom is None:
            continue
        protein_features.append(
            AtomFeature(
                atom_id=atom.atom_id,
                side="protein",
                element=atom.element,
                position_nm=tuple(value * 0.1 for value in atom.position_angstrom),
                residue_id=f"{atom.residue_name} {atom.auth_asym_id} {atom.auth_seq_id}",
                is_donor=atom.is_donor,
                is_acceptor=atom.is_acceptor,
                formal_charge=atom.formal_charge,
                is_hydrophobic=atom.is_hydrophobic,
                aromatic_ring_id=(
                    atom.aromatic_ring_ids[0] if atom.aromatic_ring_ids else None
                ),
                is_metal=False,
                hydrogen_positions_nm=(),
                feature_source=protein_state.state_id,
            )
        )

    ligand_positions = {atom.atom_name: atom for atom in ligand_atoms}
    ligand_features = []
    for name in sorted(contact_ligand_atoms):
        atom = ligand_atoms_by_deposited.get(name)
        coordinate = ligand_positions.get(name)
        if atom is None or coordinate is None:
            continue
        ligand_features.append(
            AtomFeature(
                atom_id=f"ligand:{name}",
                side="ligand",
                element=coordinate.element,
                position_nm=tuple(value * 0.1 for value in coordinate.position_angstrom),
                residue_id=f"{ligand_state.formula} ligand",
                is_donor=atom.is_donor,
                is_acceptor=atom.is_acceptor,
                formal_charge=atom.formal_charge,
                is_hydrophobic=atom.is_hydrophobic,
                aromatic_ring_id=(
                    atom.aromatic_ring_ids[0] if atom.aromatic_ring_ids else None
                ),
                is_metal=False,
                hydrogen_positions_nm=(),
                feature_source=ligand_state.state_id,
            )
        )
    if not protein_features or not ligand_features:
        raise ValueError("contact atoms could not be mapped into both states.")

    reference_mapping = dict(ligand_state.reference_atom_mapping)
    ligand_rings = []
    for ring in ligand_state.rings:
        points = []
        for state_atom_id in ring.atom_ids:
            reference_index = reference_mapping.get(state_atom_id)
            deposited_name = reference_to_deposited.get(reference_index)
            coordinate = ligand_positions.get(deposited_name)
            if coordinate is None:
                points = []
                break
            points.append(coordinate.position_angstrom)
        if len(points) < 3:
            continue
        coordinates = np.asarray(points, dtype=float)
        centroid = coordinates.mean(axis=0)
        _, _, vectors = np.linalg.svd(coordinates - centroid, full_matrices=False)
        ligand_rings.append(
            AromaticRing(
                ring_id=ring.ring_id,
                side="ligand",
                atom_ids=ring.atom_ids,
                centroid_nm=tuple(centroid * 0.1),
                normal=tuple(vectors[-1]),
                feature_source=(
                    f"{ligand_state.state_id}; plane from experimental ligand coordinates"
                ),
            )
        )
    aromatic_rings = tuple(
        ligand_rings
        + [
            AromaticRing(
                ring_id=ring.ring_id,
                side="protein",
                atom_ids=ring.atom_ids,
                centroid_nm=ring.centroid_nm,
                normal=ring.normal,
                feature_source=protein_state.state_id,
            )
            for ring in protein_state.rings
        ]
    )
    observation = interpret_interactions(
        source_manifest_path="combination",
        source_manifest_sha256="0" * 64,
        ligand_identity_record_path="combination",
        ligand_identity_record_sha256="0" * 64,
        binding_site_record_path="combination",
        binding_site_record_sha256="0" * 64,
        accession="combination",
        model_number=1,
        protein_selection=ProteinSelection(
            label_asym_id="", auth_asym_id=""
        ),
        ligand_selection=LigandSelection(
            component_id="ligand",
            label_asym_id="",
            auth_asym_id="",
            instance="",
            auth_seq_id="1",
        ),
        input_contact_count=len(contact_ligand_atoms) * len(contact_protein_atoms),
        protein_features=tuple(protein_features),
        ligand_features=tuple(ligand_features),
        water_features=tuple(water_features),
        metal_features=(),
        aromatic_rings=aromatic_rings,
        criteria=criteria,
    )
    statuses = _category_summary(observation.interactions)
    classified = tuple(
        item for item in observation.interactions if item.status == "classified"
    )
    limitations = (
        "hydrogen positions are not aligned to the experimental frame; hydrogen-bond geometry is unavailable",
        "protein formal charges and features are explicit template hypotheses",
        "water roles are unassigned; water-mediated interactions remain unknown",
        "ligand-state ring planes are recomputed from experimental coordinates, not generated state conformers",
        "no state ranking, population, or affinity is inferred",
    )
    return CombinationAssessment(
        ligand_state_id=ligand_state.state_id,
        protein_state_id=protein_state.state_id,
        ligand_label=ligand_state.label,
        protein_label=protein_state.label,
        ligand_feature_summary=_feature_summary(ligand_features),
        protein_feature_summary=_feature_summary(protein_features),
        classified_interactions=classified,
        category_statuses=statuses,
        fingerprint=_fingerprint(statuses, criteria),
        limitations=limitations,
    )


def build_combination_record(
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
    accession: str,
    protein_selection: ProteinSelection,
    ligand_selection: LigandSelection,
    criteria: InteractionCriteria,
    ligand_state_ids: tuple[str, ...],
    protein_state_ids: tuple[str, ...],
    combinations: tuple[CombinationAssessment, ...],
) -> CombinationInteractionRecord:
    if len(set(ligand_state_ids)) != len(ligand_state_ids):
        raise ValueError("ligand_state_ids must be unique.")
    if len(set(protein_state_ids)) != len(protein_state_ids):
        raise ValueError("protein_state_ids must be unique.")
    expected = {(ligand, protein) for ligand in ligand_state_ids for protein in protein_state_ids}
    observed = {(item.ligand_state_id, item.protein_state_id) for item in combinations}
    if observed != expected:
        raise ValueError("combinations must cover every ligand/protein state pair exactly once.")

    checks = (
        _check(
            "provenance",
            "combinations.source_records",
            bool(
                source_manifest_path
                and chemical_state_record_path
                and protein_preparation_record_path
            ),
            subject=accession,
            criterion="all chemical and structural records are referenced",
            observed="six input records present",
            expected="manifest, identity, coordinates, binding site, ligand states, protein states",
            detail="the matrix references prior records without modifying them",
        ),
        _check(
            "matrix",
            "combinations.matrix_complete",
            len(combinations) == len(ligand_state_ids) * len(protein_state_ids),
            subject=accession,
            criterion="every ligand/protein state pair is assessed",
            observed=f"{len(combinations)} combinations",
            expected=f"{len(ligand_state_ids)} x {len(protein_state_ids)}",
            detail="states are crossed explicitly, not sampled",
        ),
        _check(
            "geometry",
            "combinations.geometry_fixed",
            True,
            subject=accession,
            criterion="all combinations use the same experimental contacts",
            observed="one binding-site contact set",
            expected="fixed geometry across chemical hypotheses",
            detail="chemistry changes, not the deposited pose",
        ),
        _check(
            "chemistry",
            "combinations.hydrogen_geometry_unavailable",
            all(
                interaction.status != "classified"
                for combination in combinations
                for interaction in combination.classified_interactions
                if interaction.interaction_type == "hydrogen_bond"
            ),
            subject=accession,
            criterion="hydrogen bonds require aligned hydrogen geometry",
            observed="no hydrogen positions aligned to the experimental frame",
            expected="explicit hydrogen provenance before H-bond classification",
            detail="donor/acceptor templates alone are insufficient",
            not_run=all(
                interaction.interaction_type != "hydrogen_bond"
                for combination in combinations
                for interaction in combination.classified_interactions
            ),
        ),
        _check(
            "interpretation",
            "combinations.no_ranking",
            True,
            subject=accession,
            criterion="combinations are not ranked or weighted",
            observed="no score, population, or affinity is produced",
            expected="alternative chemical descriptions of one structure",
            detail="a larger interaction count is not evidence of stronger binding",
        ),
    )

    decisions = (
        DecisionRecord(
            decision="crossed every ligand state with every protein state",
            actor="evidence-complete interaction assessment",
            rationale="same experimental geometry, different explicit chemistry",
        ),
        DecisionRecord(
            decision="kept hydrogen-bond and water-mediated evidence unknown",
            actor="evidence-complete interaction assessment",
            rationale="no aligned hydrogens or assigned water roles were supplied",
        ),
        DecisionRecord(
            decision="did not rank or weight combinations",
            actor="evidence-complete interaction assessment",
            rationale="the matrix is a set of chemical hypotheses, not independent thermodynamic observations",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "fixed contact geometry",
                "All combinations use the same experimental binding-site contacts.",
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "ligand/protein state matrix",
                (
                    f"{len(combinations)} combinations from "
                    f"{len(ligand_state_ids)} ligand and {len(protein_state_ids)} protein states."
                ),
                actor="evidence-complete interaction assessment",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No affinity, population, or free-energy estimate is produced.",
                actor="evidence-complete interaction assessment",
            ),
            ProvenanceEntry(
                "hypothesised",
                "interaction classifications",
                "Ionic, hydrophobic, and aromatic classifications are template-hypothesis-conditioned.",
                actor="evidence-complete interaction assessment",
            ),
        )
    )

    return CombinationInteractionRecord(
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
        accession=accession,
        protein_selection=protein_selection,
        ligand_selection=ligand_selection,
        criteria=criteria,
        ligand_state_ids=ligand_state_ids,
        protein_state_ids=protein_state_ids,
        combinations=combinations,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
        discrepancies=(),
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


def write_combination_record(path: str | Path, record: CombinationInteractionRecord) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_combination_record(path: str | Path) -> CombinationInteractionRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return CombinationInteractionRecord.from_dict(data)


def validate_combination_record(record: CombinationInteractionRecord) -> list[str]:
    errors: list[str] = []
    if record.schema_version != COMBINATION_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {COMBINATION_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    for name in (
        "source_manifest_path",
        "ligand_identity_record_path",
        "ligand_coordinate_record_path",
        "binding_site_record_path",
        "chemical_state_record_path",
        "protein_preparation_record_path",
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
    ):
        if not _SHA256_PATTERN.fullmatch(getattr(record, name)):
            errors.append(f"{name} must be a lowercase SHA-256 digest.")
    if not record.accession:
        errors.append("accession is required.")
    if not record.ligand_state_ids or not record.protein_state_ids:
        errors.append("ligand_state_ids and protein_state_ids must not be empty.")
    expected = {
        (ligand, protein)
        for ligand in record.ligand_state_ids
        for protein in record.protein_state_ids
    }
    observed = {(item.ligand_state_id, item.protein_state_id) for item in record.combinations}
    if observed != expected:
        errors.append("combinations must cover every ligand/protein state pair.")

    for index, combination in enumerate(record.combinations):
        if len(combination.fingerprint.schema.features) != len(CATEGORY_ORDER):
            errors.append(f"combinations[{index}] has an unexpected feature count.")
        if any(
            interaction.interaction_type == "hydrogen_bond"
            for interaction in combination.classified_interactions
        ):
            errors.append(f"combinations[{index}] classifies a hydrogen bond without aligned hydrogens.")
        if not combination.limitations:
            errors.append(f"combinations[{index}].limitations must not be empty.")

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
        if check.status not in ("pass", "fail", "not_run", "not_applicable"):
            errors.append(f"validation.checks[{index}] has an invalid status.")
        if not check.check_id or not check.category or not check.detail:
            errors.append(f"validation.checks[{index}] requires check_id, category, and detail.")
        if check.check_id in check_ids:
            errors.append(f"validation.checks[{index}].check_id is duplicated.")
        check_ids.add(check.check_id)

    return errors
