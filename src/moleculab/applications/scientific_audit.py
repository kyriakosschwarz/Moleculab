"""Final capstone audit and an end-to-end statistical-mechanics demonstration.

The audit checks the provenance chain, the structural/thermodynamic boundaries,
the harmonic free-energy study, documentation links, and reproducibility
metadata.  It does not claim that passing the audit proves scientific
correctness; it records which invariants were checked.

The demonstration connects one sampled configuration to a potential energy, a
Boltzmann weight, an ensemble average, and an exact/FEP free-energy difference.
It then stops before any protein-ligand affinity claim.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..foundation import units
from ..models.potentials import HarmonicPotential
from .fep import fep_harmonic_oscillator, harmonic_free_energy_difference
from .free_energy_study import (
    FEPOverlapDiagnostics,
    TwoStateReference,
    fep_overlap_diagnostics,
    harmonic_configurational_partition_function,
    two_state_reference,
)
from .structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    ProvenanceEntry,
    ProvenanceRecord,
    sha256_file,
)

SCIENTIFIC_AUDIT_SCHEMA_VERSION = "1.0"
AUDIT_STATUSES = ("pass", "fail", "not_run")


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
class AuditCheck:
    check_id: str
    category: str
    status: str
    subject: str
    criterion: str
    observed: str
    detail: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "check_id": self.check_id,
            "category": self.category,
            "status": self.status,
            "subject": self.subject,
            "criterion": self.criterion,
            "observed": self.observed,
            "detail": self.detail,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> AuditCheck:
        data = _as_mapping(data, "audit check")
        return cls(
            check_id=_required_string(data, "check_id"),
            category=_required_string(data, "category"),
            status=_required_string(data, "status"),
            subject=_required_string(data, "subject"),
            criterion=_required_string(data, "criterion"),
            observed=_required_string(data, "observed"),
            detail=_required_string(data, "detail"),
        )


@dataclass(frozen=True)
class EndToEndDemonstration:
    temperature_kelvin: float
    sample_count: int
    mean_potential_energy_kj_per_mol: float
    analytic_mean_potential_energy_kj_per_mol: float
    proposal_effective_sample_fraction: float
    proposal_maximum_weight_fraction: float
    exact_delta_g_kj_per_mol: float
    fep_delta_g_kj_per_mol: float
    fep_standard_error_kj_per_mol: float
    overlap: FEPOverlapDiagnostics
    two_state_reference: TwoStateReference
    stop_line: tuple[str, ...]

    def to_dict(self) -> dict[str, Any]:
        return {
            "temperature_kelvin": self.temperature_kelvin,
            "sample_count": self.sample_count,
            "mean_potential_energy_kj_per_mol": (
                self.mean_potential_energy_kj_per_mol
            ),
            "analytic_mean_potential_energy_kj_per_mol": (
                self.analytic_mean_potential_energy_kj_per_mol
            ),
            "proposal_effective_sample_fraction": (
                self.proposal_effective_sample_fraction
            ),
            "proposal_maximum_weight_fraction": (
                self.proposal_maximum_weight_fraction
            ),
            "exact_delta_g_kj_per_mol": self.exact_delta_g_kj_per_mol,
            "fep_delta_g_kj_per_mol": self.fep_delta_g_kj_per_mol,
            "fep_standard_error_kj_per_mol": (
                self.fep_standard_error_kj_per_mol
            ),
            "overlap": self.overlap.to_dict(),
            "two_state_reference": self.two_state_reference.to_dict(),
            "stop_line": list(self.stop_line),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> EndToEndDemonstration:
        data = _as_mapping(data, "end-to-end demonstration")
        return cls(
            temperature_kelvin=_required_float(data, "temperature_kelvin"),
            sample_count=_required_int(data, "sample_count"),
            mean_potential_energy_kj_per_mol=_required_float(
                data, "mean_potential_energy_kj_per_mol"
            ),
            analytic_mean_potential_energy_kj_per_mol=_required_float(
                data, "analytic_mean_potential_energy_kj_per_mol"
            ),
            proposal_effective_sample_fraction=_required_float(
                data, "proposal_effective_sample_fraction"
            ),
            proposal_maximum_weight_fraction=_required_float(
                data, "proposal_maximum_weight_fraction"
            ),
            exact_delta_g_kj_per_mol=_required_float(
                data, "exact_delta_g_kj_per_mol"
            ),
            fep_delta_g_kj_per_mol=_required_float(
                data, "fep_delta_g_kj_per_mol"
            ),
            fep_standard_error_kj_per_mol=_required_float(
                data, "fep_standard_error_kj_per_mol"
            ),
            overlap=FEPOverlapDiagnostics.from_dict(
                _as_mapping(data.get("overlap"), "overlap")
            ),
            two_state_reference=TwoStateReference.from_dict(
                _as_mapping(
                    data.get("two_state_reference"),
                    "two_state_reference",
                )
            ),
            stop_line=_string_tuple(data, "stop_line"),
        )


@dataclass(frozen=True)
class ScientificAuditRecord:
    audited_artifacts: tuple[str, ...]
    checks: tuple[AuditCheck, ...]
    pass_count: int
    fail_count: int
    not_run_count: int
    status: str
    demonstration: EndToEndDemonstration
    limitations: tuple[str, ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    schema_version: str = SCIENTIFIC_AUDIT_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "audited_artifacts": list(self.audited_artifacts),
            "checks": [item.to_dict() for item in self.checks],
            "pass_count": self.pass_count,
            "fail_count": self.fail_count,
            "not_run_count": self.not_run_count,
            "status": self.status,
            "demonstration": self.demonstration.to_dict(),
            "limitations": list(self.limitations),
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> ScientificAuditRecord:
        data = _as_mapping(data, "scientific audit record")
        return cls(
            schema_version=_required_string(data, "schema_version"),
            audited_artifacts=_string_tuple(data, "audited_artifacts"),
            checks=tuple(
                AuditCheck.from_dict(_as_mapping(item, "audit check"))
                for item in _as_list(data.get("checks"), "checks")
            ),
            pass_count=_required_int(data, "pass_count"),
            fail_count=_required_int(data, "fail_count"),
            not_run_count=_required_int(data, "not_run_count"),
            status=_required_string(data, "status"),
            demonstration=EndToEndDemonstration.from_dict(
                _as_mapping(data.get("demonstration"), "demonstration")
            ),
            limitations=_string_tuple(data, "limitations"),
            decisions=tuple(
                DecisionRecord.from_dict(_as_mapping(item, "decision record"))
                for item in _as_list(data.get("decisions"), "decisions")
            ),
            provenance=ProvenanceRecord.from_dict(
                _as_mapping(data.get("provenance"), "provenance")
            ),
        )


def _check(
    check_id: str,
    category: str,
    ok: bool,
    *,
    subject: str,
    criterion: str,
    observed: str,
    detail: str,
) -> AuditCheck:
    return AuditCheck(
        check_id=check_id,
        category=category,
        status="pass" if ok else "fail",
        subject=subject,
        criterion=criterion,
        observed=observed,
        detail=detail,
    )


def _broken_relative_links(paths: Sequence[Path]) -> tuple[str, ...]:
    broken = []
    for path in paths:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        for match in re.finditer(r"\[[^\]]+\]\(([^)]+)\)", text):
            target = match.group(1)
            if target.startswith(("http://", "https://", "#", "mailto:")):
                continue
            target = target.split("#", 1)[0]
            if not target:
                continue
            if not (
                target.endswith((".md", ".json", ".svg", ".py", ".cif", ".txt"))
                or target.startswith(("../", "./"))
            ):
                continue
            if not (path.parent / target).resolve().exists():
                broken.append(f"{path}:{target}")
    return tuple(broken)


def run_end_to_end_demonstration(
    *,
    temperature_kelvin: float = 300.0,
    n_samples: int = 50_000,
    seed: int = 2024,
) -> EndToEndDemonstration:
    """Run the controlled configuration -> Boltzmann -> free-energy chain."""

    temperature = float(temperature_kelvin)
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("temperature_kelvin must be finite and positive.")
    if n_samples < 2:
        raise ValueError("n_samples must be at least 2.")

    potential0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    potential1 = HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0])
    # A uniform proposal makes the Boltzmann weights visible.  Re-weighting a
    # proposal that is already the target distribution would double-count it.
    rng = np.random.default_rng(seed)
    positions = rng.uniform(-0.75, 0.75, size=(n_samples, 3))
    energies = np.array(
        [potential0.energy(position) for position in positions],
        dtype=float,
    )
    beta = 1.0 / (units.MOLAR_GAS_CONSTANT * temperature)
    log_weights = -beta * energies
    weights = np.exp(log_weights - np.max(log_weights))
    probabilities = weights / np.sum(weights)
    mean_energy = float(np.sum(probabilities * energies))
    analytic_mean = float(1.5 * units.MOLAR_GAS_CONSTANT * temperature)
    proposal_effective_sample_size = 1.0 / float(
        np.sum(probabilities * probabilities)
    )
    proposal_effective_fraction = proposal_effective_sample_size / n_samples
    proposal_maximum_weight = float(np.max(probabilities))

    exact_delta_g = harmonic_free_energy_difference(
        potential0, potential1, temperature
    )
    fep_result = fep_harmonic_oscillator(
        potential0,
        potential1,
        temperature,
        n_samples=n_samples,
        seed=seed + 1,
    )
    overlap = fep_overlap_diagnostics(fep_result.work, temperature)
    two_state = two_state_reference(
        energy0_kj_per_mol=0.0,
        energy1_kj_per_mol=5.0,
        degeneracy0=1,
        degeneracy1=100,
        temperature=temperature,
    )

    return EndToEndDemonstration(
        temperature_kelvin=temperature,
        sample_count=n_samples,
        mean_potential_energy_kj_per_mol=mean_energy,
        analytic_mean_potential_energy_kj_per_mol=analytic_mean,
        proposal_effective_sample_fraction=proposal_effective_fraction,
        proposal_maximum_weight_fraction=proposal_maximum_weight,
        exact_delta_g_kj_per_mol=exact_delta_g,
        fep_delta_g_kj_per_mol=fep_result.delta_g,
        fep_standard_error_kj_per_mol=fep_result.standard_error,
        overlap=overlap,
        two_state_reference=two_state,
        stop_line=(
            (
                "The harmonic model has an exact Hamiltonian, partition "
                "function, and Boltzmann sampler."
            ),
            (
                "The protein-ligand record has no force-field-ready model, "
                "explicit solvent, ion environment, or protonation population."
            ),
            (
                "No conformational or alchemical sampling has been run for the "
                "complex."
            ),
            (
                "The structural evidence is not converted into a binding free "
                "energy."
            ),
        ),
    )


def run_capstone_audit(
    *,
    root: Path,
    manifest_path: Path,
    identity_path: Path,
    ligand_path: Path,
    site_path: Path,
    states_path: Path,
    protein_path: Path,
    hydrogen_path: Path,
    water_path: Path,
    combination_path: Path,
    robustness_path: Path,
    thermodynamics_path: Path,
    free_energy_study_path: Path,
    cif_path: Path,
    manifest: Any,
    identity: Any,
    ligand: Any,
    site: Any,
    states: Any,
    protein: Any,
    hydrogen: Any,
    water: Any,
    combination: Any,
    robustness: Any,
    thermodynamics: Any,
    free_energy_study: Any,
    demonstration: EndToEndDemonstration | None = None,
) -> ScientificAuditRecord:
    demonstration = (
        run_end_to_end_demonstration()
        if demonstration is None
        else demonstration
    )
    artifact_paths = (
        manifest_path,
        identity_path,
        ligand_path,
        site_path,
        states_path,
        protein_path,
        hydrogen_path,
        water_path,
        combination_path,
        robustness_path,
        thermodynamics_path,
        free_energy_study_path,
    )

    def _display_path(path: Path) -> str:
        try:
            return str(path.resolve().relative_to(root.resolve()))
        except ValueError:
            return str(path)
    required_docs = (
        root / "docs/tutorials/README.md",
        root / "docs/tutorials/physics/tutorial.md",
        root / "docs/tutorials/chemistry/README.md",
        root / "docs/tutorials/drug_discovery/README.md",
        root / "docs/tutorials/capstone_imatinib.md",
        root / "docs/reference/validation.md",
        root / "examples/README.md",
        root / "data/structures/README.md",
        root / "data/free_energy/README.md",
        root / "AGENTS.md",
    )
    documentation_paths = (
        *required_docs,
        root / "docs/development/architecture.md",
        root / "docs/tutorials/overview.md",
    )
    broken_links = _broken_relative_links(documentation_paths)

    manifest_sha = sha256_file(manifest_path)
    source_sha = str(manifest.source.sha256)
    source_size = int(manifest.source.file_size_bytes)
    source_checksum_ok = (
        sha256_file(cif_path) == source_sha and cif_path.stat().st_size == source_size
    )
    records = (
        ("identity", identity),
        ("ligand", ligand),
        ("site", site),
        ("states", states),
        ("protein", protein),
        ("hydrogen", hydrogen),
        ("water", water),
        ("combination", combination),
        ("robustness", robustness),
        ("thermodynamics", thermodynamics),
    )
    manifest_mismatches = [
        name
        for name, record in records
        if getattr(record, "source_manifest_sha256", None) != manifest_sha
    ]
    reference_chain = (
        ("ligand", ligand, "identity_record_sha256", identity_path),
        ("site", site, "ligand_coordinate_record_sha256", ligand_path),
        ("states", states, "binding_site_record_sha256", site_path),
        ("combination", combination, "chemical_state_record_sha256", states_path),
        ("combination", combination, "protein_preparation_record_sha256", protein_path),
        ("robustness", robustness, "hydrogen_evidence_record_sha256", hydrogen_path),
        ("robustness", robustness, "water_evidence_record_sha256", water_path),
        (
            "thermodynamics",
            thermodynamics,
            "combination_interaction_record_sha256",
            combination_path,
        ),
        (
            "thermodynamics",
            thermodynamics,
            "interaction_robustness_record_sha256",
            robustness_path,
        ),
    )
    reference_mismatches = [
        f"{name}.{field}"
        for name, record, field, path in reference_chain
        if getattr(record, field, None) != sha256_file(path)
    ]

    potential0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    potential1 = HarmonicPotential(stiffness=200.0, center=[0.2, 0.0, 0.0])
    temperature = 300.0
    exact = harmonic_free_energy_difference(potential0, potential1, temperature)
    expected_exact = (
        1.5
        * units.MOLAR_GAS_CONSTANT
        * temperature
        * math.log(potential1.stiffness / potential0.stiffness)
    )
    z0 = harmonic_configurational_partition_function(potential0, temperature)
    z1 = harmonic_configurational_partition_function(potential1, temperature)
    partition_delta = -units.MOLAR_GAS_CONSTANT * temperature * math.log(z1 / z0)
    gradient_position = np.array([0.31, -0.12, 0.27])
    h = 1e-6
    numeric_force = np.array(
        [
            -(
                potential0.energy(gradient_position + np.eye(3)[axis] * h)
                - potential0.energy(gradient_position - np.eye(3)[axis] * h)
            )
            / (2.0 * h)
            for axis in range(3)
        ]
    )
    force_error = float(
        np.max(np.abs(numeric_force - potential0.force(gradient_position)))
    )

    two_state = demonstration.two_state_reference
    two_state_delta_expected = (
        two_state.energy1_kj_per_mol
        - two_state.energy0_kj_per_mol
        - units.MOLAR_GAS_CONSTANT
        * demonstration.temperature_kelvin
        * math.log(two_state.degeneracy1 / two_state.degeneracy0)
    )
    case_by_id = {case.case_id: case for case in free_energy_study.cases}
    mild = case_by_id.get("mild-stiffness-change")
    poor = case_by_id.get("poor-overlap-center-shift")
    study_statuses = {item.status for item in robustness.conclusions}

    checks = (
        _check(
            "documentation.paths",
            "documentation",
            all(path.exists() for path in required_docs),
            subject="learning paths",
            criterion=(
                "all track, capstone, validation, and data documentation files "
                "exist"
            ),
            observed=(
                f"{sum(path.exists() for path in required_docs)}/"
                f"{len(required_docs)} present"
            ),
            detail=(
                "Track A, Track B, Track C, capstone, validation, examples, and "
                "data records are present"
            ),
        ),
        _check(
            "documentation.links",
            "documentation",
            not broken_links,
            subject="documentation links",
            criterion="relative links in the audited documents resolve",
            observed=f"{len(broken_links)} broken links",
            detail=(
                "checked documentation, examples, data, and AGENTS links; "
                "SMILES parentheses are ignored"
            ),
        ),
        _check(
            "documentation.capstone_batches",
            "documentation",
            all(
                marker in (root / "docs/tutorials/capstone_imatinib.md").read_text(
                    encoding="utf-8"
                )
                for marker in (
                    "4.13 Batch 11",
                    "4.14 Batch 12A",
                    "4.15 Batch 12B",
                    "4.16 Batch 13",
                    "4.17 Batch 14",
                    "4.18 Batch 15",
                    "4.19 Batch 16",
                )
            ),
            subject="capstone",
            criterion="the capstone documents the completed milestone chain",
            observed="milestone sections present",
            detail="the capstone is the integration narrative, not only a result table",
        ),
        _check(
            "physics.constants",
            "physics",
            math.isclose(
                units.MOLAR_GAS_CONSTANT,
                0.008314462618,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            and math.isclose(
                units.COULOMB_CONSTANT,
                138.935456,
                rel_tol=1e-6,
                abs_tol=1e-6,
            )
            and units.BOLTZMANN_CONSTANT == units.MOLAR_GAS_CONSTANT,
            subject="units",
            criterion="golden constants retain their documented values",
            observed=(
                f"R={units.MOLAR_GAS_CONSTANT:.12f}, "
                f"kC={units.COULOMB_CONSTANT:.6f}"
            ),
            detail="the compatibility alias BOLTZMANN_CONSTANT remains the molar gas constant",
        ),
        _check(
            "physics.harmonic_exact",
            "physics",
            math.isclose(exact, expected_exact, rel_tol=0.0, abs_tol=1e-12)
            and math.isclose(partition_delta, expected_exact, rel_tol=0.0, abs_tol=1e-12),
            subject="harmonic model",
            criterion="partition-function and analytic free-energy differences agree",
            observed=f"exact={exact:.12f}, partition={partition_delta:.12f}",
            detail="the center shift cancels from the configurational partition function",
        ),
        _check(
            "physics.force_gradient",
            "physics",
            force_error < 1e-6,
            subject="harmonic force",
            criterion="analytic force matches a central finite difference",
            observed=f"max error={force_error:.3e}",
            detail="derivative and force sign conventions are audited on a controlled potential",
        ),
        _check(
            "chemistry.state_hypotheses",
            "chemistry",
            states.population_model is None
            and states.selected_state_id is None
            and len(states.states) >= 2
            and all(state.assumptions for state in states.states),
            subject="ligand states",
            criterion="multiple ligand states remain unselected and unpopulated",
            observed=f"{len(states.states)} states",
            detail="state hypotheses are not silently converted into populations",
        ),
        _check(
            "chemistry.protein_hydrogens",
            "chemistry",
            protein.selected_state_id is None
            and len(protein.states) >= 2
            and all(
                atom.hydrogen_count == 0
                for state in protein.states
                for atom in state.atoms
            ),
            subject="protein states",
            criterion="protein states remain unselected and no hydrogens are implied",
            observed=f"{len(protein.states)} states",
            detail="protein-donor hydrogen bonds must remain unknown",
        ),
        _check(
            "chemistry.hydrogen_provenance",
            "chemistry",
            all(
                item.coordinate_source == "generated"
                and item.parent_coordinate_source == "deposited"
                for item in hydrogen.ligand_hydrogens
            )
            and hydrogen.protein_hydrogen_status == "not_run",
            subject="hydrogen evidence",
            criterion="generated ligand hydrogens remain distinct from deposited heavy atoms",
            observed=f"{len(hydrogen.ligand_hydrogens)} generated ligand hydrogens",
            detail="protein hydrogens remain not_run",
        ),
        _check(
            "chemistry.water_roles",
            "chemistry",
            all(
                item.donor_role is None
                and item.acceptor_role is None
                and not item.orientation_known
                for item in water.waters
            ),
            subject="water evidence",
            criterion="water orientation and roles remain unassigned",
            observed=f"{len(water.waters)} deposited waters",
            detail="proximity is not converted into water chemistry",
        ),
        _check(
            "drug_discovery.contacts",
            "drug_discovery",
            len(site.contacts) > 0 and len(site.waters_near_ligand) > 0,
            subject="binding site",
            criterion="the structural record contains contacts and nearby waters",
            observed=(
                f"{len(site.contacts)} contacts, "
                f"{len(site.waters_near_ligand)} nearby waters"
            ),
            detail="the interaction layer has explicit geometric observations to consume",
        ),
        _check(
            "drug_discovery.matrix",
            "drug_discovery",
            len(combination.ligand_state_ids) == 8
            and len(combination.protein_state_ids) == 2
            and len(combination.combinations) == 16
            and len(
                {len(item.classified_interactions) for item in combination.combinations}
            )
            == 1,
            subject="state matrix",
            criterion="the 8x2 state matrix is complete and count-stable",
            observed=f"{len(combination.combinations)} combinations",
            detail="equal category counts do not make specific participants identical",
        ),
        _check(
            "drug_discovery.robustness",
            "drug_discovery",
            study_statuses == {
                "invariant",
                "assumption_sensitive",
                "unsupported",
                "unknown",
            },
            subject="interaction robustness",
            criterion="all four robustness statuses remain distinct",
            observed=", ".join(sorted(study_statuses)),
            detail="no composite confidence score is present",
        ),
        _check(
            "drug_discovery.thermodynamic_boundary",
            "drug_discovery",
            thermodynamics.free_energy_status == "not_computed"
            and thermodynamics.binding_free_energy_kj_per_mol is None
            and all(term.status == "not_estimated" for term in thermodynamics.terms)
            and any(
                item.boundary_id == "counts_not_ranking"
                for item in thermodynamics.boundaries
            ),
            subject="structural thermodynamics",
            criterion="structural counts are not converted into a binding free energy",
            observed="free energy not_computed; all terms unestimated",
            detail="missing requirements remain explicit",
        ),
        _check(
            "thermodynamics.free_energy_study",
            "thermodynamics",
            len(free_energy_study.cases) == 2
            and mild is not None
            and poor is not None
            and poor.convergence[-1].overlap.effective_sample_fraction
            < mild.convergence[-1].overlap.effective_sample_fraction,
            subject="harmonic FEP study",
            criterion=(
                "good- and poor-overlap cases are both present and distinguish "
                "overlap"
            ),
            observed="mild and poor-overlap cases present",
            detail=(
                "effective sample size is a weight diagnostic, not proof of "
                "phase-space coverage"
            ),
        ),
        _check(
            "thermodynamics.two_state",
            "thermodynamics",
            math.isclose(
                two_state.population0 + two_state.population1,
                1.0,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            and math.isclose(
                two_state.delta_g_kj_per_mol,
                two_state_delta_expected,
                rel_tol=0.0,
                abs_tol=1e-12,
            ),
            subject="two-state model",
            criterion="populations and free-energy difference agree with the exact formula",
            observed=(
                f"p0+p1={two_state.population0 + two_state.population1:.12f}, "
                f"dG={two_state.delta_g_kj_per_mol:.6f}"
            ),
            detail="degeneracy is a dimensionless count in this discrete model",
        ),
        _check(
            "end_to_end.boltzmann_average",
            "thermodynamics",
            abs(
                demonstration.mean_potential_energy_kj_per_mol
                - demonstration.analytic_mean_potential_energy_kj_per_mol
            )
            / demonstration.analytic_mean_potential_energy_kj_per_mol
            < 0.05,
            subject="configuration -> ensemble",
            criterion="sampled mean potential energy agrees with equipartition",
            observed=(
                f"sampled={demonstration.mean_potential_energy_kj_per_mol:.6f}, "
                f"analytic={demonstration.analytic_mean_potential_energy_kj_per_mol:.6f}"
            ),
            detail="the harmonic sampler has a known Boltzmann distribution",
        ),
        _check(
            "end_to_end.fep_reference",
            "thermodynamics",
            abs(
                demonstration.fep_delta_g_kj_per_mol
                - demonstration.exact_delta_g_kj_per_mol
            )
            < 0.15,
            subject="ensemble -> free energy",
            criterion="finite-sample FEP agrees with the exact harmonic reference",
            observed=(
                f"FEP={demonstration.fep_delta_g_kj_per_mol:.6f}, "
                f"exact={demonstration.exact_delta_g_kj_per_mol:.6f}"
            ),
            detail="the agreement is specific to this controlled model and sample size",
        ),
        _check(
            "reproducibility.source_checksum",
            "reproducibility",
            source_checksum_ok,
            subject="deposited structure",
            criterion="the cached source matches the manifest checksum and size",
            observed=f"sha256={source_sha[:12]}..., size={source_size}",
            detail="the experimental source has not been edited in place",
        ),
        _check(
            "reproducibility.manifest_chain",
            "reproducibility",
            not manifest_mismatches,
            subject="record provenance",
            criterion="every structural record references the current manifest checksum",
            observed=f"{len(manifest_mismatches)} mismatches",
            detail="the record chain shares one source manifest",
        ),
        _check(
            "reproducibility.reference_chain",
            "reproducibility",
            not reference_mismatches,
            subject="record references",
            criterion="referenced input checksums match the current artifacts",
            observed=f"{len(reference_mismatches)} mismatches",
            detail=(
                "identity, coordinates, site, states, matrix, robustness, and "
                "thermodynamics form a checksum chain"
            ),
        ),
        _check(
            "reproducibility.study_seeds",
            "reproducibility",
            all(case.seed >= 0 for case in free_energy_study.cases)
            and len({case.seed for case in free_energy_study.cases})
            == len(free_energy_study.cases),
            subject="free-energy study",
            criterion="every study case stores a distinct deterministic seed",
            observed=f"{len(free_energy_study.cases)} seeds",
            detail="the numerical experiment is reproducible from the record",
        ),
        _check(
            "reproducibility.artifacts",
            "reproducibility",
            all(path.exists() for path in artifact_paths),
            subject="artifacts",
            criterion="every audited artifact exists",
            observed=(
                f"{sum(path.exists() for path in artifact_paths)}/"
                f"{len(artifact_paths)} present"
            ),
            detail=(
                "the audit references concrete files rather than in-memory "
                "state only"
            ),
        ),
    )

    fail_count = sum(1 for item in checks if item.status == "fail")
    not_run_count = sum(1 for item in checks if item.status == "not_run")
    pass_count = sum(1 for item in checks if item.status == "pass")
    status = "pass" if fail_count == 0 else "fail"

    limitations = (
        (
            "The audit checks recorded invariants and documentation links; it "
            "does not prove scientific correctness."
        ),
        (
            "The harmonic demonstration is a controlled model, not a "
            "protein-ligand force field or binding calculation."
        ),
        (
            "Effective sample size and maximum weight fraction are "
            "observed-weight diagnostics, not proof of complete phase-space "
            "coverage."
        ),
        (
            "A passed audit does not establish experimental affinity, population "
            "accuracy, or transferability."
        ),
        (
            "The structural chain remains dependent on the deposited structure, "
            "chemical-state hypotheses, and explicit criteria."
        ),
    )

    decisions = (
        DecisionRecord(
            decision=(
                "audited the recorded evidence chain rather than rerunning every "
                "example"
            ),
            actor="scientific audit",
            rationale=(
                "the records expose provenance, criteria, and boundaries that a "
                "summary would hide"
            ),
        ),
        DecisionRecord(
            decision=(
                "kept exact statistical mechanics, numerical estimates, and "
                "reliability diagnostics separate"
            ),
            actor="scientific audit",
            rationale="an estimate must not silently become an exact result",
        ),
        DecisionRecord(
            decision=(
                "treated overlap diagnostics as observed-weight evidence, not "
                "proof of sampling completeness"
            ),
            actor="scientific audit",
            rationale=(
                "a finite sample can look internally consistent while missing a "
                "rare region"
            ),
        ),
        DecisionRecord(
            decision="did not convert the harmonic result into a binding affinity",
            actor="scientific audit",
            rationale=(
                "the protein-ligand model lacks force field, solvent, "
                "populations, and alchemical sampling"
            ),
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "deposited structure and waters",
                (
                    "The audit checks the source checksum and the structural "
                    "record chain."
                ),
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "audit checks and harmonic demonstration",
                (
                    "The audit recomputes selected constants, gradients, "
                    "ensemble averages, and exact references."
                ),
                actor="scientific audit",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No binding affinity, population, or molecular property is predicted.",
                actor="scientific audit",
            ),
            ProvenanceEntry(
                "hypothesised",
                "chemical states and free-energy requirements",
                (
                    "State hypotheses and missing thermodynamic requirements "
                    "remain explicit."
                ),
                actor="scientific audit",
            ),
        )
    )

    return ScientificAuditRecord(
        audited_artifacts=tuple(_display_path(path) for path in artifact_paths),
        checks=checks,
        pass_count=pass_count,
        fail_count=fail_count,
        not_run_count=not_run_count,
        status=status,
        demonstration=demonstration,
        limitations=limitations,
        decisions=decisions,
        provenance=provenance,
    )


def write_scientific_audit(
    path: str | Path,
    record: ScientificAuditRecord,
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_scientific_audit(path: str | Path) -> ScientificAuditRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return ScientificAuditRecord.from_dict(data)


def validate_scientific_audit(record: ScientificAuditRecord) -> list[str]:
    errors: list[str] = []
    if record.schema_version != SCIENTIFIC_AUDIT_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {SCIENTIFIC_AUDIT_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    if not record.checks:
        errors.append("checks must not be empty.")
    check_ids = set()
    for index, check in enumerate(record.checks):
        if check.status not in AUDIT_STATUSES:
            errors.append(f"checks[{index}].status is invalid.")
        if not check.check_id or not check.category or not check.detail:
            errors.append(
                f"checks[{index}] requires check_id, category, and detail."
            )
        if check.check_id in check_ids:
            errors.append(f"checks[{index}].check_id is duplicated.")
        check_ids.add(check.check_id)
    pass_count = sum(1 for item in record.checks if item.status == "pass")
    fail_count = sum(1 for item in record.checks if item.status == "fail")
    not_run_count = sum(1 for item in record.checks if item.status == "not_run")
    if record.pass_count != pass_count:
        errors.append("pass_count does not match the checks.")
    if record.fail_count != fail_count:
        errors.append("fail_count does not match the checks.")
    if record.not_run_count != not_run_count:
        errors.append("not_run_count does not match the checks.")
    expected_status = "pass" if fail_count == 0 else "fail"
    if record.status != expected_status:
        errors.append(f"status must be {expected_status!r}.")
    if not record.limitations:
        errors.append("limitations must not be empty.")
    if not record.decisions:
        errors.append("decisions must contain at least one record.")
    if not record.provenance.records:
        errors.append("provenance.records must not be empty.")
    for index, entry in enumerate(record.provenance.records):
        if entry.label not in PROVENANCE_LABELS:
            errors.append(f"provenance.records[{index}] has an invalid label.")

    demonstration = record.demonstration
    if not math.isfinite(demonstration.exact_delta_g_kj_per_mol):
        errors.append("demonstration exact delta G must be finite.")
    if not math.isfinite(demonstration.fep_delta_g_kj_per_mol):
        errors.append("demonstration FEP delta G must be finite.")
    if not 0.0 < demonstration.overlap.effective_sample_fraction <= 1.0:
        errors.append("demonstration effective sample fraction must be in (0, 1].")
    if not 0.0 < demonstration.proposal_effective_sample_fraction <= 1.0:
        errors.append(
            "demonstration proposal effective sample fraction must be in (0, 1]."
        )
    if not 0.0 < demonstration.proposal_maximum_weight_fraction <= 1.0:
        errors.append(
            "demonstration proposal maximum weight fraction must be in (0, 1]."
        )
    if not demonstration.stop_line:
        errors.append("demonstration stop_line must not be empty.")

    return errors
