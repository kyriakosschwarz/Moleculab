# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Educational free-energy study on an exactly solvable harmonic model.

This module adds the analytical and convergence layer around the existing
harmonic FEP estimator:

* exact configurational partition functions and their free energies;
* effective sample size and weight-overlap diagnostics;
* sample-size convergence for good- and poor-overlap transformations;
* replicate uncertainty and confidence-interval coverage.

The model is a three-dimensional harmonic oscillator, not a protein-ligand
system.  The absolute configurational free energy omits the thermal-wavelength
term, so only free-energy differences are compared.
"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ...foundation.units import MOLAR_GAS_CONSTANT
from ...models.potentials import HarmonicPotential
from .fep import (
    FEPResult,
    fep_harmonic_oscillator,
    free_energy_perturbation,
    harmonic_free_energy_difference,
)
from ..structure.structure_manifest import (
    PROVENANCE_LABELS,
    DecisionRecord,
    ProvenanceEntry,
    ProvenanceRecord,
    ValidationCheck,
    ValidationRecord,
)

FREE_ENERGY_STUDY_SCHEMA_VERSION = "1.0"


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


def _required_float(data: Mapping[str, Any], key: str) -> float:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return float(value)


def _required_int(data: Mapping[str, Any], key: str) -> int:
    value = data.get(key)
    if value is None:
        raise ValueError(f"{key} is required.")
    return int(value)


def _optional_string(data: Mapping[str, Any], key: str) -> str | None:
    value = data.get(key)
    return None if value is None else str(value)


def _string_tuple(data: Mapping[str, Any], key: str) -> tuple[str, ...]:
    return tuple(str(item) for item in _as_list(data.get(key), key))


def _int_tuple(data: Mapping[str, Any], key: str) -> tuple[int, ...]:
    return tuple(int(item) for item in _as_list(data.get(key), key))


@dataclass(frozen=True)
class FEPOverlapDiagnostics:
    n_samples: int
    effective_sample_size: float
    effective_sample_fraction: float
    maximum_weight_fraction: float
    work_mean_kj_per_mol: float
    work_standard_deviation_kj_per_mol: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_samples": self.n_samples,
            "effective_sample_size": self.effective_sample_size,
            "effective_sample_fraction": self.effective_sample_fraction,
            "maximum_weight_fraction": self.maximum_weight_fraction,
            "work_mean_kj_per_mol": self.work_mean_kj_per_mol,
            "work_standard_deviation_kj_per_mol": (
                self.work_standard_deviation_kj_per_mol
            ),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> FEPOverlapDiagnostics:
        data = _as_mapping(data, "FEP overlap diagnostics")
        return cls(
            n_samples=_required_int(data, "n_samples"),
            effective_sample_size=_required_float(
                data, "effective_sample_size"
            ),
            effective_sample_fraction=_required_float(
                data, "effective_sample_fraction"
            ),
            maximum_weight_fraction=_required_float(
                data, "maximum_weight_fraction"
            ),
            work_mean_kj_per_mol=_required_float(
                data, "work_mean_kj_per_mol"
            ),
            work_standard_deviation_kj_per_mol=_required_float(
                data, "work_standard_deviation_kj_per_mol"
            ),
        )


@dataclass(frozen=True)
class TwoStateReference:
    energy0_kj_per_mol: float
    energy1_kj_per_mol: float
    degeneracy0: int
    degeneracy1: int
    partition_function: float
    population0: float
    population1: float
    delta_g_kj_per_mol: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "energy0_kj_per_mol": self.energy0_kj_per_mol,
            "energy1_kj_per_mol": self.energy1_kj_per_mol,
            "degeneracy0": self.degeneracy0,
            "degeneracy1": self.degeneracy1,
            "partition_function": self.partition_function,
            "population0": self.population0,
            "population1": self.population1,
            "delta_g_kj_per_mol": self.delta_g_kj_per_mol,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> TwoStateReference:
        data = _as_mapping(data, "two-state reference")
        return cls(
            energy0_kj_per_mol=_required_float(data, "energy0_kj_per_mol"),
            energy1_kj_per_mol=_required_float(data, "energy1_kj_per_mol"),
            degeneracy0=_required_int(data, "degeneracy0"),
            degeneracy1=_required_int(data, "degeneracy1"),
            partition_function=_required_float(data, "partition_function"),
            population0=_required_float(data, "population0"),
            population1=_required_float(data, "population1"),
            delta_g_kj_per_mol=_required_float(data, "delta_g_kj_per_mol"),
        )


@dataclass(frozen=True)
class FEPConvergencePoint:
    n_samples: int
    delta_g_kj_per_mol: float
    standard_error_kj_per_mol: float
    error_kj_per_mol: float
    confidence_interval_contains_exact: bool
    overlap: FEPOverlapDiagnostics

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_samples": self.n_samples,
            "delta_g_kj_per_mol": self.delta_g_kj_per_mol,
            "standard_error_kj_per_mol": self.standard_error_kj_per_mol,
            "error_kj_per_mol": self.error_kj_per_mol,
            "confidence_interval_contains_exact": (
                self.confidence_interval_contains_exact
            ),
            "overlap": self.overlap.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> FEPConvergencePoint:
        data = _as_mapping(data, "FEP convergence point")
        return cls(
            n_samples=_required_int(data, "n_samples"),
            delta_g_kj_per_mol=_required_float(data, "delta_g_kj_per_mol"),
            standard_error_kj_per_mol=_required_float(
                data, "standard_error_kj_per_mol"
            ),
            error_kj_per_mol=_required_float(data, "error_kj_per_mol"),
            confidence_interval_contains_exact=bool(
                data.get("confidence_interval_contains_exact")
            ),
            overlap=FEPOverlapDiagnostics.from_dict(
                _as_mapping(data.get("overlap"), "overlap")
            ),
        )


@dataclass(frozen=True)
class FEPReplicateSummary:
    n_replicates: int
    n_samples_per_replicate: int
    mean_delta_g_kj_per_mol: float
    standard_deviation_kj_per_mol: float
    minimum_delta_g_kj_per_mol: float
    maximum_delta_g_kj_per_mol: float
    mean_standard_error_kj_per_mol: float
    mean_error_kj_per_mol: float
    confidence_interval_coverage: float

    def to_dict(self) -> dict[str, Any]:
        return {
            "n_replicates": self.n_replicates,
            "n_samples_per_replicate": self.n_samples_per_replicate,
            "mean_delta_g_kj_per_mol": self.mean_delta_g_kj_per_mol,
            "standard_deviation_kj_per_mol": self.standard_deviation_kj_per_mol,
            "minimum_delta_g_kj_per_mol": self.minimum_delta_g_kj_per_mol,
            "maximum_delta_g_kj_per_mol": self.maximum_delta_g_kj_per_mol,
            "mean_standard_error_kj_per_mol": (
                self.mean_standard_error_kj_per_mol
            ),
            "mean_error_kj_per_mol": self.mean_error_kj_per_mol,
            "confidence_interval_coverage": self.confidence_interval_coverage,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> FEPReplicateSummary:
        data = _as_mapping(data, "FEP replicate summary")
        return cls(
            n_replicates=_required_int(data, "n_replicates"),
            n_samples_per_replicate=_required_int(
                data, "n_samples_per_replicate"
            ),
            mean_delta_g_kj_per_mol=_required_float(
                data, "mean_delta_g_kj_per_mol"
            ),
            standard_deviation_kj_per_mol=_required_float(
                data, "standard_deviation_kj_per_mol"
            ),
            minimum_delta_g_kj_per_mol=_required_float(
                data, "minimum_delta_g_kj_per_mol"
            ),
            maximum_delta_g_kj_per_mol=_required_float(
                data, "maximum_delta_g_kj_per_mol"
            ),
            mean_standard_error_kj_per_mol=_required_float(
                data, "mean_standard_error_kj_per_mol"
            ),
            mean_error_kj_per_mol=_required_float(
                data, "mean_error_kj_per_mol"
            ),
            confidence_interval_coverage=_required_float(
                data, "confidence_interval_coverage"
            ),
        )


@dataclass(frozen=True)
class HarmonicFEPCase:
    case_id: str
    label: str
    potential0_stiffness: float
    potential1_stiffness: float
    potential0_center_nm: tuple[float, float, float]
    potential1_center_nm: tuple[float, float, float]
    exact_partition_function0: float
    exact_partition_function1: float
    exact_configurational_free_energy0_kj_per_mol: float
    exact_configurational_free_energy1_kj_per_mol: float
    exact_delta_g_kj_per_mol: float
    convergence: tuple[FEPConvergencePoint, ...]
    replicates: FEPReplicateSummary
    seed: int
    note: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "case_id": self.case_id,
            "label": self.label,
            "potential0_stiffness": self.potential0_stiffness,
            "potential1_stiffness": self.potential1_stiffness,
            "potential0_center_nm": list(self.potential0_center_nm),
            "potential1_center_nm": list(self.potential1_center_nm),
            "exact_partition_function0": self.exact_partition_function0,
            "exact_partition_function1": self.exact_partition_function1,
            "exact_configurational_free_energy0_kj_per_mol": (
                self.exact_configurational_free_energy0_kj_per_mol
            ),
            "exact_configurational_free_energy1_kj_per_mol": (
                self.exact_configurational_free_energy1_kj_per_mol
            ),
            "exact_delta_g_kj_per_mol": self.exact_delta_g_kj_per_mol,
            "convergence": [item.to_dict() for item in self.convergence],
            "replicates": self.replicates.to_dict(),
            "seed": self.seed,
            "note": self.note,
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HarmonicFEPCase:
        data = _as_mapping(data, "harmonic FEP case")
        return cls(
            case_id=_required_string(data, "case_id"),
            label=_required_string(data, "label"),
            potential0_stiffness=_required_float(
                data, "potential0_stiffness"
            ),
            potential1_stiffness=_required_float(
                data, "potential1_stiffness"
            ),
            potential0_center_nm=tuple(
                float(item)
                for item in _as_list(
                    data.get("potential0_center_nm"), "potential0_center_nm"
                )
            ),
            potential1_center_nm=tuple(
                float(item)
                for item in _as_list(
                    data.get("potential1_center_nm"), "potential1_center_nm"
                )
            ),
            exact_partition_function0=_required_float(
                data, "exact_partition_function0"
            ),
            exact_partition_function1=_required_float(
                data, "exact_partition_function1"
            ),
            exact_configurational_free_energy0_kj_per_mol=_required_float(
                data, "exact_configurational_free_energy0_kj_per_mol"
            ),
            exact_configurational_free_energy1_kj_per_mol=_required_float(
                data, "exact_configurational_free_energy1_kj_per_mol"
            ),
            exact_delta_g_kj_per_mol=_required_float(
                data, "exact_delta_g_kj_per_mol"
            ),
            convergence=tuple(
                FEPConvergencePoint.from_dict(
                    _as_mapping(item, "FEP convergence point")
                )
                for item in _as_list(data.get("convergence"), "convergence")
            ),
            replicates=FEPReplicateSummary.from_dict(
                _as_mapping(data.get("replicates"), "replicates")
            ),
            seed=_required_int(data, "seed"),
            note=_optional_string(data, "note") or "",
        )


@dataclass(frozen=True)
class HarmonicFEPStudyRecord:
    temperature_kelvin: float
    two_state_example: TwoStateReference
    sample_sizes: tuple[int, ...]
    replicate_count: int
    replicate_sample_size: int
    cases: tuple[HarmonicFEPCase, ...]
    limitations: tuple[str, ...]
    decisions: tuple[DecisionRecord, ...]
    provenance: ProvenanceRecord
    validation: ValidationRecord
    schema_version: str = FREE_ENERGY_STUDY_SCHEMA_VERSION

    def to_dict(self) -> dict[str, Any]:
        return {
            "schema_version": self.schema_version,
            "temperature_kelvin": self.temperature_kelvin,
            "two_state_example": self.two_state_example.to_dict(),
            "sample_sizes": list(self.sample_sizes),
            "replicate_count": self.replicate_count,
            "replicate_sample_size": self.replicate_sample_size,
            "cases": [item.to_dict() for item in self.cases],
            "limitations": list(self.limitations),
            "decisions": [decision.to_dict() for decision in self.decisions],
            "provenance": self.provenance.to_dict(),
            "validation": self.validation.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: Mapping[str, Any]) -> HarmonicFEPStudyRecord:
        data = _as_mapping(data, "harmonic FEP study record")
        return cls(
            schema_version=_required_string(data, "schema_version"),
            temperature_kelvin=_required_float(data, "temperature_kelvin"),
            two_state_example=TwoStateReference.from_dict(
                _as_mapping(data.get("two_state_example"), "two_state_example")
            ),
            sample_sizes=_int_tuple(data, "sample_sizes"),
            replicate_count=_required_int(data, "replicate_count"),
            replicate_sample_size=_required_int(
                data, "replicate_sample_size"
            ),
            cases=tuple(
                HarmonicFEPCase.from_dict(_as_mapping(item, "harmonic FEP case"))
                for item in _as_list(data.get("cases"), "cases")
            ),
            limitations=_string_tuple(data, "limitations"),
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
        )


def harmonic_configurational_partition_function(
    potential: HarmonicPotential,
    temperature: float,
) -> float:
    """Return the 3D classical configurational partition function.

    The result is :math:`(2\\pi RT/k)^{3/2}` in nm^3.  It omits the thermal
    de Broglie wavelength, so absolute free energies carry an additive constant;
    free-energy differences are well defined.
    """

    if not isinstance(potential, HarmonicPotential):
        raise TypeError("potential must be a HarmonicPotential.")
    temperature_value = float(temperature)
    if not np.isfinite(temperature_value) or temperature_value <= 0.0:
        raise ValueError("temperature must be finite and positive.")
    rt = MOLAR_GAS_CONSTANT * temperature_value
    return float((2.0 * math.pi * rt / potential.stiffness) ** 1.5)


def harmonic_configurational_free_energy(
    potential: HarmonicPotential,
    temperature: float,
) -> float:
    """Return :math:`-RT\\ln Z` in kJ/mol for the configurational partition function."""

    temperature_value = float(temperature)
    partition_function = harmonic_configurational_partition_function(
        potential,
        temperature_value,
    )
    return float(
        -MOLAR_GAS_CONSTANT * temperature_value * math.log(partition_function)
    )


def two_state_partition_function(
    energy0_kj_per_mol: float,
    energy1_kj_per_mol: float,
    degeneracy0: int,
    degeneracy1: int,
    temperature: float,
) -> float:
    """Return the discrete two-state partition function."""

    temperature_value = float(temperature)
    if not np.isfinite(temperature_value) or temperature_value <= 0.0:
        raise ValueError("temperature must be finite and positive.")
    if degeneracy0 < 1 or degeneracy1 < 1:
        raise ValueError("degeneracies must be positive integers.")
    rt = MOLAR_GAS_CONSTANT * temperature_value
    return float(
        degeneracy0 * math.exp(-float(energy0_kj_per_mol) / rt)
        + degeneracy1 * math.exp(-float(energy1_kj_per_mol) / rt)
    )


def two_state_populations(
    energy0_kj_per_mol: float,
    energy1_kj_per_mol: float,
    degeneracy0: int,
    degeneracy1: int,
    temperature: float,
) -> tuple[float, float]:
    """Return the Boltzmann populations of the two discrete states."""

    partition_function = two_state_partition_function(
        energy0_kj_per_mol,
        energy1_kj_per_mol,
        degeneracy0,
        degeneracy1,
        temperature,
    )
    rt = MOLAR_GAS_CONSTANT * float(temperature)
    population0 = (
        degeneracy0 * math.exp(-float(energy0_kj_per_mol) / rt)
        / partition_function
    )
    population1 = (
        degeneracy1 * math.exp(-float(energy1_kj_per_mol) / rt)
        / partition_function
    )
    return float(population0), float(population1)


def two_state_free_energy_difference(
    energy0_kj_per_mol: float,
    energy1_kj_per_mol: float,
    degeneracy0: int,
    degeneracy1: int,
    temperature: float,
) -> float:
    """Return :math:`G_1-G_0` for two discrete states."""

    if degeneracy0 < 1 or degeneracy1 < 1:
        raise ValueError("degeneracies must be positive integers.")
    temperature_value = float(temperature)
    if not np.isfinite(temperature_value) or temperature_value <= 0.0:
        raise ValueError("temperature must be finite and positive.")
    rt = MOLAR_GAS_CONSTANT * temperature_value
    return float(
        (float(energy1_kj_per_mol) - float(energy0_kj_per_mol))
        - rt * math.log(degeneracy1 / degeneracy0)
    )


def two_state_reference(
    *,
    energy0_kj_per_mol: float,
    energy1_kj_per_mol: float,
    degeneracy0: int = 1,
    degeneracy1: int = 1,
    temperature: float,
) -> TwoStateReference:
    """Build an exact two-state reference for the study record."""

    populations = two_state_populations(
        energy0_kj_per_mol,
        energy1_kj_per_mol,
        degeneracy0,
        degeneracy1,
        temperature,
    )
    return TwoStateReference(
        energy0_kj_per_mol=float(energy0_kj_per_mol),
        energy1_kj_per_mol=float(energy1_kj_per_mol),
        degeneracy0=int(degeneracy0),
        degeneracy1=int(degeneracy1),
        partition_function=two_state_partition_function(
            energy0_kj_per_mol,
            energy1_kj_per_mol,
            degeneracy0,
            degeneracy1,
            temperature,
        ),
        population0=populations[0],
        population1=populations[1],
        delta_g_kj_per_mol=two_state_free_energy_difference(
            energy0_kj_per_mol,
            energy1_kj_per_mol,
            degeneracy0,
            degeneracy1,
            temperature,
        ),
    )


def fep_overlap_diagnostics(
    work: np.ndarray,
    temperature: float,
) -> FEPOverlapDiagnostics:
    """Return effective sample size and weight-overlap diagnostics.

    The FEP weights are :math:`w_i \\propto e^{-\\beta\\Delta U_i}`.  The
    effective sample size is :math:`1/\\sum_i p_i^2` with normalized weights
    :math:`p_i`.  A small fraction indicates that a few samples dominate the
    exponential average.
    """

    result = free_energy_perturbation(work, temperature)
    rt = MOLAR_GAS_CONSTANT * result.temperature
    log_weights = -result.work / rt
    shifted = log_weights - float(np.max(log_weights))
    weights = np.exp(shifted)
    probabilities = weights / float(np.sum(weights))
    effective_sample_size = 1.0 / float(np.sum(probabilities * probabilities))
    return FEPOverlapDiagnostics(
        n_samples=result.n_samples,
        effective_sample_size=effective_sample_size,
        effective_sample_fraction=effective_sample_size / result.n_samples,
        maximum_weight_fraction=float(np.max(probabilities)),
        work_mean_kj_per_mol=float(np.mean(result.work)),
        work_standard_deviation_kj_per_mol=float(
            np.std(result.work, ddof=1)
        ),
    )


def _convergence_point(
    result: FEPResult,
    exact_delta_g: float,
) -> FEPConvergencePoint:
    lower, upper = result.confidence_interval
    return FEPConvergencePoint(
        n_samples=result.n_samples,
        delta_g_kj_per_mol=result.delta_g,
        standard_error_kj_per_mol=result.standard_error,
        error_kj_per_mol=result.delta_g - exact_delta_g,
        confidence_interval_contains_exact=lower <= exact_delta_g <= upper,
        overlap=fep_overlap_diagnostics(result.work, result.temperature),
    )


def _replicate_summary(
    results: Sequence[FEPResult],
    exact_delta_g: float,
) -> FEPReplicateSummary:
    if not results:
        raise ValueError("results must not be empty.")
    estimates = np.array([result.delta_g for result in results], dtype=float)
    standard_errors = np.array(
        [result.standard_error for result in results], dtype=float
    )
    coverage = np.mean(
        [
            result.confidence_interval[0]
            <= exact_delta_g
            <= result.confidence_interval[1]
            for result in results
        ]
    )
    return FEPReplicateSummary(
        n_replicates=len(results),
        n_samples_per_replicate=results[0].n_samples,
        mean_delta_g_kj_per_mol=float(np.mean(estimates)),
        standard_deviation_kj_per_mol=float(np.std(estimates, ddof=1))
        if len(estimates) > 1
        else 0.0,
        minimum_delta_g_kj_per_mol=float(np.min(estimates)),
        maximum_delta_g_kj_per_mol=float(np.max(estimates)),
        mean_standard_error_kj_per_mol=float(np.mean(standard_errors)),
        mean_error_kj_per_mol=float(np.mean(estimates) - exact_delta_g),
        confidence_interval_coverage=float(coverage),
    )


def run_harmonic_fep_case(
    *,
    case_id: str,
    label: str,
    potential0: HarmonicPotential,
    potential1: HarmonicPotential,
    temperature: float,
    sample_sizes: Sequence[int],
    n_replicates: int,
    replicate_sample_size: int,
    seed: int = 2024,
    note: str = "",
) -> HarmonicFEPCase:
    """Run convergence and replicate studies for one harmonic transformation."""

    if not isinstance(potential0, HarmonicPotential) or not isinstance(
        potential1, HarmonicPotential
    ):
        raise TypeError("potential0 and potential1 must be HarmonicPotential.")
    temperature_value = float(temperature)
    if not np.isfinite(temperature_value) or temperature_value <= 0.0:
        raise ValueError("temperature must be finite and positive.")
    sample_sizes = tuple(int(size) for size in sample_sizes)
    if not sample_sizes or any(size < 2 for size in sample_sizes):
        raise ValueError("sample_sizes must contain integers of at least 2.")
    if n_replicates < 1 or replicate_sample_size < 2:
        raise ValueError(
            "n_replicates must be positive and replicate_sample_size at least 2."
        )

    exact_delta_g = harmonic_free_energy_difference(
        potential0,
        potential1,
        temperature_value,
    )
    convergence = tuple(
        _convergence_point(
            fep_harmonic_oscillator(
                potential0,
                potential1,
                temperature_value,
                n_samples=size,
                seed=seed + index,
            ),
            exact_delta_g,
        )
        for index, size in enumerate(sample_sizes)
    )
    replicate_results = tuple(
        fep_harmonic_oscillator(
            potential0,
            potential1,
            temperature_value,
            n_samples=replicate_sample_size,
            seed=seed + 10_000 + index,
        )
        for index in range(n_replicates)
    )
    return HarmonicFEPCase(
        case_id=case_id,
        label=label,
        potential0_stiffness=potential0.stiffness,
        potential1_stiffness=potential1.stiffness,
        potential0_center_nm=tuple(float(value) for value in potential0.center),
        potential1_center_nm=tuple(float(value) for value in potential1.center),
        exact_partition_function0=(
            harmonic_configurational_partition_function(
                potential0, temperature_value
            )
        ),
        exact_partition_function1=(
            harmonic_configurational_partition_function(
                potential1, temperature_value
            )
        ),
        exact_configurational_free_energy0_kj_per_mol=(
            harmonic_configurational_free_energy(
                potential0, temperature_value
            )
        ),
        exact_configurational_free_energy1_kj_per_mol=(
            harmonic_configurational_free_energy(
                potential1, temperature_value
            )
        ),
        exact_delta_g_kj_per_mol=exact_delta_g,
        convergence=convergence,
        replicates=_replicate_summary(replicate_results, exact_delta_g),
        seed=int(seed),
        note=note,
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
) -> ValidationCheck:
    return ValidationCheck(
        check_id=check_id,
        category=category,
        status="pass" if ok else "fail",
        subject=subject,
        criterion=criterion,
        observed=observed,
        expected=expected,
        detail=detail,
    )


def build_harmonic_fep_study_record(
    *,
    temperature_kelvin: float,
    two_state_example: TwoStateReference,
    sample_sizes: tuple[int, ...],
    replicate_count: int,
    replicate_sample_size: int,
    cases: tuple[HarmonicFEPCase, ...],
) -> HarmonicFEPStudyRecord:
    if not cases:
        raise ValueError("cases must not be empty.")
    if not sample_sizes:
        raise ValueError("sample_sizes must not be empty.")
    for case in cases:
        if tuple(point.n_samples for point in case.convergence) != sample_sizes:
            raise ValueError(
                f"{case.case_id} convergence points must cover sample_sizes."
            )
        if case.replicates.n_replicates != replicate_count:
            raise ValueError(
                f"{case.case_id} replicate count does not match the study."
            )
        if (
            case.replicates.n_samples_per_replicate
            != replicate_sample_size
        ):
            raise ValueError(
                f"{case.case_id} replicate sample size does not match the study."
            )

    checks = (
        _check(
            "reference",
            "free_energy.two_state_reference",
            math.isclose(
                two_state_example.population0 + two_state_example.population1,
                1.0,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            and math.isfinite(two_state_example.partition_function),
            subject="two-state model",
            criterion="the discrete populations sum to one",
            observed=(
                f"p0={two_state_example.population0:.6f}, "
                f"p1={two_state_example.population1:.6f}"
            ),
            expected="p0 + p1 = 1",
            detail=(
                "the two-state reference shows how degeneracy enters the "
                "partition function and free energy"
            ),
        ),
        _check(
            "reference",
            "free_energy.partition_functions",
            all(
                case.exact_partition_function0 > 0
                and case.exact_partition_function1 > 0
                for case in cases
            ),
            subject="harmonic model",
            criterion="exact configurational partition functions are positive",
            observed=f"{len(cases)} cases",
            expected="positive Z for both states",
            detail="absolute A omits the thermal-wavelength term; differences are meaningful",
        ),
        _check(
            "reference",
            "free_energy.exact_delta_g",
            all(
                math.isfinite(case.exact_delta_g_kj_per_mol) for case in cases
            ),
            subject="harmonic model",
            criterion="an exact free-energy difference is available for every case",
            observed=f"{len(cases)} exact values",
            expected="analytic harmonic delta G",
            detail="the exact value is the reference for convergence and bias",
        ),
        _check(
            "sampling",
            "free_energy.convergence_grid",
            all(
                len(case.convergence) == len(sample_sizes) for case in cases
            ),
            subject="harmonic model",
            criterion="every case uses the declared sample-size grid",
            observed=f"{len(sample_sizes)} sample sizes",
            expected="same convergence grid for every case",
            detail="sample-size trends are compared on the same grid",
        ),
        _check(
            "sampling",
            "free_energy.overlap_diagnostics",
            all(
                point.overlap.effective_sample_fraction > 0.0
                and point.overlap.maximum_weight_fraction <= 1.0
                for case in cases
                for point in case.convergence
            ),
            subject="harmonic model",
            criterion="effective sample size and maximum weight are recorded",
            observed="overlap diagnostics for every convergence point",
            expected="finite effective sample size and normalized weights",
            detail="a small effective-sample fraction signals poor overlap",
        ),
        _check(
            "uncertainty",
            "free_energy.replicates",
            all(case.replicates.n_replicates == replicate_count for case in cases),
            subject="harmonic model",
            criterion="replicate estimates are summarized for every case",
            observed=f"{replicate_count} replicates x {replicate_sample_size} samples",
            expected="mean, spread, and confidence-interval coverage",
            detail="repeated independent estimates expose run-to-run variability",
        ),
        _check(
            "interpretation",
            "free_energy.not_binding",
            True,
            subject="harmonic model",
            criterion="the study is not presented as a protein-ligand affinity",
            observed="3D harmonic oscillator only",
            expected="no binding, solvent, or force-field claim",
            detail="the exact reference is a toy statistical-mechanics model",
        ),
    )

    decisions = (
        DecisionRecord(
            decision="used an exactly solvable harmonic model",
            actor="free-energy study",
            rationale="an analytic reference separates estimator error from model error",
        ),
        DecisionRecord(
            decision="reported exact and estimated free energies separately",
            actor="free-energy study",
            rationale="agreement or disagreement is only meaningful against a reference",
        ),
        DecisionRecord(
            decision="included a good-overlap and a poor-overlap case",
            actor="free-energy study",
            rationale="the same estimator can behave very differently when state overlap changes",
        ),
        DecisionRecord(
            decision="reported effective sample size and maximum weight fraction",
            actor="free-energy study",
            rationale="sample count alone does not reveal whether a few configurations dominate",
        ),
        DecisionRecord(
            decision="did not present the harmonic result as an imatinib affinity",
            actor="free-energy study",
            rationale="the capstone structural record remains separate from this toy calculation",
        ),
    )

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "none",
                "No experimental data are used in this controlled harmonic study.",
                actor="free-energy study",
            ),
            ProvenanceEntry(
                "calculated",
                "exact partition functions and FEP estimates",
                (
                    "Harmonic partition functions and FEP estimates are "
                    "calculated from explicit model parameters."
                ),
                actor="free-energy study",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No binding affinity or molecular property is predicted.",
                actor="free-energy study",
            ),
            ProvenanceEntry(
                "hypothesised",
                "harmonic test model",
                (
                    "The harmonic oscillator is a deliberate teaching model, "
                    "not a molecular force field."
                ),
                actor="free-energy study",
            ),
        )
    )

    limitations = (
        (
            "The model is a three-dimensional harmonic oscillator, not a "
            "protein-ligand system."
        ),
        (
            "Absolute configurational free energies omit the thermal-wavelength "
            "term; only differences are compared."
        ),
        (
            "Forward FEP is used; no bidirectional estimator, MBAR overlap "
            "correction, or enhanced sampling is applied."
        ),
        (
            "Finite sample sizes cannot remove the bias of poor overlap, even "
            "when the reported standard error is small."
        ),
        (
            "No solvent, protonation, conformational ensemble, force field, or "
            "experimental affinity is included."
        ),
    )

    return HarmonicFEPStudyRecord(
        temperature_kelvin=float(temperature_kelvin),
        two_state_example=two_state_example,
        sample_sizes=sample_sizes,
        replicate_count=replicate_count,
        replicate_sample_size=replicate_sample_size,
        cases=cases,
        limitations=limitations,
        decisions=decisions,
        provenance=provenance,
        validation=ValidationRecord(checks=checks),
    )


def write_harmonic_fep_study(
    path: str | Path,
    record: HarmonicFEPStudyRecord,
) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(record.to_dict(), indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )
    return destination


def load_harmonic_fep_study(path: str | Path) -> HarmonicFEPStudyRecord:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    return HarmonicFEPStudyRecord.from_dict(data)


def validate_harmonic_fep_study(
    record: HarmonicFEPStudyRecord,
) -> list[str]:
    errors: list[str] = []
    if record.schema_version != FREE_ENERGY_STUDY_SCHEMA_VERSION:
        errors.append(
            f"schema_version must be {FREE_ENERGY_STUDY_SCHEMA_VERSION!r}, got "
            f"{record.schema_version!r}."
        )
    if not math.isfinite(record.temperature_kelvin) or record.temperature_kelvin <= 0:
        errors.append("temperature_kelvin must be finite and positive.")
    if not record.sample_sizes or any(size < 2 for size in record.sample_sizes):
        errors.append("sample_sizes must contain integers of at least 2.")
    if record.replicate_count < 1:
        errors.append("replicate_count must be positive.")
    if record.replicate_sample_size < 2:
        errors.append("replicate_sample_size must be at least 2.")
    if record.two_state_example.degeneracy0 < 1 or record.two_state_example.degeneracy1 < 1:
        errors.append("two-state degeneracies must be positive.")
    if not math.isclose(
        record.two_state_example.population0
        + record.two_state_example.population1,
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        errors.append("two-state populations must sum to one.")
    if record.two_state_example.partition_function <= 0:
        errors.append("two-state partition function must be positive.")
    if not record.cases:
        errors.append("cases must not be empty.")

    case_ids = set()
    for index, case in enumerate(record.cases):
        if not case.case_id or not case.label:
            errors.append(f"cases[{index}] requires case_id and label.")
        if case.case_id in case_ids:
            errors.append(f"cases[{index}].case_id is duplicated.")
        case_ids.add(case.case_id)
        if case.exact_partition_function0 <= 0 or case.exact_partition_function1 <= 0:
            errors.append(f"cases[{index}] partition functions must be positive.")
        if not math.isfinite(case.exact_delta_g_kj_per_mol):
            errors.append(f"cases[{index}] exact delta G must be finite.")
        if case.seed < 0:
            errors.append(f"cases[{index}] seed must be non-negative.")
        if tuple(point.n_samples for point in case.convergence) != record.sample_sizes:
            errors.append(f"cases[{index}] convergence grid does not match sample_sizes.")
        for point_index, point in enumerate(case.convergence):
            if point.standard_error_kj_per_mol < 0:
                errors.append(
                    f"cases[{index}].convergence[{point_index}] standard error "
                    "cannot be negative."
                )
            if not 0.0 < point.overlap.effective_sample_fraction <= 1.0:
                errors.append(
                    f"cases[{index}].convergence[{point_index}] effective sample "
                    "fraction must be in (0, 1]."
                )
            if not 0.0 < point.overlap.maximum_weight_fraction <= 1.0:
                errors.append(
                    f"cases[{index}].convergence[{point_index}] maximum weight "
                    "fraction must be in (0, 1]."
                )
        if case.replicates.n_replicates != record.replicate_count:
            errors.append(f"cases[{index}] replicate count does not match the study.")
        if (
            case.replicates.n_samples_per_replicate
            != record.replicate_sample_size
        ):
            errors.append(
                f"cases[{index}] replicate sample size does not match the study."
            )
        if not 0.0 <= case.replicates.confidence_interval_coverage <= 1.0:
            errors.append(
                f"cases[{index}] confidence-interval coverage must be in [0, 1]."
            )

    if not record.limitations:
        errors.append("limitations must not be empty.")
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

    serialized = record.to_dict()
    for forbidden in ("binding_free_energy_kj_per_mol", "affinity", "rank", "score"):
        if forbidden in serialized:
            errors.append(f"{forbidden} must not be a top-level field in this study.")

    return errors
