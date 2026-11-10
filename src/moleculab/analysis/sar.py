# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Structure–activity relationships and matched-pair comparisons.

The module is intentionally parser-independent. It consumes prepared compound
records and keeps four categories separate:

- measured assay activity;
- calculated molecular descriptors;
- computational docking scores;
- mechanistic hypotheses.

It does not detect matched molecular pairs, calculate descriptors, or predict
activity. Those operations belong to a chemistry toolkit or an external model.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from typing import Literal

import numpy as np

from ..foundation.units import MOLAR_GAS_CONSTANT

ActivityType = Literal["IC50", "Ki", "Kd"]
_ACTIVITY_TYPES = ("IC50", "Ki", "Kd")


@dataclass(frozen=True)
class AssayMeasurement:
    """One measured activity value in molar units.

    ``measurement_type`` is one of ``IC50``, ``Ki``, or ``Kd``. These quantities
    are not interchangeable without assumptions about the assay and the binding
    model.
    """

    value_m: float
    measurement_type: ActivityType
    source: str = ""
    protocol_id: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        value = float(self.value_m)
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError("value_m must be finite and positive.")
        if self.measurement_type not in _ACTIVITY_TYPES:
            raise ValueError(
                f"measurement_type must be one of {_ACTIVITY_TYPES}."
            )
        if (
            not isinstance(self.source, str)
            or not isinstance(self.protocol_id, str)
            or not isinstance(self.notes, str)
        ):
            raise TypeError("source, protocol_id, and notes must be strings.")
        object.__setattr__(self, "value_m", value)

    @property
    def p_value(self) -> float:
        """Return ``-log10(value_m)`` for this measurement."""

        return float(-np.log10(self.value_m))


@dataclass(frozen=True)
class CompoundRecord:
    """Prepared compound data with explicit provenance categories."""

    compound_id: str
    smiles: str = ""
    activity: AssayMeasurement | None = None
    descriptors: Mapping[str, float] = field(default_factory=dict)
    docking_score: float | None = None
    docking_score_unit: str | None = None
    docking_score_source: str = ""
    notes: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.compound_id, str) or not self.compound_id.strip():
            raise ValueError("compound_id must be a non-empty string.")
        if not isinstance(self.smiles, str):
            raise TypeError("smiles must be a string.")

        descriptors = dict(self.descriptors)
        for key, value in descriptors.items():
            if not isinstance(key, str):
                raise TypeError("descriptor names must be strings.")
            number = float(value)
            if not np.isfinite(number):
                raise ValueError("descriptor values must be finite.")
            descriptors[key] = number

        if self.docking_score is not None:
            score = float(self.docking_score)
            if not np.isfinite(score):
                raise ValueError("docking_score must be finite when provided.")
            if not self.docking_score_unit:
                raise ValueError(
                    "docking_score_unit is required when a score is provided."
                )
            if not self.docking_score_source:
                raise ValueError(
                    "docking_score_source is required when a score is provided."
                )
            object.__setattr__(self, "docking_score", score)

        if not isinstance(self.notes, str):
            raise TypeError("notes must be a string.")
        object.__setattr__(self, "descriptors", descriptors)


@dataclass(frozen=True)
class MatchedMolecularPair:
    """Two compounds related by one stated transformation."""

    compound_a: CompoundRecord
    compound_b: CompoundRecord
    transformation: str
    hypothesis: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.compound_a, CompoundRecord) or not isinstance(
            self.compound_b, CompoundRecord
        ):
            raise TypeError("compound_a and compound_b must be CompoundRecord objects.")
        if self.compound_a.compound_id == self.compound_b.compound_id:
            raise ValueError("matched-pair compounds must have different IDs.")
        if not isinstance(self.transformation, str) or not self.transformation.strip():
            raise ValueError("transformation must be a non-empty string.")
        if not isinstance(self.hypothesis, str):
            raise TypeError("hypothesis must be a string.")


@dataclass(frozen=True)
class MatchedPairComparison:
    """Measured, calculated, computational, and hypothesis fields separately."""

    compound_a_id: str
    compound_b_id: str
    transformation: str
    hypothesis: str
    activity_type: ActivityType | None
    p_activity_a: float | None
    p_activity_b: float | None
    delta_p_activity: float | None
    fold_change: float | None
    assay_sources: tuple[str, str]
    assay_protocols: tuple[str, str]
    assay_context_matches: bool
    descriptor_deltas: Mapping[str, float]
    docking_score_delta: float | None


def pIC50(ic50_m: float) -> float:
    """Return ``pIC50 = -log10(IC50 / 1 M)``."""

    value = float(ic50_m)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError("ic50_m must be finite and positive.")
    return float(-np.log10(value))


def ligand_efficiency(
    pic50: float,
    heavy_atom_count: int,
    *,
    temperature: float = 298.15,
) -> float:
    """Return ligand efficiency in kJ/mol per heavy atom.

    This uses the common proxy

    .. math::

        \\mathrm{LE} = \\frac{RT\\ln 10\\, \\mathrm{p}IC_{50}}{N_{\\mathrm{heavy}}}.

    It is a size-normalised potency proxy, not a measured binding free energy.
    """

    p_value = float(pic50)
    if not np.isfinite(p_value):
        raise ValueError("pic50 must be finite.")
    atoms = int(heavy_atom_count)
    if atoms <= 0:
        raise ValueError("heavy_atom_count must be positive.")
    kelvin = float(temperature)
    if not np.isfinite(kelvin) or kelvin <= 0.0:
        raise ValueError("temperature must be finite and positive.")

    rt_ln10 = MOLAR_GAS_CONSTANT * kelvin * np.log(10.0)
    return float(rt_ln10 * p_value / atoms)


def compare_matched_pair(pair: MatchedMolecularPair) -> MatchedPairComparison:
    """Compare two compounds without collapsing categories into one score.

    Descriptor sets must match. Docking scores must share unit and source when
    both are present. Missing activity or docking scores yield ``None`` deltas
    rather than zero.
    """

    if not isinstance(pair, MatchedMolecularPair):
        raise TypeError("pair must be a MatchedMolecularPair.")

    activity_a = pair.compound_a.activity
    activity_b = pair.compound_b.activity
    assay_sources = (
        "" if activity_a is None else activity_a.source,
        "" if activity_b is None else activity_b.source,
    )
    assay_protocols = (
        "" if activity_a is None else activity_a.protocol_id,
        "" if activity_b is None else activity_b.protocol_id,
    )
    assay_context_matches = bool(
        activity_a is not None
        and activity_b is not None
        and activity_a.source
        and activity_b.source
        and activity_a.source == activity_b.source
        and activity_a.protocol_id
        and activity_b.protocol_id
        and activity_a.protocol_id == activity_b.protocol_id
    )
    if activity_a is not None and activity_b is not None:
        if activity_a.measurement_type != activity_b.measurement_type:
            raise ValueError(
                "activity measurements must use the same measurement_type."
            )
        activity_type = activity_a.measurement_type
        p_a = activity_a.p_value
        p_b = activity_b.p_value
        delta = p_b - p_a
        fold = float(10.0**delta)
    else:
        activity_type = None
        p_a = None if activity_a is None else activity_a.p_value
        p_b = None if activity_b is None else activity_b.p_value
        delta = None
        fold = None

    descriptors_a = pair.compound_a.descriptors
    descriptors_b = pair.compound_b.descriptors
    if set(descriptors_a) != set(descriptors_b):
        raise ValueError("both compounds must provide the same descriptor set.")
    descriptor_deltas = {
        key: descriptors_b[key] - descriptors_a[key]
        for key in descriptors_a
    }

    score_a = pair.compound_a.docking_score
    score_b = pair.compound_b.docking_score
    if score_a is not None and score_b is not None:
        if pair.compound_a.docking_score_unit != pair.compound_b.docking_score_unit:
            raise ValueError("docking-score units must match for comparison.")
        if pair.compound_a.docking_score_source != pair.compound_b.docking_score_source:
            raise ValueError("docking-score sources must match for comparison.")
        score_delta = score_b - score_a
    else:
        score_delta = None

    return MatchedPairComparison(
        compound_a_id=pair.compound_a.compound_id,
        compound_b_id=pair.compound_b.compound_id,
        transformation=pair.transformation,
        hypothesis=pair.hypothesis,
        activity_type=activity_type,
        p_activity_a=p_a,
        p_activity_b=p_b,
        delta_p_activity=delta,
        fold_change=fold,
        assay_sources=assay_sources,
        assay_protocols=assay_protocols,
        assay_context_matches=assay_context_matches,
        descriptor_deltas=descriptor_deltas,
        docking_score_delta=score_delta,
    )
