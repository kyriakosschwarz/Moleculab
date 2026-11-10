# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for SAR matched-pair comparison and activity metrics."""

import numpy as np
import pytest

from moleculab import (
    AssayMeasurement,
    CompoundRecord,
    MatchedMolecularPair,
    compare_matched_pair,
    ligand_efficiency,
    pIC50,
)


def make_compound(
    compound_id,
    *,
    activity=None,
    descriptors=None,
    score=None,
    score_unit="arbitrary",
    score_source="example ranker",
):
    return CompoundRecord(
        compound_id=compound_id,
        activity=activity,
        descriptors={} if descriptors is None else descriptors,
        docking_score=score,
        docking_score_unit=score_unit if score is not None else None,
        docking_score_source=score_source if score is not None else "",
    )


def test_pic50_converts_molar_concentrations():
    assert pIC50(1e-6) == pytest.approx(6.0)
    assert pIC50(1e-7) == pytest.approx(7.0)
    assert pIC50(1e-5) == pytest.approx(5.0)


def test_pic50_validation():
    for value in [0.0, -1.0, np.nan, np.inf]:
        with pytest.raises(ValueError, match="ic50_m"):
            pIC50(value)


def test_assay_measurement_reports_p_value_and_type():
    measurement = AssayMeasurement(
        1e-6,
        "IC50",
        source="assay A",
        protocol_id="protocol-1",
    )

    assert measurement.p_value == pytest.approx(6.0)
    assert measurement.measurement_type == "IC50"
    assert measurement.protocol_id == "protocol-1"


def test_assay_measurement_validation():
    with pytest.raises(ValueError, match="value_m"):
        AssayMeasurement(0.0, "IC50")
    with pytest.raises(ValueError, match="measurement_type"):
        AssayMeasurement(1e-6, "unknown")


def test_ligand_efficiency_uses_pic50_and_heavy_atom_count():
    expected = 0.008314462618 * 298.15 * np.log(10.0) * 7.0 / 30.0

    assert ligand_efficiency(7.0, 30) == pytest.approx(expected)


def test_ligand_efficiency_validation():
    with pytest.raises(ValueError, match="heavy_atom_count"):
        ligand_efficiency(7.0, 0)
    with pytest.raises(ValueError, match="pic50"):
        ligand_efficiency(np.nan, 30)
    with pytest.raises(ValueError, match="temperature"):
        ligand_efficiency(7.0, 30, temperature=0.0)


def test_compound_record_validates_score_provenance_and_descriptors():
    with pytest.raises(ValueError, match="docking_score_unit"):
        make_compound("A", score=-7.0, score_unit=None)
    with pytest.raises(ValueError, match="docking_score_source"):
        make_compound("A", score=-7.0, score_source="")
    with pytest.raises(ValueError, match="finite"):
        make_compound("A", descriptors={"MW": np.nan})


def test_matched_pair_validation():
    a = make_compound("A")
    b = make_compound("B")

    with pytest.raises(ValueError, match="different IDs"):
        MatchedMolecularPair(a, a, "H -> F")
    with pytest.raises(ValueError, match="transformation"):
        MatchedMolecularPair(a, b, "")


def test_matched_pair_comparison_reports_categories_separately():
    a = make_compound(
        "A",
        activity=AssayMeasurement(1e-6, "IC50"),
        descriptors={"MW": 300.0, "cLogP": 2.0, "TPSA": 60.0},
        score=-7.2,
    )
    b = make_compound(
        "B",
        activity=AssayMeasurement(1e-7, "IC50"),
        descriptors={"MW": 320.0, "cLogP": 2.5, "TPSA": 60.0},
        score=-7.8,
    )
    pair = MatchedMolecularPair(
        a,
        b,
        "phenyl H -> F",
        hypothesis="fluorine may improve packing",
    )

    result = compare_matched_pair(pair)

    assert result.compound_a_id == "A"
    assert result.compound_b_id == "B"
    assert result.activity_type == "IC50"
    assert result.p_activity_a == pytest.approx(6.0)
    assert result.p_activity_b == pytest.approx(7.0)
    assert result.delta_p_activity == pytest.approx(1.0)
    assert result.fold_change == pytest.approx(10.0)
    assert result.descriptor_deltas == {
        "MW": 20.0,
        "cLogP": 0.5,
        "TPSA": 0.0,
    }
    assert result.docking_score_delta == pytest.approx(-0.6)
    assert result.hypothesis == "fluorine may improve packing"


def test_missing_activity_yields_missing_delta():
    a = make_compound("A", activity=AssayMeasurement(1e-6, "IC50"))
    b = make_compound("B")
    pair = MatchedMolecularPair(a, b, "H -> F")

    result = compare_matched_pair(pair)

    assert result.p_activity_a == pytest.approx(6.0)
    assert result.p_activity_b is None
    assert result.delta_p_activity is None
    assert result.fold_change is None
    assert result.activity_type is None


def test_activity_types_must_match():
    a = make_compound("A", activity=AssayMeasurement(1e-6, "IC50"))
    b = make_compound("B", activity=AssayMeasurement(1e-7, "Ki"))

    with pytest.raises(ValueError, match="measurement_type"):
        compare_matched_pair(MatchedMolecularPair(a, b, "H -> F"))


def test_assay_context_is_preserved_and_flagged():
    same_context_a = make_compound(
        "A",
        activity=AssayMeasurement(
            1e-6,
            "IC50",
            source="assay-1",
            protocol_id="protocol-1",
        ),
    )
    same_context_b = make_compound(
        "B",
        activity=AssayMeasurement(
            1e-7,
            "IC50",
            source="assay-1",
            protocol_id="protocol-1",
        ),
    )
    different_source = make_compound(
        "C",
        activity=AssayMeasurement(
            1e-7,
            "IC50",
            source="assay-2",
            protocol_id="protocol-1",
        ),
    )
    different_protocol = make_compound(
        "D",
        activity=AssayMeasurement(
            1e-7,
            "IC50",
            source="assay-1",
            protocol_id="protocol-2",
        ),
    )

    same = compare_matched_pair(
        MatchedMolecularPair(same_context_a, same_context_b, "H -> F")
    )
    source = compare_matched_pair(
        MatchedMolecularPair(same_context_a, different_source, "H -> F")
    )
    protocol = compare_matched_pair(
        MatchedMolecularPair(same_context_a, different_protocol, "H -> F")
    )

    assert same.assay_context_matches is True
    assert same.assay_sources == ("assay-1", "assay-1")
    assert same.assay_protocols == ("protocol-1", "protocol-1")
    assert source.assay_context_matches is False
    assert protocol.assay_context_matches is False
    assert same.delta_p_activity == pytest.approx(1.0)


def test_descriptor_sets_must_match():
    a = make_compound("A", descriptors={"MW": 300.0})
    b = make_compound("B", descriptors={"MW": 320.0, "cLogP": 2.5})

    with pytest.raises(ValueError, match="same descriptor set"):
        compare_matched_pair(MatchedMolecularPair(a, b, "H -> F"))


def test_docking_score_units_and_sources_must_match():
    a = make_compound("A", score=-7.0, score_unit="arbitrary", score_source="ranker-1")
    b = make_compound("B", score=-8.0, score_unit="kJ/mol", score_source="ranker-1")
    pair = MatchedMolecularPair(a, b, "H -> F")

    with pytest.raises(ValueError, match="units"):
        compare_matched_pair(pair)

    b2 = make_compound("B", score=-8.0, score_unit="arbitrary", score_source="ranker-2")
    with pytest.raises(ValueError, match="sources"):
        compare_matched_pair(MatchedMolecularPair(a, b2, "H -> F"))


def test_missing_docking_score_yields_missing_delta():
    a = make_compound("A", score=-7.0)
    b = make_compound("B")
    pair = MatchedMolecularPair(a, b, "H -> F")

    result = compare_matched_pair(pair)

    assert result.docking_score_delta is None
