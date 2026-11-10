# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for docking benchmark records and aggregate metrics."""

import pytest

from moleculab import (
    BenchmarkCase,
    DockingBenchmark,
    PoseRecovery,
    summarize_benchmark,
)


def make_recovery(passed: bool, rmsd: float | None = None) -> PoseRecovery:
    return PoseRecovery(
        rmsd=0.05 if rmsd is None and passed else (0.25 if rmsd is None else rmsd),
        threshold=0.10,
        passed=passed,
        best_mapping=((0, 0),),
        mapping_index=0,
        n_mappings=1,
        pose_score=-7.0,
        score_unit="arbitrary",
        score_source="example ranker",
    )


def make_benchmark(cases):
    return DockingBenchmark(
        name="synthetic redocking",
        protocol={"program": "example-docker", "version": "1.0"},
        cases=tuple(cases),
    )


def test_success_rate_uses_valid_cases_as_denominator():
    cases = [
        BenchmarkCase.from_pose_recovery("case-1", make_recovery(True)),
        BenchmarkCase.from_pose_recovery("case-2", make_recovery(True)),
        BenchmarkCase.from_pose_recovery("case-3", make_recovery(True)),
        BenchmarkCase.from_pose_recovery("case-4", make_recovery(False)),
        BenchmarkCase.invalid("case-5", notes="missing reference atoms"),
    ]

    result = summarize_benchmark(make_benchmark(cases))

    assert result.success_rate == pytest.approx(0.75)
    assert result.n_total == 5
    assert result.n_valid == 4
    assert result.n_recovered == 3
    assert result.n_not_recovered == 1
    assert result.n_invalid == 1
    assert result.rmsd_values == pytest.approx((0.05, 0.05, 0.05, 0.25))


def test_invalid_cases_are_retained_but_excluded_from_denominator():
    invalid = BenchmarkCase.invalid("case-2", notes="preparation failed")
    benchmark = make_benchmark(
        [
            BenchmarkCase.from_pose_recovery("case-1", make_recovery(True)),
            invalid,
        ]
    )

    result = summarize_benchmark(benchmark)

    assert result.n_invalid == 1
    assert result.invalid_cases == (invalid,)
    assert result.n_valid == 1
    assert result.success_rate == pytest.approx(1.0)


def test_no_valid_cases_gives_no_success_rate():
    benchmark = make_benchmark(
        [
            BenchmarkCase.invalid("case-1", notes="preparation failed"),
            BenchmarkCase.invalid("case-2", notes="missing atoms"),
        ]
    )

    result = summarize_benchmark(benchmark)

    assert result.n_valid == 0
    assert result.success_rate is None


def test_method_failure_enters_denominator_without_an_rmsd():
    case = BenchmarkCase.method_failure(
        "case-1",
        notes="docking run produced no pose",
    )

    result = summarize_benchmark(make_benchmark([case]))

    assert case.status == "not_recovered"
    assert case.pose_recovery is None
    assert result.n_valid == 1
    assert result.n_not_recovered == 1
    assert result.success_rate == 0.0
    assert result.rmsd_values == ()


def test_status_must_match_pose_recovery_result():
    with pytest.raises(ValueError, match="does not match"):
        BenchmarkCase(
            case_id="case-1",
            status="recovered",
            pose_recovery=make_recovery(False),
        )
    with pytest.raises(ValueError, match="does not match"):
        BenchmarkCase(
            case_id="case-1",
            status="not_recovered",
            pose_recovery=make_recovery(True),
        )


def test_invalid_case_must_not_carry_pose_recovery():
    with pytest.raises(ValueError, match="must not carry"):
        BenchmarkCase(
            case_id="case-1",
            status="invalid",
            pose_recovery=make_recovery(False),
        )


def test_duplicate_case_ids_are_rejected():
    with pytest.raises(ValueError, match="unique"):
        make_benchmark(
            [
                BenchmarkCase.from_pose_recovery("case-1", make_recovery(True)),
                BenchmarkCase.from_pose_recovery("case-1", make_recovery(False)),
            ]
        )


def test_empty_benchmark_is_rejected():
    with pytest.raises(ValueError, match="at least one"):
        make_benchmark([])


def test_protocol_and_name_are_preserved_in_summary():
    benchmark = make_benchmark(
        [BenchmarkCase.from_pose_recovery("case-1", make_recovery(True))]
    )

    result = summarize_benchmark(benchmark)

    assert result.name == "synthetic redocking"
    assert result.protocol == {"program": "example-docker", "version": "1.0"}
    assert result.cases == benchmark.cases


def test_from_pose_recovery_sets_status_and_preserves_notes():
    recovered = BenchmarkCase.from_pose_recovery(
        "case-1",
        make_recovery(True),
        notes="easy target",
    )
    failed = BenchmarkCase.from_pose_recovery(
        "case-2",
        make_recovery(False),
        notes="wrong tautomer",
    )

    assert recovered.status == "recovered"
    assert failed.status == "not_recovered"
    assert recovered.notes == "easy target"
    assert failed.notes == "wrong tautomer"


def test_invalid_factory_preserves_notes():
    case = BenchmarkCase.invalid("case-1", notes="missing reference atoms")

    assert case.status == "invalid"
    assert case.pose_recovery is None
    assert case.notes == "missing reference atoms"


def test_case_id_and_status_validation():
    with pytest.raises(ValueError, match="case_id"):
        BenchmarkCase(case_id="", status="invalid")
    with pytest.raises(ValueError, match="status"):
        BenchmarkCase(case_id="case-1", status="unknown")
    with pytest.raises(ValueError, match="PoseRecovery"):
        BenchmarkCase(case_id="case-1", status="recovered")


def test_summary_requires_a_benchmark_object():
    with pytest.raises(TypeError, match="DockingBenchmark"):
        summarize_benchmark("not a benchmark")
