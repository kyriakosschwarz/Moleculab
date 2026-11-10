"""Docking benchmark records with explicit invalid cases.

Run:

    uv run python examples/docking_benchmark.py
"""

from __future__ import annotations

from moleculab import (
    BenchmarkCase,
    DockingBenchmark,
    PoseRecovery,
    summarize_benchmark,
)

PROTOCOL = {
    "program": "example-docker",
    "version": "1.0",
    "scoring": "example ranker",
    "receptor_preparation": "synthetic",
}


def make_recovery(passed: bool, rmsd: float) -> PoseRecovery:
    return PoseRecovery(
        rmsd=rmsd,
        threshold=0.10,
        passed=passed,
        best_mapping=((0, 0),),
        mapping_index=0,
        n_mappings=1,
        pose_score=-7.0 - rmsd,
        score_unit="arbitrary",
        score_source="example ranker",
    )


def main() -> None:
    cases = [
        BenchmarkCase.from_pose_recovery(
            "case-1",
            make_recovery(True, 0.04),
            notes="clear recovery",
        ),
        BenchmarkCase.from_pose_recovery(
            "case-2",
            make_recovery(True, 0.07),
            notes="recovered with a different symmetry mapping",
        ),
        BenchmarkCase.from_pose_recovery(
            "case-3",
            make_recovery(True, 0.09),
            notes="near the threshold",
        ),
        BenchmarkCase.from_pose_recovery(
            "case-4",
            make_recovery(False, 0.24),
            notes="wrong ligand tautomer",
        ),
        BenchmarkCase.method_failure(
            "case-5",
            notes="docking run produced no pose; valid method failure",
        ),
        BenchmarkCase.invalid(
            "case-6",
            notes="missing reference atoms; excluded from denominator",
        ),
    ]
    benchmark = DockingBenchmark(
        name="synthetic redocking benchmark",
        protocol=PROTOCOL,
        cases=tuple(cases),
    )
    summary = summarize_benchmark(benchmark)

    print(f"Benchmark: {summary.name}")
    print("Protocol metadata")
    for key, value in summary.protocol.items():
        print(f"  {key}: {value}")
    print()
    print(
        f"Success rate: {summary.success_rate:.3f} "
        f"({summary.n_recovered}/{summary.n_valid} valid cases)"
    )
    print(
        f"Total: {summary.n_total}; valid: {summary.n_valid}; "
        f"recovered: {summary.n_recovered}; "
        f"not recovered: {summary.n_not_recovered}; invalid: {summary.n_invalid}"
    )
    print()
    print("Cases")
    for case in summary.cases:
        if case.pose_recovery is None:
            print(f"  {case.case_id}: {case.status} ({case.notes})")
        else:
            print(
                f"  {case.case_id}: {case.status} "
                f"RMSD={case.pose_recovery.rmsd:.3f} nm "
                f"score={case.pose_recovery.pose_score:.2f} "
                f"({case.pose_recovery.score_unit})"
            )
    print()
    print("Invalid cases are reported but do not enter the denominator.")
    print("Computational method failures remain valid cases in the denominator.")
    print("The success rate is a pose-recovery statistic, not an affinity metric.")

    assert summary.n_total == 6
    assert summary.n_valid == 5
    assert summary.n_invalid == 1
    assert summary.success_rate == 0.6


if __name__ == "__main__":
    main()
