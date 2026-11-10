"""Reproducible docking benchmarks over pose-recovery results.

A benchmark record keeps valid failures and invalid cases separate. Invalid
cases are reported but do not enter the pose-recovery success-rate denominator.
The summary does not collapse geometry, interactions, and recovery into a
composite docking-quality score.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass

from .validation import PoseRecovery

RECOVERED = "recovered"
NOT_RECOVERED = "not_recovered"
INVALID = "invalid"
_VALID_STATUSES = (RECOVERED, NOT_RECOVERED)
_STATUSES = (RECOVERED, NOT_RECOVERED, INVALID)


@dataclass(frozen=True)
class BenchmarkCase:
    """One benchmark case with an explicit validity status.

    ``recovered`` and ``not_recovered`` are valid cases and must carry a
    ``PoseRecovery`` whose ``passed`` flag agrees with the status. ``invalid``
    means the case could not be evaluated, for example because of failed
    preparation or missing reference atoms; it must not carry a pose-recovery
    result.
    """

    case_id: str
    status: str
    pose_recovery: PoseRecovery | None = None
    notes: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.case_id, str) or not self.case_id.strip():
            raise ValueError("case_id must be a non-empty string.")
        if self.status not in _STATUSES:
            raise ValueError(
                f"status must be one of {_STATUSES}, got {self.status!r}."
            )
        if not isinstance(self.notes, str):
            raise TypeError("notes must be a string.")

        if self.status in _VALID_STATUSES:
            if self.status == RECOVERED:
                if not isinstance(self.pose_recovery, PoseRecovery):
                    raise ValueError(
                        "recovered cases must carry a PoseRecovery result."
                    )
                if not self.pose_recovery.passed:
                    raise ValueError(
                        "benchmark status does not match PoseRecovery.passed."
                    )
            else:
                if self.pose_recovery is not None:
                    if not isinstance(self.pose_recovery, PoseRecovery):
                        raise ValueError(
                            "pose_recovery must be a PoseRecovery result or None."
                        )
                    if self.pose_recovery.passed:
                        raise ValueError(
                            "benchmark status does not match PoseRecovery.passed."
                        )
        elif self.pose_recovery is not None:
            raise ValueError("invalid benchmark cases must not carry PoseRecovery.")

    @classmethod
    def from_pose_recovery(
        cls,
        case_id: str,
        pose_recovery: PoseRecovery,
        notes: str = "",
    ) -> BenchmarkCase:
        """Build a valid case whose status follows the recovery result."""

        if not isinstance(pose_recovery, PoseRecovery):
            raise TypeError("pose_recovery must be a PoseRecovery result.")
        return cls(
            case_id=case_id,
            status=RECOVERED if pose_recovery.passed else NOT_RECOVERED,
            pose_recovery=pose_recovery,
            notes=notes,
        )

    @classmethod
    def invalid(cls, case_id: str, notes: str = "") -> BenchmarkCase:
        """Build an invalid case that must not enter the denominator."""

        return cls(case_id=case_id, status=INVALID, notes=notes)

    @classmethod
    def method_failure(cls, case_id: str, notes: str = "") -> BenchmarkCase:
        """Build a valid failure with no pose-recovery result.

        Use this when the protocol says a computational failure, such as a run
        that produced no pose, counts as a method failure. It enters the
        denominator but contributes no RMSD.
        """

        return cls(
            case_id=case_id,
            status=NOT_RECOVERED,
            pose_recovery=None,
            notes=notes,
        )


@dataclass(frozen=True)
class DockingBenchmark:
    """A named benchmark with one protocol and a set of cases."""

    name: str
    protocol: Mapping[str, str]
    cases: tuple[BenchmarkCase, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("name must be a non-empty string.")
        if not isinstance(self.protocol, Mapping):
            raise TypeError("protocol must be a mapping of strings.")
        protocol = dict(self.protocol)
        for key, value in protocol.items():
            if not isinstance(key, str) or not isinstance(value, str):
                raise TypeError("protocol keys and values must be strings.")
        cases = tuple(self.cases)
        if not cases:
            raise ValueError("cases must contain at least one benchmark case.")
        for case in cases:
            if not isinstance(case, BenchmarkCase):
                raise TypeError("cases must contain BenchmarkCase objects.")
        identifiers = [case.case_id for case in cases]
        if len(set(identifiers)) != len(identifiers):
            raise ValueError("case_id values must be unique within a benchmark.")
        object.__setattr__(self, "protocol", protocol)
        object.__setattr__(self, "cases", cases)


@dataclass(frozen=True)
class BenchmarkSummary:
    """Aggregate validation counts with the protocol and cases retained."""

    name: str
    protocol: Mapping[str, str]
    cases: tuple[BenchmarkCase, ...]
    n_total: int
    n_valid: int
    n_recovered: int
    n_not_recovered: int
    n_invalid: int
    success_rate: float | None

    @property
    def valid_cases(self) -> tuple[BenchmarkCase, ...]:
        """Cases that entered the denominator."""

        return tuple(case for case in self.cases if case.status in _VALID_STATUSES)

    @property
    def invalid_cases(self) -> tuple[BenchmarkCase, ...]:
        """Cases that were reported but excluded from the denominator."""

        return tuple(case for case in self.cases if case.status == INVALID)

    @property
    def rmsd_values(self) -> tuple[float, ...]:
        """Per-case RMSD values for valid cases, in benchmark order."""

        return tuple(
            case.pose_recovery.rmsd
            for case in self.valid_cases
            if case.pose_recovery is not None
        )


def summarize_benchmark(benchmark: DockingBenchmark) -> BenchmarkSummary:
    """Summarize recovered, not-recovered, and invalid benchmark cases.

    The success rate is

    .. math::

        \\text{success rate} = \\frac{n_{\\mathrm{recovered}}}
        {n_{\\mathrm{recovered}} + n_{\\mathrm{not\\ recovered}}}.

    Invalid cases are counted and retained in the summary but do not enter the
    denominator. If no valid cases exist, the success rate is ``None``.
    """

    if not isinstance(benchmark, DockingBenchmark):
        raise TypeError("benchmark must be a DockingBenchmark.")

    recovered = sum(case.status == RECOVERED for case in benchmark.cases)
    not_recovered = sum(case.status == NOT_RECOVERED for case in benchmark.cases)
    invalid = sum(case.status == INVALID for case in benchmark.cases)
    valid = recovered + not_recovered
    success_rate = None if valid == 0 else recovered / valid
    return BenchmarkSummary(
        name=benchmark.name,
        protocol=dict(benchmark.protocol),
        cases=benchmark.cases,
        n_total=len(benchmark.cases),
        n_valid=valid,
        n_recovered=recovered,
        n_not_recovered=not_recovered,
        n_invalid=invalid,
        success_rate=success_rate,
    )
