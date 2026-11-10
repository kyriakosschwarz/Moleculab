# Track C, module 07 — Docking validation and benchmarking

> **Where this fits.** [Module 05](05_docking_fundamentals.md) introduced pose
> recovery. [Module 06](06_comparing_binding_modes.md) compared binding modes.
> This module asks a different question: is the docking protocol reliable
> across a set of cases?
>
> **What you should already know.** Track C modules 04–06, especially
> `PoseRecovery`, explicit atom mapping, symmetry-aware RMSD, and success
> thresholds.
>
> **What you will be able to do afterwards.** Design a reproducible benchmark,
> distinguish valid failures from invalid cases, report success rate with the
> correct denominator, preserve protocol metadata, and explain why a pose
> benchmark is not an affinity benchmark.

## 1. Learning objectives

**Core**

- Explain redocking and cross-docking.
- Compute a pose-recovery success rate.
- Distinguish valid failures from invalid benchmark cases.
- Preserve protocol metadata with benchmark results.
- Distinguish pose-prediction benchmarks from affinity benchmarks.

**Deeper insight**

- The denominator determines what the success rate means.
- Benchmark composition and selection bias affect apparent performance.
- A pose-recovery threshold is protocol-dependent.
- A docking score is not an affinity prediction.
- Aggregate metrics must not hide per-case geometry and interaction evidence.

**Looking ahead**

- Module 08 will analyse interaction persistence over trajectories.
- Module 09 will move towards multi-ligand comparison and SAR.
- The capstone will integrate the full track.

## 2. What validation asks

Docking validation asks:

> How reliably does this protocol recover a known or plausible binding mode
> across a defined set of cases?

It does not ask:

- whether the top-ranked pose is always correct;
- whether the score equals a binding free energy;
- whether the protocol works for every protein family;
- whether a low RMSD proves the chemically correct interactions.

Validation is a property of a protocol applied to a benchmark, not a property
of a single successful pose.

## 3. Redocking

In redocking:

1. Start from a protein–ligand complex.
2. Remove the ligand.
3. Dock the ligand back into the receptor.
4. Compare the predicted pose with the experimental pose.
5. Report recovery using a specified RMSD threshold and atom mapping.

Redocking tests whether the protocol can recover a known pose when the receptor
conformation is already appropriate for that ligand.

It is a useful sanity check, but it is an easier problem than prospective
docking:

- the receptor conformation is already the bound one;
- the ligand chemistry is known;
- the binding site is known;
- the reference pose is available for threshold selection.

A high redocking success rate does not guarantee prospective performance.

## 4. Cross-docking

In cross-docking, the ligand is docked into a receptor structure from a
different complex or conformation.

Cross-docking tests:

- sensitivity to receptor conformation;
- the effect of induced fit;
- whether the protocol can handle a receptor that is not pre-organised for the
  ligand.

Common outcomes:

- the ligand recovers the site but not the exact pose;
- the ligand selects a different subpocket;
- the ligand fails because the receptor is occluded;
- the ligand succeeds only after receptor preparation or side-chain relaxation.

Cross-docking is closer to a realistic recognition problem but still differs
from prospective virtual screening.

## 5. Pose-recovery success rate

For a benchmark with $n_{\mathrm{recovered}}$ successful pose recoveries and
$n_{\mathrm{not\ recovered}}$ valid failures,

$$ \text{success rate} = \frac{n_{\mathrm{recovered}}}{n_{\mathrm{recovered}} + n_{\mathrm{not\ recovered}}}. $$

The denominator deliberately excludes invalid cases.

Why? A case that failed because of:

- a broken receptor preparation;
- missing reference atoms;
- an unusable ligand record;
- an incompatible atom mapping;
- a corrupted coordinate file;

does not tell you whether the docking protocol can recover a pose. It tells you
that the benchmark case was not evaluable.

Invalid cases must still be:

- counted;
- reported;
- retained in the record;
- explained in notes or metadata.

Silently dropping them inflates the apparent success rate.

## 6. Benchmark composition and selection bias

A benchmark is not neutral. Its composition determines what its success rate
means.

Important choices include:

- which protein families are included;
- whether only easy, well-resolved structures are used;
- whether the same ligand appears in training and test data;
- whether receptor conformations are bound or apo;
- whether waters, ions, and cofactors are retained;
- whether ligands are diverse or congeneric;
- how many cases are included;
- how the RMSD threshold is chosen.

A benchmark dominated by easy, rigid, well-resolved complexes can report a high
success rate that does not transfer to difficult targets.

### 6.1 Data leakage

If the protocol or its parameters were tuned on the same cases used for
evaluation, the benchmark is not independent. This is a form of data leakage.

Keep protocol development and benchmark evaluation separate, or report that the
benchmark is retrospective and tuned.

## 7. Protocol metadata and reproducibility

A success rate without protocol metadata is not reproducible.

Record at least:

- docking program and version;
- scoring function and version;
- receptor preparation protocol;
- ligand preparation and protonation protocol;
- retained or removed waters, ions, and cofactors;
- ligand flexibility settings;
- number of poses generated;
- RMSD threshold;
- atom-mapping and symmetry policy;
- benchmark composition and exclusions.

The `DockingBenchmark` record stores a protocol mapping as strings. The
`BenchmarkSummary` retains that mapping.

## 8. Pose prediction versus affinity benchmarks

Pose prediction and affinity prediction answer different questions.

| Benchmark | Question | Typical metric |
| --- | --- | --- |
| pose recovery | did the protocol reproduce a reference geometry? | RMSD success rate |
| interaction recovery | did it reproduce the key interactions? | fingerprint similarity or contact recovery |
| affinity prediction | did it estimate experimental binding strength? | correlation, error, ranking metrics |

A method can succeed at pose recovery and fail at affinity prediction. It can
also produce a chemically different pose that has a similar RMSD.

Do not report a pose-recovery success rate as an affinity result.

## 9. Failure modes and case statuses

A benchmark case is one of:

| Status | Meaning | Enters denominator? |
| --- | --- | --- |
| `recovered` | valid case; pose passed the recovery threshold | yes |
| `not_recovered` | valid case; pose failed or the run failed under the protocol | yes |
| `invalid` | case could not be evaluated | no |

A `not_recovered` case may carry a failed `PoseRecovery` (the pose was produced
but missed the threshold) or no `PoseRecovery` at all (a computational method
failure under the stated protocol). `BenchmarkCase.method_failure(...)` records
the latter. Do not classify a computational failure as `invalid` merely because
it produced no pose; that would let frequently failing methods appear more
successful.

Typical invalid causes:

- missing or incompatible reference atoms;
- failed receptor preparation;
- failed ligand preparation;
- ambiguous atom mapping;
- broken coordinates.

Typical valid failures:

- the search found a different pose;
- the scoring function preferred a wrong pose;
- the receptor conformation was wrong;
- the ligand state was wrong;
- a key water or cofactor was missing.

The distinction matters: a valid failure is evidence about the protocol; an
invalid case is evidence about the benchmark setup. A computational failure
that produced no pose is usually a valid failure, not an invalid case, unless
the benchmark protocol explicitly says otherwise.

## 10. Implementation and validation

The benchmark layer is in `src/moleculab/analysis/docking/benchmark.py`.

### 10.1 API

| Object or function | Purpose |
| --- | --- |
| `BenchmarkCase` | one case with explicit status and optional `PoseRecovery` |
| `BenchmarkCase.from_pose_recovery(...)` | derive recovered/not-recovered status |
| `BenchmarkCase.method_failure(...)` | record a valid computational failure without a pose |
| `BenchmarkCase.invalid(...)` | record a non-evaluable case |
| `DockingBenchmark` | named benchmark with one protocol and case list |
| `BenchmarkSummary` | counts, success rate, protocol, and retained cases |
| `summarize_benchmark(...)` | aggregate metrics with the valid-case denominator |

The summary does not combine RMSD, interaction fingerprints, score, and
occupancy into one "docking quality" number. Those metrics remain independently
visible in the case records.

### 10.2 Example

```bash
uv run python examples/drug_discovery/docking_benchmark.py
```

Observed output:

```text
Benchmark: synthetic redocking benchmark
Protocol metadata
  program: example-docker
  version: 1.0
  scoring: example ranker
  receptor_preparation: synthetic

Success rate: 0.600 (3/5 valid cases)
Total: 6; valid: 5; recovered: 3; not recovered: 2; invalid: 1

Cases
  case-1: recovered RMSD=0.040 nm score=-7.04 (arbitrary)
  case-2: recovered RMSD=0.070 nm score=-7.07 (arbitrary)
  case-3: recovered RMSD=0.090 nm score=-7.09 (arbitrary)
  case-4: not_recovered RMSD=0.240 nm score=-7.24 (arbitrary)
  case-5: not_recovered (docking run produced no pose; valid method failure)
  case-6: invalid (missing reference atoms; excluded from denominator)
```

### 10.3 Tests

```bash
uv run pytest tests/test_docking_benchmark.py -v
```

The tests cover:

- success-rate denominator;
- invalid-case exclusion and retention;
- computational method failures entering the denominator without an RMSD;
- no-valid-case behaviour;
- status/recovery consistency;
- duplicate case-ID rejection;
- empty benchmark rejection;
- protocol and name preservation;
- case factories and notes;
- input validation.

### 10.4 Acceptance criteria

A docking-benchmark module is correct when:

- success rate uses only valid cases in the denominator;
- invalid cases are counted, retained, and reported separately;
- computational method failures can enter the denominator without an RMSD;
- status is consistent with the underlying pose-recovery result;
- duplicate case IDs are rejected;
- protocol metadata is preserved;
- per-case RMSD and score provenance remain visible;
- empty benchmarks are rejected;
- the summary does not claim affinity prediction from pose recovery.

## 11. Exercises

### Reasoning

1. Why are invalid cases excluded from the success-rate denominator?
2. Why should invalid cases still be reported?
3. Why is redocking easier than cross-docking?
4. Why can a high redocking success rate overstate prospective performance?
5. Why is a pose-recovery benchmark not an affinity benchmark?
6. Why does benchmark composition affect the reported success rate?
7. Why is protocol metadata required for reproducibility?

### Calculation

1. A benchmark has 8 recovered cases, 2 valid failures, and 3 invalid cases.
   What is the success rate?
2. Repeat the calculation if the 3 invalid cases are incorrectly included in
   the denominator.
3. A benchmark has 0 valid cases. What should the summary report?
4. A protocol recovers 7 of 10 poses. What is the success rate?
5. A case has RMSD 0.12 nm and a threshold of 0.10 nm. What status should it
   receive?

### Implementation

1. Run `examples/drug_discovery/docking_benchmark.py` and reproduce the output.
2. Add a valid failure and recompute the success rate.
3. Add an invalid case and confirm the denominator is unchanged.
4. Add protocol metadata and confirm it survives summarization.
5. Attempt to create a `recovered` case with a failed `PoseRecovery` and explain
   the error.

## 12. Summary

- Docking validation measures protocol reliability across a defined benchmark.
- Redocking recovers a known pose in a bound receptor.
- Cross-docking tests receptor-conformation sensitivity.
- Success rate uses valid cases only:
  $n_{\mathrm{recovered}}/(n_{\mathrm{recovered}}+n_{\mathrm{not\ recovered}})$.
- Invalid cases are excluded from the denominator but retained and reported.
- Benchmark composition, selection bias, and data leakage affect interpretation.
- Protocol metadata is part of the result.
- Pose prediction and affinity prediction are different benchmarks.
- Aggregate success rates must not hide per-case geometry and interaction
  evidence.

**Previous:** [06 — Comparing binding modes](06_comparing_binding_modes.md).

**Next:** Trajectory contact analysis (module 08, planned).
