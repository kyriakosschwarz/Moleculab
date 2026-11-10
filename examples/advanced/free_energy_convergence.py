"""Compare exact and FEP free energies on good- and poor-overlap harmonics.

The model is a three-dimensional harmonic oscillator.  Two transformations are
studied:

* a mild stiffness change, where the two states overlap well;
* a pure center shift, where the exact free-energy difference is zero but the
  forward FEP exponential average is dominated by rare configurations.

The example reports exact partition functions, sample-size convergence,
effective sample size, weight concentration, and replicate uncertainty.  It is
not an imatinib binding calculation.

Run:

    uv run python examples/advanced/free_energy_convergence.py
"""

from __future__ import annotations

from pathlib import Path

from moleculab import (
    HarmonicPotential,
    build_harmonic_fep_study_record,
    run_harmonic_fep_case,
    two_state_reference,
    validate_harmonic_fep_study,
    write_harmonic_fep_study,
)

OUTPUT = Path("data/free_energy/harmonic-fep-study.json")
TEMPERATURE = 300.0
SAMPLE_SIZES = (1_000, 10_000, 100_000, 1_000_000)
REPLICATE_COUNT = 20
REPLICATE_SAMPLE_SIZE = 20_000


def main() -> None:
    mild0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    mild1 = HarmonicPotential(stiffness=200.0, center=[0.0, 0.0, 0.0])
    poor0 = HarmonicPotential(stiffness=100.0, center=[0.0, 0.0, 0.0])
    poor1 = HarmonicPotential(stiffness=100.0, center=[0.5, 0.0, 0.0])

    cases = (
        run_harmonic_fep_case(
            case_id="mild-stiffness-change",
            label="Good overlap: k0=100, k1=200 kJ/(mol nm^2)",
            potential0=mild0,
            potential1=mild1,
            temperature=TEMPERATURE,
            sample_sizes=SAMPLE_SIZES,
            n_replicates=REPLICATE_COUNT,
            replicate_sample_size=REPLICATE_SAMPLE_SIZE,
            seed=2024,
            note="The exact free-energy change is positive and the distributions overlap well.",
        ),
        run_harmonic_fep_case(
            case_id="poor-overlap-center-shift",
            label="Poor overlap: same stiffness, 0.5 nm center shift",
            potential0=poor0,
            potential1=poor1,
            temperature=TEMPERATURE,
            sample_sizes=SAMPLE_SIZES,
            n_replicates=REPLICATE_COUNT,
            replicate_sample_size=REPLICATE_SAMPLE_SIZE,
            seed=4048,
            note=(
                "The exact free-energy difference is zero because the center shift "
                "does not change the partition function; forward FEP still has "
                "poor overlap and high variance."
            ),
        ),
    )
    two_state = two_state_reference(
        energy0_kj_per_mol=0.0,
        energy1_kj_per_mol=5.0,
        degeneracy0=1,
        degeneracy1=100,
        temperature=TEMPERATURE,
    )
    record = build_harmonic_fep_study_record(
        temperature_kelvin=TEMPERATURE,
        two_state_example=two_state,
        sample_sizes=SAMPLE_SIZES,
        replicate_count=REPLICATE_COUNT,
        replicate_sample_size=REPLICATE_SAMPLE_SIZE,
        cases=cases,
    )
    errors = validate_harmonic_fep_study(record)
    if errors:
        for error in errors:
            print(f"free-energy study error: {error}")
        raise SystemExit(2)
    destination = write_harmonic_fep_study(OUTPUT, record)

    print("Harmonic free-energy study")
    print(f"  temperature: {record.temperature_kelvin:.1f} K")
    print(f"  record: {destination}")
    print()
    print("Exact two-state reference")
    print(
        f"  E0={record.two_state_example.energy0_kj_per_mol:.1f}, "
        f"E1={record.two_state_example.energy1_kj_per_mol:.1f} kJ/mol; "
        f"g0={record.two_state_example.degeneracy0}, "
        f"g1={record.two_state_example.degeneracy1}"
    )
    print(
        f"  Z={record.two_state_example.partition_function:.6f}; "
        f"p0={record.two_state_example.population0:.6f}, "
        f"p1={record.two_state_example.population1:.6f}; "
        f"dG={record.two_state_example.delta_g_kj_per_mol:+.6f} kJ/mol"
    )
    print()
    for case in record.cases:
        print(case.label)
        print(
            f"  exact dG: {case.exact_delta_g_kj_per_mol:+.6f} kJ/mol; "
            f"Z0={case.exact_partition_function0:.6e} nm^3; "
            f"Z1={case.exact_partition_function1:.6e} nm^3"
        )
        print(
            f"  {'samples':>10} {'dG':>12} {'std err':>10} {'error':>10} "
            f"{'ESS':>12} {'max w':>9}"
        )
        for point in case.convergence:
            print(
                f"  {point.n_samples:10d} "
                f"{point.delta_g_kj_per_mol:12.6f} "
                f"{point.standard_error_kj_per_mol:10.6f} "
                f"{point.error_kj_per_mol:10.6f} "
                f"{point.overlap.effective_sample_size:12.1f} "
                f"{point.overlap.maximum_weight_fraction:9.5f}"
            )
        replicates = case.replicates
        print(
            f"  replicates: mean={replicates.mean_delta_g_kj_per_mol:+.6f}, "
            f"std={replicates.standard_deviation_kj_per_mol:.6f}, "
            f"range=[{replicates.minimum_delta_g_kj_per_mol:+.6f}, "
            f"{replicates.maximum_delta_g_kj_per_mol:+.6f}], "
            f"95% CI coverage={replicates.confidence_interval_coverage:.2f}"
        )
        print()

    print("Lessons")
    print("  More samples reduce random error when the states overlap well.")
    print("  Poor overlap can leave a small effective sample size even at 1,000,000 samples.")
    print("  A few high-weight configurations can dominate the exponential average.")
    print("  A zero exact dG does not make a center-shift FEP calculation easy.")
    print("  The harmonic result is not an imatinib binding affinity.")


if __name__ == "__main__":
    main()
