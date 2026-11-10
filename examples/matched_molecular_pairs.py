"""Compare a matched molecular pair without collapsing data categories.

Run:

    uv run python examples/matched_molecular_pairs.py
"""

from __future__ import annotations

import numpy as np

from moleculab import (
    AssayMeasurement,
    CompoundRecord,
    MatchedMolecularPair,
    compare_matched_pair,
    ligand_efficiency,
    pIC50,
)


def main() -> None:
    compound_a = CompoundRecord(
        compound_id="A",
        smiles="example-A",
        activity=AssayMeasurement(
            1e-6,
            "IC50",
            source="assay-1",
            protocol_id="protocol-1",
        ),
        descriptors={
            "MW": 300.0,
            "cLogP": 2.0,
            "TPSA": 60.0,
            "HBD": 1.0,
            "HBA": 3.0,
            "rotatable_bonds": 4.0,
        },
        docking_score=-7.2,
        docking_score_unit="arbitrary",
        docking_score_source="example ranker",
    )
    compound_b = CompoundRecord(
        compound_id="B",
        smiles="example-B",
        activity=AssayMeasurement(
            1e-7,
            "IC50",
            source="assay-1",
            protocol_id="protocol-1",
        ),
        descriptors={
            "MW": 320.0,
            "cLogP": 2.5,
            "TPSA": 60.0,
            "HBD": 1.0,
            "HBA": 3.0,
            "rotatable_bonds": 4.0,
        },
        docking_score=-7.8,
        docking_score_unit="arbitrary",
        docking_score_source="example ranker",
    )
    pair = MatchedMolecularPair(
        compound_a,
        compound_b,
        transformation="phenyl H -> F",
        hypothesis="fluorine may improve shape complementarity, but the assay is needed",
    )
    result = compare_matched_pair(pair)

    print("Matched molecular pair")
    print(f"  transformation: {result.transformation}")
    print(f"  A: {result.compound_a_id}; B: {result.compound_b_id}")
    print()

    print("Measured activity")
    print(f"  IC50 A: {compound_a.activity.value_m * 1e9:.1f} nM")
    print(f"  IC50 B: {compound_b.activity.value_m * 1e9:.1f} nM")
    print(f"  pIC50 A: {pIC50(compound_a.activity.value_m):.3f}")
    print(f"  pIC50 B: {pIC50(compound_b.activity.value_m):.3f}")
    print(f"  delta pIC50: {result.delta_p_activity:+.3f}")
    print(f"  fold change in IC50: {result.fold_change:.3f}x")
    print(
        "  assay context matches: "
        f"{result.assay_context_matches} "
        f"({result.assay_sources[0]}/{result.assay_protocols[0]} vs "
        f"{result.assay_sources[1]}/{result.assay_protocols[1]})"
    )
    print()

    print("Calculated descriptors")
    for name, delta in result.descriptor_deltas.items():
        print(f"  delta {name}: {delta:+.2f}")
    print()

    print("Computational docking")
    print(f"  delta score: {result.docking_score_delta:+.3f} (arbitrary units)")
    print("  A docking score is a ranking quantity from one model.")
    print()

    efficiency_a = ligand_efficiency(result.p_activity_a, heavy_atom_count=20)
    efficiency_b = ligand_efficiency(result.p_activity_b, heavy_atom_count=22)
    print("Ligand efficiency (IC50 proxy)")
    print(f"  A: {efficiency_a:.3f} kJ/mol per heavy atom")
    print(f"  B: {efficiency_b:.3f} kJ/mol per heavy atom")
    print("  This is a size-normalised potency proxy, not a measured free energy.")
    print()

    print(f"Hypothesis: {result.hypothesis}")
    print("Measured activity, descriptors, docking scores, and hypotheses are separate.")
    print("IC50, Ki, and Kd are not interchangeable without explicit assumptions.")

    assert result.delta_p_activity == 1.0
    assert np.isclose(result.fold_change, 10.0)
    assert result.descriptor_deltas["MW"] == 20.0
    assert np.isclose(result.docking_score_delta, -0.6)


if __name__ == "__main__":
    main()
