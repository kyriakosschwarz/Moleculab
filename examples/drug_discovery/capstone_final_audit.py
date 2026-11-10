# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Final capstone integration, scientific audit, and controlled demonstration.

The example connects the structural records to a controlled statistical-
mechanics chain:

    configuration -> potential energy -> Boltzmann weight
    -> ensemble average -> exact/FEP free-energy difference

It then stops: the harmonic model is not a protein-ligand binding calculation.
The audit checks the provenance chain, learning-path documentation, structural
and thermodynamic boundaries, free-energy diagnostics, and reproducibility.

Run from the repository root:

    uv run python examples/drug_discovery/capstone_final_audit.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

from moleculab import (
    load_binding_site,
    load_chemical_state_enumeration,
    load_combination_record,
    load_harmonic_fep_study,
    load_hydrogen_evidence,
    load_interaction_robustness,
    load_ligand_coordinates,
    load_ligand_identity,
    load_manifest,
    load_protein_preparation,
    load_thermodynamic_interpretation,
    load_water_evidence,
    run_capstone_audit,
    run_end_to_end_demonstration,
    validate_scientific_audit,
    write_scientific_audit,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUTPUT = ROOT / "data/free_energy/capstone-scientific-audit.json"


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = _argument_parser().parse_args()
    paths = {
        "manifest": ROOT / "data/structures/2HYY.manifest.json",
        "identity": ROOT / "data/structures/2HYY.ligand-identity.json",
        "ligand": ROOT / "data/structures/2HYY.ligand-coordinates.json",
        "site": ROOT / "data/structures/2HYY.binding-site.json",
        "states": ROOT / "data/structures/2HYY.chemical-states.json",
        "protein": ROOT / "data/structures/2HYY.protein-preparation.json",
        "hydrogen": ROOT / "data/structures/2HYY.hydrogen-evidence.json",
        "water": ROOT / "data/structures/2HYY.water-evidence.json",
        "combination": ROOT / "data/structures/2HYY.combination-interactions.json",
        "robustness": ROOT / "data/structures/2HYY.interaction-robustness.json",
        "thermodynamics": (
            ROOT / "data/structures/2HYY.thermodynamic-interpretation.json"
        ),
        "free_energy_study": ROOT / "data/free_energy/harmonic-fep-study.json",
        "cif": ROOT / "data/structures/2HYY.cif",
    }
    for name, path in paths.items():
        if not path.exists():
            raise SystemExit(
                f"missing {name} record at {path}; run the earlier capstone batches first."
            )

    manifest = load_manifest(paths["manifest"])
    identity = load_ligand_identity(paths["identity"])
    ligand = load_ligand_coordinates(paths["ligand"])
    site = load_binding_site(paths["site"])
    states = load_chemical_state_enumeration(paths["states"])
    protein = load_protein_preparation(paths["protein"])
    hydrogen = load_hydrogen_evidence(paths["hydrogen"])
    water = load_water_evidence(paths["water"])
    combination = load_combination_record(paths["combination"])
    robustness = load_interaction_robustness(paths["robustness"])
    thermodynamics = load_thermodynamic_interpretation(paths["thermodynamics"])
    free_energy_study = load_harmonic_fep_study(paths["free_energy_study"])
    if len(
        {
            identity.accession,
            ligand.accession,
            site.accession,
            states.accession,
            protein.accession,
            hydrogen.accession,
            water.accession,
            combination.accession,
            robustness.accession,
            thermodynamics.accession,
            manifest.source.accession,
        }
    ) != 1:
        raise SystemExit("capstone records disagree on accession.")

    demonstration = run_end_to_end_demonstration(
        temperature_kelvin=300.0,
        n_samples=50_000,
        seed=2024,
    )
    record = run_capstone_audit(
        root=ROOT,
        manifest_path=paths["manifest"],
        identity_path=paths["identity"],
        ligand_path=paths["ligand"],
        site_path=paths["site"],
        states_path=paths["states"],
        protein_path=paths["protein"],
        hydrogen_path=paths["hydrogen"],
        water_path=paths["water"],
        combination_path=paths["combination"],
        robustness_path=paths["robustness"],
        thermodynamics_path=paths["thermodynamics"],
        free_energy_study_path=paths["free_energy_study"],
        cif_path=paths["cif"],
        manifest=manifest,
        identity=identity,
        ligand=ligand,
        site=site,
        states=states,
        protein=protein,
        hydrogen=hydrogen,
        water=water,
        combination=combination,
        robustness=robustness,
        thermodynamics=thermodynamics,
        free_energy_study=free_energy_study,
        demonstration=demonstration,
    )
    errors = validate_scientific_audit(record)
    if errors:
        for error in errors:
            print(f"scientific audit error: {error}")
        raise SystemExit(2)
    destination = write_scientific_audit(args.output, record)

    print("Controlled configuration -> ensemble -> free-energy chain")
    print("  model: 3D harmonic oscillator (not a protein-ligand system)")
    print(
        f"  sampled mean potential energy: "
        f"{record.demonstration.mean_potential_energy_kj_per_mol:.6f} kJ/mol"
    )
    print(
        f"  equipartition reference:       "
        f"{record.demonstration.analytic_mean_potential_energy_kj_per_mol:.6f} kJ/mol"
    )
    print(
        f"  proposal effective fraction:   "
        f"{record.demonstration.proposal_effective_sample_fraction:.6f}"
    )
    print(
        f"  exact delta G:                 "
        f"{record.demonstration.exact_delta_g_kj_per_mol:.6f} kJ/mol"
    )
    print(
        f"  FEP delta G:                   "
        f"{record.demonstration.fep_delta_g_kj_per_mol:.6f} "
        f"+/- {record.demonstration.fep_standard_error_kj_per_mol:.6f} kJ/mol"
    )
    print(
        f"  effective sample fraction:     "
        f"{record.demonstration.overlap.effective_sample_fraction:.6f}"
    )
    print(
        f"  two-state populations:         "
        f"p0={record.demonstration.two_state_reference.population0:.6f}, "
        f"p1={record.demonstration.two_state_reference.population1:.6f}"
    )
    print(
        f"  two-state delta G:             "
        f"{record.demonstration.two_state_reference.delta_g_kj_per_mol:.6f} kJ/mol"
    )
    print()
    print("Stop-line")
    for item in record.demonstration.stop_line:
        print(f"  - {item}")
    print()
    print(
        f"Structural boundary: free_energy_status="
        f"{thermodynamics.free_energy_status}, "
        f"binding_free_energy={thermodynamics.binding_free_energy_kj_per_mol}"
    )
    print(
        f"  {len(thermodynamics.boundaries)} boundaries, "
        f"{len(thermodynamics.terms)} missing terms, "
        f"{sum(1 for item in thermodynamics.requirements if item.status != 'met')} "
        "unmet requirements"
    )
    print()
    print(
        f"Audit: status={record.status}, pass={record.pass_count}, "
        f"fail={record.fail_count}, not_run={record.not_run_count}"
    )
    if record.fail_count:
        for check in record.checks:
            if check.status == "fail":
                print(f"  FAILED {check.check_id}: {check.observed}")
    print(f"  record: {destination}")
    print("  scope: audit and controlled demonstration; no binding affinity is computed")


if __name__ == "__main__":
    main()
