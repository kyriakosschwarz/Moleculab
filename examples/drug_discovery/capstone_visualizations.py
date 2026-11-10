# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Generate the Visualization Batch 1, 2, and 3 figures from the 2HYY records.

The figures are:
1. an RDKit 2D diagram of the reference-neutral imatinib state;
2. a PyMOL overview of imatinib in the ABL1 binding site;
3. a PyMOL binding-pocket close-up with measured contacts;
4. an RDKit comparison of selected chemical-state alternatives;
5. an interaction-evidence matrix across eight ligand and two protein states;
6. a PyMOL view of generated hydrogens, acceptor geometry, and deposited waters;
7. a multi-axis interaction-sensitivity plot;
8. an exact-versus-FEP harmonic convergence study;
9. a structural-evidence versus thermodynamic-requirement boundary.

Every figure is written with source, tool-version, camera, and limitation
metadata.  Geometry comes from the deposited structure and validated records;
no interaction or affinity is inferred.

Run from the repository root:

    uv run --with rdkit python examples/drug_discovery/capstone_visualizations.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

from moleculab import (
    FigureContext,
    draw_chemical_state_comparison,
    draw_ligand_2d,
    load_binding_site,
    load_chemical_state_enumeration,
    load_combination_record,
    load_figure_metadata,
    load_harmonic_fep_study,
    load_hydrogen_evidence,
    load_interaction_robustness,
    load_ligand_identity,
    load_manifest,
    load_thermodynamic_interpretation,
    load_water_evidence,
    plot_evidence_matrix,
    plot_free_energy_convergence,
    plot_interaction_sensitivity,
    plot_structural_thermodynamics_boundary,
    pymol_available,
    render_binding_pocket_closeup,
    render_hydrogen_water_evidence,
    render_protein_ligand_overview,
    sha256_file,
    validate_figure_metadata,
    write_figure_metadata,
)

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_MANIFEST = ROOT / "data/structures/2HYY.manifest.json"
DEFAULT_IDENTITY = ROOT / "data/structures/2HYY.ligand-identity.json"
DEFAULT_SITE = ROOT / "data/structures/2HYY.binding-site.json"
DEFAULT_STATES = ROOT / "data/structures/2HYY.chemical-states.json"
DEFAULT_STRUCTURE = ROOT / "data/structures/2HYY.cif"
DEFAULT_METADATA = Path("data/figures/capstone-visualizations.json")
FIGURE_DIRECTORY = Path("docs/tutorials/drug_discovery/visualizations")


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--site", type=Path, default=DEFAULT_SITE)
    parser.add_argument("--states", type=Path, default=DEFAULT_STATES)
    parser.add_argument("--structure", type=Path, default=DEFAULT_STRUCTURE)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--figure-directory", type=Path, default=FIGURE_DIRECTORY)
    return parser


def _existing_views(
    metadata_path: Path,
    structure_sha256: str,
) -> dict[str, tuple[float, ...]]:
    if not metadata_path.exists():
        return {}
    return {
        figure.figure_id: figure.camera_view
        for figure in load_figure_metadata(metadata_path)
        if figure.camera_view is not None
        and figure.context.structure_sha256 == structure_sha256
    }


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (
        args.manifest,
        args.identity,
        args.site,
        args.states,
        args.structure,
    ):
        if not path.exists():
            raise SystemExit(f"missing {path}; run the capstone structure batches first.")
    if not pymol_available():
        raise SystemExit(
            "PyMOL is not available. Install PyMOL to generate the 3D figures."
        )

    manifest = load_manifest(args.manifest)
    identity = load_ligand_identity(args.identity)
    binding_site = load_binding_site(args.site)
    states = load_chemical_state_enumeration(args.states)
    reference_state = next(
        state for state in states.states if state.state_id == "reference-neutral"
    )
    structure_sha256 = manifest.source.sha256
    views = _existing_views(args.metadata, structure_sha256)
    args.figure_directory.mkdir(parents=True, exist_ok=True)

    source_paths = {
        "manifest": "data/structures/2HYY.manifest.json",
        "identity": "data/structures/2HYY.ligand-identity.json",
        "site": "data/structures/2HYY.binding-site.json",
        "states": "data/structures/2HYY.chemical-states.json",
        "hydrogen": "data/structures/2HYY.hydrogen-evidence.json",
        "water": "data/structures/2HYY.water-evidence.json",
        "combination": "data/structures/2HYY.combination-interactions.json",
        "robustness": "data/structures/2HYY.interaction-robustness.json",
        "thermodynamics": (
            "data/structures/2HYY.thermodynamic-interpretation.json"
        ),
        "free_energy_study": "data/free_energy/harmonic-fep-study.json",
    }
    source_hashes = {
        name: sha256_file(ROOT / path) for name, path in source_paths.items()
    }

    def records(*names: str) -> tuple[str, ...]:
        return tuple(source_paths[name] for name in names)

    def checksums(*names: str) -> tuple[tuple[str, str], ...]:
        return tuple(
            (source_paths[name], source_hashes[name]) for name in names
        )

    hydrogen = load_hydrogen_evidence(ROOT / source_paths["hydrogen"])
    water = load_water_evidence(ROOT / source_paths["water"])
    combination = load_combination_record(ROOT / source_paths["combination"])
    robustness = load_interaction_robustness(ROOT / source_paths["robustness"])
    thermodynamics = load_thermodynamic_interpretation(
        ROOT / source_paths["thermodynamics"]
    )
    free_energy_study = load_harmonic_fep_study(
        ROOT / source_paths["free_energy_study"]
    )

    context = FigureContext(
        source_structure=manifest.source.accession,
        structure_sha256=structure_sha256,
        protein_chain=manifest.selection.protein.auth_asym_id,
        ligand_instance=manifest.selection.ligand.instance,
        model_number=manifest.selection.model_number,
        coordinate_units="angstrom",
    )

    two_d_svg = args.figure_directory / "2HYY-imatinib-2d.svg"
    two_d_png = args.figure_directory / "2HYY-imatinib-2d.png"
    two_d = draw_ligand_2d(
        ligand_identity=identity,
        chemical_state=reference_state,
        context=context,
        output_svg=two_d_svg,
        output_png=two_d_png,
        figure_id="capstone-ligand-2d",
        title="Imatinib reference-neutral chemical state",
        caption=(
            "Reference-neutral imatinib graph. Blue atoms are hydrogen-bond "
            "donors, orange atoms are acceptors, and purple atoms are both. "
            "Atom labels are reference indices; hydrogens are omitted."
        ),
        alt_text=(
            "Two-dimensional chemical structure of imatinib with donor and "
            "acceptor atoms highlighted and reference atom indices shown."
        ),
        source_records=records("identity", "states"),
        source_checksums=checksums("identity", "states"),
    )

    overview_path = args.figure_directory / "2HYY-imatinib-abl1-overview.png"
    overview = render_protein_ligand_overview(
        structure_path=args.structure,
        manifest=manifest,
        context=context,
        output_path=overview_path,
        figure_id="capstone-abl1-overview",
        title="Imatinib in the ABL1 binding site",
        caption=(
            "Imatinib (yellow sticks) in the ABL1 binding site, chain A. The "
            "pale surface marks protein residues within 4 Å of the ligand. "
            "This is one crystallographic observation, not an ensemble or an "
            "affinity estimate."
        ),
        alt_text=(
            "Protein cartoon of ABL1 with a yellow imatinib molecule in the "
            "binding pocket and a pale transparent surface around nearby residues."
        ),
        source_records=records("manifest", "site"),
        source_checksums=checksums("manifest", "site"),
        camera_view=views.get("capstone-abl1-overview"),
    )

    selected_contact_keys = {
        ("N13", "315", "OG1"),
        ("N21", "286", "OE2"),
        ("N21", "381", "O"),
    }
    pocket_contacts = [
        contact
        for contact in binding_site.contacts
        if (contact.ligand_atom, contact.auth_seq_id, contact.protein_atom)
        in selected_contact_keys
    ]
    pocket_path = args.figure_directory / "2HYY-imatinib-binding-pocket.png"
    pocket = render_binding_pocket_closeup(
        structure_path=args.structure,
        manifest=manifest,
        binding_site=binding_site,
        context=context,
        output_path=pocket_path,
        figure_id="capstone-binding-pocket",
        title="Imatinib binding-pocket contacts",
        caption=(
            "Binding-pocket close-up. Cyan residues are contact residues, "
            "yellow sticks are the ligand, and blue spheres are nearby waters. "
            "Gray dashed lines are measured donor-acceptor distances for three "
            "polar contacts; they are not classified hydrogen bonds. No "
            "energetic interpretation is implied."
        ),
        alt_text=(
            "Close-up of the imatinib binding pocket with cyan contact residues, "
            "a yellow ligand, blue waters, and gray dashed distance measurements "
            "for three polar contacts."
        ),
        source_records=records("manifest", "site", "hydrogen"),
        source_checksums=checksums("manifest", "site", "hydrogen"),
        contacts=pocket_contacts,
        max_contacts=3,
        camera_view=views.get("capstone-binding-pocket"),
    )

    comparison_ids = ("reference-neutral", "tautomer-3", "protonation-1")
    comparison_states = tuple(
        next(state for state in states.states if state.state_id == state_id)
        for state_id in comparison_ids
    )
    comparison_svg = args.figure_directory / "2HYY-imatinib-chemical-states.svg"
    comparison_png = args.figure_directory / "2HYY-imatinib-chemical-states.png"
    comparison = draw_chemical_state_comparison(
        ligand_identity=identity,
        states=comparison_states,
        context=context,
        output_svg=comparison_svg,
        output_png=comparison_png,
        figure_id="capstone-chemical-states",
        title="Selected imatinib chemical-state alternatives",
        caption=(
            "Reference-neutral, tautomer-3, and protonation-1 states shown side "
            "by side. Magenta marks formal-charge changes, purple marks "
            "donor/acceptor changes, and green marks per-atom attached-hydrogen "
            "count changes. "
            "No state is selected as dominant and no populations are shown."
        ),
        alt_text=(
            "Three side-by-side imatinib diagrams comparing a reference-neutral "
            "state, a tautomer, and a protonated state with changed atoms highlighted."
        ),
        source_records=records("identity", "states"),
        source_checksums=checksums("identity", "states"),
    )

    matrix_png = args.figure_directory / "2HYY-interaction-evidence-matrix.png"
    matrix_svg = args.figure_directory / "2HYY-interaction-evidence-matrix.svg"
    evidence_matrix = plot_evidence_matrix(
        combination_record=combination,
        context=context,
        output_png=matrix_png,
        output_svg=matrix_svg,
        figure_id="capstone-evidence-matrix",
        title="8 x 2 interaction-evidence matrix",
        caption=(
            "Six interaction categories for eight ligand states and two protein "
            "preparation states. Present, absent, unknown, and not-applicable "
            "statuses remain distinct; unknown is not zero."
        ),
        alt_text=(
            "Six-panel grid showing present, absent, unknown, and not-applicable "
            "interaction evidence for eight ligand and two protein states."
        ),
        source_records=records("combination"),
        source_checksums=checksums("combination"),
    )

    angle_text = "; ".join(
        f"{bond.donor_atom_id}-H...{bond.acceptor_atom_id} "
        f"{bond.angle_degrees:.1f} deg"
        for bond in hydrogen.bonds
    )
    supported_bonds = tuple(
        bond for bond in hydrogen.bonds if bond.status == "supported"
    )
    hydrogen_water_path = (
        args.figure_directory / "2HYY-hydrogen-water-evidence.png"
    )
    hydrogen_water = render_hydrogen_water_evidence(
        structure_path=args.structure,
        manifest=manifest,
        hydrogen=hydrogen,
        water=water,
        context=context,
        output_path=hydrogen_water_path,
        figure_id="capstone-hydrogen-water",
        title="Generated hydrogens and deposited water evidence",
        caption=(
            "Dark gray spheres are generated ligand hydrogens, cyan sticks are "
            "acceptor residues, and blue spheres are deposited waters with "
            "unassigned orientation and roles. Only the two supported "
            "candidates are drawn; the third recorded contact fails the default "
            "angle criterion. D-H...A angles from the record: "
            f"{angle_text}. Distances and angles are geometric evidence, not "
            "classified interactions."
        ),
        alt_text=(
            "Close-up with yellow ligand sticks, dark gray generated hydrogens, "
            "cyan acceptor residues, blue waters, and labeled donor-acceptor distances."
        ),
        source_records=records("manifest", "site", "hydrogen", "water"),
        source_checksums=checksums("manifest", "site", "hydrogen", "water"),
        bonds=supported_bonds,
        camera_view=views.get("capstone-hydrogen-water"),
    )

    sensitivity_png = args.figure_directory / "2HYY-interaction-sensitivity.png"
    sensitivity_svg = args.figure_directory / "2HYY-interaction-sensitivity.svg"
    sensitivity = plot_interaction_sensitivity(
        robustness_record=robustness,
        context=context,
        output_png=sensitivity_png,
        output_svg=sensitivity_svg,
        figure_id="capstone-interaction-sensitivity",
        title="Interaction sensitivity across tested assumptions",
        caption=(
            "Separate sensitivity axes for contact cutoff, hydrogen-bond angle, "
            "water cutoff, and hydrogen-placement perturbation. Conclusion "
            "colors distinguish invariant, assumption-sensitive, unsupported, "
            "and unknown statuses. No composite score is produced."
        ),
        alt_text=(
            "Four sensitivity panels plus a conclusion-status legend showing "
            "contact cutoff, hydrogen-bond angle, water cutoff, and hydrogen "
            "placement effects."
        ),
        source_records=records("robustness"),
        source_checksums=checksums("robustness"),
    )

    fe_png = args.figure_directory / "harmonic-fep-convergence.png"
    fe_svg = args.figure_directory / "harmonic-fep-convergence.svg"
    free_energy = plot_free_energy_convergence(
        study=free_energy_study,
        context=context,
        output_png=fe_png,
        output_svg=fe_svg,
        figure_id="harmonic-fep-convergence",
        title="Exact versus estimated harmonic free energy",
        caption=(
            "Controlled three-dimensional harmonic model. Dashed lines are the "
            "exact free-energy differences; points are forward-FEP estimates "
            "with standard errors as sample count grows. The right panel shows "
            "effective sample fraction for the good- and poor-overlap cases."
        ),
        alt_text=(
            "Two panels showing exact and FEP harmonic free-energy convergence "
            "with error bars and effective sample fractions."
        ),
        source_records=records("free_energy_study"),
        source_checksums=checksums("free_energy_study"),
    )

    thermo_png = (
        args.figure_directory / "2HYY-structural-vs-thermodynamics.png"
    )
    thermo_svg = (
        args.figure_directory / "2HYY-structural-vs-thermodynamics.svg"
    )
    thermo = plot_structural_thermodynamics_boundary(
        record=thermodynamics,
        context=context,
        output_png=thermo_png,
        output_svg=thermo_svg,
        figure_id="capstone-structural-vs-thermodynamics",
        title="Structural evidence is not binding free energy",
        caption=(
            "Structural counts are shown separately from the missing "
            "thermodynamic terms and preparation requirements. Every term is "
            "unestimated and the binding free energy is not computed."
        ),
        alt_text=(
            "Bar chart of structural evidence next to a table of unestimated "
            "thermodynamic terms and a large not-computed free-energy status."
        ),
        source_records=records("thermodynamics"),
        source_checksums=checksums("thermodynamics"),
    )

    figures = (
        two_d,
        overview,
        pocket,
        comparison,
        evidence_matrix,
        hydrogen_water,
        sensitivity,
        free_energy,
        thermo,
    )
    errors = validate_figure_metadata(figures)
    if errors:
        for error in errors:
            print(f"visualization metadata error: {error}")
        raise SystemExit(2)
    metadata_path = write_figure_metadata(args.metadata, figures)

    print("Visualization Batches 1, 2, and 3")
    for figure in figures:
        print(
            f"  {figure.figure_id}: {figure.image_path} "
            f"({figure.visualization_tool} {figure.tool_version})"
        )
    print(f"  metadata: {metadata_path}")
    print("  scope: structure and chemistry only; no interaction or affinity inference")


if __name__ == "__main__":
    main()
