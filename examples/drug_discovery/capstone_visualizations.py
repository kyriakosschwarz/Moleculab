"""Generate the three Visualization Batch 1 figures from the 2HYY records.

The figures are:
1. an RDKit 2D diagram of the reference-neutral imatinib state;
2. a PyMOL overview of imatinib in the ABL1 binding site;
3. a PyMOL binding-pocket close-up with measured contacts.

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
    draw_ligand_2d,
    load_binding_site,
    load_chemical_state_enumeration,
    load_figure_metadata,
    load_ligand_identity,
    load_manifest,
    pymol_available,
    render_binding_pocket_closeup,
    render_protein_ligand_overview,
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

    context = FigureContext(
        source_structure=manifest.source.accession,
        structure_sha256=structure_sha256,
        protein_chain=manifest.selection.protein.auth_asym_id,
        ligand_instance=manifest.selection.ligand.instance,
        model_number=manifest.selection.model_number,
        coordinate_units="angstrom",
        generated_from=(
            "data/structures/2HYY.manifest.json",
            "data/structures/2HYY.ligand-identity.json",
            "data/structures/2HYY.binding-site.json",
            "data/structures/2HYY.chemical-states.json",
        ),
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
        contacts=pocket_contacts,
        max_contacts=3,
        camera_view=views.get("capstone-binding-pocket"),
    )

    figures = (two_d, overview, pocket)
    errors = validate_figure_metadata(figures)
    if errors:
        for error in errors:
            print(f"visualization metadata error: {error}")
        raise SystemExit(2)
    metadata_path = write_figure_metadata(args.metadata, figures)

    print("Visualization Batch 1")
    for figure in figures:
        print(
            f"  {figure.figure_id}: {figure.image_path} "
            f"({figure.visualization_tool} {figure.tool_version})"
        )
    print(f"  metadata: {metadata_path}")
    print("  scope: structure and chemistry only; no interaction or affinity inference")


if __name__ == "__main__":
    main()
