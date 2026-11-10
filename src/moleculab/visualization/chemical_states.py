"""Chemical-state comparison figures from validated state records."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

from ..applications.chemistry.chemical_state_enumeration import ChemicalState
from ..applications.structure.ligand_identity import LigandIdentityRecord
from .metadata import FigureContext, FigureMetadata
from .molecular_2d import _require_rdkit, build_reference_molecule
from .styles import clean_svg_file

_CHARGE_COLOUR = (0.95, 0.15, 0.75, 0.45)
_ROLE_COLOUR = (0.55, 0.25, 0.85, 0.45)
_HYDROGEN_COLOUR = (0.15, 0.70, 0.35, 0.45)


def _reference_query(ligand_identity: LigandIdentityRecord) -> Any:
    Chem, _, _ = _require_rdkit()
    reference = build_reference_molecule(ligand_identity)
    parameters = Chem.AdjustQueryParameters()
    parameters.makeBondsGeneric = True
    parameters.makeAtomsGeneric = False
    return Chem.AdjustQueryProperties(reference, parameters)


def _reference_features(
    ligand_identity: LigandIdentityRecord,
    reference_state: ChemicalState,
) -> dict[int, Any]:
    deposited_to_reference = {
        deposited_atom: reference_index
        for deposited_atom, reference_index in ligand_identity.mapping.pairs
    }
    atoms = {atom.atom_id: atom for atom in reference_state.atoms}
    return {
        deposited_to_reference[atom_id]: atom
        for atom_id, atom in atoms.items()
        if atom_id in deposited_to_reference
    }


def _state_atom_to_smiles_index(
    state: ChemicalState,
    molecule: Any,
) -> dict[str, int]:
    mapping = {}
    for atom in state.atoms:
        if atom.atom_id.startswith("atom_"):
            index = int(atom.atom_id.split("_", 1)[1])
        else:
            raise ValueError(
                f"cannot map deposited atom name {atom.atom_id!r} in "
                f"{state.state_id!r}; select a state with atom_N identifiers."
            )
        if index >= molecule.GetNumAtoms():
            raise ValueError(f"atom index {index} is outside {state.state_id!r}.")
        if molecule.GetAtomWithIdx(index).GetSymbol() != atom.element:
            raise ValueError(
                f"element mismatch for {state.state_id!r} atom {atom.atom_id!r}."
            )
        mapping[atom.atom_id] = index
    return mapping


def draw_chemical_state_comparison(
    *,
    ligand_identity: LigandIdentityRecord,
    states: tuple[ChemicalState, ...],
    context: FigureContext,
    output_svg: str | Path,
    output_png: str | Path | None,
    figure_id: str,
    title: str,
    caption: str,
    alt_text: str,
    source_records: tuple[str, ...],
    source_checksums: tuple[tuple[str, str], ...],
    width: int = 1800,
    height: int = 900,
) -> FigureMetadata:
    """Compare selected chemical states without selecting a population."""

    Chem, _, __version__ = _require_rdkit()
    from rdkit.Chem.Draw import MolsToGridImage

    reference_state = next(
        (state for state in states if state.state_id == "reference-neutral"),
        states[0],
    )
    reference_features = _reference_features(ligand_identity, reference_state)
    query = _reference_query(ligand_identity)
    reference_molecule = build_reference_molecule(ligand_identity)

    molecules = []
    legends = []
    highlight_lists = []
    colour_lists = []
    for state in states:
        if state.state_id == reference_state.state_id:
            molecules.append(reference_molecule)
            highlight_lists.append([])
            colour_lists.append({})
            legends.append(
                f"{state.label}\nq={state.formal_charge:+d}; "
                f"protons={state.proton_count}"
            )
            continue

        molecule = Chem.MolFromSmiles(state.molecular_graph)
        if molecule is None:
            raise ValueError(
                f"could not parse molecular_graph for {state.state_id!r}."
            )
        match = molecule.GetSubstructMatch(query)
        if len(match) != molecule.GetNumAtoms():
            raise ValueError(
                f"could not map {state.state_id!r} onto the reference graph."
            )
        smiles_to_reference = {
            smiles_index: reference_index
            for reference_index, smiles_index in enumerate(match)
        }
        atom_to_smiles = _state_atom_to_smiles_index(state, molecule)
        state_features = {atom.atom_id: atom for atom in state.atoms}
        highlights = []
        colours = {}
        for atom_id, smiles_index in atom_to_smiles.items():
            state_atom = state_features[atom_id]
            reference_atom = reference_features.get(
                smiles_to_reference[smiles_index]
            )
            if reference_atom is None:
                continue
            notes = []
            colour = None
            if state_atom.formal_charge != reference_atom.formal_charge:
                notes.append(f"q{state_atom.formal_charge:+d}")
                colour = _CHARGE_COLOUR
            if state_atom.is_donor != reference_atom.is_donor:
                notes.append("D" if state_atom.is_donor else "not D")
                colour = colour or _ROLE_COLOUR
            if state_atom.is_acceptor != reference_atom.is_acceptor:
                notes.append("A" if state_atom.is_acceptor else "not A")
                colour = colour or _ROLE_COLOUR
            if state_atom.hydrogen_count != reference_atom.hydrogen_count:
                delta = state_atom.hydrogen_count - reference_atom.hydrogen_count
                notes.append(f"H{delta:+d}")
                colour = colour or _HYDROGEN_COLOUR
            if notes:
                molecule.GetAtomWithIdx(smiles_index).SetProp(
                    "atomNote", ";".join(notes)
                )
                highlights.append(smiles_index)
                if colour is not None:
                    colours[smiles_index] = colour
        molecules.append(molecule)
        highlight_lists.append(highlights)
        colour_lists.append(colours)
        legends.append(
            f"{state.label}\nq={state.formal_charge:+d}; "
            f"protons={state.proton_count}"
        )

    sub_width = max(300, width // min(3, len(molecules)))
    sub_height = max(300, height)
    svg = MolsToGridImage(
        molecules,
        molsPerRow=3,
        subImgSize=(sub_width, sub_height),
        legends=legends,
        highlightAtomLists=highlight_lists,
        highlightAtomColors=colour_lists,
        useSVG=True,
    )
    svg_start = svg.find("<svg")
    insertion = svg.find(">", svg_start) + 1
    svg = (
        svg[:insertion]
        + f"<title>{escape(title)}</title><desc>{escape(caption)}</desc>"
        + svg[insertion:]
    )
    output_svg_path = Path(output_svg)
    output_svg_path.parent.mkdir(parents=True, exist_ok=True)
    output_svg_path.write_text(svg, encoding="utf-8")
    clean_svg_file(output_svg_path)

    if output_png is not None:
        image = MolsToGridImage(
            molecules,
            molsPerRow=3,
            subImgSize=(sub_width, sub_height),
            legends=legends,
            highlightAtomLists=highlight_lists,
            highlightAtomColors=colour_lists,
            useSVG=False,
        )
        output_png_path = Path(output_png)
        output_png_path.parent.mkdir(parents=True, exist_ok=True)
        image.save(output_png_path)

    return FigureMetadata(
        figure_id=figure_id,
        title=title,
        image_path=str(output_svg_path),
        context=context,
        visualization_tool="RDKit",
        tool_version=__version__,
        source_records=source_records,
        source_checksums=source_checksums,
        rendering_config=(
            ("backend", "RDKit"),
            ("width", str(width)),
            ("height", str(height)),
            ("mols_per_row", "3"),
            ("charge_colour", "magenta"),
            ("role_colour", "purple"),
            ("hydrogen_colour", "green"),
            (
                "atom_notes",
                "formal charge, donor/acceptor, per-atom attached hydrogen count changes",
            ),
        ),
        rendering_config_sha256=None,
        scientific_limitations=(
            "States are compared side by side; no state is selected as dominant.",
            "No populations or free energies are shown.",
            "Atom notes mark changes relative to the reference-neutral state.",
            "Hydrogen-count changes are attached hydrogens per atom, not total molecular hydrogen count.",
            "Generated 3D hydrogen positions are not drawn in this 2D comparison.",
            "The figure does not classify interactions or imply affinity.",
        ),
        alt_text=alt_text,
        caption=caption,
        camera_view=None,
    )
