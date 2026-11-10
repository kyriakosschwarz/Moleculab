# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""RDKit 2D chemical diagrams derived from validated ligand records."""

from __future__ import annotations

from pathlib import Path
from typing import Any
from xml.sax.saxutils import escape

from ..applications.chemical_state_enumeration import ChemicalState
from ..applications.ligand_identity import LigandIdentityRecord
from .metadata import FigureContext, FigureMetadata

_DONOR_COLOUR = (0.20, 0.40, 0.95, 0.38)
_ACCEPTOR_COLOUR = (0.95, 0.42, 0.18, 0.38)
_BOTH_COLOUR = (0.62, 0.24, 0.78, 0.40)


def _require_rdkit():
    try:
        from rdkit import Chem, __version__
        from rdkit.Chem.Draw import rdMolDraw2D
    except ImportError as error:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "RDKit is required for molecular_2d figures. "
            "Run with `uv run --with rdkit ...`."
        ) from error
    return Chem, rdMolDraw2D, __version__


def build_reference_molecule(ligand_identity: LigandIdentityRecord) -> Any:
    """Build an RDKit molecule from the validated reference graph."""

    Chem, _, _ = _require_rdkit()
    reference = ligand_identity.reference
    editable = Chem.RWMol()
    for atom_record in reference.atoms:
        atom = Chem.Atom(atom_record.element)
        atom.SetFormalCharge(atom_record.formal_charge)
        editable.AddAtom(atom)
    bond_types = {
        "SINGLE": Chem.BondType.SINGLE,
        "DOUBLE": Chem.BondType.DOUBLE,
        "TRIPLE": Chem.BondType.TRIPLE,
        "AROMATIC": Chem.BondType.AROMATIC,
    }
    for bond_record in reference.bonds:
        editable.AddBond(
            bond_record.atom_index_1,
            bond_record.atom_index_2,
            bond_types[bond_record.order],
        )
    molecule = editable.GetMol()
    Chem.SanitizeMol(molecule)
    return molecule


def _feature_indices(
    chemical_state: ChemicalState,
) -> tuple[
    list[int],
    dict[int, tuple[float, float, float, float]],
    list[int],
    list[int],
]:
    state_atoms = {atom.atom_id: atom for atom in chemical_state.atoms}
    donors = []
    acceptors = []
    highlight_atoms = []
    colours = {}
    for atom_id, reference_index in chemical_state.reference_atom_mapping:
        atom = state_atoms.get(atom_id)
        if atom is None:
            continue
        is_donor = bool(atom.is_donor)
        is_acceptor = bool(atom.is_acceptor)
        if not is_donor and not is_acceptor:
            continue
        highlight_atoms.append(reference_index)
        if is_donor and is_acceptor:
            colours[reference_index] = _BOTH_COLOUR
            donors.append(reference_index)
            acceptors.append(reference_index)
        elif is_donor:
            colours[reference_index] = _DONOR_COLOUR
            donors.append(reference_index)
        else:
            colours[reference_index] = _ACCEPTOR_COLOUR
            acceptors.append(reference_index)
    return highlight_atoms, colours, donors, acceptors


def _add_svg_accessibility(svg: str, title: str, description: str) -> str:
    insertion = svg.find(">") + 1
    if insertion <= 0:
        return svg
    metadata = f"<title>{escape(title)}</title><desc>{escape(description)}</desc>"
    return svg[:insertion] + metadata + svg[insertion:]


def draw_ligand_2d(
    *,
    ligand_identity: LigandIdentityRecord,
    chemical_state: ChemicalState,
    context: FigureContext,
    output_svg: str | Path,
    output_png: str | Path | None = None,
    figure_id: str,
    title: str,
    caption: str,
    alt_text: str,
    width: int = 1400,
    height: int = 1000,
) -> FigureMetadata:
    """Draw the selected chemical state as a labelled 2D molecular diagram.

    Donor atoms are blue, acceptors orange, and atoms that are both purple.
    Atom labels are reference-graph indices.  Generated hydrogens are not drawn.
    """

    Chem, rdMolDraw2D, __version__ = _require_rdkit()
    molecule = build_reference_molecule(ligand_identity)
    highlight_atoms, colours, donors, acceptors = _feature_indices(chemical_state)
    legend = (
        f"{title}\n"
        "blue: donor | orange: acceptor | purple: both\n"
        "atom labels: reference indices; hydrogens omitted"
    )

    output_svg_path = Path(output_svg)
    output_svg_path.parent.mkdir(parents=True, exist_ok=True)
    svg_drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    svg_drawer.drawOptions().addAtomIndices = True
    svg_drawer.drawOptions().clearBackground = True
    svg_drawer.drawOptions().backgroundColour = (1.0, 1.0, 1.0, 1.0)
    svg_drawer.DrawMolecule(
        molecule,
        legend=legend,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=colours,
    )
    svg_drawer.FinishDrawing()
    svg = _add_svg_accessibility(
        svg_drawer.GetDrawingText(),
        title,
        f"{caption} Alt text: {alt_text}",
    )
    output_svg_path.write_text(svg, encoding="utf-8")

    if output_png is not None:
        output_png_path = Path(output_png)
        output_png_path.parent.mkdir(parents=True, exist_ok=True)
        png_drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        png_drawer.drawOptions().addAtomIndices = True
        png_drawer.drawOptions().clearBackground = True
        png_drawer.drawOptions().backgroundColour = (1.0, 1.0, 1.0, 1.0)
        png_drawer.DrawMolecule(
            molecule,
            legend=legend,
            highlightAtoms=highlight_atoms,
            highlightAtomColors=colours,
        )
        png_drawer.FinishDrawing()
        output_png_path.write_bytes(png_drawer.GetDrawingText())

    return FigureMetadata(
        figure_id=figure_id,
        title=title,
        image_path=str(output_svg_path),
        context=context,
        visualization_tool="RDKit",
        tool_version=__version__,
        scientific_limitations=(
            "2D depiction of the selected chemical state, not a conformer.",
            "Generated hydrogen positions and state populations are not shown.",
            "Atom labels are reference-graph indices, not deposited atom names.",
            "The drawing is not an interaction, energy, or affinity diagram.",
            (
                f"Selected state donors/acceptors: {len(donors)}/{len(acceptors)} "
                "atoms highlighted."
            ),
        ),
        alt_text=alt_text,
        caption=caption,
        camera_view=None,
    )
