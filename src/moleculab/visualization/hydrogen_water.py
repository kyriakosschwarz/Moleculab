"""Hydrogen and water evidence figures from validated evidence records."""

from __future__ import annotations

import re
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from ..applications.hydrogen_evidence import HydrogenEvidenceRecord
from ..applications.structure_manifest import StructureManifest
from ..applications.water_evidence import WaterEvidenceRecord
from .metadata import FigureContext, FigureMetadata
from .protein_ligand import _pymol_selection, _run_script


def _hydrogen_atom_name(parent_atom_id: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]", "", parent_atom_id).upper()
    return ("H" + cleaned)[:4]


def _hydrogen_pdb(hydrogen: HydrogenEvidenceRecord) -> str:
    lines = []
    for index, item in enumerate(hydrogen.ligand_hydrogens, start=1):
        x, y, z = item.position_angstrom
        name = _hydrogen_atom_name(item.parent_atom_id)
        lines.append(
            f"HETATM{index:5d} {name:<4s} HYD Z   1    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          H  "
        )
    lines.append("END")
    return "\n".join(lines) + "\n"


def _acceptor_selection(bonds: Sequence[Any]) -> str:
    selections = []
    for bond in bonds:
        chain, residue, atom = bond.acceptor_atom_id.split(":", 2)
        selections.append(f"(chain {chain} and resi {residue} and name {atom})")
    return " or ".join(selections)


def _water_selection(waters: WaterEvidenceRecord) -> str:
    if not waters.waters:
        raise ValueError("at least one deposited water is required.")
    residues = []
    for item in waters.waters:
        if item.auth_seq_id not in residues:
            residues.append(item.auth_seq_id)
    return (
        f"chain {waters.waters[0].auth_asym_id} and resn HOH and resi "
        + "+".join(residues)
    )


def build_hydrogen_water_script(
    *,
    structure_path: str | Path,
    manifest: StructureManifest,
    hydrogen: HydrogenEvidenceRecord,
    water: WaterEvidenceRecord,
    output_path: str | Path,
    bonds: Sequence[Any] | None = None,
    width: int = 1600,
    height: int = 1200,
    dpi: int = 200,
    camera_view: Sequence[float] | None = None,
) -> tuple[str, str]:
    """Return the PyMOL script and the generated-hydrogen PDB text."""

    _, _, ligand, protein = _pymol_selection(manifest)
    selected_bonds = tuple(hydrogen.bonds) if bonds is None else tuple(bonds)
    acceptors = _acceptor_selection(selected_bonds)
    water_selector = _water_selection(water)
    pdb_text = _hydrogen_pdb(hydrogen)
    distances = []
    for index, bond in enumerate(selected_bonds, start=1):
        donor = bond.donor_atom_id
        acceptor_chain, acceptor_residue, acceptor_atom = (
            bond.acceptor_atom_id.split(":", 2)
        )
        acceptor_selection = (
            f"(chain {acceptor_chain} and resi {acceptor_residue} and "
            f"name {acceptor_atom})"
        )
        distances.append(
            f"distance hb_{index}, (ligand and name {donor}), "
            f"{acceptor_selection}"
        )
    view = (
        "orient ligand\nzoom ligand, 2.0\n"
        if camera_view is None
        else "set_view ("
        + ", ".join(f"{float(value):.8f}" for value in camera_view)
        + ")\n"
    )
    image_path = repr(str(Path(output_path).resolve()))
    script = f"""load "{Path(structure_path).resolve()}", 2hyy
load "{Path(output_path).with_suffix('.h.pdb').resolve()}", generated_h
hide everything
select protein, {protein}
select ligand, {ligand}
select acceptors, {acceptors}
select waters, {water_selector}
show cartoon, protein
set cartoon_transparency, 0.70
color gray80, protein
show sticks, ligand
color yellow, ligand
show spheres, generated_h
color gray20, generated_h
set sphere_scale, 0.45
show sticks, acceptors
color cyan, acceptors
show spheres, waters
color marine, waters
set sphere_scale, 0.30
{chr(10).join(distances)}
set dash_color, gray30
set dash_gap, 0.22
set dash_length, 0.16
set dash_width, 2
set label_size, 24
set label_color, black
set label_outline_color, white
bg_color white
set antialias, 2
set ray_shadows, 0
set ambient, 0.45
set direct, 0.55
{view}ray {int(width)}, {int(height)}
python
from pymol import cmd
print("PYMOL_VERSION", cmd.get_version()[0])
print("VIEW", list(cmd.get_view()))
cmd.png({image_path}, dpi={int(dpi)})
python end
"""
    return script, pdb_text


def render_hydrogen_water_evidence(
    *,
    structure_path: str | Path,
    manifest: StructureManifest,
    hydrogen: HydrogenEvidenceRecord,
    water: WaterEvidenceRecord,
    context: FigureContext,
    output_path: str | Path,
    figure_id: str,
    bonds: Sequence[Any] | None = None,
    title: str,
    caption: str,
    alt_text: str,
    source_records: tuple[str, ...],
    source_checksums: tuple[tuple[str, str], ...],
    width: int = 1600,
    height: int = 1200,
    dpi: int = 200,
    camera_view: Sequence[float] | None = None,
) -> FigureMetadata:
    """Render generated hydrogens, acceptor geometry, and deposited waters."""

    script, pdb_text = build_hydrogen_water_script(
        structure_path=structure_path,
        manifest=manifest,
        hydrogen=hydrogen,
        water=water,
        output_path=output_path,
        bonds=bonds,
        width=width,
        height=height,
        dpi=dpi,
        camera_view=camera_view,
    )
    output_path = Path(output_path)
    pdb_path = output_path.with_suffix(".h.pdb")
    pdb_path.write_text(pdb_text, encoding="utf-8")
    try:
        version, view = _run_script(script, output_path)
    finally:
        pdb_path.unlink(missing_ok=True)
    return FigureMetadata(
        figure_id=figure_id,
        title=title,
        image_path=str(output_path),
        context=context,
        visualization_tool="PyMOL",
        tool_version=version,
        source_records=source_records,
        source_checksums=source_checksums,
        rendering_config=(
            ("backend", "PyMOL"),
            ("width", str(width)),
            ("height", str(height)),
            ("dpi", str(dpi)),
            ("ligand", "sticks yellow"),
            ("generated_hydrogens", "spheres gray20 scale 0.45"),
            ("acceptors", "sticks cyan"),
            ("waters", "spheres marine"),
            ("distance_style", "gray30 dashed width 2"),
            ("angles", "reported in caption and source record, not drawn"),
            ("shown_bonds", f"{len(bonds) if bonds is not None else len(hydrogen.bonds)}"),
            ("camera", "stored" if camera_view is not None else "auto"),
        ),
        rendering_config_sha256=None,
        scientific_limitations=(
            "Ligand hydrogens are generated positions, not experimental observations.",
            "Water orientation and donor/acceptor roles remain unassigned.",
            "Distances and angles are geometric evidence, not energies or affinities.",
            "The figure does not classify hydrogen bonds or water-mediated interactions.",
        ),
        alt_text=alt_text,
        caption=caption,
        camera_view=view,
    )
