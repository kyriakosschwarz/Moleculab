"""PyMOL renders of the selected protein-ligand structure."""

from __future__ import annotations

import importlib.util
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path
from typing import Any

from ..applications.binding_site_reconciliation import BindingSiteRecord
from ..applications.structure_manifest import StructureManifest
from .metadata import FigureContext, FigureMetadata

_PYMOL_PATH = Path("/usr/lib/python3/dist-packages/pymol")
_PYMOL_DATA = Path("/usr/share/pymol/data")
_PYMOL_SCRIPTS = Path("/usr/share/pymol/scripts")
_SYSTEM_PYTHON = Path("/usr/bin/python3")


def _launch() -> tuple[list[str], dict[str, str]]:
    environment = os.environ.copy()
    if importlib.util.find_spec("pymol") is not None:
        return [sys.executable, "-m", "pymol.__init__"], environment
    if _SYSTEM_PYTHON.exists() and _PYMOL_PATH.exists():
        environment.setdefault("PYMOL_PATH", str(_PYMOL_PATH))
        environment.setdefault("PYMOL_DATA", str(_PYMOL_DATA))
        environment.setdefault("PYMOL_SCRIPTS", str(_PYMOL_SCRIPTS))
        environment.setdefault("CHEMPY_DATA", str(_PYMOL_DATA / "chempy"))
        return [str(_SYSTEM_PYTHON), "-m", "pymol.__init__"], environment
    executable = shutil.which("pymol")
    if executable is not None:
        return [executable], environment
    raise RuntimeError(
        "PyMOL is required for protein-ligand figures. Install PyMOL or set "
        "the executable on PATH."
    )


def pymol_available() -> bool:
    """Return whether a usable PyMOL launch command can be found."""

    try:
        _launch()
    except RuntimeError:
        return False
    return True


def _pymol_selection(
    manifest: StructureManifest,
) -> tuple[str, str, str, str]:
    protein_chain = manifest.selection.protein.auth_asym_id
    ligand = manifest.selection.ligand
    component = ligand.component_id
    ligand_chain = ligand.auth_asym_id
    ligand_selection = f"chain {ligand_chain} and resn {component}"
    protein_selection = f"chain {protein_chain} and polymer"
    return protein_chain, component, ligand_selection, protein_selection


def _view_block(camera_view: Sequence[float] | None) -> str:
    if camera_view is None:
        return (
            "orient ligand or pocket\n"
            "zoom ligand or pocket, 3\n"
        )
    values = ", ".join(f"{float(value):.8f}" for value in camera_view)
    return f"set_view ({values})\n"


def _pocket_view_block(camera_view: Sequence[float] | None) -> str:
    if camera_view is None:
        return "orient ligand\nzoom ligand, 2.0\n"
    values = ", ".join(f"{float(value):.8f}" for value in camera_view)
    return f"set_view ({values})\n"


def _render_block(
    output_path: str | Path,
    *,
    width: int,
    height: int,
    dpi: int,
) -> str:
    image_path = repr(str(Path(output_path).resolve()))
    return (
        f"ray {int(width)}, {int(height)}\n"
        "python\n"
        "from pymol import cmd\n"
        "print(\"PYMOL_VERSION\", cmd.get_version()[0])\n"
        "print(\"VIEW\", list(cmd.get_view()))\n"
        f"cmd.png({image_path}, dpi={int(dpi)})\n"
        "python end\n"
    )


def build_overview_script(
    *,
    structure_path: str | Path,
    manifest: StructureManifest,
    output_path: str | Path,
    width: int = 1600,
    height: int = 1200,
    dpi: int = 200,
    camera_view: Sequence[float] | None = None,
) -> str:
    """Build a deterministic PyMOL script for the protein-ligand overview."""

    protein_chain, _, ligand, protein = _pymol_selection(manifest)
    waters = f"chain {protein_chain} and resn HOH within 4.0 of ligand"
    return f"""load "{Path(structure_path).resolve()}", 2hyy
hide everything
select protein, {protein}
select ligand, {ligand}
select pocket, byres (protein within 4.0 of ligand)
show cartoon, protein
color gray80, protein
show sticks, ligand
color yellow, ligand
show surface, pocket
set transparency, 0.85
set two_sided_lighting, 1
set surface_quality, 2
color palecyan, pocket
show spheres, {waters}
color marine, {waters}
set sphere_scale, 0.30
bg_color white
set antialias, 2
set ray_shadows, 0
set ambient, 0.45
set direct, 0.55
{_view_block(camera_view)}python
from pymol import cmd
python end
{_render_block(output_path, width=width, height=height, dpi=dpi)}"""


def build_pocket_script(
    *,
    structure_path: str | Path,
    manifest: StructureManifest,
    binding_site: BindingSiteRecord,
    output_path: str | Path,
    contacts: Sequence[Any] | None = None,
    max_contacts: int = 3,
    width: int = 1600,
    height: int = 1200,
    dpi: int = 200,
    camera_view: Sequence[float] | None = None,
) -> str:
    """Build a deterministic PyMOL script for the binding-pocket close-up."""

    protein_chain, _, ligand, protein = _pymol_selection(manifest)
    chosen = list(contacts) if contacts is not None else list(binding_site.contacts)
    chosen = sorted(chosen, key=lambda item: item.distance_angstrom)[:max_contacts]
    if not chosen:
        raise ValueError("at least one contact is required for the close-up.")

    residue_ids = []
    for contact in chosen:
        residue = (
            f"{contact.auth_seq_id}{contact.insertion_code or ''}"
            if contact.insertion_code
            else contact.auth_seq_id
        )
        if residue not in residue_ids:
            residue_ids.append(residue)
    contact_residues = (
        f"chain {protein_chain} and polymer and resi " + "+".join(residue_ids)
    )
    waters = f"chain {protein_chain} and resn HOH within 4.0 of ligand"

    distances = "\n".join(
        (
            f"distance contact_{index}, "
            f"(ligand and name {contact.ligand_atom}), "
            f"(chain {contact.auth_asym_id} and resi "
            f"{contact.auth_seq_id}{contact.insertion_code or ''} and name "
            f"{contact.protein_atom})"
        )
        for index, contact in enumerate(chosen, start=1)
    )
    return f"""load "{Path(structure_path).resolve()}", 2hyy
hide everything
select protein, {protein}
select ligand, {ligand}
select contact_residues, {contact_residues}
show cartoon, protein
set cartoon_transparency, 0.65
color gray80, protein
show sticks, contact_residues
color cyan, contact_residues
show sticks, ligand
color yellow, ligand
show spheres, {waters}
color marine, {waters}
set sphere_scale, 0.28
{distances}
set dash_color, gray30
set dash_gap, 0.22
set dash_length, 0.16
set dash_width, 2
set label_size, 28
set label_color, black
set label_outline_color, white
bg_color white
set antialias, 2
set ray_shadows, 0
set ambient, 0.45
set direct, 0.55
{_pocket_view_block(camera_view)}python
from pymol import cmd
python end
{_render_block(output_path, width=width, height=height, dpi=dpi)}"""


def _run_script(script: str, output_path: str | Path) -> tuple[str, tuple[float, ...]]:
    command, environment = _launch()
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as directory:
        script_path = Path(directory) / "render.pml"
        script_path.write_text(script, encoding="utf-8")
        result = subprocess.run(
            [*command, "-cq", str(script_path)],
            capture_output=True,
            text=True,
            env=environment,
            timeout=300,
            check=False,
        )
    if result.returncode != 0:
        raise RuntimeError(
            "PyMOL rendering failed.\n"
            f"stdout:\n{result.stdout[-3000:]}\n"
            f"stderr:\n{result.stderr[-3000:]}"
        )
    if not destination.exists():
        raise RuntimeError(f"PyMOL did not write {destination}.")
    version_match = re.search(r"PYMOL_VERSION\s+(\S+)", result.stdout)
    view_match = re.search(r"VIEW\s+\[([^\]]+)\]", result.stdout)
    if version_match is None or view_match is None:
        raise RuntimeError(
            "PyMOL output did not contain the expected version or view matrix."
        )
    view = tuple(float(value) for value in view_match.group(1).split(","))
    return version_match.group(1), view


def render_protein_ligand_overview(
    *,
    structure_path: str | Path,
    manifest: StructureManifest,
    context: FigureContext,
    output_path: str | Path,
    figure_id: str,
    title: str,
    caption: str,
    alt_text: str,
    width: int = 1600,
    height: int = 1200,
    dpi: int = 200,
    camera_view: Sequence[float] | None = None,
) -> FigureMetadata:
    """Render a protein cartoon with the ligand and a restrained pocket surface."""

    script = build_overview_script(
        structure_path=structure_path,
        manifest=manifest,
        output_path=output_path,
        width=width,
        height=height,
        dpi=dpi,
        camera_view=camera_view,
    )
    version, view = _run_script(script, output_path)
    return FigureMetadata(
        figure_id=figure_id,
        title=title,
        image_path=str(output_path),
        context=context,
        visualization_tool="PyMOL",
        tool_version=version,
        scientific_limitations=(
            "One crystallographic structure and one selected ligand copy.",
            "Cartoon and surface are representations of the deposited coordinates, not new data.",
            "The ligand orientation is from the crystal structure, not a predicted pose.",
            "No binding free energy or interaction favourability is shown.",
        ),
        alt_text=alt_text,
        caption=caption,
        camera_view=view,
    )


def render_binding_pocket_closeup(
    *,
    structure_path: str | Path,
    manifest: StructureManifest,
    binding_site: BindingSiteRecord,
    context: FigureContext,
    output_path: str | Path,
    figure_id: str,
    title: str,
    caption: str,
    alt_text: str,
    contacts: Sequence[Any] | None = None,
    max_contacts: int = 3,
    width: int = 1600,
    height: int = 1200,
    dpi: int = 200,
    camera_view: Sequence[float] | None = None,
) -> FigureMetadata:
    """Render a close-up with neutral measured-contact distance labels."""

    script = build_pocket_script(
        structure_path=structure_path,
        manifest=manifest,
        binding_site=binding_site,
        output_path=output_path,
        contacts=contacts,
        max_contacts=max_contacts,
        width=width,
        height=height,
        dpi=dpi,
        camera_view=camera_view,
    )
    version, view = _run_script(script, output_path)
    return FigureMetadata(
        figure_id=figure_id,
        title=title,
        image_path=str(output_path),
        context=context,
        visualization_tool="PyMOL",
        tool_version=version,
        scientific_limitations=(
            "Dashed lines are measured atom-atom contacts, not classified hydrogen bonds.",
            "The figure does not show hydrogen positions or angular criteria.",
            "Nearby waters are geometric observations; water roles remain unassigned.",
            "The figure does not imply energetic favourability or binding affinity.",
        ),
        alt_text=alt_text,
        caption=caption,
        camera_view=view,
    )
