"""Tests for the reproducible visualization metadata and render scripts."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from moleculab import (
    FIGURE_METADATA_SCHEMA_VERSION,
    FigureContext,
    FigureMetadata,
    build_overview_script,
    build_pocket_script,
    draw_ligand_2d,
    load_binding_site,
    load_chemical_state_enumeration,
    load_figure_metadata,
    load_ligand_identity,
    load_manifest,
    validate_figure_metadata,
    write_figure_metadata,
)

ROOT = Path(__file__).resolve().parents[1]


def context() -> FigureContext:
    manifest = load_manifest(ROOT / "data/structures/2HYY.manifest.json")
    return FigureContext(
        source_structure="2HYY",
        structure_sha256=manifest.source.sha256,
        protein_chain=manifest.selection.protein.auth_asym_id,
        ligand_instance=manifest.selection.ligand.instance,
        model_number=manifest.selection.model_number,
        generated_from=("data/structures/2HYY.manifest.json",),
    )


def metadata(tmp_path: Path, figure_id: str = "test") -> FigureMetadata:
    return FigureMetadata(
        figure_id=figure_id,
        title="Test figure",
        image_path=str(tmp_path / "figure.svg"),
        context=context(),
        visualization_tool="RDKit",
        tool_version="test",
        scientific_limitations=("test limitation",),
        alt_text="test alt text",
        caption="test caption",
    )


def test_metadata_round_trip_and_validation(tmp_path):
    figure = metadata(tmp_path)
    destination = write_figure_metadata(tmp_path / "metadata.json", (figure,))
    loaded = load_figure_metadata(destination)

    assert loaded == (figure,)
    assert FIGURE_METADATA_SCHEMA_VERSION == "1.0"
    assert validate_figure_metadata(loaded) == []


def test_duplicate_figure_ids_are_rejected(tmp_path):
    figure = metadata(tmp_path)
    errors = validate_figure_metadata((figure, figure))

    assert any("duplicated" in error for error in errors)


def test_pymol_scripts_use_record_selections_and_neutral_distances():
    manifest = load_manifest(ROOT / "data/structures/2HYY.manifest.json")
    binding_site = load_binding_site(ROOT / "data/structures/2HYY.binding-site.json")
    overview = build_overview_script(
        structure_path=ROOT / "data/structures/2HYY.cif",
        manifest=manifest,
        output_path="/tmp/overview.png",
        camera_view=(1.0,) * 18,
    )
    pocket = build_pocket_script(
        structure_path=ROOT / "data/structures/2HYY.cif",
        manifest=manifest,
        binding_site=binding_site,
        output_path="/tmp/pocket.png",
        max_contacts=3,
        camera_view=(1.0,) * 18,
    )

    assert "chain A and resn STI" in overview
    assert "show surface, pocket" in overview
    assert "set transparency, 0.85" in overview
    assert "set_view" in overview
    assert "distance contact_1" in pocket
    assert "chain A and resi" in pocket
    assert "dash_color, gray30" in pocket
    assert "not classified hydrogen bonds" not in pocket


def test_rdkit_ligand_drawing_and_accessibility(tmp_path):
    pytest.importorskip("rdkit")
    identity = load_ligand_identity(ROOT / "data/structures/2HYY.ligand-identity.json")
    states = load_chemical_state_enumeration(
        ROOT / "data/structures/2HYY.chemical-states.json"
    )
    state = next(item for item in states.states if item.state_id == "reference-neutral")
    svg_path = tmp_path / "imatinib.svg"
    png_path = tmp_path / "imatinib.png"
    figure = draw_ligand_2d(
        ligand_identity=identity,
        chemical_state=state,
        context=context(),
        output_svg=svg_path,
        output_png=png_path,
        figure_id="test-2d",
        title="Imatinib reference-neutral state",
        caption="caption",
        alt_text="alt text",
    )

    assert svg_path.exists()
    assert png_path.exists()
    assert png_path.stat().st_size > 0
    svg = svg_path.read_text(encoding="utf-8")
    assert "<title>Imatinib reference-neutral state</title>" in svg
    assert "caption Alt text: alt text" in svg
    assert figure.visualization_tool == "RDKit"
    assert figure.tool_version
    assert figure.scientific_limitations
    assert any("reference-graph indices" in item for item in figure.scientific_limitations)
