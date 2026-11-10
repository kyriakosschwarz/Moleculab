# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Tests for the reproducible visualization metadata and render scripts."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from moleculab import (
    FIGURE_METADATA_SCHEMA_VERSION,
    FigureContext,
    FigureMetadata,
    build_hydrogen_water_script,
    build_overview_script,
    build_pocket_script,
    draw_chemical_state_comparison,
    draw_ligand_2d,
    load_binding_site,
    load_chemical_state_enumeration,
    load_combination_record,
    load_figure_metadata,
    load_hydrogen_evidence,
    load_ligand_identity,
    load_manifest,
    load_water_evidence,
    plot_evidence_matrix,
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
    )


def metadata(tmp_path: Path, figure_id: str = "test") -> FigureMetadata:
    return FigureMetadata(
        figure_id=figure_id,
        title="Test figure",
        image_path=str(tmp_path / "figure.svg"),
        context=context(),
        visualization_tool="RDKit",
        tool_version="test",
        source_records=("data/structures/2HYY.manifest.json",),
        source_checksums=(("data/structures/2HYY.manifest.json", "a" * 64),),
        rendering_config=(("backend", "test"), ("width", "100")),
        rendering_config_sha256=None,
        scientific_limitations=("test limitation",),
        alt_text="test alt text",
        caption="test caption",
    )


def test_metadata_round_trip_and_validation(tmp_path):
    figure = metadata(tmp_path)
    destination = write_figure_metadata(tmp_path / "metadata.json", (figure,))
    loaded = load_figure_metadata(destination)

    assert loaded == (figure,)
    assert FIGURE_METADATA_SCHEMA_VERSION == "1.1"
    assert validate_figure_metadata(loaded) == []


def test_source_checksums_must_cover_every_source_record(tmp_path):
    figure = metadata(tmp_path)
    bad = replace(
        figure,
        source_records=(
            "data/structures/2HYY.manifest.json",
            "data/structures/2HYY.binding-site.json",
        ),
    )
    errors = validate_figure_metadata((bad,))

    assert any("source_checksums" in error for error in errors)


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
        source_records=(
            "data/structures/2HYY.ligand-identity.json",
            "data/structures/2HYY.chemical-states.json",
        ),
        source_checksums=(
            ("data/structures/2HYY.ligand-identity.json", "a" * 64),
            ("data/structures/2HYY.chemical-states.json", "b" * 64),
        ),
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


def test_chemical_state_comparison_figure(tmp_path):
    pytest.importorskip("rdkit")
    identity = load_ligand_identity(ROOT / "data/structures/2HYY.ligand-identity.json")
    states = load_chemical_state_enumeration(
        ROOT / "data/structures/2HYY.chemical-states.json"
    )
    selected = tuple(
        next(item for item in states.states if item.state_id == state_id)
        for state_id in ("reference-neutral", "tautomer-3", "protonation-1")
    )
    svg_path = tmp_path / "states.svg"
    png_path = tmp_path / "states.png"
    figure = draw_chemical_state_comparison(
        ligand_identity=identity,
        states=selected,
        context=context(),
        output_svg=svg_path,
        output_png=png_path,
        figure_id="test-states",
        title="Selected chemical states",
        caption="caption",
        alt_text="alt text",
        source_records=("identity.json", "states.json"),
        source_checksums=(("identity.json", "a" * 64), ("states.json", "b" * 64)),
    )

    assert svg_path.exists() and svg_path.stat().st_size > 0
    assert png_path.exists() and png_path.stat().st_size > 0
    assert figure.visualization_tool == "RDKit"
    assert figure.source_records == ("identity.json", "states.json")
    assert figure.rendering_config_sha256


def test_evidence_matrix_figure(tmp_path):
    pytest.importorskip("matplotlib")
    combination = load_combination_record(
        ROOT / "data/structures/2HYY.combination-interactions.json"
    )
    png_path = tmp_path / "matrix.png"
    svg_path = tmp_path / "matrix.svg"
    figure = plot_evidence_matrix(
        combination_record=combination,
        context=context(),
        output_png=png_path,
        output_svg=svg_path,
        figure_id="test-matrix",
        title="Evidence matrix",
        caption="caption",
        alt_text="alt text",
        source_records=("combination.json",),
        source_checksums=(("combination.json", "c" * 64),),
        width=8.0,
        height=8.0,
    )

    assert png_path.exists() and png_path.stat().st_size > 0
    assert svg_path.exists() and svg_path.stat().st_size > 0
    assert figure.visualization_tool == "Matplotlib"
    assert figure.image_path == str(png_path)


def test_hydrogen_water_script_keeps_generated_positions_distinct(tmp_path):
    manifest = load_manifest(ROOT / "data/structures/2HYY.manifest.json")
    hydrogen = load_hydrogen_evidence(
        ROOT / "data/structures/2HYY.hydrogen-evidence.json"
    )
    water = load_water_evidence(ROOT / "data/structures/2HYY.water-evidence.json")
    supported = tuple(
        bond for bond in hydrogen.bonds if bond.status == "supported"
    )
    script, pdb = build_hydrogen_water_script(
        structure_path=ROOT / "data/structures/2HYY.cif",
        manifest=manifest,
        hydrogen=hydrogen,
        water=water,
        output_path=tmp_path / "hydrogen-water.png",
        bonds=supported,
    )

    assert "generated_h" in script
    assert "show spheres, generated_h" in script
    assert "color gray20, generated_h" in script
    assert "distance hb_1" in script
    assert "angle_" not in script
    assert "HN13" in pdb
    assert "HN21" in pdb
