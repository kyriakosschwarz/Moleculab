"""Authoritative registry and generated gallery for capstone figures."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from .metadata import FigureMetadata

GALLERY_RELATIVE_PATH = Path(
    "docs/tutorials/capstone/visualizations/gallery.md"
)


@dataclass(frozen=True)
class FigureSpec:
    figure_id: str
    figure_number: int
    title: str
    image_path: str
    generation_function: str
    source_records: tuple[str, ...]
    tutorial_anchor: str
    teaches: str
    record_backed: bool = True


FIGURE_REGISTRY = (
    FigureSpec(
        figure_id="capstone-ligand-2d",
        figure_number=1,
        title="Imatinib chemical identity",
        image_path="docs/tutorials/capstone/visualizations/2HYY-imatinib-2d.svg",
        generation_function="visualization.molecular_2d.draw_ligand_2d",
        source_records=(
            "data/structures/2HYY.ligand-identity.json",
            "data/structures/2HYY.chemical-states.json",
        ),
        tutorial_anchor=(
            "../README.md#420-visual-learning-layer-molecular-identity-and-structure"
        ),
        teaches="molecular graph, functional groups, and reference-neutral donor/acceptor features",
    ),
    FigureSpec(
        figure_id="capstone-abl1-overview",
        figure_number=2,
        title="Imatinib in the ABL1 binding site",
        image_path=(
            "docs/tutorials/capstone/visualizations/"
            "2HYY-imatinib-abl1-overview.png"
        ),
        generation_function=(
            "visualization.protein_ligand.render_protein_ligand_overview"
        ),
        source_records=(
            "data/structures/2HYY.manifest.json",
            "data/structures/2HYY.binding-site.json",
        ),
        tutorial_anchor=(
            "../README.md#420-visual-learning-layer-molecular-identity-and-structure"
        ),
        teaches="protein–ligand context, selected chain and ligand copy, and pocket location",
    ),
    FigureSpec(
        figure_id="capstone-binding-pocket",
        figure_number=3,
        title="Imatinib binding-pocket contacts",
        image_path=(
            "docs/tutorials/capstone/visualizations/"
            "2HYY-imatinib-binding-pocket.png"
        ),
        generation_function=(
            "visualization.protein_ligand.render_binding_pocket_closeup"
        ),
        source_records=(
            "data/structures/2HYY.manifest.json",
            "data/structures/2HYY.binding-site.json",
            "data/structures/2HYY.hydrogen-evidence.json",
        ),
        tutorial_anchor=(
            "../README.md#420-visual-learning-layer-molecular-identity-and-structure"
        ),
        teaches="measured contact distances and the difference between proximity and classification",
    ),
    FigureSpec(
        figure_id="capstone-chemical-states",
        figure_number=4,
        title="Selected imatinib chemical-state alternatives",
        image_path=(
            "docs/tutorials/capstone/visualizations/"
            "2HYY-imatinib-chemical-states.svg"
        ),
        generation_function=(
            "visualization.chemical_states.draw_chemical_state_comparison"
        ),
        source_records=(
            "data/structures/2HYY.ligand-identity.json",
            "data/structures/2HYY.chemical-states.json",
        ),
        tutorial_anchor=(
            "../README.md#421-visual-learning-layer-chemical-states-and-evidence"
        ),
        teaches="charge, donor/acceptor, and per-atom attached-hydrogen changes between states",
    ),
    FigureSpec(
        figure_id="capstone-evidence-matrix",
        figure_number=5,
        title="8 × 2 interaction-evidence matrix",
        image_path=(
            "docs/tutorials/capstone/visualizations/"
            "2HYY-interaction-evidence-matrix.png"
        ),
        generation_function="visualization.interaction_evidence.plot_evidence_matrix",
        source_records=("data/structures/2HYY.combination-interactions.json",),
        tutorial_anchor=(
            "../README.md#421-visual-learning-layer-chemical-states-and-evidence"
        ),
        teaches="present, absent, unknown, and not-applicable evidence across state combinations",
    ),
    FigureSpec(
        figure_id="capstone-hydrogen-water",
        figure_number=6,
        title="Generated hydrogens and deposited water evidence",
        image_path=(
            "docs/tutorials/capstone/visualizations/"
            "2HYY-hydrogen-water-evidence.png"
        ),
        generation_function=(
            "visualization.hydrogen_water.render_hydrogen_water_evidence"
        ),
        source_records=(
            "data/structures/2HYY.manifest.json",
            "data/structures/2HYY.binding-site.json",
            "data/structures/2HYY.hydrogen-evidence.json",
            "data/structures/2HYY.water-evidence.json",
        ),
        tutorial_anchor=(
            "../README.md#421-visual-learning-layer-chemical-states-and-evidence"
        ),
        teaches="generated versus deposited positions and unassigned water chemistry",
    ),
    FigureSpec(
        figure_id="capstone-interaction-sensitivity",
        figure_number=7,
        title="Interaction sensitivity across tested assumptions",
        image_path=(
            "docs/tutorials/capstone/visualizations/"
            "2HYY-interaction-sensitivity.png"
        ),
        generation_function="visualization.robustness.plot_interaction_sensitivity",
        source_records=("data/structures/2HYY.interaction-robustness.json",),
        tutorial_anchor=(
            "../README.md#422-visual-learning-layer-sensitivity-and-thermodynamics"
        ),
        teaches="independent sensitivity axes and invariant/sensitive/unsupported/unknown statuses",
    ),
    FigureSpec(
        figure_id="harmonic-fep-convergence",
        figure_number=8,
        title="Exact versus estimated harmonic free energy",
        image_path=(
            "docs/tutorials/capstone/visualizations/"
            "harmonic-fep-convergence.png"
        ),
        generation_function="visualization.free_energy.plot_free_energy_convergence",
        source_records=("data/free_energy/harmonic-fep-study.json",),
        tutorial_anchor=(
            "../README.md#422-visual-learning-layer-sensitivity-and-thermodynamics"
        ),
        teaches="exact reference, finite-sample FEP, and effective-sample-fraction overlap",
    ),
    FigureSpec(
        figure_id="capstone-structural-vs-thermodynamics",
        figure_number=9,
        title="Structural evidence is not binding free energy",
        image_path=(
            "docs/tutorials/capstone/visualizations/"
            "2HYY-structural-vs-thermodynamics.png"
        ),
        generation_function=(
            "visualization.thermodynamics.plot_structural_thermodynamics_boundary"
        ),
        source_records=(
            "data/structures/2HYY.thermodynamic-interpretation.json",
        ),
        tutorial_anchor=(
            "../README.md#422-visual-learning-layer-sensitivity-and-thermodynamics"
        ),
        teaches="structural counts versus unestimated thermodynamic terms and requirements",
    ),
)


def figure_spec(figure_id: str) -> FigureSpec:
    try:
        return next(spec for spec in FIGURE_REGISTRY if spec.figure_id == figure_id)
    except StopIteration as error:
        raise ValueError(f"unknown figure ID {figure_id!r}.") from error


def validate_figure_registry(figures: Sequence[FigureMetadata]) -> list[str]:
    """Validate that the metadata covers the registry and its provenance contract."""

    errors: list[str] = []
    by_id = {figure.figure_id: figure for figure in figures}
    registry_ids = {spec.figure_id for spec in FIGURE_REGISTRY}
    if set(by_id) != registry_ids:
        missing = sorted(registry_ids - set(by_id))
        extra = sorted(set(by_id) - registry_ids)
        if missing:
            errors.append(f"metadata is missing registry figures: {missing}.")
        if extra:
            errors.append(f"metadata contains unregistered figures: {extra}.")
    for spec in FIGURE_REGISTRY:
        figure = by_id.get(spec.figure_id)
        if figure is None:
            continue
        if figure.image_path != spec.image_path:
            errors.append(f"{spec.figure_id}: image_path does not match the registry.")
        if tuple(sorted(figure.source_records)) != tuple(sorted(spec.source_records)):
            errors.append(f"{spec.figure_id}: source_records do not match the registry.")
        checksum_names = {key for key, _ in figure.source_checksums}
        if checksum_names != set(spec.source_records):
            errors.append(
                f"{spec.figure_id}: source_checksums do not cover every source record."
            )
        if not figure.caption or not figure.alt_text or not figure.scientific_limitations:
            errors.append(
                f"{spec.figure_id}: caption, alt text, and limitations are required."
            )
        if spec.record_backed and not spec.source_records:
            errors.append(f"{spec.figure_id}: record-backed figure has no sources.")
    return errors


def render_gallery_markdown(figures: Sequence[FigureMetadata]) -> str:
    """Generate the gallery page from the registry and metadata."""

    by_id = {figure.figure_id: figure for figure in figures}
    lines = [
        "# Capstone visualization gallery",
        "",
        "This gallery is generated from the figure registry and metadata. It is",
        "the authoritative index of the nine capstone figures, their sources,",
        "captions, and scientific limitations.",
        "",
        "| Figure | Title |",
        "| --- | --- |",
    ]
    for spec in FIGURE_REGISTRY:
        lines.append(
            f"| {spec.figure_number} | "
            f"[{spec.title}](#figure-{spec.figure_number}) |"
        )
    lines.append("")
    for spec in FIGURE_REGISTRY:
        figure = by_id.get(spec.figure_id)
        if figure is None:
            raise ValueError(f"missing metadata for {spec.figure_id!r}.")
        image_name = Path(figure.image_path).name
        lines.extend(
            [
                f"## Figure {spec.figure_number}",
                "",
                f"![{figure.alt_text}]({image_name})",
                "",
                f"**{figure.title}**",
                "",
                figure.caption,
                "",
                f"**Teaches:** {spec.teaches}",
                "",
                "**Scientific limitations:**",
                "",
                *[f"- {item}" for item in figure.scientific_limitations],
                "",
                "**Source records:**",
                "",
                *[f"- `{item}`" for item in spec.source_records],
                "",
                (
                    f"**Tool:** {figure.visualization_tool} "
                    f"`{figure.tool_version}`; rendering fingerprint "
                    f"`{figure.rendering_config_sha256}`."
                ),
                "",
                f"[Read the teaching narrative]({spec.tutorial_anchor}).",
                "",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def write_gallery(path: str | Path, figures: Sequence[FigureMetadata]) -> Path:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(render_gallery_markdown(figures), encoding="utf-8")
    return destination
