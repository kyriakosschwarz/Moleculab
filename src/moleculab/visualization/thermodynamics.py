"""Structural-evidence versus thermodynamic-requirement figures."""

from __future__ import annotations

from pathlib import Path

from ..applications.thermodynamics.thermodynamic_interpretation import (
    ThermodynamicInterpretationRecord,
)
from .metadata import FigureContext, FigureMetadata
from .styles import (
    configure_matplotlib,
    requirement_style,
    save_figure,
    style_axes,
    term_style,
)


def plot_structural_thermodynamics_boundary(
    *,
    record: ThermodynamicInterpretationRecord,
    context: FigureContext,
    output_png: str | Path,
    output_svg: str | Path | None,
    figure_id: str,
    title: str,
    caption: str,
    alt_text: str,
    source_records: tuple[str, ...],
    source_checksums: tuple[tuple[str, str], ...],
    width: float = 14.0,
    height: float = 9.0,
    dpi: int = 200,
) -> FigureMetadata:
    """Show structural counts next to the unestimated thermodynamic terms."""

    import matplotlib

    configure_matplotlib()
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    figure = plt.figure(figsize=(width, height), constrained_layout=False)
    grid = figure.add_gridspec(2, 2, height_ratios=(1.0, 0.55))
    figure.subplots_adjust(
        top=0.88,
        bottom=0.12,
        left=0.08,
        right=0.97,
        hspace=0.35,
        wspace=0.22,
    )

    structural_axis = figure.add_subplot(grid[0, 0])
    summary = record.structural_summary
    labels = ("hydrophobic\npairs", "aromatic\npairs", "supported\nH-bonds", "plausible\nwater bridges")
    values = (
        summary.hydrophobic_pair_count,
        summary.aromatic_pair_count,
        summary.supported_hydrogen_bond_count,
        summary.plausible_water_bridge_count,
    )
    colours = ("#2f855a", "#2f855a", "#805ad5", "#2b6cb0")
    structural_axis.bar(labels, values, color=colours)
    for index, value in enumerate(values):
        structural_axis.text(index, value + 0.25, str(value), ha="center", fontsize=10)
    structural_axis.set_ylim(0, max(values) + 3)
    structural_axis.set_ylabel("count in the reference state combination")
    structural_axis.set_title("Structural evidence (not energies)")
    style_axes(structural_axis)

    terms_axis = figure.add_subplot(grid[0, 1])
    terms_axis.axis("off")
    terms_axis.set_title("Missing thermodynamic terms", loc="left")
    y = 0.92
    for term in record.terms:
        style = term_style(term.status)
        terms_axis.add_patch(
            plt.Rectangle(
                (0.0, y - 0.045),
                0.035,
                0.07,
                transform=terms_axis.transAxes,
                facecolor=style.colour,
                edgecolor="none",
                clip_on=False,
            )
        )
        terms_axis.text(
            0.05,
            y,
            term.name,
            transform=terms_axis.transAxes,
            va="center",
            fontsize=9,
            color="#2d3748",
        )
        terms_axis.text(
            0.76,
            y,
            style.label,
            transform=terms_axis.transAxes,
            va="center",
            fontsize=8,
            color="#4a5568",
        )
        y -= 0.105
    terms_axis.text(
        0.0,
        0.02,
        "No term is filled with an illustrative or fitted number.",
        transform=terms_axis.transAxes,
        fontsize=8,
        color="#4a5568",
    )

    requirements_axis = figure.add_subplot(grid[1, 0])
    requirement_counts = {}
    for item in record.requirements:
        requirement_counts[item.status] = requirement_counts.get(item.status, 0) + 1
    statuses = tuple(requirement_counts)
    requirements_axis.bar(
        [requirement_style(status).label for status in statuses],
        [requirement_counts[status] for status in statuses],
        color=[requirement_style(status).colour for status in statuses],
    )
    for index, status in enumerate(statuses):
        requirements_axis.text(
            index,
            requirement_counts[status] + 0.1,
            str(requirement_counts[status]),
            ha="center",
            fontsize=10,
        )
    requirements_axis.set_ylabel("requirements")
    requirements_axis.set_title("Preparation and sampling requirements")
    requirements_axis.set_ylim(0, max(requirement_counts.values()) + 2)
    style_axes(requirements_axis)

    summary_axis = figure.add_subplot(grid[1, 1])
    summary_axis.axis("off")
    summary_axis.add_patch(
        plt.Rectangle(
            (0.02, 0.18),
            0.96,
            0.64,
            transform=summary_axis.transAxes,
            facecolor="#fff5f5",
            edgecolor="#c53030",
            linewidth=1.5,
        )
    )
    summary_axis.text(
        0.5,
        0.64,
        "Binding free energy",
        transform=summary_axis.transAxes,
        ha="center",
        fontsize=12,
        weight="bold",
        color="#742a2a",
    )
    summary_axis.text(
        0.5,
        0.44,
        record.free_energy_status,
        transform=summary_axis.transAxes,
        ha="center",
        fontsize=15,
        weight="bold",
        color="#c53030",
    )
    summary_axis.text(
        0.5,
        0.27,
        "No numeric affinity or free-energy estimate is present.",
        transform=summary_axis.transAxes,
        ha="center",
        fontsize=9,
        color="#4a5568",
    )

    figure.suptitle(title, fontsize=15, weight="bold")
    figure.text(
        0.5,
        0.012,
        "Structural counts are observations; thermodynamic terms are missing model, sampling, and uncertainty analysis.",
        ha="center",
        va="bottom",
        fontsize=9,
        color="#4a5568",
    )
    output_png_path = Path(output_png)
    output_png_path.parent.mkdir(parents=True, exist_ok=True)
    if output_svg is not None:
        Path(output_svg).parent.mkdir(parents=True, exist_ok=True)
    save_figure(figure, output_png_path, output_svg, dpi=dpi)
    plt.close(figure)

    return FigureMetadata(
        figure_id=figure_id,
        title=title,
        image_path=str(output_png_path),
        context=context,
        visualization_tool="Matplotlib",
        tool_version=matplotlib.__version__,
        source_records=source_records,
        source_checksums=source_checksums,
        rendering_config=(
            ("backend", "Matplotlib"),
            ("width_inches", f"{width:g}"),
            ("height_inches", f"{height:g}"),
            ("dpi", str(dpi)),
            ("left_panel", "structural counts"),
            ("terms_panel", "unestimated thermodynamic terms"),
            ("requirements_panel", "requirement status counts"),
            ("free_energy_status", record.free_energy_status),
        ),
        rendering_config_sha256=None,
        scientific_limitations=(
            "Structural counts are not interaction energies or free energies.",
            "Every listed thermodynamic term is unestimated in this record.",
            "No binding affinity, population, or state ranking is implied.",
            "The figure does not supply force-field, solvent, sampling, or uncertainty analysis.",
        ),
        alt_text=alt_text,
        caption=caption,
        camera_view=None,
    )
