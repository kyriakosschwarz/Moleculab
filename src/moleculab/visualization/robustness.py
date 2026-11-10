"""Interaction-sensitivity figures from robustness records."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..applications.interaction_robustness import InteractionRobustnessRecord
from .metadata import FigureContext, FigureMetadata
from .styles import (
    ROBUSTNESS_STATUSES,
    configure_matplotlib,
    robustness_style,
    save_figure,
    style_axes,
)


def plot_interaction_sensitivity(
    *,
    robustness_record: InteractionRobustnessRecord,
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
    height: float = 12.0,
    dpi: int = 200,
) -> FigureMetadata:
    """Plot contact, angular, water, and placement sensitivity as separate panels."""

    import matplotlib

    configure_matplotlib()
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    figure = plt.figure(figsize=(width, 9.5), constrained_layout=False)
    grid = figure.add_gridspec(2, 2)
    figure.subplots_adjust(
        top=0.90,
        bottom=0.14,
        left=0.09,
        right=0.98,
        hspace=0.40,
        wspace=0.26,
    )

    contact_axis = figure.add_subplot(grid[0, 0])
    cutoffs = [item.cutoff_angstrom for item in robustness_record.contact_cutoff_results]
    minimums = [
        item.minimum_classified_per_combination
        for item in robustness_record.contact_cutoff_results
    ]
    maximums = [
        item.maximum_classified_per_combination
        for item in robustness_record.contact_cutoff_results
    ]
    contact_axis.plot(cutoffs, minimums, marker="o", color="#2b6cb0", label="min")
    contact_axis.plot(cutoffs, maximums, marker="s", color="#dd6b20", label="max")
    contact_axis.fill_between(cutoffs, minimums, maximums, color="#bee3f8", alpha=0.4)
    for x, y in zip(cutoffs, maximums, strict=True):
        contact_axis.annotate(str(y), (x, y), textcoords="offset points", xytext=(0, 6), ha="center")
    contact_axis.set_title("Contact cutoff")
    contact_axis.set_xlabel("contact cutoff (Å)")
    contact_axis.set_ylabel("classified interactions\nper state combination")
    contact_axis.legend(frameon=False)
    style_axes(contact_axis)

    angle_axis = figure.add_subplot(grid[0, 1])
    thresholds = [
        item.min_angle_degrees for item in robustness_record.hydrogen_angle_results
    ]
    supported = [
        item.supported_candidate_count
        for item in robustness_record.hydrogen_angle_results
    ]
    angle_axis.step(thresholds, supported, where="mid", color="#805ad5", linewidth=2)
    angle_axis.scatter(thresholds, supported, color="#805ad5", zorder=3)
    angle_axis.set_title("Hydrogen-bond angle")
    angle_axis.set_xlabel("minimum D-H...A angle (degrees)")
    angle_axis.set_ylabel("supported candidates")
    angle_axis.set_yticks(range(0, max(supported) + 2))
    style_axes(angle_axis)

    water_axis = figure.add_subplot(grid[1, 0])
    water_cutoffs = [
        item.cutoff_angstrom for item in robustness_record.water_cutoff_results
    ]
    plausible = [
        item.plausible_count for item in robustness_record.water_cutoff_results
    ]
    water_axis.step(water_cutoffs, plausible, where="mid", color="#2f855a", linewidth=2)
    water_axis.scatter(water_cutoffs, plausible, color="#2f855a", zorder=3)
    water_axis.set_title("Water-bridge cutoff")
    water_axis.set_xlabel("water-contact cutoff (Å)")
    water_axis.set_ylabel("geometrically plausible bridges")
    water_axis.set_yticks(range(0, max(plausible) + 2))
    style_axes(water_axis)

    placement_axis = figure.add_subplot(grid[1, 1])
    candidates = [
        item.candidate_id for item in robustness_record.hydrogen_placement_results
    ]
    offsets = list(robustness_record.hydrogen_placement_offsets_degrees)
    matrix = []
    for item in robustness_record.hydrogen_placement_results:
        supported_set = set(item.supported_offsets_degrees)
        matrix.append([1 if offset in supported_set else 0 for offset in offsets])
    placement_axis.imshow(
        matrix,
        cmap=ListedColormap(("#cbd5e0", "#2f855a")),
        vmin=-0.5,
        vmax=1.5,
        aspect="auto",
    )
    placement_axis.set_xticks(range(len(offsets)))
    placement_axis.set_xticklabels([f"{offset:+g}" for offset in offsets])
    placement_axis.set_yticks(range(len(candidates)))
    placement_axis.set_yticklabels(candidates)
    placement_axis.set_xlabel("D-H orientation perturbation (degrees)")
    placement_axis.set_title("Hydrogen-placement perturbation")
    for row_index, row in enumerate(matrix):
        for column_index, value in enumerate(row):
            placement_axis.text(
                column_index,
                row_index,
                "S" if value else "-",
                ha="center",
                va="center",
                color="white" if value else "#2d3748",
                fontsize=9,
            )

    counts = {
        status: sum(
            1
            for conclusion in robustness_record.conclusions
            if conclusion.status == status
        )
        for status in ROBUSTNESS_STATUSES
    }
    handles = [
        Patch(
            facecolor=robustness_style(status).colour,
            edgecolor="none",
            label=f"{robustness_style(status).label} ({counts[status]})",
        )
        for status in ROBUSTNESS_STATUSES
    ]
    figure.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.045),
        ncol=4,
        frameon=False,
        title="Conclusion statuses (no composite score)",
        title_fontsize=10,
    )

    figure.suptitle(title, fontsize=15, weight="bold")
    figure.text(
        0.5,
        0.012,
        "Each axis is independent; invariant, sensitive, unsupported, and unknown are not probabilities.",
        ha="center",
        va="bottom",
        fontsize=9,
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
            ("panels", "contact, angle, water, placement, conclusion statuses"),
            ("invariant_colour", robustness_style("invariant").colour),
            ("assumption_sensitive_colour", robustness_style("assumption_sensitive").colour),
            ("unsupported_colour", robustness_style("unsupported").colour),
            ("unknown_colour", robustness_style("unknown").colour),
        ),
        rendering_config_sha256=None,
        scientific_limitations=(
            "Each sensitivity axis is independent; the axes are not combined into a score.",
            "Invariant means stable only across the tested assumptions.",
            "Unsupported means no tested condition supports the positive claim.",
            "Unknown means the necessary evidence is unavailable.",
            "No state ranking, population, energy, or affinity is implied.",
        ),
        alt_text=alt_text,
        caption=caption,
        camera_view=None,
    )
