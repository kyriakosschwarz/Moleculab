"""Interaction-evidence matrix figures from combination records."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..applications.combination_interactions import (
    CATEGORY_ORDER,
    CombinationInteractionRecord,
)
from .metadata import FigureContext, FigureMetadata

_STATUS_INDEX = {
    "classified": 0,
    "not_classified": 1,
    "unknown": 2,
    "not_applicable": 3,
}
_STATUS_LABEL = ("present", "absent", "unknown", "N/A")
_STATUS_COLOURS = ("#2f855a", "#cbd5e0", "#d69e2e", "#a0aec0")


def plot_evidence_matrix(
    *,
    combination_record: CombinationInteractionRecord,
    context: FigureContext,
    output_png: str | Path,
    output_svg: str | Path | None,
    figure_id: str,
    title: str,
    caption: str,
    alt_text: str,
    source_records: tuple[str, ...],
    source_checksums: tuple[tuple[str, str], ...],
    width: float = 13.0,
    height: float = 12.0,
    dpi: int = 200,
) -> FigureMetadata:
    """Plot present/absent/unknown/not-applicable evidence as a 6-panel matrix."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    ligand_ids = tuple(combination_record.ligand_state_ids)
    protein_ids = tuple(combination_record.protein_state_ids)
    by_pair = {
        (item.ligand_state_id, item.protein_state_id): item
        for item in combination_record.combinations
    }
    colormap = ListedColormap(_STATUS_COLOURS)

    figure, axes = plt.subplots(
        3,
        2,
        figsize=(width, height),
        constrained_layout=False,
    )
    figure.subplots_adjust(
        top=0.92,
        bottom=0.13,
        left=0.16,
        right=0.98,
        hspace=0.45,
        wspace=0.42,
    )
    axes = axes.ravel()
    for axis, category in zip(axes, CATEGORY_ORDER, strict=True):
        matrix = []
        for ligand_id in ligand_ids:
            row = []
            for protein_id in protein_ids:
                combination = by_pair[(ligand_id, protein_id)]
                status = dict(combination.category_statuses)[category]
                row.append(_STATUS_INDEX[status])
            matrix.append(row)
        axis.imshow(
            matrix,
            cmap=colormap,
            vmin=-0.5,
            vmax=3.5,
            aspect="auto",
        )
        axis.set_title(category.replace("_", " "))
        axis.set_xticks(range(len(protein_ids)))
        axis.set_xticklabels(protein_ids, rotation=25, ha="right")
        axis.set_yticks(range(len(ligand_ids)))
        axis.set_yticklabels(ligand_ids)
        for row_index, row in enumerate(matrix):
            for column_index, value in enumerate(row):
                axis.text(
                    column_index,
                    row_index,
                    _STATUS_LABEL[value],
                    ha="center",
                    va="center",
                    fontsize=7.5,
                    color="black",
                )
    handles = [
        Patch(facecolor=colour, edgecolor="none", label=label)
        for colour, label in zip(_STATUS_COLOURS, _STATUS_LABEL, strict=True)
    ]
    figure.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.035),
        ncol=4,
        frameon=False,
    )
    figure.suptitle(title, fontsize=15, weight="bold")
    figure.text(
        0.5,
        0.006,
        "Unknown is not zero; N/A is outside the assessed feature model.",
        ha="center",
        va="bottom",
        fontsize=9,
    )
    output_png_path = Path(output_png)
    output_png_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_png_path, dpi=dpi, facecolor="white")
    if output_svg is not None:
        output_svg_path = Path(output_svg)
        output_svg_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_svg_path, facecolor="white")
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
            ("layout", "3x2 category panels"),
            ("present_colour", _STATUS_COLOURS[0]),
            ("absent_colour", _STATUS_COLOURS[1]),
            ("unknown_colour", _STATUS_COLOURS[2]),
            ("not_applicable_colour", _STATUS_COLOURS[3]),
            ("cell_labels", "present, absent, unknown, N/A"),
        ),
        rendering_config_sha256=None,
        scientific_limitations=(
            "The matrix plots an existing evidence record; it does not reclassify interactions.",
            "Unknown is not zero and is not treated as an absent interaction.",
            "N/A means outside the assessed feature model, not a negative result.",
            "No state ranking, population, energy, or affinity is implied.",
        ),
        alt_text=alt_text,
        caption=caption,
        camera_view=None,
    )
