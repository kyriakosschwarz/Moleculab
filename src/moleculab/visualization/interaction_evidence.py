"""Interaction-evidence matrix figures from combination records."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from ..applications.interactions.combination_interactions import (
    CATEGORY_ORDER,
    CombinationInteractionRecord,
)
from .metadata import FigureContext, FigureMetadata
from .styles import (
    EVIDENCE_STYLES,
    configure_matplotlib,
    evidence_style,
    save_figure,
)

_STATUS_INDEX = {
    "classified": 0,
    "not_classified": 1,
    "unknown": 2,
    "not_applicable": 3,
}
_STATUS_NAMES = ("present", "absent", "unknown", "not_applicable")


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

    configure_matplotlib()
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from matplotlib.patches import Patch

    ligand_ids = tuple(combination_record.ligand_state_ids)
    protein_ids = tuple(combination_record.protein_state_ids)
    by_pair = {
        (item.ligand_state_id, item.protein_state_id): item
        for item in combination_record.combinations
    }
    status_styles = tuple(evidence_style(status) for status in _STATUS_NAMES)
    colormap = ListedColormap(
        [style.colour for style in status_styles]
    )

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
                    status_styles[value].label,
                    ha="center",
                    va="center",
                    fontsize=7.5,
                    color=status_styles[value].text_colour,
                )
    handles = [
        Patch(facecolor=style.colour, edgecolor="none", label=style.label)
        for style in status_styles
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
            ("layout", "3x2 category panels"),
            ("present_colour", EVIDENCE_STYLES["present"].colour),
            ("absent_colour", EVIDENCE_STYLES["absent"].colour),
            ("unknown_colour", EVIDENCE_STYLES["unknown"].colour),
            ("not_applicable_colour", EVIDENCE_STYLES["not_applicable"].colour),
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
