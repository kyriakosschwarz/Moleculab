"""Shared visual semantics for Moleculab figures.

The style module centralises evidence-state colors, robustness-state colors,
typography, and export conventions.  Scientific mapping decisions live here so
that the same status always receives the same visual treatment.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

EVIDENCE_STATUSES = ("present", "absent", "unknown", "not_applicable")
ROBUSTNESS_STATUSES = (
    "invariant",
    "assumption_sensitive",
    "unsupported",
    "unknown",
)
REQUIREMENT_STATUSES = ("met", "partial", "not_run", "not_available")
TERM_STATUSES = ("not_estimated", "not_available", "not_applicable")


@dataclass(frozen=True)
class StatusStyle:
    status: str
    label: str
    colour: str
    symbol: str


EVIDENCE_STYLES = {
    "present": StatusStyle("present", "present", "#2f855a", "P"),
    "absent": StatusStyle("absent", "absent", "#cbd5e0", "A"),
    "unknown": StatusStyle("unknown", "unknown", "#d69e2e", "U"),
    "not_applicable": StatusStyle(
        "not_applicable", "not applicable", "#a0aec0", "N/A"
    ),
}

ROBUSTNESS_STYLES = {
    "invariant": StatusStyle("invariant", "invariant", "#2f855a", "I"),
    "assumption_sensitive": StatusStyle(
        "assumption_sensitive", "assumption-sensitive", "#d69e2e", "S"
    ),
    "unsupported": StatusStyle("unsupported", "unsupported", "#c53030", "U"),
    "unknown": StatusStyle("unknown", "unknown", "#805ad5", "?"),
}

REQUIREMENT_STYLES = {
    "met": StatusStyle("met", "met", "#2f855a", "M"),
    "partial": StatusStyle("partial", "partial", "#d69e2e", "P"),
    "not_run": StatusStyle("not_run", "not run", "#cbd5e0", "N"),
    "not_available": StatusStyle(
        "not_available", "not available", "#a0aec0", "-"
    ),
}

TERM_STYLES = {
    "not_estimated": StatusStyle(
        "not_estimated", "not estimated", "#cbd5e0", "-"
    ),
    "not_available": StatusStyle(
        "not_available", "not available", "#a0aec0", "-"
    ),
    "not_applicable": StatusStyle(
        "not_applicable", "not applicable", "#a0aec0", "N/A"
    ),
}


def evidence_style(status: str) -> StatusStyle:
    try:
        return EVIDENCE_STYLES[status]
    except KeyError as error:
        raise ValueError(f"unknown evidence status {status!r}.") from error


def robustness_style(status: str) -> StatusStyle:
    try:
        return ROBUSTNESS_STYLES[status]
    except KeyError as error:
        raise ValueError(f"unknown robustness status {status!r}.") from error


def requirement_style(status: str) -> StatusStyle:
    try:
        return REQUIREMENT_STYLES[status]
    except KeyError as error:
        raise ValueError(f"unknown requirement status {status!r}.") from error


def term_style(status: str) -> StatusStyle:
    try:
        return TERM_STYLES[status]
    except KeyError as error:
        raise ValueError(f"unknown thermodynamic term status {status!r}.") from error


def configure_matplotlib() -> None:
    """Apply the shared plotting conventions."""

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "svg.hashsalt": "moleculab",
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.edgecolor": "#4a5568",
            "axes.labelcolor": "#2d3748",
            "xtick.color": "#2d3748",
            "ytick.color": "#2d3748",
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def clean_svg_file(path: str | Path) -> None:
    """Remove trailing whitespace from generated SVG text."""

    destination = Path(path)
    text = destination.read_text(encoding="utf-8")
    text = re.sub(r"\s*<dc:date>[^<]*</dc:date>", "", text)
    cleaned = "\n".join(line.rstrip() for line in text.splitlines()) + "\n"
    if cleaned != text:
        destination.write_text(cleaned, encoding="utf-8")


def save_figure(
    figure: Any,
    output_png: str | Path,
    output_svg: str | Path | None,
    *,
    dpi: int = 200,
) -> None:
    figure.savefig(Path(output_png), dpi=dpi, facecolor="white")
    if output_svg is not None:
        figure.savefig(Path(output_svg), facecolor="white")
        clean_svg_file(output_svg)


def style_axes(axis: Any) -> None:
    axis.grid(axis="y", color="#e2e8f0", linewidth=0.8)
    axis.set_axisbelow(True)
