"""Generate the architecture diagram used by docs/development/architecture.md.

Run from the repository root:

    uv run python docs/generate_architecture_visuals.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ASSETS = Path(__file__).resolve().parent / "assets"
ASSETS.mkdir(parents=True, exist_ok=True)

BLUE = "#2b6cb0"
ORANGE = "#dd6b20"
GREEN = "#2f855a"
RED = "#c53030"
GRAY = "#4a5568"
LIGHT = "#ebf4ff"

plt.rcParams.update(
    {
        "figure.dpi": 140,
        "savefig.dpi": 220,
        "font.size": 10,
        "font.family": "DejaVu Sans",
    }
)


def _box(ax, x, y, width, height, text, color, fontsize=9.0):
    rect = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.04,rounding_size=0.08",
        linewidth=1.3,
        edgecolor=color,
        facecolor="white",
        zorder=2,
    )
    ax.add_patch(rect)
    ax.text(
        x + width / 2,
        y + height / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        color=GRAY,
        zorder=3,
    )


def architecture_layers():
    fig, ax = plt.subplots(figsize=(12.4, 7.8))
    ax.set_xlim(0, 12.4)
    ax.set_ylim(0, 7.8)
    ax.axis("off")

    bands = [
        (
            "Applications\nand analysis",
            5.55,
            1.00,
            [
                "docking",
                "free energy\nFEP / TI / MBAR",
                "solvent and\nwater boxes",
                "MSD / RDF",
                "QM/MM",
                "structural records,\nevidence, state matrices",
            ],
            RED,
        ),
        (
            "Algorithms\nand protocols",
            4.15,
            1.00,
            [
                "Velocity Verlet",
                "steepest descent",
                "thermostat",
                "barostat",
                "constraints\nSHAKE / RATTLE",
            ],
            ORANGE,
        ),
        (
            "Composition",
            2.85,
            1.00,
            [
                "System",
                "Topology",
                "NonbondedModel",
                "Box",
                "neighbour\nlist",
                "Ewald",
            ],
            GREEN,
        ),
        (
            "Models and\ncontracts",
            1.55,
            1.00,
            [
                "Particle",
                "pair\npotentials",
                "bonded\nterms",
                "protocols",
            ],
            BLUE,
        ),
        (
            "Foundation",
            0.55,
            0.75,
            [
                "units",
                "validation",
                "backend\n(optional)",
            ],
            GRAY,
        ),
    ]

    left = 2.05
    right = 12.15
    gap = 0.10
    for title, y, height, items, color in bands:
        ax.text(
            0.15,
            y + height / 2,
            title,
            ha="left",
            va="center",
            fontsize=10.5,
            weight="bold",
            color=GRAY,
        )
        width = (right - left - gap * (len(items) - 1)) / len(items)
        for index, item in enumerate(items):
            x = left + index * (width + gap)
            _box(ax, x, y, width, height, item, color, fontsize=8.6)

    ax.text(
        6.1,
        7.45,
        "Moleculab layered architecture",
        ha="center",
        fontsize=15,
        weight="bold",
        color=GRAY,
    )
    ax.text(
        6.1,
        0.18,
        "Dependencies point downward: lower layers do not import higher layers.",
        ha="center",
        fontsize=10,
        color=GRAY,
    )

    for y_top, y_bottom in [(5.52, 5.18), (4.12, 3.88), (2.82, 2.58), (1.52, 1.33)]:
        ax.add_patch(
            FancyArrowPatch(
                (6.1, y_top),
                (6.1, y_bottom),
                arrowstyle="-|>",
                mutation_scale=14,
                linewidth=1.5,
                color=GRAY,
            )
        )

    fig.savefig(ASSETS / "architecture_layers.svg", format="svg", bbox_inches="tight")
    plt.close(fig)


def main():
    architecture_layers()
    for path in ASSETS.glob("*.svg"):
        text = path.read_text()
        cleaned = "\n".join(line.rstrip() for line in text.splitlines()) + "\n"
        if cleaned != text:
            path.write_text(cleaned)
    print(f"Wrote architecture visuals to {ASSETS}")


if __name__ == "__main__":
    main()
