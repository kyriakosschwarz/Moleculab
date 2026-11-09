"""Generate the figures and animations used by the individual tutorials.

Run from the repository root:

    uv run python docs/tutorials/individual/generate_visuals.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

from moleculab import (
    HarmonicPotential,
    Particle,
    kinetic_energy,
    simulate_velocity_verlet,
)

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
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.22,
        "font.family": "DejaVu Sans",
        "legend.frameon": False,
    }
)


def save_svg(fig, name):
    fig.savefig(ASSETS / name, format="svg", bbox_inches="tight")
    plt.close(fig)


def save_gif(animation, name, fps=25):
    animation.save(ASSETS / name, writer=PillowWriter(fps=fps))
    plt.close(animation._fig)


def strip_svg_trailing_whitespace():
    """Keep generated SVG diffs free of trailing spaces."""

    for path in ASSETS.glob("*.svg"):
        text = path.read_text()
        cleaned = "\n".join(line.rstrip() for line in text.splitlines()) + "\n"
        if cleaned != text:
            path.write_text(cleaned)


def unit_scale():
    fig, ax = plt.subplots(figsize=(10.5, 3.6))
    objects = [
        (0.1, "atom", "0.1 nm"),
        (0.15, "bond", "0.15 nm"),
        (1.0, "small molecule", "1 nm"),
        (6.0, "protein", "6 nm"),
        (100.0, "virus", "100 nm"),
        (10_000.0, "cell", "10 µm"),
        (1_700_000_000.0, "human", "1.7 m"),
    ]
    ax.set_xscale("log")
    ax.set_xlim(0.05, 5e9)
    ax.set_ylim(-1.2, 1.2)
    ax.set_yticks([])
    ax.axhline(0, color=GRAY, lw=2, alpha=0.5)
    for index, (value, label, text) in enumerate(objects):
        y = 0.35 if index % 2 == 0 else -0.35
        ax.plot([value], [0], "o", ms=11, color=BLUE if index % 2 == 0 else ORANGE)
        ax.plot([value, value], [0, y * 0.75], color=GRAY, lw=1)
        ax.text(
            value,
            y,
            f"{label}\n{text}",
            ha="center",
            va="center" if y > 0 else "center",
            fontsize=10,
            bbox=dict(boxstyle="round,pad=0.25", fc=LIGHT, ec="none", alpha=0.9),
        )
    ax.set_xlabel("length (nm, logarithmic scale)")
    ax.set_title("The molecular world spans many orders of magnitude")
    ax.grid(False)
    save_svg(fig, "unit_scale.svg")


def unit_bridge():
    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 6)
    ax.axis("off")

    def box(x, y, w, h, text, color):
        rect = FancyBboxPatch(
            (x, y),
            w,
            h,
            boxstyle="round,pad=0.08,rounding_size=0.12",
            linewidth=1.4,
            edgecolor=color,
            facecolor="white",
        )
        ax.add_patch(rect)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=10)

    def arrow(x1, y1, x2, y2, label):
        ax.add_patch(
            FancyArrowPatch(
                (x1, y1),
                (x2, y2),
                arrowstyle="-|>",
                mutation_scale=15,
                linewidth=1.4,
                color=GRAY,
            )
        )
        ax.text((x1 + x2) / 2, (y1 + y2) / 2 + 0.18, label, ha="center", fontsize=9, color=GRAY)

    ax.text(5, 5.55, "The internal unit system makes a = F/m numerically direct", ha="center", fontsize=14, weight="bold")
    box(0.6, 3.8, 2.7, 1.15, "force\nkJ/(mol nm)", BLUE)
    box(0.6, 2.25, 2.7, 1.15, "mass\nDa", ORANGE)
    box(5.75, 3.03, 3.25, 1.15, "acceleration\nnm/ps²", GREEN)
    arrow(3.35, 4.37, 5.65, 3.95, "F = 10")
    arrow(3.35, 2.82, 5.65, 3.28, "m = 12")
    ax.text(6.95, 2.48, "a = F/m = 10/12 = 0.8333 nm/ps²", ha="center", fontsize=12.5, color=GREEN, weight="bold")
    ax.text(
        5.0,
        1.82,
        "Why no extra constant?  The SI conversion factors for kJ/mol, nm, ps and Da cancel exactly.",
        ha="center",
        fontsize=10.5,
        color=GRAY,
    )
    box(0.4, 0.55, 4.1, 0.9, "SI force per particle: kJ → J, divide by N_A", BLUE)
    box(5.0, 0.55, 4.6, 0.9, "SI mass: Da → kg, then convert m/s² → nm/ps²", ORANGE)
    ax.text(5.0, 0.16, "Mixing Å, kcal/mol, fs or kg without converting breaks this shortcut.", ha="center", fontsize=9.5, color=RED)
    save_svg(fig, "unit_bridge.svg")


def learning_ladder():
    fig, ax = plt.subplots(figsize=(11.0, 4.6))
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 4.6)
    ax.axis("off")

    steps = [
        ("01", "Units", "nm, ps, kJ/mol, Da\nmake numbers physical"),
        ("02", "Particle", "mass, position, velocity\nkinetic energy"),
        ("03", "Landscape", "U(x) and F = -dU/dx\nharmonic model"),
        ("04", "Motion", "Velocity Verlet\nenergy stays bounded"),
    ]
    width = 2.35
    gap = 0.34
    start = 0.28
    for index, (number, title, body) in enumerate(steps):
        x = start + index * (width + gap)
        color = [BLUE, ORANGE, GREEN, RED][index]
        rect = FancyBboxPatch(
            (x, 2.55),
            width,
            1.35,
            boxstyle="round,pad=0.06,rounding_size=0.12",
            linewidth=1.6,
            edgecolor=color,
            facecolor="white",
        )
        ax.add_patch(rect)
        ax.text(x + 0.18, 3.62, number, fontsize=18, weight="bold", color=color)
        ax.text(x + width / 2, 3.42, title, ha="center", fontsize=11.5, weight="bold")
        ax.text(x + width / 2, 2.92, body, ha="center", va="center", fontsize=8.4, color=GRAY)
        if index < len(steps) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (x + width + 0.03, 3.22),
                    (x + width + gap - 0.03, 3.22),
                    arrowstyle="-|>",
                    mutation_scale=14,
                    linewidth=1.4,
                    color=GRAY,
                )
            )

    ax.text(5.5, 1.95, "One particle, reused at every step", ha="center", fontsize=12, weight="bold", color=GRAY)
    details = [
        ("mass", "12 Da"),
        ("position", "(0.1, 0, 0) nm"),
        ("velocity", "(0.1, 0, 0) nm/ps"),
        ("speed", "0.10 nm/ps"),
        ("kinetic energy", "0.06 kJ/mol"),
    ]
    for index, (label, value) in enumerate(details):
        x = 0.55 + index * 2.05
        ax.text(x, 1.55, label, fontsize=9, color=GRAY)
        ax.text(x, 1.22, value, fontsize=10, weight="bold", color=BLUE)
        ax.plot([x, x + 1.75], [1.05, 1.05], color=LIGHT, lw=4)
    ax.text(5.5, 0.55, "Tutorial 01 chooses the language  →  02 describes the state  →  03 adds a landscape  →  04 makes it move",
            ha="center", fontsize=9.5, color=GRAY)
    save_svg(fig, "learning_ladder.svg")


def stage_flow():
    fig, ax = plt.subplots(figsize=(10.8, 2.9))
    ax.set_xlim(0, 10.8)
    ax.set_ylim(0, 2.9)
    ax.axis("off")
    labels = [
        ("Physical\nmeaning", BLUE),
        ("Mathematical\ndefinition", ORANGE),
        ("Units and\ndimensions", GREEN),
        ("Python\nimplementation", RED),
        ("Numerical\nverification", GRAY),
    ]
    width = 1.85
    gap = 0.24
    start = 0.22
    for index, (label, color) in enumerate(labels):
        x = start + index * (width + gap)
        rect = FancyBboxPatch(
            (x, 0.95),
            width,
            1.05,
            boxstyle="round,pad=0.06,rounding_size=0.12",
            linewidth=1.5,
            edgecolor=color,
            facecolor="white",
        )
        ax.add_patch(rect)
        ax.text(x + width / 2, 1.47, label, ha="center", va="center", fontsize=10.5)
        if index < len(labels) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (x + width + 0.02, 1.47),
                    (x + width + gap - 0.02, 1.47),
                    arrowstyle="-|>",
                    mutation_scale=13,
                    linewidth=1.3,
                    color=GRAY,
                )
            )
    ax.text(5.4, 2.5, "The pattern behind every section of this course", ha="center", fontsize=12.5, weight="bold")
    ax.text(5.4, 0.45, "A model is not finished until the equation, the units, the code, and the test all agree.",
            ha="center", fontsize=9.5, color=GRAY)
    save_svg(fig, "stage_flow.svg")


def particle_state():
    fig, ax = plt.subplots(figsize=(10.4, 4.8))
    ax.set_xlim(0, 10.4)
    ax.set_ylim(0, 4.8)
    ax.axis("off")

    stored = [
        ("mass", "scalar, Da", "particle parameter"),
        ("position", "vector, nm", "current state"),
        ("velocity", "vector, nm/ps", "current state"),
    ]
    derived = [
        ("speed", "|v|", "from velocity"),
        ("kinetic energy", "½ m |v|²", "from mass and velocity"),
        ("force", "-∇U(r)", "from the potential"),
        ("acceleration", "F / m", "from force and mass"),
    ]

    ax.text(1.75, 4.35, "Stored state", ha="center", fontsize=12.5, weight="bold", color=BLUE)
    ax.text(7.0, 4.35, "Derived quantities", ha="center", fontsize=12.5, weight="bold", color=GREEN)

    for index, (name, symbol, reason) in enumerate(stored):
        y = 3.35 - index * 1.02
        rect = FancyBboxPatch(
            (0.35, y - 0.34),
            2.8,
            0.68,
            boxstyle="round,pad=0.05,rounding_size=0.1",
            linewidth=1.4,
            edgecolor=BLUE,
            facecolor="white",
        )
        ax.add_patch(rect)
        ax.text(0.55, y + 0.1, name, fontsize=10.5, weight="bold")
        ax.text(0.55, y - 0.16, f"{symbol}  ·  {reason}", fontsize=8.6, color=GRAY)

    for index, (name, symbol, reason) in enumerate(derived):
        y = 3.65 - index * 0.92
        rect = FancyBboxPatch(
            (5.55, y - 0.3),
            4.45,
            0.6,
            boxstyle="round,pad=0.05,rounding_size=0.1",
            linewidth=1.3,
            edgecolor=GREEN,
            facecolor="white",
        )
        ax.add_patch(rect)
        ax.text(5.75, y + 0.09, name, fontsize=10.2, weight="bold")
        ax.text(7.55, y + 0.09, symbol, fontsize=9.4, color=GREEN)
        ax.text(5.75, y - 0.17, reason, fontsize=8.4, color=GRAY)

    ax.add_patch(
        FancyArrowPatch(
            (3.25, 2.55),
            (5.45, 2.55),
            arrowstyle="-|>",
            mutation_scale=16,
            linewidth=1.5,
            color=GRAY,
            connectionstyle="arc3,rad=-0.15",
        )
    )
    ax.text(4.35, 2.78, "always recomputed\nfrom the stored state", ha="center", fontsize=8.8, color=GRAY)
    ax.text(5.2, 0.35, "One source of truth: the state. Everything else is a calculation.",
            ha="center", fontsize=9.5, color=GRAY)
    save_svg(fig, "particle_state.svg")


def verlet_flow():
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    ax.set_xlim(0, 9.2)
    ax.set_ylim(0, 5.0)
    ax.axis("off")

    steps = [
        ("1", r"$\mathbf{F}_n = -\nabla U(\mathbf{r}_n)$", "force at the current position"),
        ("2", r"$\mathbf{a}_n = \mathbf{F}_n / m$", "current acceleration"),
        ("3", r"$\mathbf{r}_{n+1} = \mathbf{r}_n + \mathbf{v}_n\Delta t + \frac{1}{2}\mathbf{a}_n\Delta t^2$", "move the position"),
        ("4", r"$\mathbf{F}_{n+1} = -\nabla U(\mathbf{r}_{n+1})$", "force at the new position"),
        ("5", r"$\mathbf{v}_{n+1} = \mathbf{v}_n + \frac{1}{2}(\mathbf{a}_n + \mathbf{a}_{n+1})\Delta t$", "average the acceleration"),
    ]
    y_positions = [4.35, 3.55, 2.72, 1.9, 1.08]
    for index, (number, equation, note) in enumerate(steps):
        y = y_positions[index]
        color = GREEN if index < 2 else (BLUE if index == 2 else (ORANGE if index == 3 else RED))
        circle = plt.Circle((0.7, y), 0.22, color=color, alpha=0.9)
        ax.add_patch(circle)
        ax.text(0.7, y, number, ha="center", va="center", fontsize=10, color="white", weight="bold")
        ax.text(1.15, y + 0.11, equation, fontsize=12.5, va="center")
        ax.text(1.15, y - 0.2, note, fontsize=9, va="center", color=GRAY)
        if index < len(steps) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (1.15, y - 0.34),
                    (1.15, y_positions[index + 1] + 0.25),
                    arrowstyle="-|>",
                    mutation_scale=12,
                    linewidth=1.2,
                    color=GRAY,
                )
            )
    ax.text(4.6, 0.35, "Repeat: $\\mathbf{F}_{n+1}$ and $\\mathbf{a}_{n+1}$ become the next step's current values.",
            ha="center", fontsize=9.5, color=GRAY)
    save_svg(fig, "verlet_flow.svg")


def particle_velocity():
    fig, ax = plt.subplots(figsize=(6.4, 5.0))
    vx, vy = 0.08, 0.06
    ax.quiver(0, 0, vx, vy, angles="xy", scale_units="xy", scale=1, color=BLUE, width=0.012)
    ax.plot([0, vx], [0, 0], color=GRAY, ls="--", lw=1)
    ax.plot([vx, vx], [0, vy], color=GRAY, ls="--", lw=1)
    ax.plot([0], [0], "o", color=ORANGE, ms=13)
    ax.text(vx / 2, -0.012, "$v_x = 0.08$ nm/ps", ha="center", color=GRAY)
    ax.text(vx + 0.006, vy / 2, "$v_y = 0.06$ nm/ps", va="center", color=GRAY)
    ax.text(vx * 0.58, vy * 0.58 + 0.012, r"$v$", color=BLUE, fontsize=15)
    ax.text(0.0, 0.105, r"$|v|=\sqrt{v_x^2+v_y^2}=0.10$ nm/ps", fontsize=12)
    ax.set_xlim(-0.02, 0.13)
    ax.set_ylim(-0.02, 0.12)
    ax.set_aspect("equal")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_title("Velocity is a vector: components and speed")
    ax.grid(False)
    save_svg(fig, "particle_velocity.svg")


def kinetic_energy_curves():
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    speed = np.linspace(0, 0.5, 400)
    for mass, color in [(1.0, ORANGE), (12.0, BLUE), (39.948, GREEN)]:
        energy = 0.5 * mass * speed**2
        ax.plot(speed, energy, color=color, lw=2.3, label=f"m = {mass:g} Da")
    ax.plot([0.1], [0.06], "o", color=BLUE, ms=8)
    ax.annotate(
        "v = 0.1 nm/ps, m = 12 Da\nK = 0.06 kJ/mol",
        xy=(0.1, 0.06),
        xytext=(0.20, 0.22),
        arrowprops=dict(arrowstyle="->", color=GRAY),
        color=GRAY,
    )
    ax.set_xlabel("speed |v| (nm/ps)")
    ax.set_ylabel("kinetic energy K (kJ/mol)")
    ax.set_title("Kinetic energy grows with mass and with speed squared")
    ax.legend()
    save_svg(fig, "kinetic_energy_curves.svg")


def maxwell_speed():
    fig, ax = plt.subplots(figsize=(8.2, 4.6))
    R = 0.008314462618
    mass = 12.0
    speed = np.linspace(0, 1.6, 500)
    for temperature, color in [(150, ORANGE), (300, BLUE), (600, GREEN)]:
        exponent = -mass * speed**2 / (2.0 * R * temperature)
        density = speed**2 * np.exp(exponent)
        density /= np.trapezoid(density, speed)
        ax.plot(speed, density, lw=2.3, color=color, label=f"{temperature} K")
        mean = np.trapezoid(speed * density, speed)
        ax.axvline(mean, color=color, ls=":", lw=1.2, alpha=0.8)
    ax.axvline(0.1, color=GRAY, ls="--", lw=1.6, label="example velocity 0.1 nm/ps")
    ax.set_xlabel("speed |v| (nm/ps)")
    ax.set_ylabel("probability density")
    ax.set_title("Thermal speeds for a 12 Da particle, with the example marked")
    ax.legend(title="temperature / example", loc="upper right")
    save_svg(fig, "maxwell_speed.svg")


def particle_motion_gif():
    mass = 12.0
    velocity = np.array([0.08, 0.05])
    times = np.linspace(0, 2.0, 100)

    fig, ax = plt.subplots(figsize=(6.0, 4.4))
    ax.set_xlim(0, 0.2)
    ax.set_ylim(0, 0.16)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_title("Constant-velocity motion: no force, no change in speed")
    ax.set_aspect("equal")
    point, = ax.plot([], [], "o", color=BLUE, ms=14)
    trail, = ax.plot([], [], color=BLUE, lw=1.4, alpha=0.5)
    velocity_arrow = ax.quiver([0], [0], [0], [0], angles="xy", scale_units="xy", scale=1, color=ORANGE, width=0.012)
    energy_text = ax.text(0.01, 0.145, "", fontsize=10)
    speed = np.linalg.norm(velocity)
    kinetic = 0.5 * mass * speed**2

    def update(index):
        position = velocity * times[index]
        point.set_data([position[0]], [position[1]])
        trail.set_data(velocity[0] * times[: index + 1], velocity[1] * times[: index + 1])
        velocity_arrow.set_offsets([position])
        velocity_arrow.set_UVC(velocity[0], velocity[1])
        energy_text.set_text(f"time = {times[index]:.2f} ps\nK = {kinetic:.3f} kJ/mol")
        return point, trail, velocity_arrow, energy_text

    animation = FuncAnimation(fig, update, frames=len(times), interval=50, blit=False)
    save_gif(animation, "particle_motion.gif", fps=25)


def harmonic_energy_force():
    fig, (ax_energy, ax_force) = plt.subplots(1, 2, figsize=(10.4, 4.4), sharex=True)
    k = 100.0
    x = np.linspace(-0.2, 0.2, 400)
    energy = 0.5 * k * x**2
    force = -k * x
    ax_energy.plot(x, energy, color=BLUE, lw=2.4)
    ax_energy.axvline(0, color=GRAY, lw=0.8)
    ax_energy.plot([0.1], [0.5], "o", color=BLUE, ms=8)
    ax_energy.annotate("x = 0.1 nm\nU = 0.5 kJ/mol", xy=(0.1, 0.5), xytext=(0.06, 1.45), arrowprops=dict(arrowstyle="->", color=GRAY))
    ax_energy.set_ylabel("potential energy U (kJ/mol)")
    ax_energy.set_title("Energy is a bowl")
    ax_energy.set_xlabel("x (nm)")
    ax_force.plot(x, force, color=RED, lw=2.4)
    ax_force.axhline(0, color=GRAY, lw=0.8)
    ax_force.axvline(0, color=GRAY, lw=0.8)
    ax_force.plot([0.1], [-10], "o", color=RED, ms=8)
    ax_force.annotate("x = 0.1 nm\nF = -10 kJ/(mol nm)", xy=(0.1, -10), xytext=(-0.17, -8), arrowprops=dict(arrowstyle="->", color=GRAY))
    ax_force.set_ylabel("force F (kJ/(mol nm))")
    ax_force.set_title("Force is the downhill slope")
    ax_force.set_xlabel("x (nm)")
    fig.suptitle("A harmonic potential and its force are two views of the same model", y=1.02)
    save_svg(fig, "harmonic_energy_force.svg")


def harmonic_quiver():
    fig, ax = plt.subplots(figsize=(6.2, 5.2))
    k = 100.0
    x = np.linspace(-0.2, 0.2, 240)
    y = np.linspace(-0.2, 0.2, 240)
    X, Y = np.meshgrid(x, y)
    U = 0.5 * k * (X**2 + Y**2)
    contour = ax.contourf(X, Y, U, levels=18, cmap="Blues", alpha=0.85)
    fig.colorbar(contour, ax=ax, label="U (kJ/mol)")
    qx = np.linspace(-0.18, 0.18, 13)
    qy = np.linspace(-0.18, 0.18, 13)
    QX, QY = np.meshgrid(qx, qy)
    FX, FY = -k * QX, -k * QY
    magnitude = np.hypot(FX, FY)
    ax.quiver(QX, QY, FX, FY, magnitude, cmap="Oranges", pivot="middle", scale=3800, width=0.004)
    ax.plot([0], [0], "o", color=RED, ms=10)
    ax.text(0.005, 0.012, "equilibrium", color=RED, fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("y (nm)")
    ax.set_title("Force vectors point toward the energy minimum")
    ax.grid(False)
    save_svg(fig, "harmonic_quiver.svg")


def force_derivative():
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    k = 100.0
    x = np.linspace(0.0, 0.2, 400)
    u = 0.5 * k * x**2
    x0 = 0.1
    h = 0.025
    u0 = 0.5 * k * x0**2
    slope = -k * x0
    tangent_x = np.array([x0 - 0.055, x0 + 0.055])
    tangent_y = u0 + slope * (tangent_x - x0)
    secant_x = np.array([x0 - h, x0 + h])
    secant_y = 0.5 * k * secant_x**2
    ax.plot(x, u, color=BLUE, lw=2.4, label="U(x)")
    ax.plot(secant_x, secant_y, "o", color=ORANGE, ms=7, label="U(x-h), U(x+h)")
    ax.plot(secant_x, secant_y, color=ORANGE, lw=1.5, ls="--", label="secant")
    ax.plot(tangent_x, tangent_y, color=RED, lw=2.0, label="tangent at x")
    ax.plot([x0], [u0], "o", color=RED, ms=8)
    ax.annotate(
        "slope = dU/dx = -F\nF = -slope",
        xy=(x0, u0),
        xytext=(0.125, 0.35),
        arrowprops=dict(arrowstyle="->", color=GRAY),
        color=GRAY,
    )
    ax.set_xlim(0, 0.2)
    ax.set_ylim(0, 1.25)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("U (kJ/mol)")
    ax.set_title("A finite-difference test verifies force = -dU/dx")
    ax.legend(loc="upper left")
    save_svg(fig, "force_derivative.svg")


def oscillator_data():
    mass = 12.0
    stiffness = 12.0 * np.pi**2
    particle = Particle(mass=mass, position=[0.1, 0.0, 0.0], velocity=[0.0, 0.0, 0.0])
    potential = HarmonicPotential(stiffness=stiffness, center=[0.0, 0.0, 0.0])
    trajectory = simulate_velocity_verlet(particle, potential, dt=0.001, steps=2000)
    return particle, potential, trajectory


def oscillator_position_velocity(trajectory):
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    ax.plot(trajectory.time, trajectory.position[:, 0], color=BLUE, lw=2.2, label="x(t)")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("position x (nm)", color=BLUE)
    ax.tick_params(axis="y", labelcolor=BLUE)
    ax2 = ax.twinx()
    ax2.plot(trajectory.time, trajectory.velocity[:, 0], color=ORANGE, lw=1.8, ls="--", label="v(t)")
    ax2.set_ylabel("velocity v (nm/ps)", color=ORANGE)
    ax2.tick_params(axis="y", labelcolor=ORANGE)
    ax2.grid(False)
    ax.set_title("Harmonic oscillator: position and velocity over one period")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], loc="upper right")
    save_svg(fig, "oscillator_position_velocity.svg")


def oscillator_energy(trajectory):
    fig, ax = plt.subplots(figsize=(8.6, 4.8))
    ax.plot(trajectory.time, trajectory.kinetic_energy, color=ORANGE, lw=2.0, label="kinetic K")
    ax.plot(trajectory.time, trajectory.potential_energy, color=BLUE, lw=2.0, label="potential U")
    ax.plot(trajectory.time, trajectory.total_energy, color=GREEN, lw=2.6, label="total E")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("energy (kJ/mol)")
    ax.set_title("Energy sloshes between K and U while total E stays nearly constant")
    ax.legend(loc="center right")
    save_svg(fig, "oscillator_energy.svg")


def energy_scale():
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    labels = [
        "thermal RT at 300 K",
        "typical LJ well",
        "hydrogen bond",
        "point-charge Coulomb\n(±e at 0.3 nm)",
        "covalent bond",
    ]
    values = [2.5, 1.0, 20.0, 463.0, 400.0]
    colors = [GREEN, BLUE, ORANGE, RED, GRAY]
    y = np.arange(len(labels))
    ax.barh(y, values, color=colors, alpha=0.85)
    ax.set_xscale("log")
    ax.set_xlim(0.5, 1200)
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.set_xlabel("energy scale (kJ/mol, logarithmic)")
    ax.set_title("Not all energies are on the same scale")
    for yi, value in zip(y, values):
        ax.text(value * 1.1, yi, f"{value:g}", va="center", fontsize=9)
    ax.grid(axis="x", which="both", alpha=0.2)
    save_svg(fig, "energy_scale.svg")


def work_energy():
    fig, (ax_force, ax_energy) = plt.subplots(1, 2, figsize=(10.4, 4.4), sharex=True)
    k = 100.0
    x = np.linspace(0.0, 0.2, 400)
    force = -k * x
    energy = 0.5 * k * x**2
    ax_force.plot(x, force, color=RED, lw=2.3)
    ax_force.fill_between(x, force, 0.0, where=(x <= 0.1), color=RED, alpha=0.18)
    ax_force.axhline(0, color=GRAY, lw=0.8)
    ax_force.set_xlabel("x (nm)")
    ax_force.set_ylabel("F (kJ/(mol nm))")
    ax_force.set_title("Force does work along the path")
    ax_force.annotate("work done\n= shaded area", xy=(0.05, -5), xytext=(0.10, -13), arrowprops=dict(arrowstyle="->", color=GRAY))
    ax_energy.plot(x, energy, color=BLUE, lw=2.3)
    ax_energy.fill_between(x, energy, 0.0, where=(x <= 0.1), color=BLUE, alpha=0.18)
    ax_energy.set_xlabel("x (nm)")
    ax_energy.set_ylabel("U (kJ/mol)")
    ax_energy.set_title("The work becomes potential energy")
    save_svg(fig, "work_energy.svg")


def harmonic_curvature():
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    x = np.linspace(-0.2, 0.2, 400)
    for stiffness, color in [(50.0, ORANGE), (100.0, BLUE), (200.0, GREEN)]:
        ax.plot(x, 0.5 * stiffness * x**2, lw=2.2, color=color, label=f"k = {stiffness:g}")
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("U (kJ/mol)")
    ax.set_title("Stiffer springs have steeper walls")
    ax.legend()
    save_svg(fig, "harmonic_curvature.svg")


def finite_difference_error():
    k = 100.0
    x0 = 0.1
    exact = -k * x0
    h = np.logspace(-8, -1, 90)
    error = np.empty_like(h)
    for index, step in enumerate(h):
        numerical = -((0.5 * k * (x0 + step) ** 2) - (0.5 * k * (x0 - step) ** 2)) / (2.0 * step)
        error[index] = abs(numerical - exact)
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    ax.loglog(h, error, color=BLUE, lw=2.0)
    ax.set_xlabel("finite-difference step h (nm)")
    ax.set_ylabel("absolute force error (kJ/(mol nm))")
    ax.set_title("Finite differences: truncation vs floating-point error")
    ax.grid(True, which="both", alpha=0.25)
    save_svg(fig, "finite_difference_error.svg")


def dt_stability():
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    particle = Particle(mass=12.0, position=[0.1, 0.0, 0.0], velocity=[0.0, 0.0, 0.0])
    potential = HarmonicPotential(stiffness=12.0 * np.pi**2, center=[0.0, 0.0, 0.0])
    exact_time = np.linspace(0.0, 2.0, 400)
    ax.plot(exact_time, 0.1 * np.cos(np.pi * exact_time), color=GRAY, ls="--", lw=1.2, label="exact")
    for dt, color in [(0.001, GREEN), (0.01, BLUE), (0.02, ORANGE), (0.05, RED)]:
        steps = int(round(2.0 / dt))
        trajectory = simulate_velocity_verlet(particle, potential, dt=dt, steps=steps)
        ax.plot(trajectory.time, trajectory.position[:, 0], color=color, lw=1.6, label=f"dt = {dt:g} ps")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("x (nm)")
    ax.set_title("The same oscillator at several timesteps")
    ax.legend(ncol=2)
    save_svg(fig, "dt_stability.svg")


def oscillator_phase(trajectory):
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ax.plot(trajectory.position[:, 0], trajectory.velocity[:, 0], color=BLUE, lw=2.0)
    ax.plot(trajectory.position[0, 0], trajectory.velocity[0, 0], "o", color=ORANGE, ms=9, label="start")
    ax.axhline(0, color=GRAY, lw=0.8)
    ax.axvline(0, color=GRAY, lw=0.8)
    ax.set_xlabel("position x (nm)")
    ax.set_ylabel("velocity v (nm/ps)")
    ax.set_title("Phase-space orbit of the oscillator")
    ax.set_aspect("equal")
    ax.legend()
    save_svg(fig, "oscillator_phase.svg")


def oscillator_gif(trajectory):
    fig, (ax_spring, ax_trace) = plt.subplots(1, 2, figsize=(9.2, 4.4), gridspec_kw={"width_ratios": [1.4, 1.0]})
    ax_spring.set_xlim(-0.22, 0.22)
    ax_spring.set_ylim(-0.7, 0.7)
    ax_spring.set_aspect("equal")
    ax_spring.axis("off")
    ax_spring.set_title("Mass on a spring")
    ax_spring.plot([-0.21, -0.21], [-0.35, 0.35], color=GRAY, lw=3)
    ax_spring.plot([-0.2, 0.0], [0, 0], color=GRAY, lw=1)
    mass_block = Rectangle((-0.045, -0.16), 0.09, 0.32, fc=BLUE, ec="none", alpha=0.9)
    ax_spring.add_patch(mass_block)
    spring_x = np.linspace(-0.2, 0.0, 20)
    spring_y = 0.16 * np.sin(np.linspace(0, 18 * np.pi, len(spring_x)))
    spring, = ax_spring.plot(spring_x, spring_y, color=ORANGE, lw=1.4)
    force_arrow = ax_spring.quiver([0], [0], [0], [0], angles="xy", scale_units="xy", scale=1, color=RED, width=0.014)
    ax_spring.text(-0.2, -0.55, "equilibrium at x = 0", color=GRAY, fontsize=9)

    ax_trace.set_xlim(0, 2.0)
    ax_trace.set_ylim(-0.12, 0.12)
    ax_trace.set_xlabel("time (ps)")
    ax_trace.set_ylabel("x (nm)")
    ax_trace.set_title("Position trace")
    trace, = ax_trace.plot([], [], color=BLUE, lw=1.8)

    frames = 240
    indices = np.linspace(0, len(trajectory.time) - 1, frames).astype(int)

    def update(frame):
        index = indices[frame]
        x = float(trajectory.position[index, 0])
        v = float(trajectory.velocity[index, 0])
        mass_block.set_x(x - 0.045)
        spring.set_xdata(np.linspace(-0.2, x - 0.045, 20))
        spring.set_ydata(0.16 * np.sin(np.linspace(0, 18 * np.pi, 20)))
        force_arrow.set_offsets([[x, 0.0]])
        force_arrow.set_UVC(-100.0 * x * 0.001, 0.0)
        trace.set_data(trajectory.time[: index + 1], trajectory.position[: index + 1, 0])
        return mass_block, spring, force_arrow, trace

    animation = FuncAnimation(fig, update, frames=frames, interval=40, blit=False)
    save_gif(animation, "oscillator.gif", fps=25)


def energy_drift():
    particle = Particle(mass=12.0, position=[0.1, 0.0, 0.0], velocity=[0.0, 0.0, 0.0])
    potential = HarmonicPotential(stiffness=12.0 * np.pi**2, center=[0.0, 0.0, 0.0])
    dts = np.array([0.01, 0.005, 0.002, 0.001, 0.0005, 0.0002])
    drifts = []
    for dt in dts:
        steps = int(round(2.0 / dt))
        trajectory = simulate_velocity_verlet(particle, potential, dt=dt, steps=steps)
        initial = trajectory.total_energy[0]
        drifts.append(np.max(np.abs(trajectory.total_energy - initial)) / abs(initial))
    fig, ax = plt.subplots(figsize=(7.4, 4.8))
    ax.loglog(dts, drifts, "o-", color=BLUE, lw=2.0, ms=7)
    ax.set_xlabel("time step Δt (ps)")
    ax.set_ylabel("maximum relative energy drift")
    ax.set_title("Velocity Verlet energy error scales roughly as Δt²")
    ax.invert_xaxis()
    ax.grid(True, which="both", alpha=0.25)
    save_svg(fig, "energy_drift.svg")


def euler_vs_verlet():
    mass = 12.0
    stiffness = 12.0 * np.pi**2
    dt = 0.005
    steps = 400
    times = np.arange(steps + 1) * dt
    x = 0.1
    v = 0.0
    euler_energy = [0.5 * mass * v**2 + 0.5 * stiffness * x**2]
    for _ in range(steps):
        acceleration = -stiffness * x / mass
        x_new = x + v * dt
        v_new = v + acceleration * dt
        x, v = x_new, v_new
        euler_energy.append(0.5 * mass * v**2 + 0.5 * stiffness * x**2)
    particle = Particle(mass=mass, position=[0.1, 0.0, 0.0], velocity=[0.0, 0.0, 0.0])
    potential = HarmonicPotential(stiffness=stiffness, center=[0.0, 0.0, 0.0])
    verlet = simulate_velocity_verlet(particle, potential, dt=dt, steps=steps)
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    ax.plot(times, np.array(euler_energy) / euler_energy[0], color=RED, lw=2.0, label="explicit Euler")
    ax.plot(verlet.time, verlet.total_energy / verlet.total_energy[0], color=GREEN, lw=2.2, label="velocity Verlet")
    ax.axhline(1.0, color=GRAY, lw=0.8, ls=":")
    ax.set_xlabel("time (ps)")
    ax.set_ylabel("total energy / initial energy")
    ax.set_title("Euler drifts; Velocity Verlet stays bounded (Δt = 0.005 ps)")
    ax.legend()
    save_svg(fig, "euler_vs_verlet.svg")


def main():
    learning_ladder()
    stage_flow()
    unit_scale()
    unit_bridge()
    particle_state()
    particle_velocity()
    kinetic_energy_curves()
    maxwell_speed()
    particle_motion_gif()
    harmonic_energy_force()
    harmonic_quiver()
    force_derivative()
    energy_scale()
    work_energy()
    harmonic_curvature()
    finite_difference_error()
    particle, potential, trajectory = oscillator_data()
    oscillator_position_velocity(trajectory)
    oscillator_energy(trajectory)
    oscillator_phase(trajectory)
    oscillator_gif(trajectory)
    energy_drift()
    euler_vs_verlet()
    dt_stability()
    verlet_flow()
    strip_svg_trailing_whitespace()
    print(f"Wrote visualizations to {ASSETS}")


if __name__ == "__main__":
    main()
