"""Exact-versus-estimated harmonic free-energy convergence figures."""

from __future__ import annotations

from pathlib import Path

from ..applications.free_energy_study import HarmonicFEPStudyRecord
from .metadata import FigureContext, FigureMetadata
from .styles import configure_matplotlib, save_figure, style_axes

_CASE_COLOURS = {
    0: "#2b6cb0",
    1: "#dd6b20",
    2: "#2f855a",
}


def plot_free_energy_convergence(
    *,
    study: HarmonicFEPStudyRecord,
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
    height: float = 6.5,
    dpi: int = 200,
) -> FigureMetadata:
    """Plot exact references, sample-size convergence, and overlap diagnostics."""

    import matplotlib

    configure_matplotlib()
    import matplotlib.pyplot as plt

    figure, (energy_axis, overlap_axis) = plt.subplots(
        1,
        2,
        figsize=(width, height),
        constrained_layout=False,
    )
    figure.subplots_adjust(
        top=0.86,
        bottom=0.24,
        left=0.08,
        right=0.98,
        wspace=0.28,
    )

    for index, case in enumerate(study.cases):
        colour = _CASE_COLOURS.get(index, "#4a5568")
        samples = [point.n_samples for point in case.convergence]
        estimates = [point.delta_g_kj_per_mol for point in case.convergence]
        errors = [point.standard_error_kj_per_mol for point in case.convergence]
        energy_axis.errorbar(
            samples,
            estimates,
            yerr=errors,
            marker="o",
            capsize=4,
            linewidth=1.8,
            color=colour,
            label=case.label,
        )
        energy_axis.axhline(
            case.exact_delta_g_kj_per_mol,
            color=colour,
            linestyle="--",
            linewidth=1.2,
            alpha=0.8,
        )
        fractions = [
            point.overlap.effective_sample_fraction
            for point in case.convergence
        ]
        overlap_axis.plot(
            samples,
            fractions,
            marker="o",
            linewidth=1.8,
            color=colour,
            label=case.label,
        )
        replicate = case.replicates
        short_label = "mild" if index == 0 else "poor"
        figure.text(
            0.01,
            0.115 - 0.045 * index,
            f"{short_label}: replicate mean {replicate.mean_delta_g_kj_per_mol:+.3f} "
            f"+/- {replicate.standard_deviation_kj_per_mol:.3f} kJ/mol; "
            f"coverage {replicate.confidence_interval_coverage:.2f}",
            fontsize=8,
            color=colour,
            ha="left",
        )

    energy_axis.set_xscale("log")
    energy_axis.set_xlabel("samples")
    energy_axis.set_ylabel("ΔG (kJ/mol)")
    energy_axis.set_title("Exact reference and FEP estimate")
    energy_axis.legend(frameon=False, fontsize=8, loc="best")
    energy_axis.text(
        0.02,
        0.98,
        "dashed: exact harmonic ΔG",
        transform=energy_axis.transAxes,
        va="top",
        fontsize=8,
        color="#4a5568",
    )
    style_axes(energy_axis)

    overlap_axis.set_xscale("log")
    overlap_axis.set_yscale("log")
    overlap_axis.set_xlabel("samples")
    overlap_axis.set_ylabel("effective sample fraction")
    overlap_axis.set_title("Overlap diagnostic")
    overlap_axis.legend(frameon=False, fontsize=8, loc="best")
    style_axes(overlap_axis)

    figure.suptitle(title, fontsize=15, weight="bold")
    figure.text(
        0.5,
        0.012,
        "Controlled 3D harmonic model; not a protein-ligand or imatinib affinity calculation.",
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
            ("x_scale", "log samples"),
            ("left_panel", "exact reference and FEP error bars"),
            ("right_panel", "effective_sample_fraction"),
            ("replicate_summary", "mean, standard_deviation, CI_coverage"),
        ),
        rendering_config_sha256=None,
        scientific_limitations=(
            "The model is a three-dimensional harmonic oscillator, not a protein-ligand system.",
            "Only forward FEP is shown; no bidirectional estimator or MBAR overlap correction is applied.",
            "Finite-sample bias under poor overlap can remain even when the standard error is small.",
            "Effective sample size is an observed-weight diagnostic, not proof of complete phase-space coverage.",
            "No imatinib binding free energy is computed.",
        ),
        alt_text=alt_text,
        caption=caption,
        camera_view=None,
    )
