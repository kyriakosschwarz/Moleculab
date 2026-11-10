# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Reproducible visual outputs derived from validated Moleculab records.

The visualization layer is a consumer: it does not assign protonation states,
infer interactions, or modify scientific geometry.  Each figure is accompanied
by metadata describing the source structure, tool versions, camera view, and
scientific limitations.
"""

from .chemical_states import draw_chemical_state_comparison
from .free_energy import plot_free_energy_convergence
from .hydrogen_water import (
    build_hydrogen_water_script,
    render_hydrogen_water_evidence,
)
from .interaction_evidence import plot_evidence_matrix
from .metadata import (
    FIGURE_METADATA_SCHEMA_VERSION,
    FigureContext,
    FigureMetadata,
    load_figure_metadata,
    rendering_config_sha256,
    validate_figure_metadata,
    write_figure_metadata,
)
from .molecular_2d import draw_ligand_2d
from .protein_ligand import (
    build_overview_script,
    build_pocket_script,
    pymol_available,
    render_binding_pocket_closeup,
    render_protein_ligand_overview,
)
from .registry import (
    FIGURE_REGISTRY,
    GALLERY_RELATIVE_PATH,
    FigureSpec,
    figure_spec,
    render_gallery_markdown,
    validate_figure_registry,
    write_gallery,
)
from .robustness import plot_interaction_sensitivity
from .styles import (
    EVIDENCE_STYLES,
    REQUIREMENT_STYLES,
    ROBUSTNESS_STYLES,
    TERM_STYLES,
    StatusStyle,
    clean_svg_file,
    configure_matplotlib,
    contrast_ratio,
    evidence_style,
    requirement_style,
    robustness_style,
    save_figure,
    style_axes,
    term_style,
)
from .thermodynamics import plot_structural_thermodynamics_boundary

__all__ = [
    "FIGURE_METADATA_SCHEMA_VERSION",
    "EVIDENCE_STYLES",
    "FIGURE_REGISTRY",
    "GALLERY_RELATIVE_PATH",
    "FigureContext",
    "FigureMetadata",
    "FigureSpec",
    "REQUIREMENT_STYLES",
    "ROBUSTNESS_STYLES",
    "StatusStyle",
    "TERM_STYLES",
    "build_hydrogen_water_script",
    "build_overview_script",
    "build_pocket_script",
    "clean_svg_file",
    "configure_matplotlib",
    "contrast_ratio",
    "draw_chemical_state_comparison",
    "draw_ligand_2d",
    "evidence_style",
    "figure_spec",
    "load_figure_metadata",
    "plot_evidence_matrix",
    "plot_free_energy_convergence",
    "plot_interaction_sensitivity",
    "plot_structural_thermodynamics_boundary",
    "pymol_available",
    "render_binding_pocket_closeup",
    "render_hydrogen_water_evidence",
    "render_protein_ligand_overview",
    "rendering_config_sha256",
    "render_gallery_markdown",
    "requirement_style",
    "robustness_style",
    "save_figure",
    "style_axes",
    "term_style",
    "validate_figure_metadata",
    "validate_figure_registry",
    "write_figure_metadata",
    "write_gallery",
]
