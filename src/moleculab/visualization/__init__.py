# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Reproducible visual outputs derived from validated Moleculab records.

The visualization layer is a consumer: it does not assign protonation states,
infer interactions, or modify scientific geometry.  Each figure is accompanied
by metadata describing the source structure, tool versions, camera view, and
scientific limitations.
"""

from .chemical_states import draw_chemical_state_comparison
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

__all__ = [
    "FIGURE_METADATA_SCHEMA_VERSION",
    "FigureContext",
    "FigureMetadata",
    "build_hydrogen_water_script",
    "build_overview_script",
    "build_pocket_script",
    "draw_chemical_state_comparison",
    "draw_ligand_2d",
    "load_figure_metadata",
    "plot_evidence_matrix",
    "pymol_available",
    "render_binding_pocket_closeup",
    "render_hydrogen_water_evidence",
    "render_protein_ligand_overview",
    "rendering_config_sha256",
    "validate_figure_metadata",
    "write_figure_metadata",
]
