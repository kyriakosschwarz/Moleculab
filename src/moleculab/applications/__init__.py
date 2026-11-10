# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Applications grouped by scientific domain.

Applications answer scientific questions using the lower layers. They may
combine models, algorithms, and protocols, but they do not redefine the core
physics. The subpackages are:

``structure``
    experimental-structure manifests, ligand identity and coordinate records,
    and binding-site reconciliation;
``chemistry``
    acid-base bookkeeping, chemical-state preparation and enumeration, and
    protein-pocket template chemistry;
``interactions``
    evidence-gated interaction interpretation, three-state fingerprints,
    hydrogen and water evidence, state-aware assessment, and robustness records;
``thermodynamics``
    FEP/TI/MBAR, harmonic and soft-core alchemy studies, and explicit
    structural-to-thermodynamic boundaries;
``docking``
    rigid and flexible receptor-ligand docking;
``solvents``
    coarse-grained solvent and three-site water builders;
``qmmm``
    the finite-difference QM/MM model;
``audit``
    the final scientific audit and integration records.

None of the record layers parses coordinates, chooses a structure, or invents
an atom mapping on the caller's behalf.
"""
