# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Prepare explicit chemical states for the selected protein pocket.

This is Batch 10 of the capstone.  It records residue-template formal charges,
donor/acceptor/hydrophobic/aromatic features, and explicit protonation
hypotheses for the pocket atoms selected by Batch 4.  Experimental coordinates
are preserved; no protein hydrogens are generated; water and hetero roles remain
unassigned; no state is selected.

Run with gemmi available:

    uv run --with gemmi python examples/drug_discovery/capstone_protein_preparation.py
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import gemmi
import numpy as np

from moleculab import (
    PreparedProteinAtom,
    PreparedProteinRing,
    PreparedProteinState,
    PreparedWaterRole,
    ProteinHeteroRole,
    ProteinResidueState,
    build_protein_preparation_record,
    residue_atom_features,
    load_binding_site,
    load_manifest,
    load_protein_preparation,
    sha256_file,
    validate_protein_preparation,
    write_protein_preparation,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_SITE = Path("data/structures/2HYY.binding-site.json")
DEFAULT_SOURCE = Path("data/structures/2HYY.cif")
DEFAULT_OUTPUT = Path("data/structures/2HYY.protein-preparation.json")

_AROMATIC_RINGS = {
    "PHE": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TYR": ("CG", "CD1", "CD2", "CE1", "CE2", "CZ"),
    "TRP": ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"),
    "HIS": ("CG", "ND1", "CD2", "CE1", "NE2"),
}


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--site", type=Path, default=DEFAULT_SITE)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def _unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        return value[1:-1]
    return value


def _none_if_missing(value: str) -> str | None:
    return None if value in ("?", ".") else value


def _atom_rows(block: gemmi.cif.Block):
    return [
        [_unquote(value) for value in row]
        for row in block.find(
            "_atom_site.",
            [
                "group_PDB",
                "label_comp_id",
                "label_asym_id",
                "auth_asym_id",
                "label_seq_id",
                "auth_seq_id",
                "pdbx_PDB_ins_code",
                "label_atom_id",
                "type_symbol",
                "Cartn_x",
                "Cartn_y",
                "Cartn_z",
                "label_alt_id",
                "pdbx_PDB_model_num",
            ],
        )
    ]


def _residue_key(site_residue) -> tuple[str, ...]:
    return (
        site_residue.label_asym_id,
        site_residue.auth_asym_id,
        site_residue.label_seq_id or "",
        site_residue.auth_seq_id,
        site_residue.insertion_code or "",
        site_residue.residue_name,
    )


def _pocket_atoms(block: gemmi.cif.Block, site, selection) -> list[list[str]]:
    keys = {_residue_key(residue) for residue in site.residues}
    atoms = []
    for values in _atom_rows(block):
        if values[0] != "ATOM":
            continue
        if values[3] != selection.auth_asym_id:
            continue
        if values[13] not in ("?", ".") and int(values[13]) != 1:
            continue
        key = (
            values[2],
            values[3],
            _none_if_missing(values[4]) or "",
            values[5],
            _none_if_missing(values[6]) or "",
            values[1],
        )
        if key in keys:
            atoms.append(values)
    return atoms


def _ring_geometry(atoms: list[list[str]], ring_names: tuple[str, ...]):
    positions = {
        values[7]: (float(values[9]), float(values[10]), float(values[11]))
        for values in atoms
        if values[7] in ring_names
    }
    if len(positions) != len(ring_names):
        return None
    points = np.asarray([positions[name] for name in ring_names], dtype=float)
    centroid = points.mean(axis=0)
    _, _, vectors = np.linalg.svd(points - centroid, full_matrices=False)
    return tuple(centroid * 0.1), tuple(vectors[-1])


def _residue_protonation(residue_name: str, his_state: str) -> str:
    if residue_name == "ASP":
        return "deprotonated"
    if residue_name == "GLU":
        return "deprotonated"
    if residue_name == "LYS":
        return "protonated"
    if residue_name == "ARG":
        return "protonated"
    if residue_name == "HIS":
        return "protonated HIP" if his_state == "HIP" else "neutral HIE"
    return "standard template"


def _state(site, atom_rows: list[list[str]], selection, state_id: str, his_state: str) -> PreparedProteinState:
    rings = []
    ring_ids_by_atom: dict[str, list[str]] = {}
    for residue_name, ring_names in _AROMATIC_RINGS.items():
        for values in atom_rows:
            if values[1] != residue_name:
                continue
            residue_id = f"{values[3]}:{values[5]}"
            ring_id = f"{residue_id}:{residue_name}:ring"
            if ring_id in {ring.ring_id for ring in rings}:
                continue
            geometry = _ring_geometry(
                [row for row in atom_rows if row[1] == residue_name and row[5] == values[5]],
                ring_names,
            )
            if geometry is None:
                continue
            centroid, normal = geometry
            rings.append(
                PreparedProteinRing(
                    ring_id=ring_id,
                    residue_id=residue_id,
                    atom_ids=tuple(
                        f"{values[3]}:{values[5]}:{name}" for name in ring_names
                    ),
                    centroid_nm=centroid,
                    normal=normal,
                    source="residue-template ring definition; SVD plane from experimental coordinates",
                )
            )
            for name in ring_names:
                ring_ids_by_atom.setdefault(f"{values[3]}:{values[5]}:{name}", []).append(ring_id)

    atoms = []
    for values in atom_rows:
        residue_name = values[1]
        atom_name = values[7]
        charge, donor, acceptor, hydrophobic = residue_atom_features(
            residue_name, atom_name, his_state
        )
        atom_id = f"{values[3]}:{values[5]}:{atom_name}"
        atoms.append(
            PreparedProteinAtom(
                atom_id=atom_id,
                atom_name=atom_name,
                element=values[8],
                residue_name=residue_name,
                label_asym_id=values[2],
                auth_asym_id=values[3],
                label_seq_id=_none_if_missing(values[4]),
                auth_seq_id=values[5],
                insertion_code=_none_if_missing(values[6]),
                alternate_location=_none_if_missing(values[12]),
                position_angstrom=(float(values[9]), float(values[10]), float(values[11])),
                formal_charge=charge,
                is_donor=donor,
                is_acceptor=acceptor,
                is_hydrophobic=hydrophobic,
                aromatic_ring_ids=tuple(ring_ids_by_atom.get(atom_id, ())),
                hydrogen_count=0,
                hydrogen_positions_nm=(),
                coordinate_source="deposited",
                hydrogen_source=None,
                feature_source="explicit residue/atom template table",
            )
        )

    residue_states = []
    for residue in site.residues:
        residue_atoms = [
            values for values in atom_rows
            if (values[2], values[3], _none_if_missing(values[4]) or "", values[5], _none_if_missing(values[6]) or "", values[1])
            == _residue_key(residue)
        ]
        total_charge = sum(
            residue_atom_features(values[1], values[7], his_state)[0]
            for values in residue_atoms
        )
        residue_states.append(
            ProteinResidueState(
                residue_id=f"{residue.auth_asym_id}:{residue.auth_seq_id}",
                residue_name=residue.residue_name,
                label_asym_id=residue.label_asym_id,
                auth_asym_id=residue.auth_asym_id,
                label_seq_id=residue.label_seq_id,
                auth_seq_id=residue.auth_seq_id,
                insertion_code=residue.insertion_code,
                protonation_state=_residue_protonation(residue.residue_name, his_state),
                formal_charge=total_charge,
                source="explicit residue template, not experimentally validated",
            )
        )

    water_roles = tuple(
        PreparedWaterRole(
            water_id=f"{water.auth_asym_id}:{water.auth_seq_id}:{water.atom_name}",
            role="unassigned",
            source="experimental water; no role assigned",
        )
        for water in site.waters_near_ligand
    )
    hetero_roles = tuple(
        ProteinHeteroRole(
            hetero_id=f"{residue.auth_asym_id}:{residue.auth_seq_id}:{residue.residue_name}",
            residue_name=residue.residue_name,
            role="unknown",
            source="preserved hetero residue; no ion/cofactor role assigned",
        )
        for residue in site.hetero_residues
    )
    return PreparedProteinState(
        state_id=state_id,
        label=(
            "Standard residue templates" if his_state != "HIP"
            else "Standard templates with protonated histidine"
        ),
        protonation_model=(
            "standard residue templates; HIS neutral HIE" if his_state != "HIP"
            else "standard residue templates; HIS protonated HIP"
        ),
        hydrogen_method=None,
        assumptions=(
            "Formal charges are localized on template atoms for educational eligibility checks.",
            "Donor/acceptor and hydrophobic assignments are template hypotheses, not experimental observations.",
            "No protein hydrogens were generated; angular hydrogen-bond geometry is unavailable.",
            "Water and hetero roles remain unassigned.",
            "No pH-dependent population model was supplied.",
        ),
        atoms=tuple(atoms),
        residue_states=tuple(residue_states),
        rings=tuple(rings),
        water_roles=water_roles,
        hetero_roles=hetero_roles,
    )


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (args.manifest, args.site, args.source):
        if not path.exists():
            raise SystemExit(f"missing {path}; run the earlier capstone batches first.")
    manifest = load_manifest(args.manifest)
    site = load_binding_site(args.site)
    if site.accession != manifest.source.accession:
        raise SystemExit("manifest and binding-site records disagree on accession.")

    block = gemmi.cif.read_file(str(args.source)).sole_block()
    selection = manifest.selection.protein
    atom_rows = _pocket_atoms(block, site, selection)
    if not atom_rows:
        raise SystemExit("no pocket atoms matched the binding-site residue set.")

    states = (
        _state(site, atom_rows, selection, "standard-templates", "HIE"),
        _state(site, atom_rows, selection, "histidine-protonated", "HIP"),
    )
    record = build_protein_preparation_record(
        source_manifest_path=str(args.manifest),
        source_manifest_sha256=sha256_file(args.manifest),
        binding_site_record_path=str(args.site),
        binding_site_record_sha256=sha256_file(args.site),
        accession=manifest.source.accession,
        selection=selection,
        scope=(
            f"all atoms of {len(site.residues)} pocket residues selected by the "
            "binding-site record"
        ),
        toolkit="Moleculab residue templates",
        toolkit_version="1.0",
        parameters=(
            "explicit residue/atom template table",
            "no protein hydrogens generated",
            "no pH-dependent population model",
            "water and hetero roles left unassigned",
        ),
        states=states,
        selected_state_id=None,
    )
    errors = validate_protein_preparation(record)
    if errors:
        for error in errors:
            print(f"protein preparation error: {error}")
        raise SystemExit(2)

    destination = write_protein_preparation(args.output, record)
    print(f"{record.accession} protein pocket preparation")
    print(f"  scope: {record.scope}")
    print(f"  states: {len(record.states)}; selected_state_id: {record.selected_state_id!r}")
    for state in record.states:
        print(
            f"    {state.state_id}: {len(state.atoms)} atoms, "
            f"{len(state.residue_states)} residues, {len(state.rings)} rings, "
            f"{sum(1 for a in state.atoms if a.is_donor)} donors, "
            f"{sum(1 for a in state.atoms if a.is_acceptor)} acceptors, "
            f"{sum(1 for a in state.atoms if a.is_hydrophobic)} hydrophobic, "
            f"{sum(1 for a in state.atoms if a.formal_charge)} charged atoms"
        )
    print(
        f"  water roles: {len(record.states[0].water_roles)} unassigned; "
        f"hetero roles: {len(record.states[0].hetero_roles)} unknown"
    )
    print("  hydrogens: none generated")
    print(f"  record: {destination}")
    print("  scope note: template chemistry only; no experimental protonation or affinity")


if __name__ == "__main__":
    main()
