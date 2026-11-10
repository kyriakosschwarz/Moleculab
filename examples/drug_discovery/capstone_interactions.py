# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Interpret ligand-protein contacts with explicit chemical evidence.

This is Batch 5 of the capstone.  It consumes the Batch 4 contact record and
classifies an interaction only when donor/acceptor, charge, hydrophobic, ring,
water, or metal evidence is explicitly supplied.  For the deposited 2HYY
structure, most of that chemistry is unavailable, so the expected result is
zero confidently classified interactions rather than a speculative diagram.

It does not infer donor/acceptor identity from elements, assign formal charges,
or convert interactions into binding affinity.

Run with gemmi available:

    uv run --with gemmi python examples/drug_discovery/capstone_interactions.py
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import gemmi

from moleculab import (
    AtomFeature,
    ProteinAtomRecord,
    load_binding_site,
    load_ligand_coordinates,
    load_ligand_identity,
    load_manifest,
    interpret_interactions,
    sha256_file,
    validate_interactions,
    write_interactions,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_IDENTITY = Path("data/structures/2HYY.ligand-identity.json")
DEFAULT_LIGAND = Path("data/structures/2HYY.ligand-coordinates.json")
DEFAULT_SITE = Path("data/structures/2HYY.binding-site.json")
DEFAULT_SOURCE = Path("data/structures/2HYY.cif")
DEFAULT_OUTPUT = Path("data/structures/2HYY.interactions.json")


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--ligand", type=Path, default=DEFAULT_LIGAND)
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
                "pdbx_PDB_model_num",
            ],
        )
    ]


def _protein_lookup(
    block: gemmi.cif.Block,
    auth_chain: str,
    model_number: int,
) -> dict[tuple[str, ...], ProteinAtomRecord]:
    lookup = {}
    for values in _atom_rows(block):
        if values[0] != "ATOM" or values[3] != auth_chain:
            continue
        if values[12] not in ("?", ".") and int(values[12]) != model_number:
            continue
        record = ProteinAtomRecord(
            atom_name=values[7],
            element=values[8],
            residue_name=values[1],
            label_asym_id=values[2],
            auth_asym_id=values[3],
            label_seq_id=_none_if_missing(values[4]),
            auth_seq_id=values[5],
            insertion_code=_none_if_missing(values[6]),
            alternate_location=None,
            model_number=int(values[12]) if values[12] not in ("?", ".") else 1,
            position_angstrom=(float(values[9]), float(values[10]), float(values[11])),
            occupancy=None,
            b_iso_or_equiv=None,
        )
        key = (
            record.label_asym_id,
            record.auth_asym_id,
            record.label_seq_id,
            record.auth_seq_id,
            record.insertion_code,
            record.residue_name,
            record.atom_name,
        )
        lookup[key] = record
    return lookup


def _water_lookup(
    block: gemmi.cif.Block,
    auth_chain: str,
    model_number: int,
) -> dict[tuple[str, str], tuple[float, float, float]]:
    lookup = {}
    for values in _atom_rows(block):
        if values[0] != "HETATM" or values[1] != "HOH" or values[3] != auth_chain:
            continue
        if values[12] not in ("?", ".") and int(values[12]) != model_number:
            continue
        lookup[(values[5], values[7])] = (
            float(values[9]),
            float(values[10]),
            float(values[11]),
        )
    return lookup


def _protein_feature(atom: ProteinAtomRecord) -> AtomFeature:
    return AtomFeature(
        atom_id=f"protein:{atom.auth_asym_id}:{atom.auth_seq_id}:{atom.atom_name}",
        side="protein",
        element=atom.element,
        position_nm=tuple(value * 0.1 for value in atom.position_angstrom),
        residue_id=(
            f"{atom.residue_name} {atom.auth_asym_id} {atom.auth_seq_id}"
        ),
        is_donor=None,
        is_acceptor=None,
        formal_charge=None,
        is_hydrophobic=None,
        aromatic_ring_id=None,
        is_metal=None,
        hydrogen_positions_nm=(),
        feature_source=None,
    )


def _ligand_feature(atom, component_id: str, label_asym: str, auth_seq: str) -> AtomFeature:
    return AtomFeature(
        atom_id=f"ligand:{atom.atom_name}",
        side="ligand",
        element=atom.element,
        position_nm=tuple(value * 0.1 for value in atom.position_angstrom),
        residue_id=f"{component_id} {label_asym} {auth_seq}",
        is_donor=None,
        is_acceptor=None,
        formal_charge=None,
        is_hydrophobic=None,
        aromatic_ring_id=None,
        is_metal=None,
        hydrogen_positions_nm=(),
        feature_source=None,
    )


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (
        args.manifest,
        args.identity,
        args.ligand,
        args.site,
        args.source,
    ):
        if not path.exists():
            raise SystemExit(
                f"missing {path}; run the earlier capstone batches first."
            )

    manifest = load_manifest(args.manifest)
    identity = load_ligand_identity(args.identity)
    ligand = load_ligand_coordinates(args.ligand)
    site = load_binding_site(args.site)
    if ligand.accession != manifest.source.accession or site.accession != manifest.source.accession:
        raise SystemExit("capstone records disagree on accession.")
    if site.ligand_selection != identity.selection:
        raise SystemExit("binding-site and identity records disagree on ligand selection.")

    block = gemmi.cif.read_file(str(args.source)).sole_block()
    protein_lookup = _protein_lookup(
        block,
        manifest.selection.protein.auth_asym_id,
        manifest.selection.model_number,
    )
    water_lookup = _water_lookup(
        block,
        manifest.selection.protein.auth_asym_id,
        manifest.selection.model_number,
    )

    protein_features = {}
    ligand_features = {}
    for contact in site.contacts:
        key = (
            contact.label_asym_id,
            contact.auth_asym_id,
            contact.label_seq_id,
            contact.auth_seq_id,
            contact.insertion_code,
            contact.residue_name,
            contact.protein_atom,
        )
        if key not in protein_lookup:
            raise SystemExit(f"cannot resolve contact protein atom {key}.")
        protein = protein_lookup[key]
        protein_features[key] = _protein_feature(protein)
        ligand_atom = next(
            (atom for atom in ligand.atoms if atom.atom_name == contact.ligand_atom),
            None,
        )
        if ligand_atom is None:
            raise SystemExit(f"cannot resolve contact ligand atom {contact.ligand_atom}.")
        ligand_features[ligand_atom.atom_name] = _ligand_feature(
            ligand_atom,
            identity.selection.component_id,
            identity.selection.label_asym_id or "",
            identity.selection.auth_seq_id,
        )

    water_features = []
    for water in site.waters_near_ligand:
        position = water_lookup.get((water.auth_seq_id, water.atom_name))
        if position is None:
            raise SystemExit(f"cannot resolve nearby water {water.auth_seq_id}.")
        water_features.append(
            AtomFeature(
                atom_id=f"water:{water.auth_seq_id}:{water.atom_name}",
                side="water",
                element="O",
                position_nm=tuple(value * 0.1 for value in position),
                residue_id=f"HOH {water.auth_asym_id} {water.auth_seq_id}",
                is_donor=None,
                is_acceptor=None,
                formal_charge=None,
                is_hydrophobic=None,
                aromatic_ring_id=None,
                is_metal=None,
                hydrogen_positions_nm=(),
                feature_source=None,
            )
        )

    record = interpret_interactions(
        source_manifest_path=str(args.manifest),
        source_manifest_sha256=sha256_file(args.manifest),
        ligand_identity_record_path=str(args.identity),
        ligand_identity_record_sha256=sha256_file(args.identity),
        binding_site_record_path=str(args.site),
        binding_site_record_sha256=sha256_file(args.site),
        accession=manifest.source.accession,
        model_number=manifest.selection.model_number,
        protein_selection=manifest.selection.protein,
        ligand_selection=identity.selection,
        input_contact_count=len(site.contacts),
        protein_features=tuple(protein_features.values()),
        ligand_features=tuple(ligand_features.values()),
        water_features=tuple(water_features),
        metal_features=(),
        aromatic_rings=(),
    )
    errors = validate_interactions(record)
    if errors:
        for error in errors:
            print(f"interaction record error: {error}")
        raise SystemExit(2)

    destination = write_interactions(args.output, record)
    statuses = Counter(item.status for item in record.interactions)

    print(f"{record.accession} interaction interpretation")
    print(f"  input contacts: {record.input_contact_count}")
    print(
        f"  features: {len(protein_features)} protein, {len(ligand_features)} ligand, "
        f"{len(water_features)} water"
    )
    print(
        f"  interactions: {statuses['classified']} classified, "
        f"{statuses['not_classified']} not classified, "
        f"{statuses['unknown']} unknown, "
        f"{statuses['not_applicable']} not applicable"
    )
    print("  classified types: none")
    print("  unknown reasons: donor/acceptor, formal charge, hydrophobic, ring, and water features unavailable")
    print(f"  record: {destination}")
    print("  scope: no chemical features or charges were inferred; no affinity was computed")

    failures = [check for check in record.validation.checks if check.status == "fail"]
    if failures:
        for check in failures:
            print(f"  FAILED {check.check_id}: {check.detail}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
