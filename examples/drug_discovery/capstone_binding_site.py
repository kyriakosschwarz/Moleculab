"""Reconcile the validated ligand with its protein binding site.

This is Batch 4 of the capstone.  It selects protein atoms around the ligand
using an explicit cutoff, records residue identities in both mmCIF namespaces,
and reports atom-level ligand-protein contacts, nearby waters, other hetero
residues, missing residue numbers, and alternate locations.

It does not assign interactions or affinity from proximity, add hydrogens,
prepare a force field, minimise, dock, or estimate a binding free energy.

Run with gemmi available:

    uv run --with gemmi python examples/drug_discovery/capstone_binding_site.py
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path

import gemmi

from moleculab import (
    HeteroResidueRecord,
    ProteinAtomRecord,
    WaterRecord,
    load_ligand_coordinates,
    load_manifest,
    reconcile_binding_site,
    sha256_file,
    validate_binding_site,
    write_binding_site,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_LIGAND = Path("data/structures/2HYY.ligand-coordinates.json")
DEFAULT_SOURCE = Path("data/structures/2HYY.cif")
DEFAULT_OUTPUT = Path("data/structures/2HYY.binding-site.json")


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--ligand", type=Path, default=DEFAULT_LIGAND)
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
                "occupancy",
                "B_iso_or_equiv",
                "label_alt_id",
                "pdbx_PDB_model_num",
            ],
        )
    ]


def _protein_atoms(
    block: gemmi.cif.Block,
    auth_chain: str,
    model_number: int,
) -> tuple[ProteinAtomRecord, ...]:
    atoms = []
    for values in _atom_rows(block):
        if values[0] != "ATOM":
            continue
        if values[3] != auth_chain:
            continue
        if values[15] not in ("?", ".") and int(values[15]) != model_number:
            continue
        atoms.append(
            ProteinAtomRecord(
                atom_name=values[7],
                element=values[8],
                residue_name=values[1],
                label_asym_id=values[2],
                auth_asym_id=values[3],
                label_seq_id=_none_if_missing(values[4]),
                auth_seq_id=values[5],
                insertion_code=_none_if_missing(values[6]),
                alternate_location=_none_if_missing(values[14]),
                model_number=(
                    int(values[15]) if values[15] not in ("?", ".") else 1
                ),
                position_angstrom=(float(values[9]), float(values[10]), float(values[11])),
                occupancy=(
                    None if values[12] in ("?", ".") else float(values[12])
                ),
                b_iso_or_equiv=(
                    None if values[13] in ("?", ".") else float(values[13])
                ),
            )
        )
    if not atoms:
        raise ValueError(f"no protein atoms found for auth chain {auth_chain!r}.")
    return tuple(atoms)


def _water_atoms(
    block: gemmi.cif.Block,
    auth_chain: str,
    model_number: int,
) -> tuple[WaterRecord, ...]:
    waters = []
    for values in _atom_rows(block):
        if values[0] != "HETATM" or values[1] != "HOH":
            continue
        if values[3] != auth_chain:
            continue
        if values[15] not in ("?", ".") and int(values[15]) != model_number:
            continue
        waters.append(
            WaterRecord(
                atom_name=values[7],
                residue_name=values[1],
                label_asym_id=values[2],
                auth_asym_id=values[3],
                label_seq_id=_none_if_missing(values[4]),
                auth_seq_id=values[5],
                insertion_code=_none_if_missing(values[6]),
                position_angstrom=(float(values[9]), float(values[10]), float(values[11])),
                occupancy=(
                    None if values[12] in ("?", ".") else float(values[12])
                ),
                b_iso_or_equiv=(
                    None if values[13] in ("?", ".") else float(values[13])
                ),
            )
        )
    return tuple(waters)


def _hetero_residues(
    block: gemmi.cif.Block,
    auth_chain: str,
    model_number: int,
    ligand_component: str,
) -> tuple[HeteroResidueRecord, ...]:
    grouped: dict[tuple[str, ...], int] = defaultdict(int)
    for values in _atom_rows(block):
        if values[0] != "HETATM":
            continue
        if values[1] == "HOH" or values[1] == ligand_component:
            continue
        if values[3] != auth_chain:
            continue
        if values[15] not in ("?", ".") and int(values[15]) != model_number:
            continue
        key = (
            values[1],
            values[2],
            values[3],
            values[4],
            values[5],
            values[6],
        )
        grouped[key] += 1
    return tuple(
        HeteroResidueRecord(
            residue_name=key[0],
            label_asym_id=key[1],
            auth_asym_id=key[2],
            label_seq_id=_none_if_missing(key[3]),
            auth_seq_id=key[4],
            insertion_code=_none_if_missing(key[5]),
            atom_count=count,
        )
        for key, count in sorted(grouped.items())
    )


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (args.manifest, args.ligand, args.source):
        if not path.exists():
            raise SystemExit(
                f"missing {path}; run the earlier capstone batches first."
            )

    manifest = load_manifest(args.manifest)
    ligand_record = load_ligand_coordinates(args.ligand)
    if ligand_record.accession != manifest.source.accession:
        raise SystemExit("manifest and ligand coordinate records disagree on accession.")
    if ligand_record.selection != manifest.selection.ligand:
        raise SystemExit("manifest and ligand coordinate records disagree on selection.")

    block = gemmi.cif.read_file(str(args.source)).sole_block()
    protein_selection = manifest.selection.protein
    ligand_selection = manifest.selection.ligand
    model_number = manifest.selection.model_number
    protein_atoms = _protein_atoms(block, protein_selection.auth_asym_id, model_number)
    water_atoms = _water_atoms(block, protein_selection.auth_asym_id, model_number)
    hetero_residues = _hetero_residues(
        block,
        protein_selection.auth_asym_id,
        model_number,
        ligand_selection.component_id,
    )

    record = reconcile_binding_site(
        source_manifest_path=str(args.manifest),
        source_manifest_sha256=sha256_file(args.manifest),
        ligand_coordinate_record_path=str(args.ligand),
        ligand_coordinate_record_sha256=sha256_file(args.ligand),
        accession=manifest.source.accession,
        model_number=model_number,
        protein_selection=protein_selection,
        ligand_selection=ligand_selection,
        protein_atoms=protein_atoms,
        ligand_atoms=ligand_record.atoms,
        water_atoms=water_atoms,
        hetero_residues=hetero_residues,
    )
    errors = validate_binding_site(record)
    if errors:
        for error in errors:
            print(f"binding-site record error: {error}")
        raise SystemExit(2)

    destination = write_binding_site(args.output, record)
    statuses = Counter(check.status for check in record.validation.checks)
    discrepancies = Counter(item.status for item in record.discrepancies)
    failures = [check for check in record.validation.checks if check.status == "fail"]

    print(f"{record.accession} binding-site reconciliation")
    print(
        f"  selection: protein {record.protein_selection.auth_asym_id} + "
        f"ligand {record.ligand_selection.component_id} "
        f"{record.ligand_selection.label_asym_id}"
    )
    print(
        f"  pocket: {record.pocket_atom_count} atoms in "
        f"{record.pocket_residue_count} residues within "
        f"{record.criteria.pocket_cutoff_angstrom} A"
    )
    print(
        f"  contacts: {len(record.contacts)} ligand-protein atom pairs; "
        f"closest {record.contacts[0].distance_angstrom:.3f} A"
        if record.contacts
        else "  contacts: none"
    )
    print(
        f"  environment: {record.water_count_in_selected_chain} waters in chain, "
        f"{len(record.waters_near_ligand)} within "
        f"{record.criteria.water_cutoff_angstrom} A; "
        f"{len(record.hetero_residues)} other hetero residues"
    )
    print(
        f"  gaps: {', '.join(record.residue_number_gaps) if record.residue_number_gaps else 'none'}; "
        f"alternate-location atoms: {record.alternate_location_atom_count}"
    )
    print(
        f"  validation: {statuses['pass']} pass, {statuses['fail']} fail, "
        f"{statuses['not_run']} not run, {statuses['not_applicable']} not applicable"
    )
    print(
        f"  discrepancies: {discrepancies['resolved']} resolved, "
        f"{discrepancies['open']} open, "
        f"{discrepancies['not_applicable']} not applicable"
    )
    print(f"  record: {destination}")
    print("  scope: proximity only; no interaction or affinity claim")

    if failures:
        for check in failures:
            print(f"  FAILED {check.check_id}: {check.detail}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
