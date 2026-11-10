# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Reconcile the selected 2HYY ligand with the imatinib reference identity.

This is Batch 2 of the capstone.  It reads the ligand selected by the Batch 1
manifest, preserves the deposited atoms and heavy-atom bonds, builds an explicit
imatinib reference graph with RDKit, compares the two graphs, and writes an
identity record with explicit atom mapping and discrepancies.

It does not add hydrogens, assign formal charges, change protonation, repair
coordinates, dock a ligand, or estimate affinity.

Run with gemmi and RDKit available:

    uv run --with gemmi --with rdkit python examples/drug_discovery/capstone/capstone_ligand_identity.py
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import FindMolChiralCenters, rdMolDescriptors

import gemmi

from moleculab import (
    DepositedAtom,
    DepositedBond,
    DepositedLigand,
    NormalizationRecord,
    ReferenceAtom,
    ReferenceBond,
    ReferenceLigand,
    load_manifest,
    reconcile_ligand_identity,
    sha256_file,
    validate_ligand_identity,
    write_ligand_identity,
)

IMATINIB_SMILES = (
    "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)"
    "NC4=NC=CC(=N4)C5=CN=CC=C5"
)
REFERENCE_SOURCE = (
    "imatinib reference SMILES used by the capstone chemical-record example; "
    "canonicalized and kekulized by RDKit"
)
DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_SOURCE = Path("data/structures/2HYY.cif")
DEFAULT_OUTPUT = Path("data/structures/2HYY.ligand-identity.json")

_ORDER_MAP = {
    "sing": "SINGLE",
    "doub": "DOUBLE",
    "trip": "TRIPLE",
    "arom": "AROMATIC",
}
_RDKIT_ORDER_MAP = {
    "SINGLE": Chem.BondType.SINGLE,
    "DOUBLE": Chem.BondType.DOUBLE,
    "TRIPLE": Chem.BondType.TRIPLE,
    "AROMATIC": Chem.BondType.AROMATIC,
}


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def _unquote(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in ("'", '"'):
        return value[1:-1]
    return value


def _value(block: gemmi.cif.Block, tag: str) -> str | None:
    value = block.find_value(tag)
    if value is None or value in ("?", "."):
        return None
    return _unquote(value)


def _chem_comp(block: gemmi.cif.Block, component_id: str) -> tuple[str, str]:
    for row in block.find("_chem_comp.", ["id", "name", "formula"]):
        values = [_unquote(value) for value in row]
        if values[0] == component_id:
            return values[1], values[2]
    raise ValueError(f"_chem_comp has no record for component {component_id!r}.")


def _component_atoms(
    block: gemmi.cif.Block,
    component_id: str,
) -> dict[str, tuple[str, str | None]]:
    atoms = {}
    for row in block.find(
        "_chem_comp_atom.",
        ["comp_id", "atom_id", "type_symbol", "pdbx_stereo_config"],
    ):
        values = [_unquote(value) for value in row]
        if values[0] != component_id:
            continue
        stereo = None if values[3] in ("N", "?", ".") else values[3]
        atoms[values[1]] = (values[2], stereo)
    if not atoms:
        raise ValueError(
            f"_chem_comp_atom has no records for component {component_id!r}."
        )
    return atoms


def _deposited_ligand(
    block: gemmi.cif.Block,
    selection,
    model_number: int,
) -> DepositedLigand:
    atom_rows = []
    for row in block.find(
        "_atom_site.",
        [
            "label_comp_id",
            "label_asym_id",
            "auth_asym_id",
            "auth_seq_id",
            "label_seq_id",
            "pdbx_PDB_ins_code",
            "label_atom_id",
            "type_symbol",
            "pdbx_formal_charge",
            "pdbx_PDB_model_num",
        ],
    ):
        values = [_unquote(value) for value in row]
        if values[0] != selection.component_id:
            continue
        if selection.label_asym_id is not None and values[1] != selection.label_asym_id:
            continue
        if values[2] != selection.auth_asym_id:
            continue
        if values[3] != selection.auth_seq_id:
            continue
        if values[9] not in ("?", ".") and int(values[9]) != model_number:
            continue
        atom_rows.append(values)
    if not atom_rows:
        raise ValueError(
            "the Batch 1 ligand selection did not match any deposited atom records."
        )

    component_atoms = _component_atoms(block, selection.component_id)
    atoms = tuple(
        DepositedAtom(
            atom_name=values[6],
            element=values[7],
            formal_charge=(
                None if values[8] in ("?", ".") else int(values[8])
            ),
        )
        for values in atom_rows
    )
    atom_names = {atom.atom_name for atom in atoms}
    bond_rows = []
    omitted_bond_records = 0
    omitted_atom_names: set[str] = set()
    for row in block.find(
        "_chem_comp_bond.",
        ["comp_id", "atom_id_1", "atom_id_2", "value_order", "pdbx_aromatic_flag"],
    ):
        values = [_unquote(value) for value in row]
        if values[0] != selection.component_id:
            continue
        if values[1] in atom_names and values[2] in atom_names:
            bond_rows.append(values)
        else:
            omitted_bond_records += 1
            omitted_atom_names.update(
                name for name in (values[1], values[2]) if name not in atom_names
            )
    omitted_bond_atom_elements = tuple(
        sorted(
            {
                component_atoms[name][0]
                for name in omitted_atom_names
                if name in component_atoms
            }
        )
    )

    bonds = tuple(
        DepositedBond(
            atom_name_1=values[1],
            atom_name_2=values[2],
            order=_ORDER_MAP[values[3]],
            aromatic=values[4] == "Y",
        )
        for values in bond_rows
    )
    _component_name, formula = _chem_comp(block, selection.component_id)
    return DepositedLigand(
        component_id=selection.component_id,
        label_asym_id=selection.label_asym_id or "",
        auth_asym_id=selection.auth_asym_id,
        auth_seq_id=selection.auth_seq_id,
        label_seq_id=selection.label_seq_id,
        insertion_code=selection.insertion_code,
        formula=formula,
        atoms=atoms,
        bonds=bonds,
        atoms_are_heavy=all(atom.element.upper() != "H" for atom in atoms),
        formal_charges_present=any(atom.formal_charge is not None for atom in atoms),
        hydrogen_atoms_present=any(atom.element.upper() == "H" for atom in atoms),
        omitted_bond_records=omitted_bond_records,
        omitted_bond_atom_elements=omitted_bond_atom_elements,
        stereochemistry=tuple(
            sorted(
                f"{name}:{config}"
                for name, (_element, config) in component_atoms.items()
                if name in atom_names and config is not None
            )
        ),
    )


def _reference_ligand() -> tuple[ReferenceLigand, Chem.Mol]:
    molecule = Chem.MolFromSmiles(IMATINIB_SMILES)
    if molecule is None:
        raise ValueError("the imatinib reference SMILES failed to parse.")
    Chem.AssignStereochemistry(molecule, cleanIt=True, force=True)
    aromatic = Chem.Mol(molecule)
    kekule = Chem.Mol(molecule)
    Chem.Kekulize(kekule, clearAromaticFlags=True)

    atoms = tuple(
        ReferenceAtom(
            index=atom.GetIdx(),
            element=atom.GetSymbol(),
            formal_charge=atom.GetFormalCharge(),
        )
        for atom in molecule.GetAtoms()
    )
    bonds = []
    for bond in molecule.GetBonds():
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        kekule_bond = kekule.GetBondBetweenAtoms(begin, end)
        aromatic_bond = aromatic.GetBondBetweenAtoms(begin, end)
        if kekule_bond is None or aromatic_bond is None:
            raise ValueError("reference bond copies are inconsistent.")
        bonds.append(
            ReferenceBond(
                atom_index_1=begin,
                atom_index_2=end,
                order=str(kekule_bond.GetBondType()),
                aromatic=bool(aromatic_bond.GetIsAromatic()),
            )
        )
    stereocenters = tuple(
        f"{index}:{label}"
        for index, label in FindMolChiralCenters(
            molecule,
            includeUnassigned=True,
            useLegacyImplementation=False,
        )
    )
    reference = ReferenceLigand(
        name="imatinib",
        source=REFERENCE_SOURCE,
        smiles=Chem.MolToSmiles(molecule),
        formula=rdMolDescriptors.CalcMolFormula(molecule),
        formal_charge=Chem.GetFormalCharge(molecule),
        hydrogen_count=sum(atom.GetTotalNumHs() for atom in molecule.GetAtoms()),
        atoms=atoms,
        bonds=tuple(bonds),
        atoms_are_heavy=True,
        stereocenters=stereocenters,
    )
    return reference, molecule


def _deposited_molecule(deposited: DepositedLigand) -> Chem.Mol:
    molecule = Chem.RWMol()
    index_by_name = {}
    for atom in deposited.atoms:
        index_by_name[atom.atom_name] = molecule.AddAtom(Chem.Atom(atom.element))
    for bond in deposited.bonds:
        molecule.AddBond(
            index_by_name[bond.atom_name_1],
            index_by_name[bond.atom_name_2],
            _RDKIT_ORDER_MAP[bond.order],
        )
    result = molecule.GetMol()
    Chem.SanitizeMol(result)
    return result


def _generic_graph(molecule: Chem.Mol) -> Chem.Mol:
    graph = Chem.RWMol(molecule)
    for bond in graph.GetBonds():
        bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    result = graph.GetMol()
    for atom in result.GetAtoms():
        atom.SetIsAromatic(False)
    return result


def _candidate_mappings(
    deposited: DepositedLigand,
    reference_molecule: Chem.Mol,
) -> list[dict[str, int]]:
    deposited_molecule = _deposited_molecule(deposited)
    target = _generic_graph(deposited_molecule)
    query = _generic_graph(reference_molecule)
    matches = target.GetSubstructMatches(query, uniquify=False)
    candidates = []
    for match in matches:
        candidates.append(
            {
                deposited.atoms[deposited_index].atom_name: reference_index
                for reference_index, deposited_index in enumerate(match)
            }
        )
    return candidates


def main() -> None:
    args = _argument_parser().parse_args()
    if not args.manifest.exists():
        raise SystemExit(f"missing {args.manifest}; run capstone_structure.py first.")
    if not args.source.exists():
        raise SystemExit(f"missing {args.source}; run capstone_structure.py first.")

    manifest = load_manifest(args.manifest)
    selection = manifest.selection.ligand
    structure = gemmi.read_structure(str(args.source))
    if len(structure) == 0:
        raise SystemExit(f"{args.source} contains no models.")
    block = gemmi.cif.read_file(str(args.source)).sole_block()
    deposited = _deposited_ligand(block, selection, manifest.selection.model_number)
    reference, reference_molecule = _reference_ligand()
    candidates = _candidate_mappings(deposited, reference_molecule)

    record = reconcile_ligand_identity(
        source_manifest_path=str(args.manifest),
        source_manifest_sha256=sha256_file(args.manifest),
        accession=manifest.source.accession,
        selection=selection,
        deposited=deposited,
        reference=reference,
        candidate_mappings=candidates,
        normalization=NormalizationRecord(
            bond_order_rule=(
                "raw Kekule orders are compared first; aromatic bonds are "
                "treated as equivalent only in the separate aromatic-normalized "
                "comparison"
            ),
            aromaticity_method=(
                "RDKit sanitization and built-in aromaticity perception"
            ),
            toolkit="RDKit",
            toolkit_version=rdBase.rdkitVersion,
        ),
    )
    errors = validate_ligand_identity(record)
    if errors:
        for error in errors:
            print(f"identity record error: {error}")
        raise SystemExit(2)

    destination = write_ligand_identity(args.output, record)
    statuses = Counter(check.status for check in record.validation.checks)
    discrepancy_statuses = Counter(item.status for item in record.discrepancies)
    failures = [check for check in record.validation.checks if check.status == "fail"]

    print(f"{record.accession} ligand identity reconciliation")
    print(
        f"  selected: {record.selection.component_id} label "
        f"{record.selection.label_asym_id} auth {record.selection.auth_asym_id} "
        f"seq {record.selection.auth_seq_id}"
    )
    omitted_elements = ", ".join(deposited.omitted_bond_atom_elements) or "none"
    print(
        f"  deposited: {len(deposited.atoms)} atoms, {len(deposited.bonds)} heavy-atom "
        f"bonds, {deposited.omitted_bond_records} component bond records omitted "
        f"(atoms absent from coordinates: {omitted_elements})"
    )
    print(f"  reference: {reference.name}, {reference.formula}, charge {reference.formal_charge}")
    print(
        f"  mapping: {record.mapping.candidate_count} candidate(s), "
        f"ambiguous={record.mapping.ambiguous}, "
        f"ordering={record.mapping.ordering_basis}"
    )
    print(
        f"  normalization: {record.normalization.toolkit} "
        f"{record.normalization.toolkit_version} "
        f"({record.normalization.aromaticity_method})"
    )
    print(
        f"  validation: {statuses['pass']} pass, {statuses['fail']} fail, "
        f"{statuses['not_run']} not run"
    )
    print(
        f"  discrepancies: {discrepancy_statuses['resolved']} resolved, "
        f"{discrepancy_statuses['open']} open, "
        f"{discrepancy_statuses['not_applicable']} not applicable"
    )
    print(f"  record: {destination}")
    print("  scope: no coordinates, hydrogens, protonation, or formal charges were changed")

    if failures:
        for check in failures:
            print(f"  FAILED {check.check_id}: {check.detail}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
