# SPDX-FileCopyrightText: 2026 Moleculab Contributors
# SPDX-License-Identifier: MIT OR Apache-2.0

"""Prepare explicit chemical states without altering experimental records.

This is Batch 7 of the capstone.  It creates a separate preparation record with
two coexisting states:

- ``deposited-only``: no protonation, tautomer, hydrogen, charge, or feature
  assignment;
- ``reference-neutral``: explicit neutral imatinib reference state with formal
  charges from the validated reference graph, generated hydrogens, and chemical
  features from a documented toolkit configuration.

The deposited coordinates, identity record, coordinate record, and binding-site
record are untouched.  No state is selected automatically.

Run with RDKit available:

    uv run --with rdkit --with gemmi python examples/drug_discovery/capstone/capstone_chemical_preparation.py
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem import ChemicalFeatures

from moleculab import (
    PreparedAtom,
    PreparedRing,
    PreparedState,
    PreparedWaterRole,
    build_chemical_preparation_record,
    load_binding_site,
    load_chemical_preparation,
    load_ligand_coordinates,
    load_ligand_identity,
    load_manifest,
    sha256_file,
    validate_chemical_preparation,
    write_chemical_preparation,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_IDENTITY = Path("data/structures/2HYY.ligand-identity.json")
DEFAULT_LIGAND = Path("data/structures/2HYY.ligand-coordinates.json")
DEFAULT_SITE = Path("data/structures/2HYY.binding-site.json")
DEFAULT_OUTPUT = Path("data/structures/2HYY.chemical-preparation.json")

_ORDER_MAP = {
    "SINGLE": Chem.BondType.SINGLE,
    "DOUBLE": Chem.BondType.DOUBLE,
    "TRIPLE": Chem.BondType.TRIPLE,
    "AROMATIC": Chem.BondType.AROMATIC,
}


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--ligand", type=Path, default=DEFAULT_LIGAND)
    parser.add_argument("--site", type=Path, default=DEFAULT_SITE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def _plane(positions: dict[str, tuple[float, float, float]], names: tuple[str, ...]):
    points = np.asarray([positions[name] for name in names], dtype=float)
    centroid = points.mean(axis=0)
    _, _, vectors = np.linalg.svd(points - centroid, full_matrices=False)
    normal = vectors[-1]
    return tuple(centroid * 0.1), tuple(normal)


def _water_roles(site) -> tuple[PreparedWaterRole, ...]:
    return tuple(
        PreparedWaterRole(
            water_id=f"{water.auth_asym_id}:{water.auth_seq_id}:{water.atom_name}",
            role="unassigned",
            source="experimental water; no interaction role assigned",
        )
        for water in site.waters_near_ligand
    )


def _deposited_only_state(identity, ligand, site) -> PreparedState:
    atoms = tuple(
        PreparedAtom(
            atom_id=atom.atom_name,
            element=atom.element,
            formal_charge=None,
            is_donor=None,
            is_acceptor=None,
            is_hydrophobic=None,
            aromatic_ring_ids=(),
            hydrogen_count=0,
            hydrogen_positions_nm=(),
            coordinate_source="deposited",
            hydrogen_source=None,
            feature_source=None,
        )
        for atom in ligand.atoms
    )
    return PreparedState(
        state_id="deposited-only",
        label="Deposited heavy-atom state without chemical assignment",
        protonation_state="unknown",
        tautomer="unknown",
        formal_charge=None,
        hydrogen_method=None,
        feature_source=None,
        assumptions=(
            "No protonation, tautomer, hydrogen, formal-charge, donor/acceptor, "
            "hydrophobic, or aromatic feature assignment was made.",
            "Only deposited heavy-atom coordinates are represented.",
        ),
        atoms=atoms,
        rings=(),
        water_roles=_water_roles(site),
    )


def _reference_state(identity, ligand, site) -> PreparedState:
    indices: dict[str, int] = {}
    editable = Chem.RWMol()
    for atom in identity.deposited.atoms:
        indices[atom.atom_name] = editable.AddAtom(Chem.Atom(atom.element))
    for bond in identity.deposited.bonds:
        editable.AddBond(
            indices[bond.atom_name_1],
            indices[bond.atom_name_2],
            _ORDER_MAP[bond.order],
        )
    molecule = editable.GetMol()
    for pair in identity.mapping.pairs:
        reference_atom = identity.reference.atoms[pair[1]]
        molecule.GetAtomWithIdx(indices[pair[0]]).SetFormalCharge(
            reference_atom.formal_charge
        )
    conformer = Chem.Conformer(molecule.GetNumAtoms())
    positions = {atom.atom_name: atom.position_angstrom for atom in ligand.atoms}
    for name, index in indices.items():
        conformer.SetAtomPosition(index, positions[name])
    molecule.AddConformer(conformer)
    Chem.SanitizeMol(molecule)

    prepared = Chem.AddHs(molecule, addCoords=True)
    Chem.SanitizeMol(prepared)

    feature_factory = ChemicalFeatures.BuildFeatureFactory(
        str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef")
    )
    features = feature_factory.GetFeaturesForMol(prepared)
    donors: set[int] = set()
    acceptors: set[int] = set()
    hydrophobic: set[int] = set()
    for feature in features:
        family = feature.GetFamily()
        if family == "Donor":
            donors.update(feature.GetAtomIds())
        elif family == "Acceptor":
            acceptors.update(feature.GetAtomIds())
        elif family in ("Hydrophobe", "LumpedHydrophobe"):
            hydrophobic.update(feature.GetAtomIds())

    heavy_names = tuple(atom.atom_name for atom in identity.deposited.atoms)
    aromatic_rings: list[PreparedRing] = []
    for ring_index, ring in enumerate(Chem.GetSymmSSSR(prepared)):
        if not all(prepared.GetAtomWithIdx(index).GetIsAromatic() for index in ring):
            continue
        names = tuple(heavy_names[index] for index in ring if index < len(heavy_names))
        if len(names) < 3:
            continue
        centroid, normal = _plane(positions, names)
        aromatic_rings.append(
            PreparedRing(
                ring_id=f"aromatic-{ring_index + 1}",
                atom_ids=names,
                centroid_nm=centroid,
                normal=normal,
                source="RDKit aromaticity perception; SVD plane from prepared heavy atoms",
            )
        )
    ring_ids_by_atom: dict[str, list[str]] = {}
    for ring in aromatic_rings:
        for name in ring.atom_ids:
            ring_ids_by_atom.setdefault(name, []).append(ring.ring_id)

    prepared_atoms = []
    for index, atom in enumerate(identity.deposited.atoms):
        rd_atom = prepared.GetAtomWithIdx(index)
        hydrogen_positions = []
        for neighbour in rd_atom.GetNeighbors():
            if neighbour.GetAtomicNum() == 1:
                position = prepared.GetConformer().GetAtomPosition(neighbour.GetIdx())
                hydrogen_positions.append(
                    (position.x * 0.1, position.y * 0.1, position.z * 0.1)
                )
        prepared_atoms.append(
            PreparedAtom(
                atom_id=atom.atom_name,
                element=atom.element,
                formal_charge=rd_atom.GetFormalCharge(),
                is_donor=index in donors,
                is_acceptor=index in acceptors,
                is_hydrophobic=index in hydrophobic,
                aromatic_ring_ids=tuple(ring_ids_by_atom.get(atom.atom_name, ())),
                hydrogen_count=len(hydrogen_positions),
                hydrogen_positions_nm=tuple(hydrogen_positions),
                coordinate_source="deposited",
                hydrogen_source=(
                    "RDKit AddHs(addCoords=True)" if hydrogen_positions else None
                ),
                feature_source=(
                    "RDKit BaseFeatures.fdef; formal charges from validated reference graph"
                ),
            )
        )

    return PreparedState(
        state_id="reference-neutral",
        label="Neutral imatinib reference state on the deposited heavy-atom conformer",
        protonation_state="reference neutral; no pH-specific protonation selection",
        tautomer="reference tautomer from validated SMILES; alternatives not enumerated",
        formal_charge=0,
        hydrogen_method="RDKit AddHs(addCoords=True)",
        feature_source=(
            "RDKit BaseFeatures.fdef, aromaticity perception, and validated reference formal charges"
        ),
        assumptions=(
            "Formal charges are taken from the validated imatinib reference graph.",
            "Hydrogen coordinates are generated, not experimentally observed.",
            "No pH-specific protonation or alternative tautomer enumeration was performed.",
            "Water roles remain unassigned; experimental waters are preserved separately.",
        ),
        atoms=tuple(prepared_atoms),
        rings=tuple(aromatic_rings),
        water_roles=_water_roles(site),
    )


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (args.manifest, args.identity, args.ligand, args.site):
        if not path.exists():
            raise SystemExit(f"missing {path}; run the earlier capstone batches first.")

    manifest = load_manifest(args.manifest)
    identity = load_ligand_identity(args.identity)
    ligand = load_ligand_coordinates(args.ligand)
    site = load_binding_site(args.site)
    if not (
        identity.accession
        == ligand.accession
        == site.accession
        == manifest.source.accession
    ):
        raise SystemExit("capstone records disagree on accession.")

    states = (
        _deposited_only_state(identity, ligand, site),
        _reference_state(identity, ligand, site),
    )
    record = build_chemical_preparation_record(
        source_manifest_path=str(args.manifest),
        source_manifest_sha256=sha256_file(args.manifest),
        ligand_identity_record_path=str(args.identity),
        ligand_identity_record_sha256=sha256_file(args.identity),
        ligand_coordinate_record_path=str(args.ligand),
        ligand_coordinate_record_sha256=sha256_file(args.ligand),
        binding_site_record_path=str(args.site),
        binding_site_record_sha256=sha256_file(args.site),
        accession=manifest.source.accession,
        selection=identity.selection,
        toolkit="RDKit",
        toolkit_version=rdBase.rdkitVersion,
        parameters=(
            "formal charges transferred from the validated reference graph",
            "hydrogens generated with RDKit AddHs(addCoords=True)",
            "donor/acceptor/hydrophobic features from RDKit BaseFeatures.fdef",
            "aromatic rings from RDKit aromaticity perception",
            "ring centroids and normals from prepared heavy-atom coordinates",
            "no pH-specific protonation or tautomer enumeration",
        ),
        states=states,
        selected_state_id=None,
    )
    errors = validate_chemical_preparation(record)
    if errors:
        for error in errors:
            print(f"preparation record error: {error}")
        raise SystemExit(2)

    destination = write_chemical_preparation(args.output, record)
    reference = next(state for state in record.states if state.state_id == "reference-neutral")
    families = Counter(
        "donor" if atom.is_donor else "acceptor" if atom.is_acceptor else "other"
        for atom in reference.atoms
    )
    print(f"{record.accession} chemical-state preparation")
    print(f"  states: {len(record.states)}; selected_state_id: {record.selected_state_id!r}")
    for state in record.states:
        print(
            f"    {state.state_id}: {len(state.atoms)} atoms, "
            f"{len(state.rings)} rings, charge={state.formal_charge}"
        )
    print(
        f"  reference state: {sum(a.hydrogen_count for a in reference.atoms)} generated hydrogens, "
        f"{len(reference.rings)} aromatic rings, "
        f"{sum(1 for a in reference.atoms if a.is_donor)} donors, "
        f"{sum(1 for a in reference.atoms if a.is_acceptor)} acceptors"
    )
    print(
        f"  water roles: {len(reference.water_roles)} preserved, all unassigned"
    )
    print(f"  toolkit: {record.toolkit} {record.toolkit_version}")
    print(f"  record: {destination}")
    print("  scope: separate preparation record; experimental records unchanged")


if __name__ == "__main__":
    main()
