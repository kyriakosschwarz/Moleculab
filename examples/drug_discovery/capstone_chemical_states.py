"""Enumerate protonation and tautomer states without selecting or weighting one.

This is Batch 8 of the capstone.  It reuses the Batch 7 reference-neutral state
as the reference, enumerates RDKit tautomers, and creates hypothetical
protonation alternatives at explicit ionizable sites.  Every state carries its
own graph, proton count, formal charge, features, rings, and coordinate
provenance.

No pKa model, pH-dependent population, or affinity estimate is supplied.

Run with RDKit available:

    uv run --with rdkit --with gemmi python examples/drug_discovery/capstone_chemical_states.py
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem, RDConfig, rdBase
from rdkit.Chem import AllChem, ChemicalFeatures, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize

from moleculab import (
    ChemicalState,
    IonizableSite,
    PreparedAtom,
    PreparedRing,
    build_chemical_state_enumeration_record,
    load_chemical_preparation,
    load_chemical_state_enumeration,
    load_ligand_identity,
    load_manifest,
    sha256_file,
    validate_chemical_state_enumeration,
    write_chemical_state_enumeration,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_IDENTITY = Path("data/structures/2HYY.ligand-identity.json")
DEFAULT_LIGAND = Path("data/structures/2HYY.ligand-coordinates.json")
DEFAULT_SITE = Path("data/structures/2HYY.binding-site.json")
DEFAULT_PREPARATION = Path("data/structures/2HYY.chemical-preparation.json")
DEFAULT_OUTPUT = Path("data/structures/2HYY.chemical-states.json")


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--ligand", type=Path, default=DEFAULT_LIGAND)
    parser.add_argument("--site", type=Path, default=DEFAULT_SITE)
    parser.add_argument("--preparation", type=Path, default=DEFAULT_PREPARATION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def _feature_sets(molecule: Chem.Mol):
    factory = ChemicalFeatures.BuildFeatureFactory(
        str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef")
    )
    donors: set[int] = set()
    acceptors: set[int] = set()
    hydrophobic: set[int] = set()
    for feature in factory.GetFeaturesForMol(molecule):
        family = feature.GetFamily()
        if family == "Donor":
            donors.update(feature.GetAtomIds())
        elif family == "Acceptor":
            acceptors.update(feature.GetAtomIds())
        elif family in ("Hydrophobe", "LumpedHydrophobe"):
            hydrophobic.update(feature.GetAtomIds())
    return donors, acceptors, hydrophobic, factory


def _plane(conformer, indices):
    points = np.asarray(
        [
            (
                conformer.GetAtomPosition(index).x,
                conformer.GetAtomPosition(index).y,
                conformer.GetAtomPosition(index).z,
            )
            for index in indices
        ],
        dtype=float,
    )
    centroid = points.mean(axis=0)
    _, _, vectors = np.linalg.svd(points - centroid, full_matrices=False)
    return tuple(centroid * 0.1), tuple(vectors[-1])


def _prepared_state_from_molecule(
    molecule: Chem.Mol,
    *,
    state_id: str,
    label: str,
    kind: str,
    parent_state_id: str,
    tautomer_of: str | None,
    pka_sites: tuple[IonizableSite, ...],
    enumeration_method: str,
    enumeration_parameters: tuple[str, ...],
    assumptions: tuple[str, ...],
) -> ChemicalState:
    heavy_count = molecule.GetNumAtoms()
    prepared = Chem.AddHs(molecule)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = 0xF00D
    if AllChem.EmbedMolecule(prepared, parameters) != 0:
        raise RuntimeError(f"failed to embed {state_id}")
    Chem.SanitizeMol(prepared)
    conformer = prepared.GetConformer()
    donors, acceptors, hydrophobic, _factory = _feature_sets(prepared)

    ring_ids_by_atom: dict[int, list[str]] = {}
    rings = []
    for ring_index, ring in enumerate(Chem.GetSymmSSSR(prepared)):
        if not all(prepared.GetAtomWithIdx(index).GetIsAromatic() for index in ring):
            continue
        indices = tuple(int(index) for index in ring)
        names = tuple(f"atom_{index}" for index in indices)
        centroid, normal = _plane(conformer, indices)
        rings.append(
            PreparedRing(
                ring_id=f"aromatic-{ring_index + 1}",
                atom_ids=names,
                centroid_nm=centroid,
                normal=normal,
                source="RDKit aromaticity perception; ETKDGv3 conformer",
            )
        )
        for index, name in zip(indices, names, strict=True):
            ring_ids_by_atom.setdefault(index, []).append(f"aromatic-{ring_index + 1}")

    atoms = []
    for index in range(heavy_count):
        rd_atom = prepared.GetAtomWithIdx(index)
        hydrogens = []
        for neighbour in rd_atom.GetNeighbors():
            if neighbour.GetAtomicNum() == 1:
                position = conformer.GetAtomPosition(neighbour.GetIdx())
                hydrogens.append(
                    (position.x * 0.1, position.y * 0.1, position.z * 0.1)
                )
        atoms.append(
            PreparedAtom(
                atom_id=f"atom_{index}",
                element=rd_atom.GetSymbol(),
                formal_charge=rd_atom.GetFormalCharge(),
                is_donor=index in donors,
                is_acceptor=index in acceptors,
                is_hydrophobic=index in hydrophobic,
                aromatic_ring_ids=tuple(ring_ids_by_atom.get(index, ())),
                hydrogen_count=len(hydrogens),
                hydrogen_positions_nm=tuple(hydrogens),
                coordinate_source="generated",
                hydrogen_source="RDKit AddHs + ETKDGv3",
                feature_source="RDKit BaseFeatures.fdef",
            )
        )

    return ChemicalState(
        state_id=state_id,
        kind=kind,
        label=label,
        parent_state_id=parent_state_id,
        molecular_graph=Chem.MolToSmiles(molecule, canonical=True),
        formula=rdMolDescriptors.CalcMolFormula(molecule),
        proton_count=sum(
            1 for atom in prepared.GetAtoms() if atom.GetAtomicNum() == 1
        ),
        formal_charge=Chem.GetFormalCharge(molecule),
        tautomer_of=tautomer_of,
        ph_assumption=None,
        pka_sites=pka_sites,
        enumeration_method=enumeration_method,
        enumeration_parameters=enumeration_parameters,
        assumptions=assumptions,
        atoms=tuple(atoms),
        reference_atom_mapping=tuple(
            (f"atom_{index}", index) for index in range(heavy_count)
        ),
        rings=tuple(rings),
        water_roles=(),
    )


def _reference_state(preparation, identity) -> ChemicalState:
    prepared = next(
        state for state in preparation.states if state.state_id == "reference-neutral"
    )
    reference_indices = dict(identity.mapping.pairs)
    return ChemicalState(
        state_id="reference-neutral",
        kind="reference",
        label="Batch 7 reference-neutral prepared state",
        parent_state_id=None,
        molecular_graph=identity.reference.smiles,
        formula=identity.reference.formula,
        proton_count=31,
        formal_charge=0,
        tautomer_of=None,
        ph_assumption=None,
        pka_sites=(),
        enumeration_method="Batch 7 reference-neutral state",
        enumeration_parameters=(),
        assumptions=prepared.assumptions,
        atoms=prepared.atoms,
        reference_atom_mapping=tuple(
            (atom.atom_id, reference_indices[atom.atom_id])
            for atom in prepared.atoms
        ),
        rings=prepared.rings,
        water_roles=prepared.water_roles,
    )


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (
        args.manifest,
        args.identity,
        args.ligand,
        args.site,
        args.preparation,
    ):
        if not path.exists():
            raise SystemExit(f"missing {path}; run the earlier capstone batches first.")

    manifest = load_manifest(args.manifest)
    identity = load_ligand_identity(args.identity)
    preparation = load_chemical_preparation(args.preparation)
    reference = _reference_state(preparation, identity)

    smiles = identity.reference.smiles
    molecule = Chem.MolFromSmiles(smiles)
    if molecule is None:
        raise SystemExit("reference SMILES failed to parse.")

    enumerator = rdMolStandardize.TautomerEnumerator()
    tautomers = []
    seen = {Chem.MolToSmiles(molecule, canonical=True)}
    for index, tautomer in enumerate(enumerator.Enumerate(molecule)):
        canonical = Chem.MolToSmiles(tautomer, canonical=True)
        if canonical in seen:
            continue
        seen.add(canonical)
        tautomers.append(
            _prepared_state_from_molecule(
                tautomer,
                state_id=f"tautomer-{len(tautomers) + 1}",
                label=f"RDKit tautomer {len(tautomers) + 1}",
                kind="tautomer",
                parent_state_id="reference-neutral",
                tautomer_of="reference-neutral",
                pka_sites=(),
                enumeration_method="RDKit TautomerEnumerator",
                enumeration_parameters=(
                    "enumerate all returned tautomers",
                    "ETKDGv3 randomSeed=0xF00D",
                ),
                assumptions=(
                    "Same formula and proton count as the reference state.",
                    "Tautomer enumeration is not a population estimate.",
                    "Generated conformer is not aligned to experimental coordinates.",
                ),
            )
        )

    factory = ChemicalFeatures.BuildFeatureFactory(
        str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef")
    )
    ionizable = [
        feature.GetAtomIds()[0]
        for feature in factory.GetFeaturesForMol(molecule)
        if feature.GetFamily() == "PosIonizable"
    ]
    protonated = []
    for index, atom_index in enumerate(ionizable):
        editable = Chem.RWMol(molecule)
        atom = editable.GetAtomWithIdx(atom_index)
        atom.SetFormalCharge(atom.GetFormalCharge() + 1)
        atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
        protonated_molecule = editable.GetMol()
        Chem.SanitizeMol(protonated_molecule)
        protonated.append(
            _prepared_state_from_molecule(
                protonated_molecule,
                state_id=f"protonation-{index + 1}",
                label=f"Hypothetical protonation at reference atom {atom_index}",
                kind="protonation",
                parent_state_id="reference-neutral",
                tautomer_of=None,
                pka_sites=(
                    IonizableSite(
                        site_id=f"site-{index + 1}",
                        atom_id=f"reference-atom-{atom_index}",
                        site_type="base",
                        pka=None,
                        pka_source="not supplied; hypothetical enumeration",
                        pka_uncertainty=None,
                        note="No measured or modelled pKa supplied; no population estimated.",
                    ),
                ),
                enumeration_method="explicit protonation of one RDKit PosIonizable site",
                enumeration_parameters=(
                    "BaseFeatures.fdef PosIonizable site selection",
                    "ETKDGv3 randomSeed=0xF00D",
                ),
                assumptions=(
                    "Hypothetical protonation alternative, not a dominant species.",
                    "No pKa value, pH condition, or population model was supplied.",
                    "Generated conformer is not aligned to experimental coordinates.",
                ),
            )
        )

    states = (reference, *tautomers, *protonated)
    record = build_chemical_state_enumeration_record(
        source_manifest_path=str(args.manifest),
        source_manifest_sha256=sha256_file(args.manifest),
        ligand_identity_record_path=str(args.identity),
        ligand_identity_record_sha256=sha256_file(args.identity),
        ligand_coordinate_record_path=str(args.ligand),
        ligand_coordinate_record_sha256=sha256_file(args.ligand),
        binding_site_record_path=str(args.site),
        binding_site_record_sha256=sha256_file(args.site),
        chemical_preparation_record_path=str(args.preparation),
        chemical_preparation_record_sha256=sha256_file(args.preparation),
        accession=manifest.source.accession,
        selection=identity.selection,
        reference_state_id=reference.state_id,
        states=states,
        toolkit="RDKit",
        toolkit_version=rdBase.rdkitVersion,
        parameters=(
            "RDKit TautomerEnumerator for tautomer alternatives",
            "RDKit BaseFeatures.fdef PosIonizable sites for hypothetical protonation",
            "ETKDGv3 conformer generation with fixed random seed",
            "no pKa model, pH population model, or affinity estimate",
        ),
        selected_state_id=None,
        ph_assumption=None,
        population_model=None,
    )
    errors = validate_chemical_state_enumeration(record)
    if errors:
        for error in errors:
            print(f"state record error: {error}")
        raise SystemExit(2)

    destination = write_chemical_state_enumeration(args.output, record)
    kinds = Counter(state.kind for state in record.states)
    print(f"{record.accession} chemical-state enumeration")
    print(f"  states: {len(record.states)} ({dict(kinds)})")
    print(f"  selected_state_id: {record.selected_state_id!r}")
    print(f"  pH assumption: {record.ph_assumption!r}")
    print(f"  population model: {record.population_model!r}")
    for state in record.states:
        print(
            f"    {state.state_id}: {state.kind}, charge={state.formal_charge:+d}, "
            f"protons={state.proton_count}, graph={state.molecular_graph[:60]}"
        )
    print(f"  toolkit: {record.toolkit} {record.toolkit_version}")
    print(f"  record: {destination}")
    print("  scope: alternatives coexist; no state or population was selected")


if __name__ == "__main__":
    main()
