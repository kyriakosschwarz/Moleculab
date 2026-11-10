"""Retrieve 2HYY, preserve the deposited file, and emit a provenance manifest.

This is batch 1 of the capstone: experimental-structure retrieval and
identification.  It records the exact file, checksum, revision, protein chain,
and ligand instance before any binding-site or interaction analysis is run.

It does not dock a ligand, predict affinity, or estimate a binding free energy.

Run with gemmi available:

    uv run --with gemmi python examples/drug_discovery/capstone_structure.py
"""

from __future__ import annotations

import argparse
import datetime as dt
import urllib.request
from pathlib import Path

from moleculab import (
    AnalysisRecord,
    ChemicalIdentityRecord,
    DecisionRecord,
    LigandSelection,
    ProcessingRecord,
    ProteinSelection,
    ProvenanceEntry,
    ProvenanceRecord,
    StructureManifest,
    StructureSelection,
    StructureSource,
    ValidationCheck,
    ValidationRecord,
    load_manifest,
    sha256_file,
    validate_manifest,
    write_manifest,
)

try:
    import gemmi
except ImportError as exc:  # pragma: no cover - depends on the environment
    raise SystemExit(
        "gemmi is required for this example. Run it with: "
        "uv run --with gemmi python examples/drug_discovery/capstone_structure.py"
    ) from exc

DEFAULT_ACCESSION = "2HYY"
DEFAULT_SOURCE = Path("data/structures/2HYY.cif")
DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_PROTEIN_AUTH_CHAIN = "A"
DEFAULT_PROTEIN_LABEL_CHAIN = "A"
DEFAULT_LIGAND_COMPONENT = "STI"
DEFAULT_LIGAND_LABEL_CHAIN = "E"
REFERENCE = "imatinib"
REFERENCE_FORMULA = "C29H31N7O"


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--accession", default=DEFAULT_ACCESSION)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument(
        "--download-url",
        default=None,
        help="Defaults to the RCSB mmCIF download URL for --accession.",
    )
    parser.add_argument("--protein-auth-chain", default=DEFAULT_PROTEIN_AUTH_CHAIN)
    parser.add_argument("--protein-label-chain", default=DEFAULT_PROTEIN_LABEL_CHAIN)
    parser.add_argument("--ligand-component", default=DEFAULT_LIGAND_COMPONENT)
    parser.add_argument("--ligand-label-chain", default=DEFAULT_LIGAND_LABEL_CHAIN)
    return parser


def _retrieve(source: Path, download_url: str) -> str:
    """Download *source* once; never overwrite an existing deposited file."""

    if source.exists():
        return "reused"

    source.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        download_url,
        headers={"User-Agent": "moleculab-capstone/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        payload = response.read()
    if not payload.startswith(b"data_"):
        raise ValueError(f"{download_url} did not return an mmCIF data block.")

    try:
        with source.open("xb") as handle:
            handle.write(payload)
    except FileExistsError:
        return "reused"
    return "downloaded"


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


def _formula_key(formula: str) -> str:
    return "".join(_unquote(formula).replace(" ", "").split())


def _selected_ligands(model: gemmi.Model, component: str) -> list[tuple[str, str, gemmi.Residue]]:
    ligands: list[tuple[str, str, gemmi.Residue]] = []
    for chain in model:
        for residue in chain:
            if residue.name == component:
                ligands.append((chain.name, residue.subchain, residue))
    return ligands


def _atom_site_selection(
    block: gemmi.cif.Block,
    component_id: str,
    label_asym_id: str,
) -> list[str]:
    """Return the atom-site fields that identify the selected residue."""

    rows = block.find(
        "_atom_site.",
        [
            "label_comp_id",
            "label_asym_id",
            "label_seq_id",
            "pdbx_PDB_ins_code",
            "pdbx_PDB_model_num",
        ],
    )
    for row in rows:
        values = [_unquote(value) for value in row]
        if values[0] == component_id and values[1] == label_asym_id:
            return values
    raise ValueError(
        f"no _atom_site record for component {component_id!r} label chain "
        f"{label_asym_id!r}."
    )


def _altloc_count(residues: list[gemmi.Residue]) -> int:
    count = 0
    for residue in residues:
        for atom in residue:
            if atom.altloc not in ("", "\x00"):
                count += 1
    return count


def _check(
    category: str,
    check_id: str,
    ok: bool,
    *,
    subject: str,
    criterion: str,
    observed: str,
    expected: str,
    detail: str,
    not_run: bool = False,
) -> ValidationCheck:
    status = "not_run" if not_run else ("pass" if ok else "fail")
    return ValidationCheck(
        check_id=check_id,
        category=category,
        status=status,
        subject=subject,
        criterion=criterion,
        observed=observed,
        expected=expected,
        detail=detail,
    )


def build_manifest(args: argparse.Namespace) -> tuple[StructureManifest, str]:
    download_url = args.download_url or (
        f"https://files.rcsb.org/download/{args.accession}.cif"
    )
    args.download_url = download_url
    retrieval = _retrieve(args.source, download_url)
    checksum = sha256_file(args.source)
    file_size = args.source.stat().st_size
    if args.manifest.exists():
        previous = load_manifest(args.manifest)
        if previous.source.sha256 != checksum:
            raise ValueError(
                f"{args.manifest} records SHA-256 {previous.source.sha256}, but "
                f"{args.source} now has {checksum}. Refusing to replace a "
                "provenance record with a different source file."
            )
    retrieval_date = (
        dt.datetime.fromtimestamp(args.source.stat().st_mtime, tz=dt.timezone.utc)
        .date()
        .isoformat()
    )

    document = gemmi.cif.read_file(str(args.source))
    block = document.sole_block()
    structure = gemmi.read_structure(str(args.source))
    model = structure[0]
    parsed_format = structure.input_format.name
    coordinate_format = "mmCIF" if "mmcif" in parsed_format.lower() else parsed_format

    entry_id = _value(block, "_entry.id")
    method = _value(block, "_exptl.method")
    resolution_text = _value(block, "_refine.ls_d_res_high")
    resolution = float(resolution_text) if resolution_text is not None else None

    revision_rows = [
        [value for value in row]
        for row in block.find(
            "_pdbx_audit_revision_history.",
            ["major_revision", "minor_revision", "revision_date"],
        )
    ]
    if revision_rows:
        major, minor, revision_date = revision_rows[-1]
        revision = f"{major}.{minor}"
    else:
        revision = None
        revision_date = None

    chains = {chain.name: chain for chain in model}
    protein_chain = chains.get(args.protein_auth_chain)
    if protein_chain is None:
        raise ValueError(
            f"protein auth chain {args.protein_auth_chain!r} is absent from {args.source}."
        )
    polymer_residues = [residue for residue in protein_chain if residue.het_flag == "A"]
    if not polymer_residues:
        raise ValueError(
            f"protein chain {args.protein_auth_chain!r} contains no polymer residues."
        )
    protein_label = polymer_residues[0].subchain or args.protein_label_chain

    ligands = _selected_ligands(model, args.ligand_component)
    if not ligands:
        raise ValueError(
            f"ligand component {args.ligand_component!r} is absent from {args.source}."
        )
    selected = [item for item in ligands if item[1] == args.ligand_label_chain]
    if len(selected) != 1:
        raise ValueError(
            f"ligand instance {args.ligand_label_chain!r} matched {len(selected)} "
            "deposited copies; select exactly one."
        )
    ligand_auth_chain, ligand_label, ligand_residue = selected[0]
    available = sorted({label for _, label, _ in ligands})
    atom_site = _atom_site_selection(block, args.ligand_component, ligand_label)
    label_seq_id = None if atom_site[2] in ("?", ".") else atom_site[2]
    insertion_code = None if atom_site[3] in ("?", ".") else atom_site[3]
    alternate_location = None
    model_number = int(atom_site[4]) if atom_site[4] not in ("?", ".") else 1
    auth_seq_id = str(ligand_residue.seqid.num)

    chem_rows = [
        [value for value in row]
        for row in block.find(
            "_chem_comp.",
            ["id", "name", "formula", "formula_weight"],
        )
    ]
    chem_row = next((row for row in chem_rows if row[0] == args.ligand_component), None)
    if chem_row is None:
        raise ValueError(
            f"_chem_comp has no record for component {args.ligand_component!r}."
        )
    component_id, component_name, formula_text, _formula_weight = (
        _unquote(value) for value in chem_row
    )
    formula = _formula_key(formula_text)

    source = StructureSource(
        database="RCSB PDB",
        accession=entry_id or args.accession,
        record_url=f"https://www.rcsb.org/structure/{args.accession}",
        download_url=args.download_url,
        retrieval_date=retrieval_date,
        revision=revision,
        revision_date=revision_date,
        experimental_method=method,
        resolution_angstrom=resolution,
        coordinate_format=coordinate_format,
        file_name=args.source.name,
        sha256=checksum,
        file_size_bytes=file_size,
        deposited_asymmetric_unit=True,
        biological_assembly=None,
    )
    selection = StructureSelection(
        model_number=model_number,
        protein=ProteinSelection(
            label_asym_id=protein_label,
            auth_asym_id=args.protein_auth_chain,
        ),
        ligand=LigandSelection(
            component_id=component_id,
            label_asym_id=ligand_label,
            auth_asym_id=ligand_auth_chain,
            instance=ligand_label,
            auth_seq_id=auth_seq_id,
            label_seq_id=label_seq_id,
            insertion_code=insertion_code,
            alternate_location=alternate_location,
            available_instances=tuple(available),
        ),
    )
    processing = ProcessingRecord(
        coordinate_transforms=(),
        removed_waters=(),
        removed_ions=(),
        removed_cofactors=(),
        alternate_location_policy=(
            "No selection applied; deposited alternate locations remain unchanged "
            "in the cached source file."
        ),
        hydrogen_addition="None; deposited atoms retained as-is.",
        waters="retained",
        ions_cofactors="retained",
        missing_atoms_or_residues=(
            "Not assessed in batch 1; no atoms or residues were rebuilt or removed."
        ),
        assembly_id=None,
        assembly_transformations_applied=False,
        original_coordinate_frame_preserved=True,
        translation_applied=None,
        rotation_applied=None,
    )
    decisions = (
        DecisionRecord(
            decision="selected one protein chain and one ligand copy",
            actor="capstone_structure.py command configuration",
            rationale=(
                f"protein auth chain {args.protein_auth_chain} and ligand "
                f"component {component_id} label chain {ligand_label}; the other "
                "deposited copies remain in the source file but are excluded from "
                "this selection record."
            ),
        ),
        DecisionRecord(
            decision="retained the deposited asymmetric unit without transforming coordinates",
            actor="capstone_structure.py",
            rationale=(
                "no biological assembly expansion, translation, or rotation was "
                "applied; the deposited coordinate frame is preserved."
            ),
        ),
        DecisionRecord(
            decision="deferred formal charge, stereochemistry, and atom mapping",
            actor="capstone_structure.py",
            rationale=(
                "these chemical-identity decisions require the ligand "
                "reconciliation workflow planned for batch 2."
            ),
        ),
    )
    chemical_identity = ChemicalIdentityRecord(
        reference=REFERENCE,
        component_id=component_id,
        name=component_name,
        formula=formula,
        formal_charge=None,
        stereochemistry=None,
        atom_mapping=None,
        notes=(
            "The deposited _chem_comp category does not provide a formal charge, "
            "and atom_site formal charges are not populated. No charge or "
            "stereochemistry assignment was applied in batch 1."
        ),
    )
    analysis = AnalysisRecord(
        coordinate_unit="angstrom",
        internal_distance_unit="nm",
        force_field=None,
        charge_model=None,
        interaction_criteria=None,
    )

    ligand_altlocs = _altloc_count([ligand_residue])
    protein_altlocs = _altloc_count(polymer_residues)
    checks = [
        _check(
            "source_integrity",
            "source.accession_matches_deposition",
            entry_id == args.accession,
            subject=args.source.name,
            criterion="_entry.id equals the requested accession",
            observed=entry_id or "missing",
            expected=args.accession,
            detail=f"requested {args.accession}; file _entry.id is {entry_id}",
        ),
        _check(
            "source_integrity",
            "source.checksum_recorded",
            len(checksum) == 64,
            subject=args.source.name,
            criterion="SHA-256 digest is a lowercase 64-character hex string",
            observed=checksum,
            expected="64 lowercase hexadecimal characters",
            detail="the digest is recorded for later file-identity checks",
        ),
        _check(
            "source_integrity",
            "source.coordinate_format",
            coordinate_format == "mmCIF",
            subject=args.source.name,
            criterion="the deposited file is parsed as mmCIF",
            observed=coordinate_format,
            expected="mmCIF",
            detail=f"parsed as {coordinate_format}",
        ),
        _check(
            "source_integrity",
            "source.revision_recorded",
            revision is not None and revision_date is not None,
            subject=args.source.name,
            criterion="entry revision and revision date are present",
            observed=f"{revision} ({revision_date})",
            expected="revision with ISO date",
            detail=f"revision {revision} dated {revision_date}",
        ),
        _check(
            "source_integrity",
            "source.experimental_method",
            method is not None,
            subject=args.source.name,
            criterion="experimental method is present",
            observed=method or "missing",
            expected="not null",
            detail=f"method {method}",
        ),
        _check(
            "source_integrity",
            "source.resolution_recorded",
            resolution is not None and resolution > 0,
            subject=args.source.name,
            criterion="deposited resolution is positive",
            observed=f"{resolution} angstrom" if resolution is not None else "missing",
            expected="> 0 angstrom",
            detail=f"resolution {resolution} angstrom",
        ),
        _check(
            "source_integrity",
            "source.structure_scope",
            True,
            subject=args.source.name,
            criterion="coordinate scope is stated explicitly",
            observed="deposited asymmetric unit; no transformation",
            expected="deposited asymmetric unit",
            detail="no biological assembly expansion, translation, or rotation was applied",
        ),
        _check(
            "structural_integrity",
            "selection.protein_chain_present",
            True,
            subject=f"{args.protein_auth_chain}:{protein_label}",
            criterion="selected chain exists and contains polymer residues",
            observed=f"{len(polymer_residues)} polymer residues",
            expected="at least one polymer residue",
            detail=(
                f"auth chain {args.protein_auth_chain}, label chain "
                f"{protein_label}, {len(polymer_residues)} polymer residues"
            ),
        ),
        _check(
            "structural_integrity",
            "selection.ligand_instance_present",
            True,
            subject=f"{component_id}:{ligand_label}",
            criterion="selected ligand instance exists exactly once",
            observed=(
                f"{len(selected)} copy; auth chain {ligand_auth_chain}, "
                f"auth seq {auth_seq_id}"
            ),
            expected="exactly one copy",
            detail=(
                f"component {component_id} label chain {ligand_label}, "
                f"auth chain {ligand_auth_chain}, auth seq {auth_seq_id}"
            ),
        ),
        _check(
            "structural_integrity",
            "selection.single_ligand_copy_selected",
            len(selected) == 1,
            subject=f"{component_id}:{ligand_label}",
            criterion="the selection names one deposited ligand copy",
            observed=f"{len(selected)} matching copy",
            expected="1",
            detail=(
                f"available label chains {', '.join(available)}; selected "
                f"{ligand_label}"
            ),
        ),
        _check(
            "structural_integrity",
            "selection.ligand_associated_with_protein_chain",
            ligand_auth_chain == args.protein_auth_chain,
            subject=f"{component_id}:{ligand_label}",
            criterion="ligand author chain matches the selected protein chain",
            observed=ligand_auth_chain,
            expected=args.protein_auth_chain,
            detail=(
                f"ligand auth chain {ligand_auth_chain}, protein auth chain "
                f"{args.protein_auth_chain}"
            ),
        ),
        _check(
            "structural_integrity",
            "selection.ligand_atoms_present",
            len(ligand_residue) > 0,
            subject=f"{component_id}:{ligand_label}",
            criterion="selected ligand has deposited atoms",
            observed=str(len(ligand_residue)),
            expected="> 0",
            detail=f"{len(ligand_residue)} deposited ligand atoms",
        ),
        _check(
            "structural_integrity",
            "selection.alternate_locations_recorded",
            True,
            subject=f"{component_id}:{ligand_label}",
            criterion="alternate locations are recorded and not discarded",
            observed=f"ligand {ligand_altlocs}, protein {protein_altlocs} altloc atoms",
            expected="recorded; none discarded",
            detail=(
                f"ligand altloc atoms {ligand_altlocs}, protein altloc atoms "
                f"{protein_altlocs}; no alternate location selected"
            ),
        ),
        _check(
            "chemical_identity",
            "chemical.formula_matches_reference",
            formula == REFERENCE_FORMULA,
            subject=f"{component_id}:{ligand_label}",
            criterion="deposited formula equals the intended imatinib formula",
            observed=formula,
            expected=REFERENCE_FORMULA,
            detail=(
                f"deposited {formula_text} normalized to {formula}; reference "
                f"{REFERENCE_FORMULA}"
            ),
        ),
        _check(
            "chemical_identity",
            "chemical.formal_charge_assessed",
            False,
            subject=f"{component_id}:{ligand_label}",
            criterion="formal charge is compared with the reference",
            observed="not provided",
            expected="reference formal charge or documented assignment method",
            detail=(
                "no formal charge is provided in _chem_comp and atom_site "
                "charges are unpopulated"
            ),
            not_run=True,
        ),
        _check(
            "chemical_identity",
            "chemical.stereochemistry_assessed",
            False,
            subject=f"{component_id}:{ligand_label}",
            criterion="stereochemistry and atom mapping are compared with the reference",
            observed="not assessed",
            expected="reference stereochemistry comparison",
            detail="stereochemistry and atom mapping are deferred to batch 2",
            not_run=True,
        ),
        _check(
            "analysis_integrity",
            "analysis.coordinate_conversion_boundary",
            True,
            subject=args.source.name,
            criterion="deposited coordinates remain in angstrom until analysis",
            observed="source angstrom; internal nm",
            expected="convert only at the analysis boundary",
            detail=(
                "deposited coordinates remain in angstrom; conversion to nm is "
                "deferred to the analysis boundary"
            ),
        ),
        _check(
            "analysis_integrity",
            "analysis.interaction_criteria_defined",
            False,
            subject="interaction analysis",
            criterion="distance and angle criteria are defined",
            observed="not defined",
            expected="explicit distance and angle criteria",
            detail="no interaction criteria are defined in the retrieval batch",
            not_run=True,
        ),
        _check(
            "analysis_integrity",
            "analysis.docking_or_affinity_computed",
            False,
            subject="docking and affinity",
            criterion="docking or binding free energy is computed",
            observed="not computed",
            expected="out of scope for this milestone",
            detail="not computed; out of scope for this milestone",
            not_run=True,
        ),
    ]

    provenance = ProvenanceRecord(
        records=(
            ProvenanceEntry(
                "experimental",
                "coordinate file",
                f"RCSB PDB {source.accession} revision {revision} dated "
                f"{revision_date}; downloaded from {args.download_url}.",
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "experimental",
                "ligand component description",
                f"Component {component_id} with deposited formula {formula_text} "
                f"and name {component_name}.",
                actor="RCSB PDB deposition",
            ),
            ProvenanceEntry(
                "calculated",
                "file integrity",
                f"SHA-256 {checksum}; {file_size} bytes; retrieval date "
                f"{retrieval_date}.",
                actor="capstone_structure.py",
            ),
            ProvenanceEntry(
                "predicted",
                "none",
                "No docking pose, affinity, or free-energy estimate is produced "
                "in this batch.",
                actor="capstone_structure.py",
            ),
            ProvenanceEntry(
                "hypothesised",
                "none",
                "No mechanistic interpretation is made from the selected "
                "coordinates in this batch.",
                actor="capstone_structure.py",
            ),
        )
    )
    manifest = StructureManifest(
        source=source,
        selection=selection,
        processing=processing,
        decisions=decisions,
        chemical_identity=chemical_identity,
        analysis=analysis,
        provenance=provenance,
        validation=ValidationRecord(checks=tuple(checks)),
    )
    return manifest, retrieval


def main() -> None:
    args = _argument_parser().parse_args()
    manifest, retrieval = build_manifest(args)

    errors = validate_manifest(manifest)
    if errors:
        for error in errors:
            print(f"manifest error: {error}")
        raise SystemExit(2)

    destination = write_manifest(args.manifest, manifest)
    statuses = [check.status for check in manifest.validation.checks]
    failures = [check for check in manifest.validation.checks if check.status == "fail"]

    print(f"{manifest.source.accession} source manifest")
    print(f"  retrieval: {retrieval} {args.source}")
    print(f"  sha256: {manifest.source.sha256}")
    print(f"  revision: {manifest.source.revision} ({manifest.source.revision_date})")
    print(
        "  selection: protein auth {} / label {}; ligand {} auth {} / label {} "
        "instance {}, auth seq {}".format(
            manifest.selection.protein.auth_asym_id,
            manifest.selection.protein.label_asym_id,
            manifest.selection.ligand.component_id,
            manifest.selection.ligand.auth_asym_id,
            manifest.selection.ligand.label_asym_id,
            manifest.selection.ligand.instance,
            manifest.selection.ligand.auth_seq_id,
        )
    )
    print(
        f"  validation: {statuses.count('pass')} pass, "
        f"{statuses.count('not_run')} not run, {statuses.count('fail')} fail"
    )
    print(f"  manifest: {destination}")
    print("  scope: no docking score, affinity, or binding free energy is computed")

    if failures:
        for check in failures:
            print(f"  FAILED {check.check_id}: {check.detail}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
