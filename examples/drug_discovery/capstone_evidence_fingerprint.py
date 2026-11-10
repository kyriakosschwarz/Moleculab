"""Build a three-state interaction-evidence fingerprint.

This is Batch 6 of the capstone.  It consumes the Batch 5 interaction
interpretation and records category-level evidence as present, absent, or
unknown.  Unknown categories are excluded from assessable denominators rather
than being encoded as zero.

It does not compute a similarity score, interaction strength, or affinity.

Run:

    uv run python examples/drug_discovery/capstone_evidence_fingerprint.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

from moleculab import (
    build_interaction_evidence,
    load_interactions,
    sha256_file,
    validate_interaction_evidence,
    write_interaction_evidence,
)

DEFAULT_MANIFEST = Path("data/structures/2HYY.manifest.json")
DEFAULT_IDENTITY = Path("data/structures/2HYY.ligand-identity.json")
DEFAULT_SITE = Path("data/structures/2HYY.binding-site.json")
DEFAULT_INTERACTIONS = Path("data/structures/2HYY.interactions.json")
DEFAULT_OUTPUT = Path("data/structures/2HYY.evidence-fingerprint.json")


def _argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--identity", type=Path, default=DEFAULT_IDENTITY)
    parser.add_argument("--site", type=Path, default=DEFAULT_SITE)
    parser.add_argument("--interactions", type=Path, default=DEFAULT_INTERACTIONS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser


def main() -> None:
    args = _argument_parser().parse_args()
    for path in (args.manifest, args.identity, args.site, args.interactions):
        if not path.exists():
            raise SystemExit(
                f"missing {path}; run the earlier capstone batches first."
            )

    interactions = load_interactions(args.interactions)
    record = build_interaction_evidence(
        interaction_record=interactions,
        source_manifest_path=str(args.manifest),
        source_manifest_sha256=sha256_file(args.manifest),
        ligand_identity_record_path=str(args.identity),
        ligand_identity_record_sha256=sha256_file(args.identity),
        binding_site_record_path=str(args.site),
        binding_site_record_sha256=sha256_file(args.site),
        interaction_record_path=str(args.interactions),
        interaction_record_sha256=sha256_file(args.interactions),
    )
    errors = validate_interaction_evidence(record)
    if errors:
        for error in errors:
            print(f"evidence record error: {error}")
        raise SystemExit(2)

    destination = write_interaction_evidence(args.output, record)
    fingerprint = record.fingerprint
    print(f"{record.accession} interaction evidence fingerprint")
    print(
        f"  features: {len(fingerprint.schema.features)} category features "
        f"({fingerprint.present_count} present, {fingerprint.absent_count} absent, "
        f"{fingerprint.unknown_count} unknown)"
    )
    print(
        f"  assessable denominator: {fingerprint.assessable_count} of "
        f"{len(fingerprint.schema.features)} categories"
    )
    print("  category statuses:")
    for interaction_type, status in record.category_statuses:
        print(f"    {interaction_type}: {status}")
    print(f"  schema: {record.schema.schema_id}")
    print(f"  criteria fingerprint: {record.schema.criteria_fingerprint}")
    print(f"  record: {destination}")
    print("  scope: unknown evidence is not zero; no affinity or similarity score is computed")


if __name__ == "__main__":
    main()
