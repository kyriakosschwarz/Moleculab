# Experimental-structure cache

This directory holds a local, unmodified copy of a deposited coordinate file
and its provenance manifest.

- `2HYY.cif` is downloaded from RCSB PDB and is **not** committed. The manifest
  records its SHA-256, size, retrieval date, and entry revision.
- `2HYY.manifest.json` is committed and is the authoritative selection record
  for the capstone.
- `2HYY.ligand-identity.json` is the Batch 2 reconciliation record. It preserves
  the deposited heavy-atom graph, the reference graph, the explicit atom
  mapping, and every discrepancy.
- `2HYY.ligand-coordinates.json` is the Batch 3 coordinate record. It preserves
  coordinates, occupancy, displacement parameters, alternate locations, bond
  lengths, angles, ring planarity, and close contacts.
- `2HYY.binding-site.json` is the Batch 4 binding-site record. It preserves the
  pocket cutoff, residue identities, atom-level contacts, nearby waters, other
  hetero residues, missing residue numbers, and alternate locations.
- `2HYY.interactions.json` is the Batch 5 interaction interpretation. It records
  classified interactions only where donor/acceptor, charge, hydrophobic, ring,
  water, or metal evidence is explicit; unknown chemistry remains unknown.
- `2HYY.evidence-fingerprint.json` is the Batch 6 evidence fingerprint. It keeps
  present, absent, and unknown distinct and uses only assessable categories in
  occupancy denominators.
- `2HYY.chemical-preparation.json` is the Batch 7 preparation record. It holds
  coexisting deposited-only and reference-neutral states without silently
  selecting one or modifying the earlier experimental records.
- Do not edit a downloaded file in place. If the checksum no longer matches the
  manifest, the workflow stops instead of silently replacing the provenance
  record.

Regenerate the records in order:

```bash
uv run --with gemmi python examples/drug_discovery/capstone_structure.py
uv run --with gemmi --with rdkit python examples/drug_discovery/capstone_ligand_identity.py
uv run --with gemmi python examples/drug_discovery/capstone_ligand_coordinates.py
uv run --with gemmi python examples/drug_discovery/capstone_binding_site.py
uv run --with gemmi python examples/drug_discovery/capstone_interactions.py
uv run python examples/drug_discovery/capstone_evidence_fingerprint.py
uv run --with rdkit --with gemmi python examples/drug_discovery/capstone_chemical_preparation.py
```
