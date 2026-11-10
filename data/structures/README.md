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
- `2HYY.chemical-states.json` is the Batch 8 enumeration. It holds reference,
  tautomer, and hypothetical protonation alternatives without selecting a state
  or estimating populations.
- `2HYY.state-interactions.json` is the Batch 9 state-aware assessment. It maps
  each state to the same experimental contacts, records state-specific feature
  differences, and leaves interaction classification unknown where protein
  chemistry is unavailable.
- `2HYY.protein-preparation.json` is the Batch 10 protein-pocket preparation. It
  holds explicit residue-template states with localized formal charges,
  donor/acceptor/hydrophobic/aromatic features, unassigned water/hetero roles,
  and no generated protein hydrogens.
- `2HYY.combination-interactions.json` is the Batch 11 matrix of all 16
  ligand/protein state combinations on the same experimental contacts, with
  hypothesis-conditioned interaction classification and no ranking.
- `2HYY.hydrogen-evidence.json` is the Batch 12A hydrogen record. It separates
  generated ligand hydrogens from deposited heavy atoms, validates donor–hydrogen
  connectivity, and reports D–H···A geometry. Protein hydrogens are `not_run`
  because Batch 10 generated none, so protein-donor hydrogen bonds remain
  unknown.
- `2HYY.water-evidence.json` is the Batch 12B water record. It preserves
  deposited water identity, coordinates, occupancy, and B-factors, leaves
  orientation and donor/acceptor roles unassigned, and reports geometric bridge
  plausibility without classifying a water-mediated interaction.
- `2HYY.interaction-robustness.json` is the Batch 13 sensitivity record. It
  varies contact cutoffs, hydrogen-bond angles, generated-hydrogen orientation,
  ligand/protein states, and water cutoffs independently, then labels each
  conclusion invariant, assumption-sensitive, unsupported, or unknown. It
  contains no composite confidence score, ranking, or affinity estimate.
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
uv run --with rdkit --with gemmi python examples/drug_discovery/capstone_chemical_states.py
uv run python examples/drug_discovery/capstone_state_interactions.py
uv run --with gemmi python examples/drug_discovery/capstone_protein_preparation.py
uv run --with gemmi python examples/drug_discovery/capstone_combination_interactions.py
uv run python examples/drug_discovery/capstone_hydrogen_evidence.py
uv run --with gemmi python examples/drug_discovery/capstone_water_evidence.py
uv run python examples/drug_discovery/capstone_interaction_robustness.py
```
