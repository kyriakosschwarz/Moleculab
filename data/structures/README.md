# Experimental-structure cache

This directory holds a local, unmodified copy of a deposited coordinate file
and its provenance manifest.

- `2HYY.cif` is downloaded from RCSB PDB and is **not** committed. The manifest
  records its SHA-256, size, retrieval date, and entry revision.
- `2HYY.manifest.json` is committed and is the authoritative selection record
  for the capstone.
- Do not edit a downloaded file in place. If the checksum no longer matches the
  manifest, the workflow stops instead of silently replacing the provenance
  record.

Regenerate the manifest with:

```bash
uv run --with gemmi python examples/drug_discovery/capstone_structure.py
```
