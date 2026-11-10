# Figure metadata

`capstone-visualizations.json` describes the Visualization Batch 1–3 figures
stored in `docs/tutorials/drug_discovery/visualizations/`.

Each entry records:

- figure ID, title, and image path;
- source structure, source checksum, protein chain, ligand instance, model, and
  coordinate units;
- visualization tool and exact tool version;
- per-record source paths and checksums;
- the rendering configuration and its SHA-256 fingerprint;
- scientific limitations, caption, and alt text;
- the PyMOL camera view when applicable.

The camera and rendering fingerprint make the visual settings inspectable.
Scientific reproducibility means the same atoms, coordinates, and measurements
are represented; pixel identity is not guaranteed across rendering
environments.

The code-level registry in `src/moleculab/visualization/registry.py` defines
the authoritative figure IDs, paths, generation functions, source records, and
teaching anchors. The generated gallery at
`docs/tutorials/drug_discovery/visualizations/gallery.md` is built from that
registry and the metadata.

The metadata and gallery are rewritten by:

```bash
uv run --with rdkit python examples/drug_discovery/capstone_visualizations.py
```

The figures are derived from the deposited structure and validated chemical
records. They do not add interaction classification, energetics, or affinity.
