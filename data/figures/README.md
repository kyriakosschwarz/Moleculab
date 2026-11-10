# Figure metadata

`capstone-visualizations.json` describes the Visualization Batch 1 figures
stored in `docs/tutorials/drug_discovery/visualizations/`.

Each entry records:

- figure ID, title, and image path;
- source structure, source checksum, protein chain, ligand instance, model, and
  coordinate units;
- visualization tool and exact tool version;
- generated-from records;
- scientific limitations, caption, and alt text;
- the PyMOL camera view when applicable.

The metadata is rewritten by:

```bash
uv run --with rdkit python examples/drug_discovery/capstone_visualizations.py
```

The figures are derived from the deposited structure and validated chemical
records. They do not add interaction classification, energetics, or affinity.
