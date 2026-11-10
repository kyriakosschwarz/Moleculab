# Capstone visualizations

These figures are generated directly from the 2HYY provenance chain and the
validated chemical-state records. They are not AI-generated molecular images
and they do not perform new chemical interpretation.

The companion metadata is
[`data/figures/capstone-visualizations.json`](../../../../data/figures/capstone-visualizations.json).
It records the source structure checksum, protein chain, ligand instance, model
number, tool versions, camera view, generated-from records, alt text, and
scientific limitations.

## Figure 1 — Imatinib chemical identity

![Reference-neutral imatinib with donor and acceptor atoms highlighted](2HYY-imatinib-2d.svg)

The RDKit diagram shows the reference-neutral imatinib graph:

- blue: hydrogen-bond donor atoms;
- orange: hydrogen-bond acceptor atoms;
- purple: atoms that are both;
- atom labels: reference-graph indices;
- hydrogens: omitted.

**Teaches:** molecular graph, functional groups, donor/acceptor features, and
chemical-state identity.

**Does not establish:** a conformer, a protonation population, an interaction
pattern, an energy, or an affinity. The deposited atom names are different from
the reference indices shown here.

## Figure 2 — Imatinib in the ABL1 binding site

![Imatinib in the ABL1 binding site](2HYY-imatinib-abl1-overview.png)

The PyMOL overview shows the selected protein chain as a cartoon, imatinib as
yellow sticks, protein residues within 4 Å of the ligand as a pale surface, and
nearby waters as blue spheres.

**Teaches:** protein structure, binding-site location, and the relationship
between ligand and pocket.

**Does not establish:** binding affinity, binding-mode uniqueness, an
equilibrium ensemble, or a validated pose. The figure uses one crystallographic
model and one selected ligand copy.

## Figure 3 — Binding-pocket contacts

![Imatinib binding-pocket contacts](2HYY-imatinib-binding-pocket.png)

The close-up shows cyan contact residues, the yellow ligand, nearby waters, and
gray dashed lines for three measured donor–acceptor distances. The selected
contacts correspond to the N13···THR315, N21···GLU286, and N21···ASP381
geometries used in the hydrogen-evidence record.

**Teaches:** measured geometry, contact pairs, and the difference between
proximity and a classified interaction.

**Does not establish:** a hydrogen bond, a favourable interaction energy, a
water-mediated interaction, or an affinity. Dashed lines are distances, not
interaction classifications. Generated hydrogens and angular criteria are not
shown in this figure.

## Evidence and visual language

The first visualization batch uses a restrained legend:

| Visual cue | Meaning |
| --- | --- |
| gray cartoon and pale surface | deposited protein structure |
| yellow sticks | deposited ligand geometry |
| cyan sticks | selected contact residues |
| blue spheres | deposited nearby waters |
| gray dashed lines | measured atom–atom distances |
| blue/orange/purple atom highlights | donor/acceptor/both features from the selected chemical state |

No color means “energetically favourable.” No nearby contact is drawn as a
hydrogen bond merely because it is close.

## Regeneration

Run from the repository root:

```bash
uv run --with rdkit python examples/drug_discovery/capstone_visualizations.py
```

RDKit produces the 2D diagram. PyMOL produces the two 3D renders. The script is
deterministic: when metadata for the same structure checksum already exists, the
stored PyMOL camera view is reused. The metadata file is rewritten with the tool
versions and camera views actually used.

If PyMOL is unavailable, the example stops with a clear error rather than
silently producing a different kind of figure.

## Scope

Visualization Batch 1 covers molecular identity and structure. Later batches can
add:

- chemical-state comparisons and population-independent highlighting;
- interaction-evidence matrices and uncertainty views;
- robustness and sensitivity plots;
- exact-versus-sampled free-energy figures;
- optional interactive 3D exploration.

All later figures should follow the same rule: consume validated records, state
what the figure teaches, and disclose what it does not establish.
