# Develop Moleculab

This is the contributor entry point. It covers how the software is structured,
how it is validated, and where to start when changing it.

## Start here

| Material | Use it for |
| --- | --- |
| [Architecture](architecture.md) | layers, dependency rules, data contracts, and decision records |
| [Validation](../reference/validation.md) | what the tests establish and how to run them |
| [Performance](performance.md) | profiling, backend choices, and acceleration boundaries |
| [Acceleration backends](acceleration.md) | CPU, ROCm/PyTorch, and NPU limitations |
| [Roadmap](roadmap.md) | long-term scope and deferred work |
| [AGENTS.md](../../AGENTS.md) | repository conventions for contributors and agents |
| [`tests/`](../../tests/) | executable validation suite |

## Working rules

1. Keep physical models, numerical algorithms, protocols, analysis, and
   applications in their existing layers.
2. Preserve `pair_potentials`, `topology`, `box`, `nonbonded_model`,
   `constraints`, `indexed_pair_potentials`, particle `atom_type`, and particle
   `charge` through every `System` rebuild.
3. Add an independent test for every physical relationship, not merely a test
   that code runs.
4. Add a runnable example under [`examples/`](../../examples/README.md) and link it from
   the relevant tutorial or examples index.
5. Keep parser and chemistry-toolkit adapters outside the physics layer.
6. Run `uv run pytest -q` before finishing.

## Package layers

```text
foundation → models → composition
            → algorithms + protocols
            → analysis + applications
```

The full dependency map and rationale are in
[`architecture.md`](architecture.md).

## Documentation rules

- There is one project overview: [`README.md`](../../README.md).
- There is one learning entry point: [`tutorials/README.md`](../tutorials/README.md).
- There is one development entry point: this page.
- Every concept has one authoritative explanation; other pages link to it.
- Every example is linked from a tutorial or the examples index.
- Experimental, calculated, predicted, and hypothesised results remain
  distinguishable.
