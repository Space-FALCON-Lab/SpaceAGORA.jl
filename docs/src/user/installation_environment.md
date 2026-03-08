# Installation and Environment

Use this page when you need to set up a machine for local development, examples,
or docs builds.

This page is for users and contributors who need to understand the repository
environments before running anything heavier than a quick smoke path.

Shortest successful command:

```bash
julia --project=.AGORA -e 'using Pkg; Pkg.instantiate()'
```

What to read next:

- [Quickstart](quickstart.md)
- [Assets & Modes](../assets.md)
- [Distributed and HPC](../distributed_hpc.md)

## The repository environments

### `/.AGORA`

`/.AGORA` is the committed execution environment used by examples, tests,
benchmarks, and most local commands in this repository.

Use it for:

- example scripts
- local verification runs
- benchmark launchers
- CLI commands

### `docs/`

The docs build uses its own project:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

## Local output paths

Generated reports and builds stay local. Common ignored output roots are:

- `output/`
- `docs/build/`
- `docs/site/`
- `docs/src/generated/`

You can regenerate common local artifacts with:

```bash
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl runtime-analysis quick
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl telemetry quick
julia --project=.AGORA scripts/regenerate_ignored_outputs.jl docs
```

## When you need more than setup

Environment setup only tells you that the Julia side is ready. If you need to
know whether the machine also has the right data and licensed assets, go to
[Assets & Modes](../assets.md).
