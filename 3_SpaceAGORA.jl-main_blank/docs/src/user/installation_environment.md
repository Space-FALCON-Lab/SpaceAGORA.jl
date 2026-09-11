# Installation and Environment

Use this page when you need to set up a machine for local development, examples,
or docs builds.

This page is for users and contributors who need to understand the repository
environments before running anything heavier than a quick smoke path.

Shortest successful command:

```text
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

What to read next:

- [Quickstart](quickstart.md)
- [Assets & Modes](../assets.md)
- [Distributed and HPC](../distributed_hpc.md)

## Julia version

The root `Project.toml` pins `julia = "1.12"`. Install that Julia release
before instantiating either project.

## The repository environments

### Root project (`--project=.`)

The root project at the repository root is the canonical execution environment
used by examples, tests, benchmarks, and most local commands in this
repository.

Use it for:

- example scripts
- local verification runs
- benchmark launchers
- CLI commands

### `docs/`

The docs build uses its own project:

```text
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. docs/make.jl
```

## Local output paths

Generated reports and builds stay local. Common ignored output roots are:

- `output/`
- `docs/build/`
- `docs/site/`
- `docs/src/generated/`

You can regenerate common local artifacts with:

```text
julia --project=. scripts/regenerate_ignored_outputs.jl runtime-analysis quick
julia --project=. scripts/regenerate_ignored_outputs.jl telemetry quick
julia --project=. scripts/regenerate_ignored_outputs.jl docs
```

## When you need more than setup

Environment setup only tells you that the Julia side is ready. If you need to
know whether the machine also has the right data and licensed assets, go to
[Assets & Modes](../assets.md).
