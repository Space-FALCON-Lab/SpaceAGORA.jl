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

### `examples/odyssey_surrogate_env/`

The Odyssey surrogate example has its own environment. It adds the public
`GRAMSuite` Julia wrapper, which provides the native-free grid atmosphere, and
uses the current SpaceAGORA checkout. Its setup script instantiates it, and the
resolved manifest stays local:

```text
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate_env/setup.jl
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl
```

See the [Odyssey walkthrough](../tutorials/odyssey_surrogate.md).

### Optional MuJoCo companion

[SpaceAGORAMuJoCo.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl)
provides standalone proximity scenes and externally propagated spacecraft inside
`run_simulation`. Its repository owns the package, native binding, tests and
Basilisk comparison harness. Julia 1.12 is required; native scenes currently
support Linux x86_64 and use the pinned MuJoCo 3.11.0 artifact.

Build an independent environment from the companion checkout:

```text
git clone https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl.git
cd SpaceAGORAMuJoCo.jl
julia --startup-file=no scripts/setup_env.jl /path/to/mujoco-env
julia --startup-file=no --project=/path/to/mujoco-env
```

Then load the packages explicitly:

```julia
using SpaceAGORA, SpaceAGORAMuJoCo
```

The companion's setup installs the exact core revision recorded in its
`Project.toml` and creates its own manifest. It does not use a parent checkout,
copy a core manifest, or call a HYPR installer. The package version-family
constraint alone is not a compatibility guarantee; the source pair and
`ExternalPropagationServices` contract version must agree. Loading core alone
does not load this adapter or download MuJoCo.

For development, `SPACEAGORA_CORE_PATH` selects a local SpaceAGORA checkout
instead of that source pin. Core CI uses this override with an immutable
companion revision so changes to core are tested against the actual candidate.
See the companion's [installation guide](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl#installation)
for overrides. Start with [Orbital Scenes with MuJoCo](mujoco_scenes.md) for
the frame conventions, ownership rules and complete examples. The supported connection is documented in
[External propagation services](../maintainer/external_propagation_services.md)
and the shared [optional adapter standard](../maintainer/adapter_contract.md).

### `docs/`

The docs build uses its own project:

```text
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. docs/make.jl
```

## Local output paths

Generated reports and builds stay local. Common ignored output roots are:

- `output/`
- `odyssey_surrogate_results/` and `odyssey_surrogate_<cap>/`, the Odyssey
  surrogate example's defaults when run from the repository root
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
