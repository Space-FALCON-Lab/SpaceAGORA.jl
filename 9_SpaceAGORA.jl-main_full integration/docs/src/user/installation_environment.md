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

### Optional package commands

From a SpaceAGORA Git clone, a shared helper manages the optional repositories
under `packages/robotics/` and `packages/gnc/`. It uses Julia's standard library and command-line Git
for source inspection and checkout; it does not require loading SpaceAGORA first.

```sh
julia --startup-file=no scripts/optional_packages.jl status
julia --startup-file=no scripts/optional_packages.jl checkout hypr
julia --startup-file=no scripts/optional_packages.jl check hypr
julia --startup-file=no scripts/optional_packages.jl setup hypr /path/to/hypr-env
```

The selectors are `hypr`, `mujoco` and `edg`. `status` reports the checkouts; `checkout`
acquires the selected source at its recorded Git revision; `check` verifies its
revision, package identity, version and clean state. It does not run scientific
or native tests. `setup` also performs checkout and verification, then calls the
existing installer with the selected local package and current core checkout.
Use a separate environment outside the source checkouts. Installation can
download dependencies and perform the installer’s normal preparation work.

Only the selected submodule is initialized. This route does not recursively
fetch other packages or change the existing GRAMSuite setup under
`data/GRAMSuite.jl`. A dirty or mismatched package checkout is refused rather
than reset. The environment refers to the source paths, so keep those checkouts
in place while using it.

GitHub source ZIPs omit submodule contents and the Git metadata needed to verify
the recorded revisions. Clone SpaceAGORA to use these commands, or follow a
companion's standalone installation guide. Merely seeing the folder in a clone
does not mean the package is installed in the active Julia environment.

### Optional HYPR planner

HYPR lives at `packages/robotics/HYPR.jl` and supports RPO and robot-arm planning.
The existing `packages/SpaceAGORAHYPR/` compatibility shim loads its SpaceAGORA
adapter. After `setup hypr`, start Julia with `--project=/path/to/hypr-env` and
load `SpaceAGORA` and `SpaceAGORAHYPR`. The [planner guide](rpo_planner_pilot.md#Installing-optional-HYPR)
also documents the existing `scripts/setup_hypr.jl` route, which acquires HYPR
from its external repository without using this submodule layout.

### Optional MuJoCo companion

[SpaceAGORAMuJoCo.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl)
provides standalone proximity scenes and externally propagated spacecraft inside
`run_simulation`. Its repository owns the package, native binding, tests and
Basilisk comparison harness. Julia 1.12 is required; native scenes currently
support Linux x86_64 and use the pinned MuJoCo 3.11.0 artifact.

The pinned source is available at `packages/robotics/SpaceAGORAMuJoCo.jl`.
From the SpaceAGORA repository root:

```sh
julia --startup-file=no scripts/optional_packages.jl setup mujoco /path/to/mujoco-env
julia --startup-file=no --project=/path/to/mujoco-env
```

Then load the packages explicitly:

```julia
using SpaceAGORA, SpaceAGORAMuJoCo
```

The shared helper passes this core checkout to the companion's installer through
`SPACEAGORA_CORE_PATH`. The companion creates its own environment and manifest.
Loading core alone does not load this adapter or download MuJoCo; the native
artifact is obtained when first needed by a native operation.

The standalone route remains available: clone the companion separately and run
its `scripts/setup_env.jl`. Without an explicit core override, that installer
fetches the exact core revision recorded in the companion's `Project.toml`.
The standalone default and the current-checkout route are distinct source pairs.
For either route, a package version-family constraint alone is insufficient;
the source pair and `ExternalPropagationServices` contract must agree.
See the companion's [installation guide](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl#installation)
for standalone commands and overrides. Start with [Orbital Scenes with MuJoCo](mujoco_scenes.md) for
the frame conventions, ownership rules and complete examples. The supported connection is documented in
[External propagation services](../maintainer/external_propagation_services.md)
and the shared [optional adapter standard](../maintainer/adapter_contract.md).

### Optional energy-depletion guidance

SpaceAGORAEDG contains the experimental energy-depletion planning implementation.
Its linked checkout belongs in `packages/gnc/SpaceAGORAEDG.jl`. SpaceAGORA owns the
simulation state, callback schedule, accumulated heat and application of panel
commands; the companion computes decisions through `EDGServices` contract 1.0.0.

For a core revision that records a published EDG companion pin:

```sh
julia --startup-file=no scripts/optional_packages.jl setup edg /path/to/edg-env
julia --startup-file=no --project=/path/to/edg-env
```

Load both packages before creating an EDG mission:

```julia
using SpaceAGORA, SpaceAGORAEDG
SpaceAGORA.EDGServices.require_provider()
```

The CLI preserves the prepared EDG environment for the built-in EDG examples.
Start it with `--project=/path/to/edg-env`, or set
`SPACEAGORA_EDG_ENV=/path/to/edg-env`. With `--visualize`, the selected visualizer
environment must also declare SpaceAGORAEDG; this avoids launching a mission
with only one of its required companions.

The installer checks the recorded revision and package identity, then creates a
separate environment using this core checkout. An absent pin, incompatible
service contract or missing provider raises an error. It does not select a
different guidance method or install a provider merely because core was loaded.
The standalone companion installer instead uses its recorded immutable core
revision. See the companion's `docs/installation.md` for that route and explicit
local development overrides.

Installation and load-order checks establish that the packages connect. They do
not establish flight readiness or full reproduction of the EDG study. The
approved experimental cases retain the original 3 Hz guidance/control rates and
original heat and structural limits. Known sampled constraint violations remain
reported as violations. Native validation still needs the identified GRAM and
SPICE assets. The native-free Odyssey surrogate exercise is a separate example
with its own documented configuration, not a substitute for that campaign.

The existing numerical EDG suites continue to run in both normal and coverage
CI, with the companion required. The parent coverage gate also collects the
companion kernels at the existing thresholds and rejects missing measurements;
extracting a file does not remove it from that gate. A separate core-only check verifies that
loading SpaceAGORA without EDG remains possible. To reproduce those tests in an
isolated development environment:

```sh
julia --startup-file=no scripts/setup_edg.jl /path/to/edg-test-env
julia --startup-file=no --project=/path/to/edg-test-env test/runtests.jl
```

The test setup retains the core's dependency versions and leaves its tracked
project files unchanged. `SPACEAGORA_EDG_PATH=/absolute/path/to/SpaceAGORAEDG.jl`
is an explicit development override for this test helper; it is not a release
pin. Full suites may require the other documented optional packages and assets.

### `docs/`

The docs build uses its own project:

```text
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. docs/make.jl
```

## Your own project

Keep your script and generated results in a project outside the SpaceAGORA
checkout. First select a recorded commit in your existing checkout and use
absolute paths. The setup below explicitly installs that checkout into the new
project, resolves its dependencies without copying a repository manifest, and
declares the packages imported by the
[custom-force walkthrough](custom_effector.md). It does not add private modules
to Julia's load path.

Save this as `setup_user_project.jl` outside the checkout:

```julia
using Pkg
length(ARGS) == 2 || error("Usage: setup_user_project.jl CHECKOUT NEW_PROJECT")
checkout, project = abspath.(ARGS)
isdir(project) && !isempty(readdir(project)) &&
    error("Choose a new empty project directory")
mkpath(project)
Pkg.activate(project)
Pkg.develop(PackageSpec(path=checkout); preserve=Pkg.PRESERVE_ALL)
Pkg.add([PackageSpec(name="StaticArrays"), PackageSpec(name="CSV"),
         PackageSpec(name="DataFrames")]; preserve=Pkg.PRESERVE_ALL)
Pkg.instantiate()
```

Run setup, then save `my_effector.jl` from the walkthrough in the new project:

```text
julia --startup-file=no setup_user_project.jl /absolute/SpaceAGORA.jl /absolute/my-project
JULIA_LOAD_PATH="@:@stdlib" julia --startup-file=no --project=/absolute/my-project /absolute/my-project/my_effector.jl
```

The load-path setting prevents unrelated packages in your default global project
from hiding a missing direct dependency. On Windows, Julia uses `;` as the load
path separator, so set `JULIA_LOAD_PATH` to `@;@stdlib` before launching Julia.
The script writes results beside itself, independently of the terminal's working
directory. Inspect `SimulationResults` and the three written files as described
on [Simulation Outputs](outputs.md).

Setup uses your existing Julia depot and can download missing packages. An
explicit cached/offline check may instead set `JULIA_PKG_OFFLINE=true` and
`JULIA_PKG_PRECOMPILE_AUTO=0`; a missing cached dependency then remains a setup
failure, not permission to substitute another version. Success with a populated
depot is not a clean-depot installation result. Do not modify the source checkout
while a separate project is using its development path; record the checkout
commit and both generated environment files with any result you keep.

The optional planner companion has its own
[documented installation](rpo_planner_pilot.md#Installing-optional-HYPR). Installing
the core project alone does not install or qualify HYPR, native GRAM, SPICE data,
or a full mission workflow.

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
