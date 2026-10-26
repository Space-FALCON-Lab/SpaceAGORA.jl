# Installation and Environment

This page explains how to set up a machine for SpaceAGORA: the Julia version
it needs, the Julia environments in the repository and what each is for, the
optional companion packages, and a project of your own outside the checkout.

The first command installs SpaceAGORA and the package versions its committed
manifest records. Terminal, from the repository root (the directory that
contains `Project.toml`):

```text
julia --project=. -e "using Pkg; Pkg.instantiate(julia_version_strict=true)"
```

Command blocks on this page say where they are typed: "Terminal" is the shell
prompt (PowerShell or Command Prompt on Windows, a terminal on Linux/macOS),
and "Julia REPL" is the `julia>` prompt inside a running Julia session.

What to read next:

- [Quickstart](quickstart.md)
- [Assets & Modes](../assets.md)
- [Distributed and HPC](../distributed_hpc.md)

## Julia version

SpaceAGORA requires Julia 1.12; any 1.12.x patch release works. The root
`Project.toml` and the other environments in the repository (examples, test,
docs and the Odyssey example) declare `julia = "~1.12"`, which excludes Julia
1.11 and earlier as well as 1.13 and later. Install a 1.12 release before
instantiating any of them (with juliaup: `juliaup add 1.12` and
`juliaup default 1.12`) and confirm it with `julia --version`.

An unsupported release fails in one of two ways:

- Resolving an environment, which the Odyssey example setup and every
  `Pkg.add`, `Pkg.update` or `Pkg.resolve` do, stops with an error saying that
  the `julia = ~1.12` compat entry does not include your Julia version. Julia
  1.12 words it as `julia version requirement for package at ... not
  satisfied`; other releases word it differently.
- `Pkg.instantiate()` of the root project installs the committed manifest
  without resolving, so it never reads that compat entry. With
  `julia_version_strict=true`, as in the commands on this site, it stops at
  once with `The active manifest file has dependencies that were resolved with
  a different julia version (1.12.3)`. Without the option it only warns, and
  the command fails later on missing standard-library sources or during
  precompilation. On Julia 1.11 or older the option itself is unknown and the
  command fails with `type Context has no field julia_version_strict`; the
  remedy is the same.

[Troubleshooting](troubleshooting.md#Wrong-Julia-version) lists these symptoms
with their fixes.

## Check the installation

A working setup passes this check. Terminal, from the repository root:

```text
julia --version
julia --project=. -e "using SpaceAGORA; println(pkgdir(SpaceAGORA))"
```

The first command prints `julia version 1.12.` followed by the patch number.
The second prints the absolute path of your checkout. The first time it runs
it also precompiles SpaceAGORA and its dependencies, which takes a while and
shows precompilation progress; nothing else should be printed. To make the same
check in a Julia REPL, start Julia from the repository root with
`julia --project=.`:

```julia
julia> Base.active_project()
"/absolute/path/to/SpaceAGORA.jl/Project.toml"

julia> using SpaceAGORA

julia>
```

`Base.active_project()` must return the `Project.toml` of your checkout, which
ends in `SpaceAGORA.jl/Project.toml` for a clone under its default name. Any
other path means Julia was started from another directory or without
`--project=.`, and `using SpaceAGORA` then fails with `Package SpaceAGORA not
found in current path`. See
[Troubleshooting](troubleshooting.md#Wrong-working-directory).

## The repository environments

![Which Julia environment to start: the root project for SpaceAGORA itself and the command-line tool; the examples environment for the examples and plotting scripts; the Odyssey example environment for the frozen-atmosphere tutorial; the viewer environment for the 3D viewer; one environment per optional companion for HYPR, MuJoCo and EDG; and your own project for scenarios kept outside the checkout. Each row shows how to start Julia and the one-time setup command.](../assets/figures/environments_map.svg)

Each environment in the figure has a section below.

### Root project (`--project=.`)

The root project, at the repository root, is the environment SpaceAGORA itself
runs in. Instantiating it installs SpaceAGORA, and its committed manifest pins
the version of every package SpaceAGORA uses; the docs, test and examples
environments are seeded from that manifest when they are set up. Use the root
project to load SpaceAGORA in the Julia REPL or in your own scripts, to run the
command-line tool (which starts its examples in the examples environment), and
to run the installation check above and the scripts under `scripts/` that need
only SpaceAGORA.

### Examples environment

The examples, the plotting scripts and the benchmark studies run in the
examples environment, `examples/Project.toml`. It takes SpaceAGORA from this
checkout and adds the packages the examples use beyond SpaceAGORA itself, such
as Plots and Roots. Prepare it once after the root instantiate, and again after
a pull that changes the root manifest. Terminal, from the repository root:

```text
julia --startup-file=no examples/setup.jl
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

The setup copies the root manifest to the ignored `examples/Manifest.toml` and
resolves from that copy, so the packages SpaceAGORA uses keep the versions the
root manifest records. Only the records of the libraries bundled with Julia,
such as OpenSSL_jll, follow the installed Julia.

Most examples include `examples/common.jl`. Started with `--project=.`, such an
example switches to the examples environment by itself once that environment
is prepared, and otherwise stops and prints the setup command. The cloth-panel
and solar-panel demos do not include it and load PlotlyJS directly, so run
them, like every example, with `--project=examples`. Plotting packages are not
dependencies of SpaceAGORA, so a project of your own that makes plots adds them
itself.

### `examples/odyssey_surrogate_env/`

The Odyssey surrogate example has its own environment. It adds the public
`GRAMSuite` Julia wrapper, which provides the native-free grid atmosphere, and
takes SpaceAGORA from this checkout. Its setup script instantiates it; the
resolved manifest is ignored by Git and stays on your machine. Terminal, from
the repository root:

```text
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate_env/setup.jl
julia --project=examples/odyssey_surrogate_env examples/odyssey_surrogate.jl
```

The setup also installs the EDG guidance companion. It initializes the
companion's Git submodule at `packages/gnc/SpaceAGORAEDG.jl`, which fetches the
recorded revision from GitHub, and adds it to this environment. That rewrites
the tracked `examples/odyssey_surrogate_env/Project.toml`, so `git status` then
shows the file as modified; this is expected. The environment resolves its own
manifest, SpaceAGORA's dependencies included, so it does not depend on the root
instantiate, although the installation order on this site runs that first. The
first setup needs network access. On an empty Julia depot it downloads on the
order of a gigabyte of packages and binary artifacts, and with compiled caches
and the example's assets the depot grows by a few gigabytes. See the
[Odyssey walkthrough](../tutorials/odyssey_surrogate.md).

### Optional package commands

In a Git clone of SpaceAGORA, one helper script manages the optional companion
repositories under `packages/robotics/` and `packages/gnc/`. It needs only
Julia's standard library and the `git` command, and it does not load
SpaceAGORA. Terminal, from the repository root:

```sh
julia --startup-file=no scripts/optional_packages.jl status
julia --startup-file=no scripts/optional_packages.jl checkout hypr
julia --startup-file=no scripts/optional_packages.jl check hypr
julia --startup-file=no scripts/optional_packages.jl setup hypr /path/to/hypr-env
```

The selectors are `hypr`, `mujoco` and `edg`. `status` reports each checkout
without fetching or writing anything. `checkout` initializes the selected
submodule at its recorded Git revision. `check` verifies the checkout's
revision, origin, package name, UUID and version, and that it has no local
changes; it does not run scientific or native tests. `setup` performs
`checkout` and `check`, then runs the package's installer
(`scripts/setup_hypr.jl` for HYPR, the companion's own `scripts/setup_env.jl`
for MuJoCo and EDG) with that local checkout and this core checkout. The
installer creates the environment in the directory you name, which must be
outside the SpaceAGORA checkout. Setup can download dependencies and runs the
installer's usual preparation steps.

Only the selected submodule is initialized: the helper does not fetch other
submodules and leaves the GRAMSuite setup under `data/GRAMSuite.jl` alone. A
package checkout with local changes or at another revision is refused, never
reset. The environment refers to the source paths, so keep those checkouts in
place while you use it.

GitHub source ZIPs leave out submodule contents and the Git metadata needed to
verify the recorded revisions. Clone SpaceAGORA to use these commands, or
follow a companion's standalone installation guide. A package folder present
in a clone does not mean that the package is installed in the Julia
environment you have active.

### 0.7 companion sources

This checkout pins the companion source revisions below. Their package
versions are those of the companions' latest release tags, but the tags do not
declare compatibility with SpaceAGORA 0.7, so the pins name other commits that
do. Use the setup commands above, from this checkout, so that its recorded pins
select these revisions.

| Companion | Package version | Source revision |
| --- | --- | --- |
| HYPR | 0.5.0 | [`2816df11`](https://github.com/Space-FALCON-Lab/HYPR.jl/commit/2816df11c86206cc348ec2c4f82732f958f18ad2) |
| SpaceAGORAEDG | 0.5.1 | [`424d8aa8`](https://github.com/Space-FALCON-Lab/SpaceAGORAEDG.jl/commit/424d8aa822e4f2d0a710231b86c28a1bc1dfd636) |
| SpaceAGORAMuJoCo | 0.5.0 | [`7b5a0031`](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl/commit/7b5a00318f6883d7cdc7c7e7c0c18830f96ce57f) |

Each companion's own default core pin, which its standalone installer uses,
still names a SpaceAGORA 0.6 release. SpaceAGORA's CI jobs for the optional
packages test the revisions above against this checkout. A companion's broad
compatibility range (for example `SpaceAGORA = "0.6, 0.7"`) does not mean that
every revision in that range has been tested with this core.

### Optional HYPR planner

HYPR lives at `packages/robotics/HYPR.jl` and provides RPO and robot-arm
planning. Loading both packages activates HYPR's SpaceAGORA adapter. After
`setup hypr`, start Julia with `--project=/path/to/hypr-env` and load
`SpaceAGORA` and `HYPR`. The
[planner guide](rpo_planner_pilot.md#Installing-optional-HYPR) also describes
the `scripts/setup_hypr.jl` route, which fetches HYPR from its own repository
at the recorded revision instead of using the submodule.

### Optional MuJoCo companion

[SpaceAGORAMuJoCo.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl)
provides standalone proximity scenes and externally propagated spacecraft
inside `run_simulation`. Its repository holds the package, the native binding,
its tests and the Basilisk comparison harness. It requires Julia 1.12; native
scenes support Linux x86_64 only and use the pinned MuJoCo 3.11.0 artifact.

The pinned source is at `packages/robotics/SpaceAGORAMuJoCo.jl`. Terminal,
from the SpaceAGORA repository root:

```sh
julia --startup-file=no scripts/optional_packages.jl setup mujoco /path/to/mujoco-env
julia --startup-file=no --project=/path/to/mujoco-env
```

Then, in the Julia REPL that the second command started, load both packages:

```julia
using SpaceAGORA, SpaceAGORAMuJoCo
```

The helper passes this core checkout to the companion's installer through
`SPACEAGORA_CORE_PATH`, and the installer creates the environment and its
manifest. Loading SpaceAGORA alone loads neither the adapter nor MuJoCo; the
native artifact is downloaded the first time a native operation needs it.

The companion can also be installed on its own: clone it and run its
`scripts/setup_env.jl`. Without an explicit core override, that installer
fetches the exact SpaceAGORA revision the companion records in its
`SpaceAGORASource.toml` (in `Project.toml` on the v0.5.0 release tag), so the
two routes pair the companion with different core sources. On either route the
package versions alone do not establish compatibility: the source pair and the
`ExternalPropagationServices` contract version (2.0.0 in this checkout) must
agree. The companion's
[installation guide](https://github.com/Space-FALCON-Lab/SpaceAGORAMuJoCo.jl#installation)
gives the standalone commands and overrides. Start with
[Orbital Scenes with MuJoCo](mujoco_scenes.md) for the frame conventions,
ownership rules and complete examples. The connection itself is documented in
[External propagation services](../reference/external_propagation_services.md)
and the shared [optional adapter standard](../reference/adapter_contract.md).

### Optional energy-depletion guidance

SpaceAGORAEDG holds the experimental energy-depletion guidance (EDG)
implementation. Its checkout belongs at `packages/gnc/SpaceAGORAEDG.jl`.
SpaceAGORA owns the simulation state, the callback schedule, the accumulated
heat and the application of panel commands; the companion computes the guidance
and control decisions through the `EDGServices` contract, version 3.0.0.

The helper sets it up from the companion revision this checkout records.
Terminal, from the repository root:

```sh
julia --startup-file=no scripts/optional_packages.jl setup edg /path/to/edg-env
julia --startup-file=no --project=/path/to/edg-env
```

Load both packages before creating an EDG mission. Julia REPL:

```julia
using SpaceAGORA, SpaceAGORAEDG
SpaceAGORA.EDGServices.require_provider()
```

The command-line tool runs the two built-in EDG examples
(`AGORA_Odyssey_Energy_Depletion_Control_Test.jl` and `odyssey_surrogate.jl`)
in the prepared EDG environment rather than the examples environment. Start the
tool with `--project=/path/to/edg-env`, or set
`SPACEAGORA_EDG_ENV=/path/to/edg-env`. With `--visualize`, the selected
visualizer environment must also declare SpaceAGORAEDG, so that a mission never
starts with only one of the two companions it needs.

The installer checks the recorded revision and the package identity, then
creates a separate environment that uses this core checkout. A missing pin, an
incompatible service contract or a missing provider raises an error.
SpaceAGORA does not fall back to another guidance method, and loading core
never installs a provider. The companion's standalone installer instead uses
the core revision that the companion records; its `docs/installation.md`
describes that route and the overrides for local development.

Installation and load-order checks show that the packages connect. They do not
show flight readiness or a full reproduction of the EDG study. The
experimental EDG examples keep the original 3 Hz guidance and control rates and
the original heat and structural limits, and known constraint violations in the
sampled output are still reported as violations. Validation with native GRAM
needs the GRAM and SPICE assets the study used. The native-free Odyssey
surrogate example is a separate exercise with its own documented
configuration, not a substitute for that campaign.

CI runs the numerical EDG test suites with the companion installed, in both
normal and coverage runs. The coverage gate holds the companion's source files
to its thresholds along with core's, and fails when any of them has no coverage
measurement. A separate check confirms that SpaceAGORA still loads without EDG.
To reproduce those tests in an isolated development environment, terminal,
from the repository root:

```sh
julia --startup-file=no scripts/setup_edg.jl /path/to/edg-test-env
julia --startup-file=no test/setup.jl
julia --startup-file=no examples/setup.jl
SPACEAGORA_EDG_TEST_ENV=/path/to/edg-test-env julia --startup-file=no --project=test test/runtests.jl
```

The tests run in the repository's test environment, which holds the test
tools, and start their example and command-line child processes in the
examples environment; the EDG environment is added to the load path as the
provider. Both setup scripts keep the core's dependency versions and leave its
tracked project files unchanged. Setting
`SPACEAGORA_EDG_PATH=/absolute/path/to/SpaceAGORAEDG.jl` substitutes a local
EDG checkout, which must be version 0.5.1, for development with this test
helper; it is not a release pin. The full suites can also need the other
optional packages and assets described on this site.

## Your own project

Keep your scripts and results in a project outside the SpaceAGORA checkout.
First check out the SpaceAGORA commit you want to use, and use absolute paths.
The setup script below installs that checkout into the new project as a
development dependency, resolves the dependencies without copying a repository
manifest, and adds the packages that the
[custom-force walkthrough](custom_effector.md) imports. It does not add private
modules to Julia's load path.

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

Run the setup, then save `my_effector.jl` from the walkthrough in the new
project. Terminal, from any directory, since both commands use absolute paths:

```text
julia --startup-file=no setup_user_project.jl /absolute/SpaceAGORA.jl /absolute/my-project
JULIA_LOAD_PATH="@:@stdlib" julia --startup-file=no --project=/absolute/my-project /absolute/my-project/my_effector.jl
```

The `JULIA_LOAD_PATH` setting keeps packages in your default global
environment from hiding a missing direct dependency. On Windows, Julia
separates load-path entries with `;`, so set `JULIA_LOAD_PATH` to `@;@stdlib`
before starting Julia. The walkthrough script writes its results next to
itself, whatever the terminal's working directory. Inspect the
`SimulationResults` and the three written files as described on
[Simulation Outputs](outputs.md).

Setup uses your existing Julia depot and downloads the packages it is missing.
To check that setup works from the depot alone, set `JULIA_PKG_OFFLINE=true`
and `JULIA_PKG_PRECOMPILE_AUTO=0`; a dependency missing from the depot then
makes setup fail instead of being replaced by another version. A setup that
succeeds with a well-stocked depot does not show that it would succeed on a
clean one. Do not change the source checkout while a separate project uses it
as a development path, and keep the checkout's commit and both generated
environment files with any result you keep.

The planner companion has its own
[installation instructions](rpo_planner_pilot.md#Installing-optional-HYPR).
Installing the core project alone does not install HYPR, native GRAM or SPICE
data, and does not set up a full mission workflow.

## Local output paths

Generated results stay on your machine. Git ignores `output/`, and also
`odyssey_surrogate_results/` and `odyssey_surrogate_<cap>/`, where the Odyssey
surrogate example writes by default when it runs from the repository root.

## When you need more than setup

A working environment means the Julia side is ready. Whether the machine also
has the data and licensed assets a run needs is the subject of
[Assets & Modes](../assets.md).
