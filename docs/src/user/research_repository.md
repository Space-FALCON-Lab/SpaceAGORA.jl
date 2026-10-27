# Building Your Own Research Repository

Research software built on SpaceAGORA can live in a repository of your own: a
study, a package with a model or an algorithm, or a companion package. That
repository is tested against a recorded SpaceAGORA commit with the same checks
the core runs on its built-in models.

The kit has three parts:

- starters, the repository layouts under `templates/starters/`, with a
  generator that creates a repository from one;
- conformance checks, `SpaceAGORA.Conformance`, which run the engine's own
  machinery on your models and algorithms;
- a reusable CI workflow, `.github/workflows/contributor-ci.yml`, which your
  repository calls at the same commit it records.

## Who owns what

You own the scientific method, its parameters, its validation and the claims
made from it. SpaceAGORA owns scheduling (callbacks, rates and solver stops),
state ownership (run state in the model, deep copies, binding), reproducible
randomness (supplied and derived random streams), the connections between
components (links, families and install targets), lifecycle handling (pre-run
checks and resets) and the testing machinery.

A conformance pass checks interface behavior on the supplied fixtures. It does
not validate the scientific model, the method or the results, and the
starters' READMEs and the CI summary say so.

## Choose a starter

| Starter | For | Its example | Manifest |
|---|---|---|---|
| `project` | a research study on the public API | a force model, `LowThrustPush` (a `wrench` method), and a short case in `scripts/run_case.jl` | committed, so every result has its exact package versions |
| `package` | a feedback model, or an algorithm for an existing family | an RPO reference planner, run in the RPO pilot by `examples/run_in_simulation.jl` | not committed |
| `companion` | a companion that serves a core service module (also called a facade) | an external propagator on `ExternalPropagationServices` | not committed |

A study or a package should use only names exported by `SpaceAGORA` and
`SpaceAGORA.Conformance`. The workflow scans its Julia files and refuses
references to `SpaceAGORA.SimulationModel`, `SpaceAGORA.SimulationEngine` and
underscore-prefixed names. A package can implement a feedback model directly
or add a planner to an existing family, but it never adds a family of its own.
A new family or a new service module is trusted code owned by the core, so it
starts as a SpaceAGORA pull request (the repository's contributor guide,
`docs/contributing/contributor_kit_integrations.md`, lists what that needs),
and the companion starter is then pointed at it.

For a controller package, follow
[Your First Controller Package](feedback_package.md). It adapts the existing
`package` starter and shows which planner files to replace; there is no
separate controller starter and no new core family to create.

## Create the repository

From a clean checkout of the SpaceAGORA commit you want to build on:

```bash
julia --startup-file=no scripts/new_contributor_repository.jl package ../MyPlanner --name MyPlanner
cd ../MyPlanner
```

Run the remaining commands from `MyPlanner`. They assume that the core
checkout beside it is named `SpaceAGORA.jl`; adjust the path if yours has
another name.

The generator copies the starter as it is in that commit, gives the new
repository a fresh UUID and a `[compat]` entry for the core's release series,
and writes the commit (`HEAD`) in three places: the pin file
`SpaceAGORASource.toml`, the workflow call in `.github/workflows/ci.yml`, and
the core checkout command in the generated `README.md`. A `project` also gets
the commit's Manifest as its committed Manifest, and a `companion` requires the
service module's current contract version. A package or a companion gets the
MIT license with this year and the copyright holder you pass as
`--holder "<name>"`; without it, the placeholder `Your name or organization`
stays for you to replace. The generator refuses a checkout with tracked
changes, a target directory that exists and is not empty, a directory inside
the checkout, and a name that is not a package name or that contains the
starters' placeholder `StarterName`.

## The pin

`SpaceAGORASource.toml` records the core your repository is tested against:

```toml
url = "https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git"
rev = "<the full 40-character commit>"
uuid = "afbfb69f-5c0b-4832-b760-43725dff8540"
version = "0.7.1"          # optional: the core's version at rev
companions = ["mujoco"]    # optional: hypr, mujoco or edg, at that commit's own pins
```

The pin lives in its own file, not under `[sources]` in `Project.toml`. Pkg
cannot install SpaceAGORA from a URL and revision because the core has
submodules, so the kit tests against a local clone of the core at the pinned
commit. The setup script refuses a `SpaceAGORA` entry in `[sources]`, which
`Pkg.test` would use in place of the tested core.

To move to a newer core, check that core out beside your repository and run
its `repin`. It changes `rev` and `version` in the pin file and the commit in
the workflow call together, and rewrites the `SpaceAGORA` entry in `[compat]`
when it does not admit the new core's version:

```bash
julia --startup-file=no --project=@stdlib ../SpaceAGORA.jl/scripts/contributor_environment.jl repin package .
```

`repin` does not change the README's checkout command, so replace the commit
there with the new `rev` by hand. Then rebuild the environment, and for a
project, commit the new Manifest. Record the commit with every result you
keep.

## Set up and test locally

There are two routes; use one of them, not both. Each builds the environment
in a directory that must be new or empty, so a second build into the same
directory is refused:

- The generated README's commands are self-contained: they clone a new
  `../MyPlanner-core` at the pinned commit and build `../MyPlanner-env` from it.
- The commands below reuse a sibling `../SpaceAGORA.jl` checkout instead.

If you generated `MyPlanner` yourself, use the clean SpaceAGORA checkout that
generated it, which is already at the pinned commit. If you instead cloned an
existing contributor repository and do not have that core checkout yet,
create it from the contributor repository's root:

```bash
git clone --filter=blob:none -c lfs.fetchexclude="*" https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git ../SpaceAGORA.jl
julia --startup-file=no --project=@stdlib
```

At the Julia prompt, read the recorded revision and check it out in the new
clone:

```julia
using TOML
pin = TOML.parsefile("SpaceAGORASource.toml")
run(`git -C ../SpaceAGORA.jl checkout $(pin["rev"])`)
exit()
```

With the sibling checkout at the recorded revision, build the separate
environment from the contributor repository's root. For a repository other
than this `MyPlanner` package, put its kind and name in place of `package` and
`MyPlanner`:

```bash
julia --startup-file=no ../SpaceAGORA.jl/scripts/contributor_environment.jl package . ../MyPlanner-env --manifest core
```

`--filter=blob:none` skips the contents of old file versions: Git downloads the
files of a commit when you check it out. `-c lfs.fetchexclude="*"` keeps this
clone, its checkouts and later pulls from downloading SpaceAGORA's large data
files (Git LFS), which a contributor repository does not use. Without Git LFS
installed the setting does nothing, and nothing is downloaded either. The
commands on this page run as written in bash or zsh, PowerShell and Command
Prompt.

The first argument of `contributor_environment.jl` is the starter's kind, and
the environment directory must lie outside both the repository and the
SpaceAGORA checkout. `--manifest core` seeds the environment from the core's
Manifest and keeps its versions; a project uses `--manifest committed`, its
own Manifest. The script refuses a pin file whose `rev` is not the checked-out
commit, so the commit you record is the commit you test. Then run the tests in
that environment: `Pkg.test` for a package or a companion, `test/runtests.jl`
for a project. Each starter's README gives the exact commands.

## The conformance checks

Each `check_*` function returns a `ConformanceReport` and needs no test
framework. With `using Test`, the matching `test_*` function records the same
checks as tests.

| What you wrote | Check | It runs |
|---|---|---|
| a force or torque model | `test_force_torque_model(model)` | the pre-run hook check, the supported sample signature, repeated and concurrent runs, a deep copy, a reseeded global random stream, pure output hooks, the declared environment |
| a sensor kind | `test_sensor(sensor)` | its run-state fields and requirements, the measurement shape, the run-level checks, seeds that repeat, the reset of its own state, a run resumed from a checkpoint against an uninterrupted one |
| a guidance, navigation or control model | `test_guidance_model`, `test_navigation_model(model; sensors)`, `test_control_model` | the pre-run hook check, the run-level checks, the reset at the run's start (`initialize_*!`), fresh state from `bind_spacecraft`, a run resumed from a checkpoint against an uninterrupted one; for a control model also pure output hooks, positive propellant consumption and zero output for spacecraft outside its scope |
| an algorithm of a family | `test_algorithm(algorithm, requests; trusted, simulation)` | the family's representative requests (`rpo_conformance_requests()`, `burn_plan_conformance_requests()`, `joint_reference_conformance_requests()`), a deterministic walk, an untouched stream when declared, the request driver, a short simulation |
| a family or a service module (a core pull request) | `test_algorithm_family`, `test_service_facade(facade; version)` | the family's obligations and rejection policy; the service module's shape at its contract version |

The run-level checks use a short Earth orbit that needs no native GRAM, SPICE
kernels or private data. A model that needs more is checked on a fixture of
its own: `conformance_fixture(; configuration=...)` with a configuration that
already contains the models under test, or the keywords of
`conformance_fixture` for the models, sensors and rates. Run the checks with
two or more threads, so that concurrent runs are compared with sequential
ones. [Writing a GNC Algorithm](gnc_algorithm.md#Test-your-algorithm) shows
them on a controller, an estimator and a planner.

Each check records the release that added it. A check's status is `:pass`,
`:fail`, `:skip` (the check could not apply, and its detail says why) or
`:advisory_fail`: a check added in the current minor release reports a
failure without failing, and becomes required in the next minor release. The
checks of SpaceAGORA 0.3.0 have been required since 0.3.0.
`checkpoint_round_trip`, added in 0.4.0, was advisory in 0.4 and is required
from 0.5.0: the fixture run, checkpointed halfway, interrupted after the
checkpoint and resumed from it, must end like an uninterrupted run, so the
model's `checkpoint_state` must carry the run state it keeps (see
[Run state in checkpoint files](../extensibility.md#Run-state-in-checkpoint-files)).
`outputs_in_scope`, the control-model check added in 0.7.0, is advisory in
0.7 and required from 0.8.0.

## Continuous integration

The starter's `.github/workflows/ci.yml` calls the reusable workflow at the
pinned commit:

```yaml
jobs:
  spaceagora:
    uses: Space-FALCON-Lab/SpaceAGORA.jl/.github/workflows/contributor-ci.yml@<sha>
    with:
      spaceagora-rev: "<sha>"
      kind: package
      examples: examples/run_in_simulation.jl
```

The workflow checks out both repositories and refuses anything but a full
commit that matches the pin file. For a project or a package it runs the
internal-name scan. It then builds the environment as above, precompiles
strictly, and runs the tests and the `examples` scripts under
`--depwarn=error` on two threads. For a companion it checks that an
incompatible core is refused and, with `load-orders: true` (as the companion
starter sets), loads the companion and the core in both orders. It ends by
checking that no tracked file changed, and writes a summary table with one
row per check. Its other inputs are `manifest` (`core`, `committed`, or
`resolve`, an advisory fresh resolution), `depwarn`, `threads`, `directory`
(the package inside a larger repository), `test-entrypoint` and
`julia-version`.

A standalone companion may declare SpaceAGORA in `[weakdeps]` instead of
`[deps]`. The shared workflow adds the pinned core only to its separate
integration environment and keeps the companion's compatibility bound. Keep
the standalone tests in a job of their own. If `test/runtests.jl` deliberately
requires SpaceAGORA to be absent, set `test-entrypoint: test/contributor_ci.jl`
to select the integration suite. This path is relative to `directory`, must
resolve to a file inside it, and runs with `test/` as its working directory;
the default is `test/runtests.jl`. The load-order and incompatible-core probes
still run in separate processes.

SpaceAGORA's own CI generates each starter from the commit under test and runs
it through this workflow, so the starters stay in step with the core. It also
runs the kit's demonstration repository,
[SpaceAGORAKitDemo.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAKitDemo.jl),
a package built from the package starter: every night against `main`, and in
each release pull request against the release candidate, with its pin moved to
the commit under test.

## Where next

- [Writing a GNC Algorithm](gnc_algorithm.md): controllers, estimators and
  planners, with the checks.
- [Extensibility](../extensibility.md): every extension point and its template.
- [Extension Conventions](../reference/extension_conventions.md): the rules every
  hook follows.
- [Optional Adapter Standard](../reference/adapter_contract.md): the rules for a
  companion's service module.
- [Your own project](installation_environment.md#Your-own-project): a quick
  separate project without the kit, for a script.
