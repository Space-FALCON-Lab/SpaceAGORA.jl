# Building Your Own Research Repository

Use this page to keep research software built on SpaceAGORA in a repository of
your own: a study, an algorithm package or a companion package, tested against
a recorded SpaceAGORA commit by the same checks the core runs on its built-in
models.

The kit has three parts:

- **starters**, the repository layouts under `templates/starters/`, and a
  generator that creates a repository from one;
- **conformance checks**, `SpaceAGORA.Conformance`, which run the engine's own
  machinery on your models and algorithms;
- **a reusable CI workflow**, `.github/workflows/contributor-ci.yml`, which your
  repository calls with the same commit it records.

## Who owns what

You own the scientific method, its parameters, its validation and the claims
made from it. SpaceAGORA owns scheduling (callbacks, rates and solver stops),
state ownership (run state in the model, deep copies, binding), reproducible
randomness (supplied and derived streams), connections between components
(links, families and install targets), lifecycle handling (pre-run checks and
resets) and the testing machinery.

A conformance pass certifies interface behaviour only. It does not validate a
model, a method or a result; the starters' READMEs and the CI summary say so.

## Choose a starter

| Starter | For | Its example | Manifest |
|---|---|---|---|
| `project` | a research study on the public API | a force model (`wrench`) and a short case in `scripts/run_case.jl` | committed, so every result has its exact package versions |
| `package` | an algorithm package for an existing family | an RPO reference planner, run in the RPO pilot by `examples/run_in_simulation.jl` | not committed |
| `companion` | a companion that serves a core service facade | an external propagator on `ExternalPropagationServices` | not committed |

A study or a package uses only names exported by `SpaceAGORA` and
`SpaceAGORA.Conformance`; the workflow scans it for internal names. A package
adds algorithms to a family; it never adds a family of its own. A new family or
a new service facade is trusted, core-owned code, so it starts as a SpaceAGORA
pull request (the repository's contributor guide,
`docs/contributing/contributor_kit_integrations.md`, lists what it needs), and
the companion starter is then pointed at it.

## Create the repository

From a clean checkout of the SpaceAGORA commit you want to build on:

```bash
julia --startup-file=no scripts/new_contributor_repository.jl package ../MyPlanner --name MyPlanner
```

The generator copies the starter from that commit, gives the repository a fresh
UUID and a `[compat]` entry for the core's release series, and pins the commit
(`HEAD`) in two places: the pin file `SpaceAGORASource.toml` and the workflow
call in `.github/workflows/ci.yml`. A `project` also gets the commit's Manifest
as its committed Manifest; a `companion` requires the facade's current contract
version. It refuses a checkout with tracked changes, a non-empty directory, a
directory inside the checkout and a name that is not a package name.

## The pin

`SpaceAGORASource.toml` records the core your repository is tested against:

```toml
url = "https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git"
rev = "<the full 40-character commit>"
uuid = "afbfb69f-5c0b-4832-b760-43725dff8540"
version = "0.5.0"          # optional: the core's version at rev
companions = ["mujoco"]    # optional: hypr, mujoco or edg, at that commit's own pins
```

The pin lives in its own file, not in `[sources]` of `Project.toml`: Pkg cannot
install SpaceAGORA from a URL and revision, and `Pkg.test` would read such an
entry instead of the tested environment. To move to a newer core, check that
core out beside your repository and run its `repin`, which changes `rev`, the
`version` and the workflow call together, and the `SpaceAGORA` entry in
`[compat]` when the release series changes:

```bash
julia --startup-file=no --project=@stdlib ../SpaceAGORA.jl/scripts/contributor_environment.jl repin package .
```

Then rebuild the environment, and for a project commit the new Manifest.
Record the commit with every result you keep.

## Set up and test locally

Clone SpaceAGORA beside your repository at the pinned commit, then build a
separate environment that develops that checkout by path:

```bash
git clone -c lfs.fetchexclude='*' https://github.com/Space-FALCON-Lab/SpaceAGORA.jl.git ../SpaceAGORA.jl
git -C ../SpaceAGORA.jl checkout <the rev in SpaceAGORASource.toml>
julia --startup-file=no ../SpaceAGORA.jl/scripts/contributor_environment.jl package . ../MyPlanner-env --manifest core
```

`-c lfs.fetchexclude='*'` keeps this clone, its checkouts and later pulls from
downloading SpaceAGORA's large data files (Git LFS), which a contributor
repository does not use. Without Git LFS installed the setting does nothing,
and nothing is downloaded either. The first argument is the starter's kind.
`--manifest core` seeds the environment from the core's Manifest and keeps its
versions; a project uses `--manifest committed`, its own Manifest. The script
refuses a pin file whose `rev` is not the checked-out commit, so the commit you
record is the commit you test. Then run the tests in that environment:
`Pkg.test` for a package or a companion, `test/runtests.jl` for a project. Each
starter's README gives the exact commands.

## The conformance checks

Each `check_*` function returns a `ConformanceReport` and needs no test
framework; with `using Test`, the matching `test_*` function records the same
checks as tests.

| What you wrote | Check | It runs |
|---|---|---|
| a force or torque model | `test_force_torque_model(model)` | the pre-run hook check, repeated and concurrent runs, a deep copy, a reseeded global stream, pure output hooks, declared environment |
| a sensor kind | `test_sensor(sensor)` | its run-state fields and requirements, measurement shape, seeds that repeat, reset of its own state, a run resumed from a checkpoint against an uninterrupted one |
| a guidance, navigation or control model | `test_guidance_model`, `test_navigation_model(model; sensors)`, `test_control_model` | the pre-run hook check, the run-level checks, the reset at the run's start (`initialize_*!`), pure output hooks and positive consumption (control), fresh state from `bind_spacecraft`, a run resumed from a checkpoint against an uninterrupted one |
| an algorithm of a family | `test_algorithm(algorithm, requests; trusted, simulation)` | the family's representative requests (`rpo_conformance_requests()`, `burn_plan_conformance_requests()`, `joint_reference_conformance_requests()`), a deterministic walk, an untouched stream when declared, the request driver, a short simulation |
| a family or a facade (a core pull request) | `test_algorithm_family`, `test_service_facade(facade; version)` | the family's obligations and rejection policy; the facade's shape at its contract version |

The run-level checks use a short Earth orbit that needs no native GRAM, SPICE
kernels or private data. A model that needs more is checked on a fixture of its
own: `conformance_fixture(; configuration=...)`, or the keywords of
`conformance_fixture` for the models, sensors and rates. Run the checks with two
or more threads so that concurrent runs are compared with sequential ones.
[Writing a GNC Algorithm](gnc_algorithm.md#Test-your-algorithm) shows them on a
controller, an estimator and a planner.

Each check records the release that added it. A check is `:pass`, `:fail`,
`:skip` (it could not apply, and its detail says why) or `:advisory_fail`: a
check added in the current minor release reports a failure without failing,
and becomes required in the next minor release. The checks of SpaceAGORA 0.3.0
are required from 0.3.0. `checkpoint_round_trip`, added in 0.4.0, is advisory
in 0.4 and required from 0.5.0: a run of the fixture checkpointed halfway,
interrupted after the checkpoint and resumed from it must end like an
uninterrupted run, so the model's `checkpoint_state` carries the run state it
keeps (see [Run state in checkpoint files](../extensibility.md#Run-state-in-checkpoint-files)).

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

The workflow checks out both repositories, refuses anything but a full commit
that matches the pin file, builds the environment as above, precompiles
strictly, runs the tests and the `examples` scripts under `--depwarn=error` on
two threads, and for a companion refuses an incompatible core and loads it in
both orders. It ends by checking that no tracked file changed, and writes one
summary line per check. Other inputs: `manifest` (`core`, `committed`, or
`resolve`, an advisory fresh resolution), `depwarn`, `threads`, `directory` (the
package inside a larger repository) and `julia-version`.

SpaceAGORA's own CI generates each starter from the commit under test and runs
it through this workflow, so the starters stay in step with the core. It also
runs the kit's demonstration repository,
[SpaceAGORAKitDemo.jl](https://github.com/Space-FALCON-Lab/SpaceAGORAKitDemo.jl),
a package built from the package starter, every night against `main` and in each
release pull request against the release candidate, with its pin moved to the
commit under test.

## Where next

- [Writing a GNC Algorithm](gnc_algorithm.md): controllers, estimators and
  planners, with the checks.
- [Extensibility](../extensibility.md): every extension point and its template.
- [Extension Conventions](../reference/extension_conventions.md): the rules every
  hook follows.
- [Optional Adapter Standard](../reference/adapter_contract.md): the rules for a
  companion's facade.
- [Your own project](installation_environment.md#Your-own-project): a quick
  separate project without the kit, for a script.
