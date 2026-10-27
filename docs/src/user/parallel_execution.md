# Parallel Execution

Use this page when you want a simulation, a Monte Carlo campaign or a
constellation ensemble to use more than one core.

Shortest successful command:

```text
julia --project=examples --threads=auto examples/AGORA_Earth_MonteCarlo.jl
```

By default this example runs its samples on a fixed number of worker tasks,
one per Julia thread up to the number of samples; with
`SPACEAGORA_MC_PARALLEL=1` it uses the parallel flag described below instead.

What to read next:

- [Distributed and HPC](../distributed_hpc.md)
- [Configuring a Simulation](simulation_configuration.md)

## Turning on parallel execution

Parallel execution is one switch in the solver setup:
`SolverConfig(parallel=true)`. With it set, SpaceAGORA decides how to use the
threads Julia was started with. A single run threads its callbacks, effectors
and right-hand side where its adaptive policy expects threading to help. A
campaign lets the predictive planner choose serial, threaded or process-worker
execution for the campaign's shape before the first sample runs. Nothing else
needs configuring.

For a single simulation, set the flag in the configuration; in this fragment
the other fields are elided:

```julia
using SpaceAGORA

args = SimulationConfiguration(
    # ... the usual fields ...
    solver_config=SolverConfig(parallel=true),
)
run_simulation(args)
```

A constellation ensemble reads the flag from the configuration it is given.
`run_monte_carlo` only sees your function, so it takes the flag as a keyword.
In this fragment, `args` is a multi-spacecraft configuration with the flag set
and `make_config_for_seed` is your own function:

```julia
result = run_constellation_ensemble(args; return_solution=true)

result = run_monte_carlo(1:100; parallel=true) do seed
    run_simulation(make_config_for_seed(seed))
end
```

Requirements and behavior:

- Start Julia with threads, for example `julia --threads=auto`. With a single
  thread the flag can only use process workers for campaigns, and says so
  once.
- The first parallel run on a machine measures the machine's cost constants
  and stores them in
  `output/parallel_policy_state/cost_constants_<fingerprint>.toml` under the
  working directory (`SPACEAGORA_COST_CONSTANTS_PATH` sets another path). How
  long that takes depends on the machine, the thread count and what is
  already compiled. Later runs reuse the constants while they are valid. A
  file measured with fewer threads than the current session is measured
  again, because the constants only reach the thread count of the session
  that measured them, and so is a file that cannot be used. Processes that
  start together in one working directory measure once: the others wait for
  the file, and a `.lock` file beside it marks the measurement in progress. To
  redo the measurement, delete the file or run
  `julia --project=. --threads=<T> scripts/calibration/calibrate_machine.jl`.
- The flag's settings apply only for the duration of the call. Nothing is
  left in the process environment afterwards, also when the run throws or when
  several flagged runs overlap on different tasks. While any flagged run is
  active, the settings are visible to every task in the process, because they
  live in the process environment.
- `run_monte_carlo(...; parallel=true)` does not send samples to process
  workers unless you describe the workload with `route_features` (see
  [Internal and benchmark controls](#Internal-and-benchmark-controls)). A
  process worker is a separate Julia process: a function defined in your
  script, anything it reaches in your script's `Main`, and SPICE kernels you
  furnished yourself exist there only if you load them there, so pass
  `route_features` only once the workers are prepared. A constellation
  ensemble derives its features from its configuration and can use process
  workers. A configuration that holds types or functions defined in your
  script needs the same preparation, or `parallel=false` with a fixed
  `threads=N`.
- The flag picks the worker count itself, so `parallel=true` together with an
  integer `threads=` is an `ArgumentError`. Every member of a campaign must
  carry the same flag; members that disagree are an `ArgumentError` too.
- With `parallel=false`, the default, nothing changes: a single run behaves as
  it does without the flag, and campaigns are serial unless you give a fixed
  `threads=N` (see below).
- A run or campaign started inside another campaign's worker keeps the
  enclosing campaign's thread budget rather than claiming the whole machine
  again.

Parallel runs are not bit-for-bit reproducible across thread counts in every
configuration (see
[Internal and benchmark controls](#Internal-and-benchmark-controls) for the
route that is). Validate a study against a run with the flag off before
relying on it.

## Monte Carlo campaigns

Use `run_monte_carlo` to run many independent simulations from an example
script or a study without writing your own orchestration. The runner
parallelizes across samples: each worker task calls your function with one
seed. `parallel=true` (above) lets the planner choose; this section describes
a fixed worker count instead.

Start Julia with the number of threads you want available:

```bash
julia --project=examples --threads=8 examples/AGORA_Earth_MonteCarlo.jl
```

Then cap the campaign runner at that thread count or fewer. As above,
`make_config_for_seed` stands for your own function:

```julia
using SpaceAGORA

result = run_monte_carlo(1:100; threads=8) do seed
    args = make_config_for_seed(seed)
    run_simulation(args; return_solution=true)
end

println("successful samples: $(length(result.successful))")
println("failed samples: $(length(result.failed))")
```

`threads` does not create Julia threads at runtime. If the script asks for
more Monte Carlo threads than Julia was started with, the runner throws an
`ArgumentError` that says how many threads are available.

For reusable scripts, build a spec explicitly:

```julia
spec = MonteCarloSpec(seeds=1001:1100, threads=8, fail_fast=false)

result = run_monte_carlo(spec) do seed
    args = make_config_for_seed(seed)
    run_simulation(args)
end
```

By default, failed samples are collected in `result.failed` and do not stop
the rest of the campaign. Set `fail_fast=true` to have the runner rethrow the
first failure instead.

### Process-backend outer parallelism

When the campaign planner picks `:process`, or when you want to control it
directly, campaigns dispatch to `SpaceAGORA.ParallelProcess`, a small worker
pool built on `Distributed`:

```julia
pool = campaign_process_pool()
worker_ids = ensure_process_workers!(pool, 8)
```

- `campaign_process_pool()` returns the process-global `ProcessPool` that
  parallel campaigns share. Reusing it means each worker is started and
  bootstrapped once per session rather than once per call.
- `ensure_process_workers!(pool, n; warmup_fn=nothing)` grows the pool to at
  least `n` workers, starting only the shortfall with `addprocs`. It loads
  `SpaceAGORA`, `GRAMSuite` (when it can) and the default SPICE kernel set on
  each new worker, and returns the current worker ids. `warmup_fn` is a
  zero-argument closure that mirrors the real per-sample call; the pool runs
  it on each new worker before the first real dispatch, so that compilation
  and initialization for that call happen there. How much it covers depends on
  the workload, the argument types, the Julia environment and the cache state,
  so later calls can still compile.
- `shutdown_process_pool!(pool)` removes every worker with `rmprocs` and
  clears the pool. It is mainly useful in tests: campaign code leaves the
  process-global pool running between calls on purpose.

Each process worker is started with `--threads=1`, so it does not share the
coordinator's Julia thread pool, and thread-based parallelism inside a
worker's own `run_simulation` call does not depend on how many process
workers are active. A campaign that uses a non-default SPICE kernel
directory, or a mission around another body whose kernels are not in the
default set, must furnish its own kernels on the pool's workers before
dispatching, for example with `remotecall_wait` on the ids that
`ensure_process_workers!` returns.

### GRAM atmosphere models in threaded campaigns

By default SpaceAGORA serializes native GRAM calls through one process-wide
lock, so threaded Monte Carlo samples that all query GRAM wait on that one lock
however many threads are available. When every sample builds its own
`GRAMAtmosphereModel` (or receives its own `deepcopy`), set:

```bash
export SPACEAGORA_GRAM_LOCK_SCOPE=model
```

SpaceAGORA then serializes calls on each model instance against that instance
only, so samples with separate instances are not held behind one lock. This
rests on the same premise as the isolated-pool batch path
(`SPACEAGORA_GRAM_ISOLATED_POOL`): distinct GRAM model instances may run
concurrently as long as each single instance is serialized. Do not enable it
when several threads share one model instance; the default `global` scope is
always safe. Process-based campaigns (separate workers started with
`addprocs`) do not need this setting, since each process has its own lock.

### Automatic native GRAM pooling and wind histories

The density callback's isolated pool, which spreads native GRAM queries across
separate model instances, is selected automatically only when Julia has more
than one thread and one callback batch holds at least
`SPACEAGORA_GRAM_ISOLATED_POOL_THRESHOLD` native-GRAM queries (1024 by
default); it uses up to `SPACEAGORA_GRAM_ISOLATED_POOL_MAX_WORKERS` instances
(by default 4, or fewer with fewer threads). When a batch requests winds,
automatic selection also requires nominal winds. The pinned GRAM wrapper's
default wind mode (`auto`) selects perturbed winds, whose values depend on
each model instance's query history, so such batches keep the locked route,
including look-ahead cache sampling. Density-only queries (`wind=false`,
including every query of a run with `EnvironmentModel.wind = false`) stay
eligible.

Setting `SPACEAGORA_GRAM_ISOLATED_POOL=on` explicitly allows separate
stochastic histories. With perturbed winds it can change wind diagnostics and
trajectories as the pool width or the thread count changes. Routing never
switches the requested wind mode or seed to run faster.

### Real GRAM without the vacuum-predicted cache

`SPACEAGORA_VACUUM_GRAM_CACHE`, off by default, answers density queries from
splines of log density and temperature fitted along a predicted drag-free
(two-body plus J2) trajectory, and rebuilds them when the actual trajectory
strays from the prediction. Its knot count, look-ahead and deviation threshold
are `SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS` (default 20),
`SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S` (default 600 s) and
`SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M` (default 5000 m). It supports
per-satellite GRAM queries at constellation scale. Without it, and with the
default per-stage sampling, the dynamics query native GRAM at every stage, and
an adaptive solver that reacts to query-dependent perturbations makes many
such queries. Freezing the sample for a step is an explicit alternative, whose
accuracy must be checked for the application.

`SPACEAGORA_DENSITY_FREEZE_PER_STEP` defaults to `off`. An unset or empty
value, `0` or `off` keeps per-stage sampling. Set `auto` explicitly to freeze
only when the run's density model has history-dependent winds and
`environment_model.wind` is true. Set `1` or `on` to freeze every density
model, including nominal-wind GRAM. Neither choice is ever selected
automatically.

Under explicit `auto`, a `GRAMAtmosphereModelSurrogate` counts as native GRAM
only when it has a configured `point_fallback_below_m`, below which it queries
native GRAM point by point. The default Earth, Mars and Venus surrogates have
none, so a surrogate run answered from its grid keeps per-stage sampling and a
threaded density callback. A query outside the grid, which falls back to
native GRAM on a path that warns once, does not change that.

```bash
# Explicitly choose frozen sampling for history-dependent winds.
export SPACEAGORA_DENSITY_FREEZE_PER_STEP=auto
```

`run_simulation` runs a density callback once per accepted solver step. With
freezing enabled, the right-hand side uses that sample at every stage of the
step instead of requesting atmosphere values at each stage time.

Holding an atmosphere sample over a step adds a numerical approximation. The
adaptive solver's local error estimate does not measure that sampling error,
so tightening the tolerances alone does not establish accuracy. Check the
outputs your application needs while reducing the step-size limits.

`dt_max_atmosphere` applies from the start whenever any active spacecraft is
inside `EnvironmentModel.EI`, and stays in force until every active spacecraft
is outside. Outside that boundary, `dt_max_orbit` sets the maximum step. A
frozen sample is still held outside the boundary, so refining only the
atmospheric cap can leave an error contribution from the orbital cap. Check
the actual step sizes and refine both limits when needed.

With nominal winds, a per-stage reference and step refinement can separate
the hold approximation from integration error. With perturbed winds, changing
the number or order of native queries changes the random-wind sequence, even
at the same seed, so a same-seed comparison of frozen and per-stage sampling
includes both effects, and the spread between a few seeds is not an accuracy
acceptance limit. Establish application-specific limits on position,
velocity, apsides and thermal outputs before accepting a sampling policy; a
check made for one mission configuration or phase policy does not establish
limits for another.

The default keeps per-stage sampling. Enabling freezing can change
GRAM-backed trajectories and thermal outputs, even at the same seed. Set
`SPACEAGORA_DENSITY_FREEZE_PER_STEP=0` to return to per-stage sampling.

A freeze in effect takes precedence over the vacuum-predicted cache
(`SPACEAGORA_VACUUM_GRAM_CACHE=1`): the right-hand side then reads the
once-per-step sample rather than interpolating the look-ahead spline at every
stage, and only the per-step density callback can consult or rebuild the
spline. This is deliberate. A cache rebuild is a burst of native queries, and
rebuilds issued from threaded right-hand-side stages would advance GRAM's
perturbed-wind walk in thread-scheduling order. Under explicit `auto` this
applies only to runs whose winds are history-dependent (below); with nominal
winds, or with `0`, the right-hand side interpolates the spline at each stage.

### Reproducible perturbed winds

GRAM's perturbed winds are a correlated random walk that advances on every
native query of an atmosphere instance, so the winds a query returns depend on
every earlier query on that instance. Initialization and the explicit
sampling choice affect reproducibility:

- A run with history-dependent winds reseeds its GRAM models (the configured
  model, per-satellite instances and isolated-pool instances) with their
  configured `seed` immediately before the solve, after every pre-solve probe
  (right-hand-side calibration, contention probe, callback-width
  calibration), and clears the GRAM track and look-ahead caches. Reseeding
  puts native GRAM back into its first-update state, whose one-time
  initialization is not safe to run concurrently on several instances, so
  each reseeded model is then warmed with one serial query before the solve
  starts. This removes the dependence on earlier probes and previous runs, but
  does not by itself fix a query order that changes with thread scheduling.
  Runs without history-dependent winds are not reseeded.
- With history-dependent winds, the per-step density callback evaluates the
  spacecraft serially, in index order, whatever
  `SPACEAGORA_DENSITY_CALLBACK_PARALLEL` says. Explicitly selecting `auto` or
  `1` also makes that callback supply the dynamics' native atmosphere
  samples, so per-stage right-hand-side queries, which would otherwise advance
  the walk in thread-scheduling order and in a number that depends on the
  solver route, read that sample instead. The results are then independent of
  the thread count.
- Native Earth-GRAM loads part of its wind-variability tables lazily and reads
  them before loading them on the first atmosphere update in a process. The
  GRAMSuite extension therefore loads them once per process before any model
  is used.

Per-stage sampling, the default, keeps the dependence on query history and
thread scheduling for perturbed winds. The repository's thread-count
reproducibility tests cover only the explicit frozen route. Custom guidance
or control code that issues additional native queries must preserve their
order too; thread safety alone does not make the winds reproducible.

## Constellation ensembles

For multi-satellite configurations whose members do not interact (no
inter-satellite links, no coordinated GNC), `run_constellation_ensemble`
splits the constellation into independent single-satellite propagations and
runs them in parallel across satellites, as a Monte Carlo campaign runs its
samples:

```julia
# Planner-chosen, from the configuration's SolverConfig(parallel=true):
result = run_constellation_ensemble(args; return_solution=true)
# Or a fixed worker count, with the flag off:
result = run_constellation_ensemble(args; threads=8, return_solution=true)
solutions = [s.value for s in result.successful]
```

Compared with propagating the constellation as one coupled state vector, each
satellite goes to a worker once for its whole propagation, instead of the
satellites being dispatched across threads at every time step, and each
satellite keeps its own adaptive step size instead of following the smallest
step of the whole constellation.

On the `:process` outer route each member runs in a worker process with its
own environment variables, GRAM lock and thread pool, so members share no
process-global state. With `SolverConfig(parallel=true)` the planner can pick
`:process` itself and start a `Distributed` worker pool with
`ensure_process_workers!`; no manual `addprocs` call or cluster setup is
needed for this path. See
[Process-backend outer parallelism](#Process-backend-outer-parallelism)
above. Multi-node and scheduler launches for benchmark and study scripts are
a separate path with explicit `addprocs`; see
[Distributed and HPC](../distributed_hpc.md).

The runner refuses a configuration with interlinks, with externally
propagated spacecraft, or with configuration-level guidance, navigation or
control effectors, because effectors that coordinate satellites cannot act
across ensemble members. If every configuration-level effector acts on each
satellite on its own (`spacecraft_scope` `AllSpacecraft()`), opt in with
`allow_gnc_effectors=true`; GNC for one satellite is declared on its
`SpacecraftModel` and needs no flag. A dynamic effector bound to one
spacecraft by its run position, such as `ReactionWheelMomentumModel`, is
refused as well, because every member runs its spacecraft alone at position 1.
Keep the single `run_simulation` path for coupled constellations (RPO,
formation control).

## Optional native GRAM density workers

`SPACEAGORA_GRAM_PROCESS_POOL` defaults to `off`. For a large batch of
satellites using native GRAM, `auto` starts a separate pool of density worker
processes once a batch reaches `SPACEAGORA_GRAM_PROCESS_POOL_THRESHOLD`
satellites (default 64); `on` uses it for smaller batches too.
`SPACEAGORA_GRAM_PROCESS_POOL_WORKERS` defaults to half the logical CPU count,
with a minimum of one; set it explicitly to control memory use and the number
of processes. This service does not apply to surrogate models, and it declines
work inside an outer parallel run.

Each worker process owns its own native GRAM instance, so the wind-history
rule of the in-process isolated pool
([Automatic native GRAM pooling and wind histories](#Automatic-native-GRAM-pooling-and-wind-histories))
applies here too. Under `auto`, a run that requests winds while they depend on
query history (the pinned wrapper's default perturbed mode) keeps local
evaluation; runs with `EnvironmentModel.wind = false` or
`SPACEAGORA_GRAM_WIND_MODE=nominal` stay eligible. `on` is an explicit opt-in
to separate stochastic histories: with perturbed winds it changes the winds,
and so the trajectory, relative to local evaluation, and the result can
depend on the worker count.

Workers rebuild the model from its saved constructor settings, including the
epoch and the resolved data paths, and reuse a model only while that recipe
matches. Each batch carries the recipe, so overlapping requests use the right
configuration. Workers whose model setup fails are excluded, and the remaining
workers can still serve the batch. If no worker is usable or a query fails,
the batch falls back to local evaluation. A model wrapped directly from a
native core has no saved recipe and uses local evaluation without starting
density workers. A failed replacement clears the worker's cache, so the next
valid request rebuilds its native model.

GRAM construction does not transfer the coordinator's SPICE kernel pool. The
process bootstrap loads the default Earth kernels; missions around other
bodies, or custom kernel sets, need the caller to furnish the kernels on the
density workers before dispatch, and unavailable ephemerides cause local
fallback. A recipe whose worker setup or batch fails is remembered for the
rest of the session, and later batches with the same settings fall back
locally at once, without rebuilding native models on every worker. Recovery
is explicit: change the configuration (a different recipe is tried again),
call `SpaceAGORA.ParallelProcess.clear_density_service_failures!()` after
correcting the worker setup, or restart the pool with
`SpaceAGORA.ParallelProcess.shutdown_density_workers!()`, which also forgets
remembered failures. `SpaceAGORA.ParallelProcess.density_service_failures()`
lists them.

The workers reproduce the construction settings, not a random stream that has
already advanced or manual changes to a native handle. Runtime environment
settings such as `SPACEAGORA_GRAM_WIND_MODE` are not part of the recipe: set
them before the workers start. Identical recipes can share worker state
across requests; this pool does not isolate stochastic streams by
simulation. Compare your results against local evaluation before using it
for a study. Per-satellite model instances and stateful interpolation caches
keep their existing exclusions.

## When not to use parallelism

- **Debugging**: leave the flag off (and `threads` unset) to rule out
  concurrency as a source of non-determinism.
- **Single-satellite no-GRAM runs**: a simple orbit-only scenario has little
  work for inner parallelism to split. Start serial and compare.
- **Correctness validation**: run the same scenario with the flag off and on
  before relying on the parallel result.

## HPC and process worker setup

For multi-node or process-worker execution, see
[Distributed and HPC](../distributed_hpc.md) for the `SPACEAGORA_PERF_PROCS`,
`SPACEAGORA_PERF_WORKER_PROJECT` and `SPACEAGORA_PERF_MACHINE_LABEL`
environment variables.

## Internal and benchmark controls

The route, profile and environment-variable controls behind
`SolverConfig(parallel=true)` serve benchmarking and performance diagnosis.
They are not a supported interface and may change without notice. They are
documented for contributors in the repository, in
`docs/contributing/parallel_internals.md`.
