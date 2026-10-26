# Parallel Execution

Use this page when you want a simulation, a Monte Carlo campaign or a
constellation ensemble to use more than one core.

Shortest successful command:

```text
julia --project=examples --threads=auto examples/AGORA_Earth_MonteCarlo.jl
```

What to read next:

- [Distributed and HPC](../distributed_hpc.md)
- [Simulation Configuration](simulation_configuration.md)

## Turning on parallel execution

Parallel execution is one switch in the solver setup:
`SolverConfig(parallel=true)`. With it set, SpaceAGORA decides how to use
the threads Julia was started with. A single run threads its callbacks,
effectors and right-hand side where its adaptive policy predicts a gain; a
campaign lets the predictive planner choose serial, threaded or
process-worker execution for the campaign's shape before the first sample
runs. Nothing else needs configuring.

A single simulation:

```julia
using SpaceAGORA

args = SimulationConfiguration(
    # ... the usual fields ...
    solver_config=SolverConfig(parallel=true),
)
run_simulation(args)
```

A campaign. A constellation ensemble reads the flag from the configuration
it is given; `run_monte_carlo` takes it as a keyword, since it only sees your
function:

```julia
result = run_constellation_ensemble(args; return_solution=true)

result = run_monte_carlo(1:100; parallel=true) do seed
    run_simulation(make_config_for_seed(seed))
end
```

Requirements and behavior:

- Start Julia with threads, for example `julia --threads=auto`. With a
  single thread the flag can only use process workers for campaigns, and
  says so once.
- The first parallel run on a machine measures the machine's cost
  constants and stores them in
  `output/parallel_policy_state/cost_constants_<fingerprint>.toml` under
  the working directory. Calibration duration depends on the machine, thread
  count and compilation state. Later runs reuse valid constants; a wider Julia
  thread count or an unusable cache can trigger calibration, since measured
  speedups apply only up to the recorded thread count. Processes that start together in one working directory measure
  once: the others wait for the file (a `.lock` file beside it marks the
  measurement in progress). To redo the measurement, delete the file or run
  `julia --project=. --threads=<T> scripts/calibration/calibrate_machine.jl`.
- The flag's settings apply for the duration of the call only. Nothing is
  left in the process environment afterwards, including when the run
  throws or when several flagged runs overlap on different tasks. While any
  flagged run is active, the flag's settings are visible to every task in
  the process, since they live in the process environment.
- `run_monte_carlo(...; parallel=true)` only sees your function, so it
  does not send samples to process workers unless you describe the
  workload with `route_features` (see
  [Internal and benchmark controls](#Internal-and-benchmark-controls)).
  A process worker is a separate Julia process: a function defined in your
  script, anything it reaches in your script's `Main`, and SPICE kernels
  you furnished yourself exist there only if you load them there. Pass
  `route_features` only once the workers are prepared. A constellation
  ensemble derives its features from its configuration and can use
  process workers; a configuration that holds types or functions defined
  in your script needs the same preparation, or `parallel=false` with a
  fixed `threads=N`.
- The flag picks the worker count itself, so `parallel=true` together with
  an integer `threads=` is an `ArgumentError`. Every member of a campaign
  must carry the same flag.
- With `parallel=false` (the default) nothing changes: a single run behaves
  exactly as before, and campaigns are serial unless you give a fixed
  `threads=N` (see below).
- A run or campaign started inside another campaign's worker keeps the
  enclosing campaign's thread budget rather than claiming the machine
  again.

Parallel runs are not bit-for-bit reproducible across thread counts in
every configuration (see
[Internal and benchmark controls](#Internal-and-benchmark-controls) for the
route that is). Validate a study against a run with the flag off before
relying on it.

## Monte Carlo campaigns

Use `run_monte_carlo` when you want to run many independent simulations from
an example script or study without copying benchmark-specific orchestration
code. The runner applies parallelism across Monte Carlo samples: each worker
task calls your function with one seed. `parallel=true` (above) lets the
planner choose; this section describes a fixed worker count instead.

Start Julia with the number of threads you want available:

```bash
julia --project=examples --threads=8 examples/AGORA_Earth_MonteCarlo.jl
```

Then cap the campaign runner to that same thread count, or a smaller count:

```julia
using SpaceAGORA

result = run_monte_carlo(1:100; threads=8) do seed
    args = make_config_for_seed(seed)
    run_simulation(args; return_solution=true)
end

println("successful samples: $(length(result.successful))")
println("failed samples: $(length(result.failed))")
```

`threads` does not create new Julia threads at runtime. If the script requests
more Monte Carlo threads than Julia was launched with, the runner throws a
clear `ArgumentError`.

For reusable scripts, build a spec explicitly:

```julia
spec = MonteCarloSpec(seeds=1001:1100, threads=8, fail_fast=false)

result = run_monte_carlo(spec) do seed
    args = make_config_for_seed(seed)
    run_simulation(args)
end
```

By default, failed samples are captured in `result.failed` and do not stop the
rest of the campaign. Set `fail_fast=true` when you want the runner to rethrow
on the first failed sample instead.

### Process-backend outer parallelism

When the campaign planner picks `:process` (or when you want to control it
directly), campaigns dispatch to `SpaceAGORA.ParallelProcess`, a small
`Distributed`-based worker pool built for this purpose:

```julia
pool = campaign_process_pool()
worker_ids = ensure_process_workers!(pool, 8)
```

- `campaign_process_pool()` returns the process-global `ProcessPool` that
  parallel campaigns already share; reusing it avoids repaying the
  one-time `SpaceAGORA`/`GRAMSuite` worker precompilation cost on every call.
- `ensure_process_workers!(pool, n; warmup_fn=nothing)` grows the pool to at
  least `n` workers (spawning only the shortfall via `addprocs`), bootstraps
  each new worker with `SpaceAGORA`, `GRAMSuite` (best-effort), and the
  default SPICE kernel set, and returns the current worker ids. Pass
  `warmup_fn`, a zero-argument closure that mirrors the real per-sample call,
  to attempt representative compilation and initialization on each new worker
  before the first timed dispatch. Its cost and coverage depend on the workload,
  argument types, Julia environment and cache state; later calls can still
  require additional compilation.
- `shutdown_process_pool!(pool)` removes every worker via `rmprocs` and
  clears the pool; mainly useful for tests, since campaign code otherwise
  leaves the process-global pool warm across calls by design.

Each process worker is started with `--threads=1`, so it does not share the
coordinator's Julia thread pool: inner thread-based parallelism inside a
worker's own `run_simulation` call is unaffected by how many process workers
are active. A campaign built on a non-default SPICE kernel directory, or a
non-Earth-primary mission whose kernels aren't in the shared default set,
must still furnish its own kernels on the pool's workers (for example via
`remotecall_wait` on `ensure_process_workers!`'s return value) before
dispatching.

### GRAM atmosphere models in threaded campaigns

Native GRAM calls are serialized through a single process-wide lock by default,
so threaded Monte Carlo samples that all query GRAM contend on one lock no
matter how many threads are available. When every sample builds its own
`GRAMAtmosphereModel` (or receives its own `deepcopy`), set:

```bash
export SPACEAGORA_GRAM_LOCK_SCOPE=model
```

so each model instance serializes only against itself and samples evaluate
concurrently. This relies on the same instance-isolation premise as the
isolated-pool batch path (`SPACEAGORA_GRAM_ISOLATED_POOL`): distinct GRAM model
instances may run concurrently as long as any single instance is serialized. Do
not enable it if several threads share one model instance and you have not
measured the workload — the default `global` scope is always safe. Process-based
campaigns (separate workers via `addprocs`) do not need this: each process has
its own lock already.

### Automatic native GRAM pooling and wind histories

The isolated pool is selected automatically only above its native-query
threshold. When a batch requests winds, automatic selection also requires
nominal winds: the pinned wrapper's default (`auto`) selects perturbed winds,
whose values depend on each model instance's query history. Such batches keep
the locked route, including look-ahead cache sampling. Density-only queries
(`wind=false`, including every query of a run with `EnvironmentModel.wind =
false`) keep their existing eligibility.

Setting `SPACEAGORA_GRAM_ISOLATED_POOL=on` explicitly allows separate stochastic
histories. With perturbed winds it can change wind diagnostics and trajectories
as pool width or thread count changes. Routing never switches the requested
wind mode or seed to obtain a speedup.

### Real GRAM without the vacuum-predicted cache

`SPACEAGORA_VACUUM_GRAM_CACHE` (the drag-free trajectory spline described
above) supports per-satellite GRAM queries at constellation scale. Direct,
uncached native queries can be expensive when an adaptive solver reacts to
query-dependent perturbations. Frozen sampling is an explicit alternative;
its accuracy must be checked for the application.

`SPACEAGORA_DENSITY_FREEZE_PER_STEP` defaults to `off`. An unset or empty
value, `0`, or `off` keeps per-stage sampling. Set `auto` explicitly to freeze
only when the run's density model has history-dependent winds and
`environment_model.wind` is true. Set `1` or `on` to freeze every density model,
including nominal-wind GRAM. Neither opt-in is selected automatically.

Under explicit `auto`, a `GRAMAtmosphereModelSurrogate` counts as native
GRAM only when it has a configured `point_fallback_below_m`, below which it
queries native GRAM point by point. The default Earth, Mars and Venus
surrogates have none, so a surrogate run answered from its grid keeps
per-stage sampling and a threaded density callback. A query outside the grid,
which falls back to native GRAM on a warn-once path, does not change that.

```bash
# Explicitly choose frozen sampling for history-dependent winds.
export SPACEAGORA_DENSITY_FREEZE_PER_STEP=auto
```

`run_simulation` fires a density callback once per accepted solver step.
With freezing enabled, the RHS uses that sample at every stage of the step
instead of requesting atmosphere values at each stage time.

Holding an atmosphere sample over a step introduces an additional numerical
approximation. The adaptive solver's local error estimate does not measure
that sampling error, so tightening its tolerances alone does not establish
accuracy. Check the outputs needed by the application while reducing the
sampling timestep limits.

`dt_max_atmosphere` applies from startup when any active spacecraft is inside
`EnvironmentModel.EI`, and remains active until every active spacecraft is
outside. Outside that boundary, `dt_max_orbit` controls the maximum step.
A frozen sample is still held outside the boundary, so refining only the
atmospheric cap may leave an error contribution from the orbital cap.
Check the actual step sizes and refine both limits when needed.

With nominal winds, a per-stage reference and step refinement can separate the
hold approximation from integration error. With perturbed winds, changing the
number or order of native queries changes the random-wind sequence, even at
the same seed. A same-seed frozen/per-stage comparison therefore includes both
effects; the spread between a few seeds is not an accuracy acceptance limit.
Establish application-specific position, velocity, apsis and thermal limits
before accepting a sampling policy. Earlier one-off measurements do not
establish those limits for a different phase policy or mission configuration.

The default retains per-stage sampling. Explicitly enabling freezing can
change GRAM-backed trajectories and thermal outputs, even at the same seed.
Set `SPACEAGORA_DENSITY_FREEZE_PER_STEP=0` to select per-stage sampling again.

A freeze in effect takes precedence over the vacuum-predicted cache
(`SPACEAGORA_VACUUM_GRAM_CACHE=1`): the RHS then reads the once-per-step
sample rather than interpolating the look-ahead spline at every stage, and
only the per-step density callback can consult or rebuild the spline. This is
deliberate. A cache rebuild is a burst of native queries, and rebuilds issued
from threaded RHS stages would advance GRAM's perturbed-wind walk in
thread-scheduling order. Under explicit `auto` this applies only to runs whose winds are
history-dependent (below); with nominal winds, or with `0`, the RHS
interpolates the spline per stage as before.

### Reproducible perturbed winds

GRAM's perturbed winds are a correlated random walk that advances on every
native query of an atmosphere instance, so the winds a query returns depend on
every earlier query on that instance. Initialization and the explicit sampling
choice affect reproducibility:

- Each such run reseeds its GRAM models (the configured model, per-satellite
  instances and isolated-pool instances) with their configured `seed`
  immediately before the solve, after every pre-solve probe (RHS calibration,
  contention probe, callback-width calibration), and clears the GRAM track and
  look-ahead caches. Reseeding puts native GRAM back into its first-update
  state, whose one-time initialization is not safe to run concurrently on
  several instances, so each reseeded model is then warmed with one serial
  query before the solve starts. This removes dependence on earlier probes
  and previous runs, but does not by itself fix a query order that changes
  with thread scheduling. Runs without history-dependent winds are not
  reseeded.
- With perturbed winds, explicitly selecting `auto` or `1` makes the
  per-step density callback supply the dynamics' native atmosphere samples, and
  that callback then evaluates spacecraft serially in index order whatever
  `SPACEAGORA_DENSITY_CALLBACK_PARALLEL` says. Per-stage RHS queries, which
  would otherwise advance the walk in thread-scheduling order and in a number
  that depends on the solver route, read that sample instead. Results are then
  independent of the thread count.
- Native Earth-GRAM loads part of its wind-variability tables lazily and reads
  them before loading them on the first atmosphere update in a process; the
  GRAMSuite extension loads them once per process before a model is used.

Per-stage sampling, including the default, retains query-history and
thread-scheduling dependence for perturbed winds and can be slow. The
one/four-thread regression guarantees cover the explicit frozen dynamics
route. Custom guidance or control code that issues additional native queries
must preserve their order too; an assumption of thread safety does not
establish stochastic reproducibility.

## Constellation ensembles

For multi-satellite configurations whose members do not interact (no
inter-satellite links, no coordinated GNC), `run_constellation_ensemble` splits
the constellation into independent single-satellite propagations and applies
Monte Carlo-style outer parallelism across satellites:

```julia
# Planner-chosen, from the configuration's SolverConfig(parallel=true):
result = run_constellation_ensemble(args; return_solution=true)
# Or a fixed worker count, with the flag off:
result = run_constellation_ensemble(args; threads=8, return_solution=true)
solutions = [s.value for s in result.successful]
```

Compared to propagating the constellation as one coupled state vector, this
dispatches each satellite to a worker once for its entire propagation (instead
of paying per-timestep thread dispatch across satellites) and lets each
satellite keep its own adaptive step size (instead of forcing every satellite
to the global minimum step).

Current limitation: with in-process worker threads (the `:threads` outer
route), per-step environment-variable configuration reads in the RHS and
callback plumbing serialize concurrent members (Julia `ENV` access is
process-global), which can erase the outer-parallel gain for light dynamics.

The `:process` outer route avoids this: each satellite is dispatched to its
own OS process with its own `ENV`, GRAM lock, and thread pool, so members
never contend on process-global state. With `SolverConfig(parallel=true)` the
planner picks `:process` itself for workload shapes where it wins,
auto-bootstrapping a `Distributed` worker pool via `ensure_process_workers!` —
no manual `addprocs` call or cluster setup is required for this
library-level path. See
[Process-backend outer parallelism](#Process-backend-outer-parallelism)
below. Multi-node/scheduler launches for benchmark and study scripts are a
separate, explicit-`addprocs` path — see
[Distributed and HPC](../distributed_hpc.md).

The runner refuses configurations with guidance, navigation, or control
effectors, because effectors that coordinate satellites cannot act across
ensemble members. If every configuration-level effector acts on each satellite on
its own (`spacecraft_scope` `AllSpacecraft()`), opt in with
`allow_gnc_effectors=true`; GNC for one satellite is declared on its
`SpacecraftModel` and needs no flag. Keep the monolithic `run_simulation` path
for genuinely coupled constellations (RPO, formation control).

## Optional native GRAM density workers

`SPACEAGORA_GRAM_PROCESS_POOL` defaults to `off`. For a large batch of
satellites using native GRAM, `auto` enables a separate density-worker pool
at or above `SPACEAGORA_GRAM_PROCESS_POOL_THRESHOLD` (default 64); `on` also
enables it for smaller batches. `SPACEAGORA_GRAM_PROCESS_POOL_WORKERS` defaults
to half the logical CPU count, with a minimum of one. Set it explicitly to
control memory and process cost. This service does not apply to
surrogate models, and it declines work inside an outer parallel run.

Each worker process owns its own native GRAM instance, so the wind-history rule
of the in-process isolated pool (Automatic native GRAM pooling and wind
histories, above) applies here too. Under `auto`, a run that requests winds while they depend on
query history (the pinned wrapper's default perturbed mode) keeps local
evaluation; runs with `EnvironmentModel.wind = false` or
`SPACEAGORA_GRAM_WIND_MODE=nominal` stay eligible. `on` is an explicit opt-in to
separate stochastic histories: with perturbed winds it changes the winds, and so
the trajectory, relative to local evaluation, and the result can depend on the
worker count.

Workers reconstruct the model from its saved constructor settings, including
the epoch and resolved data paths. They reuse a model only while that recipe
matches. Each batch carries the recipe so overlapping requests use the right
configuration. Workers whose model setup fails are excluded; the remaining
workers can still serve the batch. If no worker is usable or a query fails, the
batch falls back to local evaluation. A model wrapped directly from a native core has no saved recipe
and uses local evaluation without starting density workers. A failed replacement
clears the worker cache so the next valid request rebuilds its native model.

GRAM construction paths do not transfer the coordinator's SPICE kernel pool.
The existing process bootstrap loads default Earth kernels. Non-Earth missions
or custom kernel sets still require the caller to furnish the needed kernels on
the density workers before dispatch; unavailable ephemerides cause local fallback.
A recipe whose worker setup or batch fails is remembered for the rest of the
session: later batches with the same settings fall back locally at once, without
rebuilding native models on every worker. Recovery is explicit. Change the
configuration (a different recipe is tried once more), call
`SpaceAGORA.ParallelProcess.clear_density_service_failures!()` after correcting
the worker setup, or restart the pool with `shutdown_density_workers!()`, which
also forgets remembered failures. `density_service_failures()` lists them.

This preserves construction settings, not an already advanced random stream or
manual changes to a native handle. Runtime environment policy, such as
`SPACEAGORA_GRAM_WIND_MODE`, is not part of the recipe: set it before workers
start. Identical recipes can share worker state
across requests; this pool does not isolate stochastic streams by simulation.
Compare your results against local evaluation before using it for a study.
Per-satellite model instances and stateful interpolation caches retain their
existing exclusions. No speedup is guaranteed for small batches.

## When not to use parallelism

- **Debugging**: leave the flag off (and `threads` unset) to eliminate
  concurrency as a source of non-determinism.
- **Single-satellite no-GRAM runs**: inner parallelism often adds overhead that
  exceeds the benefit for simple orbit-only scenarios. Start serial and
  measure.
- **Correctness validation**: run the same scenario with the flag off and on
  before relying on the parallel result.

## HPC and process worker setup

For multi-node or process-worker execution, see [Distributed and HPC](../distributed_hpc.md)
for guidance on the `SPACEAGORA_PERF_PROCS`, `SPACEAGORA_PERF_WORKER_PROJECT`,
and `SPACEAGORA_PERF_MACHINE_LABEL` environment variables.

## Internal and benchmark controls

The route, profile and environment-variable controls behind
`SolverConfig(parallel=true)` serve the paper benchmark harness and
performance diagnosis. They are not a supported interface and may change
without notice. They are documented for contributors in the repository, in
`docs/contributing/parallel_internals.md`.
