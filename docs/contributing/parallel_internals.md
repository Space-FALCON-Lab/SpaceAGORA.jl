# Parallel runtime internals

This page documents the machinery that `SolverConfig(parallel=true)` drives:
the routing and policy settings, the named profiles and the environment
variables behind them. They are kept for the paper benchmark harness, which
needs fixed routes and older profiles to reproduce its measurements, and for
diagnosing performance. None of it is a supported interface, and it may
change without notice; users should start from
[Parallel Execution](../src/user/parallel_execution.md). The profile functions
are not exported; reach them as `SpaceAGORA.ParallelProfiles.<name>`.

## Runtime controls

The parallel runtime is controlled along these axes:

| Axis | Values | Set with | Purpose |
|---|---|---|---|
| Outer route | `none`, `threads`, `process`, `auto` | the campaign runner's `threads` keyword (an integer or `:auto`) or `parallel=true`; `SPACEAGORA_PERF_PARALLEL_BACKEND` for the benchmark harness | Whether campaign or constellation work runs serially, on Julia threads, or on process workers |
| Callback and effector modes | `off`, `auto`, `on` | `SPACEAGORA_DENSITY_CALLBACK_PARALLEL`, `SPACEAGORA_CONTROL_CALLBACK_PARALLEL`, `SPACEAGORA_THERMAL_CALLBACK_PARALLEL`, `SPACEAGORA_MULTIBODY_PARALLEL`, `SPACEAGORA_EFFECTOR_PARALLEL` | Inner parallelism of the density, control and thermal callbacks, the multibody work and the dynamic effectors |
| RHS execution mode | `auto`, `serial`, `satellite`, `per_satellite`, `flat` | `SPACEAGORA_RHS_EXECUTION_MODE` | The dominant threading strategy of the right-hand side |
| Inner scheduler | `static`, `dynamic` | `SPACEAGORA_PARALLEL_POLICY_INNER_SCHEDULER` | Fixed worker assignment or dynamic work stealing |
| Adaptive policy | `0`, `1` | `SPACEAGORA_PARALLEL_POLICY_ADAPTIVE` | Lets measured run times adjust worker counts and route choices |
| Persistent hints | `0`, `1` | `SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS` | Reuses measured policy choices across repeated runs when state persistence is on |

`SPACEAGORA_PERF_PARALLEL_BACKEND` and `SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE`
are read only by the benchmark harness
(`benchmarks/studies/performance_runtime_analysis/`); a profile sets them for
it. The package itself takes the outer route from the campaign runners'
keywords.

On a machine with many cores running a multi-satellite campaign, start with
`auto` outer routing, `auto` inner modes and the default `auto` RHS execution
mode. Turn on the adaptive policy and persistent hints for repeated
performance studies only after checking the results against a forced serial
run.

## Parallel profiles

A profile is a named bundle of the settings above, defined in
`src/orchestration/parallel/routing/profile_definitions.jl`.
`SolverConfig(parallel=true)` runs under `R7`. Apply a profile explicitly with
`SpaceAGORA.ParallelProfiles.with_parallel_profile(f, "R6")`, or with
`ParallelConfig(profile="R6")` in a `SimulationEngineConfig`; both expand the
whole bundle. Setting `SPACEAGORA_PARALLEL_PROFILE` on its own does not apply
a profile. Alone it only names the persisted policy state files and, for
`R0`, forces a serial right-hand side.

| Profile | Outer route | Inner modes | What it adds |
|---|---|---|---|
| `R0` | serial | off | The true serial baseline: nothing is threaded, at any thread count |
| `R1_a` | threads | off | Outer parallelism only, on Julia threads |
| `R1_b` | process | off | Outer parallelism only, on worker processes |
| `R2` | serial | auto | Inner parallelism only |
| `R3` | auto | auto | Heuristic outer route with static inner widths |
| `R4` | auto, adaptive | auto | The outer-route bandit and the adaptive inner policy |
| `R5` | auto, adaptive | auto | R4 plus persistent hints and the measured-reward width chooser |
| `R6` | auto, adaptive | auto | R5 plus `SPACEAGORA_PARALLEL_POLICY_V2` |
| `R7` | auto, adaptive | auto | R6 plus `SPACEAGORA_CAMPAIGN_PLANNER=predictive` |

Pre-solve RHS plan calibration is not part of any profile. It runs, unless
`SPACEAGORA_RHS_CALIBRATE=off`, on every run with more than one inner thread,
some dynamic effectors and at least two active spacecraft
(`_calibrate_rhs_plan_if_needed!` in `engine/rhs_calibration.jl`), and reuses
a cached result when one exists.

### R6: the revised policy

Everything `R6` changes sits behind one switch,
`SPACEAGORA_PARALLEL_POLICY_V2`, so that with the switch off it is exactly
`R5`, and a paired comparison of the two measures the switch alone. Under the
switch:

- A Monte Carlo campaign goes to the outer route that finishes it in fewer
  rounds of concurrent samples. The process route counts its affordable
  workers plus the coordinator's spare threads, which it also fills with
  samples (mixed dispatch; `SPACEAGORA_PARALLEL_MIXED_DISPATCH=0` turns that
  off); the threads route counts the thread budget. `SPACEAGORA_PERF_PROCS`
  caps the pool. A tie goes to the process route, and the router then spends
  two campaigns measuring the other route; that is the only trial it forces.
  Otherwise it does not spend campaigns trying the other route and switches
  only on history it already holds. With mixed dispatch off, the process route
  wins when its workers match or exceed the thread budget. A campaign too small
  to amortize a dispatch round stays on threads, and native GRAM point density
  takes the process route, as it does without the switch.
- The process pool is sized by memory as well as cores. Each worker is
  charged the coordinator's own resident set (never under 2 GiB) plus a
  per-spacecraft term for native GRAM (2 MiB by default, nothing for the
  offline surrogate), against the machine's memory or its cgroup limit, less
  a reserve and less what the coordinator already holds. The process route is
  offered only when at least two workers fit, and the split ladder stops at
  the number that fit. `SPACEAGORA_MEMORY_BUDGET_GB`,
  `SPACEAGORA_PERF_WORKER_MEMORY_GB` and `SPACEAGORA_GRAM_SAT_MEMORY_MB`
  override the three terms. The per-spacecraft default comes from short
  missions that did not return their solutions
  ([GRAM memory footprint](../architecture/gram_memory_footprint.md)); for
  longer missions, denser saving or samples that return their solutions it
  can be too low, so on a tight memory budget measure a worker and set the
  overrides.
- On the threads route, the outer split width is learned by racing the
  candidate widths inside the first campaign of a new workload, with at least
  three samples per worker per width (`_SPLIT_RACE_ROUNDS`). A campaign too
  small to race, and every process or mixed campaign, takes the widest split.
  Widths are never explored one per campaign.
- A cached "the heuristic won" verdict from RHS plan calibration is honoured
  on long solves instead of re-running the sweep in every process.
- Callback and effector widths take the static `min(items, budget)`. The
  density callback also gets a short pre-solve width sweep, and a narrower
  width is pinned only if it beats the static one by the calibration margin.
- The coordinator collects garbage before dispatching a campaign, so a
  threaded campaign's heap does not stall the next campaign's process
  dispatch.

Whether `R6` is faster than `R5` depends on the workload, the launch
configuration and the machine; the lab's recorded comparisons include slower
cases as well as faster ones. Compare the routes for your workload and
machine with `scripts/calibration/paired_profile_probe.jl` and
`scripts/calibration/paired_campaign_probe.jl --src-runner --profile=full`,
and keep the source revision, configuration, calibration state, the cold or
warm distinction and the run records with any timing you report. Use paired
measurements and their observed spread rather than a single number.

### R7: the predictive planner

`R7` is `R6` with the campaign route chosen by a predictive planner instead of
the bandit. The planner scores the candidate routes for the campaign in front
of it (its sample count, spacecraft count, force model and mission length)
against per-machine calibration constants, and dispatches the route it prices
as cheapest. The plan is chosen before the first sample. The dispatch guard
can then narrow a mixed plan: it closes consumers before their next job,
without interrupting a running sample or starting a second dispatch. The
[planner and guard contract](../architecture/predictive_routing_r7.md#the-guard)
explains those conditions and the route switch, which is off by default.
`SPACEAGORA_CAMPAIGN_PLANNER` (`bandit`, the default, or `predictive`) is the
only difference between `R6` and `R7`, so a paired run of the two measures the
planner and nothing else.

The calibration constants are created automatically on the first
`parallel=true` run on a machine, or explicitly with

```bash
julia --project=. --threads=<T> scripts/calibration/calibrate_machine.jl
```

run once per machine, at the thread count you intend to simulate at. It writes
`output/parallel_policy_state/cost_constants_<fingerprint>.toml`, which every
later run on that machine reads instead of measuring anything itself. `R7`
applied through a profile, without the flag, does not calibrate; it then
prices the routes with no model of contention, so run the script before
comparing `R7` with `R6` on a new machine.

## Environment variables

For CLI runs or scripted batches, set the controls before Julia starts. The
examples run in the examples environment; prepare it once with
`julia --startup-file=no examples/setup.jl`. In a POSIX shell:

```bash
export SPACEAGORA_DENSITY_CALLBACK_PARALLEL=auto
export SPACEAGORA_CONTROL_CALLBACK_PARALLEL=auto
export SPACEAGORA_THERMAL_CALLBACK_PARALLEL=auto
export SPACEAGORA_MULTIBODY_PARALLEL=auto
export SPACEAGORA_EFFECTOR_PARALLEL=auto
export SPACEAGORA_RHS_BATCH_PARALLEL=auto
export SPACEAGORA_RHS_EXECUTION_MODE=auto
export SPACEAGORA_RHS_CALIBRATE=auto
export SPACEAGORA_PARALLEL_POLICY_ADAPTIVE=1
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

In PowerShell write each line as `$env:SPACEAGORA_EFFECTOR_PARALLEL = "auto"`,
and in `cmd.exe` as `set SPACEAGORA_EFFECTOR_PARALLEL=auto`.

To force a serial baseline for validation:

```bash
export SPACEAGORA_DENSITY_CALLBACK_PARALLEL=off
export SPACEAGORA_CONTROL_CALLBACK_PARALLEL=off
export SPACEAGORA_THERMAL_CALLBACK_PARALLEL=off
export SPACEAGORA_MULTIBODY_PARALLEL=off
export SPACEAGORA_EFFECTOR_PARALLEL=off
export SPACEAGORA_RHS_BATCH_PARALLEL=off
export SPACEAGORA_RHS_EXECUTION_MODE=serial
export SPACEAGORA_RHS_CALIBRATE=off
export SPACEAGORA_PARALLEL_POLICY_ADAPTIVE=0
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

In automatic mode the effector loop runs serially while at most one satellite
is active, which avoids per-step task overhead.
`SPACEAGORA_EFFECTOR_PARALLEL=on` forces threading when the remaining
eligibility checks allow it. Threaded and serial evaluation give bit-identical
results either way.

That guarantee holds within one RHS execution route. The flat constellation
route (`SPACEAGORA_RHS_EXECUTION_MODE=flat`, also one of the candidates the
start-up calibration sweep may pick for two or more satellites) evaluates the
batched gravity terms through SIMD kernels and sums the contributions in a
different order, so its trajectories agree with the serial and per-satellite
routes to rounding, not bit for bit. For runs that must be reproducible across
machines and thread counts, pin the route and skip the sweep:

```bash
export SPACEAGORA_RHS_CALIBRATE=off
export SPACEAGORA_RHS_EXECUTION_MODE=per_satellite   # or serial
```

To run several scenarios with different settings in one process, scope the
variables with `withenv`. The parallelism and cache settings are read once,
when `run_simulation` starts, into a run-scoped snapshot (together with any
active `SimulationEngineConfig` overrides), so the right-hand side and the
callbacks never read `ENV` during integration. Wrap the whole `run_simulation`
call in `withenv`; changing a variable while a run is in flight does not
affect that run.

Per-step solution storage is kept only when something reads the trajectory:
`return_solution=true`, or the gravity-backbone split mode. Otherwise the
solver runs with `save_everystep=false` and `save_on=false` and keeps only the
endpoints. Writing results does not need it, because the results come from
the run's saving callback. Skipping the storage reduces memory use; how much
depends on the workload. Explicit `SPACEAGORA_SOLVER_SAVE_EVERYSTEP` and
`SPACEAGORA_SOLVER_SAVE_ON` values override this default in either
direction. With storage on, a run holds one copy of the state (plus the
interpolant stages) per accepted step, and the before and after saves of
continuous events; the per-step housekeeping callbacks (planet frame, density
and thermal samples, quaternion projection) add no saves of their own.

## Adaptive campaign routing (`threads=:auto`)

`parallel=true` passes `threads=:auto` to the campaign runners under `R7`.
Passed directly, without the flag, `threads=:auto` uses the campaign planner
the environment selects: `SPACEAGORA_CAMPAIGN_PLANNER`, by default the
outer-route bandit `SpaceAGORA.ParallelProfiles.select_outer_route!`, which
keeps runtime statistics per workload signature. In the fragment below,
`make_config_for_seed` and `args` stand for your own configuration code:

```julia
result = run_monte_carlo(1:100; threads=:auto,
                         route_features=SpaceAGORA.SimulationCampaigns.campaign_route_features(
                             samples=100, n_sats=1,
                             density_family="exponential",
                             mission_time_s=5400.0
                         )) do seed
    run_simulation(make_config_for_seed(seed); return_solution=true)
end

# Constellation ensembles derive their features from the configuration itself:
result = run_constellation_ensemble(args; threads=:auto, return_solution=true)
```

`SpaceAGORA.SimulationCampaigns.campaign_route_features` describes the
campaign shape (sample count, satellites per sample, density-model family,
mission length); its `SimulationConfiguration` method derives those fields
from the configuration. After every campaign the runner records per-sample
success and amortized wall-clock time with
`SpaceAGORA.ParallelProfiles.record_outer_route_feedback!`, so later campaigns
of the same shape can adapt their choice. What they choose depends on the
workload, host load and available history; nothing guarantees convergence to
the fastest route. The history accumulates in the process-global
`SpaceAGORA.SimulationCampaigns.campaign_outer_route_state()`. Inspect it with
`SpaceAGORA.ParallelProfiles.outer_route_stats_snapshot`, reset it with
`SpaceAGORA.ParallelProfiles.reset_outer_route_state!`, or pass an isolated
`SpaceAGORA.ParallelProfiles.OuterRouteState()` through `route_state` (useful
for tests and one-off studies). These names are internal: they are not
exported and may change without notice.

While the adaptive route runs threaded workers, the runner sets
`SPACEAGORA_OUTER_PARALLEL_ACTIVE=1` and, unless you exported one yourself, an
`SPACEAGORA_INNER_THREAD_BUDGET` of `nthreads() ÷ workers`, so per-sample inner
threading and the outer campaign split the thread pool instead of
oversubscribing it. A nested adaptive campaign (an `:auto` campaign running
inside another campaign's worker, where `SPACEAGORA_OUTER_PARALLEL_ACTIVE` is
already set) yields to the enclosing split: it runs serially and records no
feedback, so contended timings never enter the shared route statistics.

The spherical-harmonics gravity SIMD batch runs on every right-hand-side
evaluation rather than once per sample, and its own nested worker batch
contends with outer workers. It therefore runs serially whenever outer
parallelism is active. `SPACEAGORA_HARMONICS_BATCH_ALLOW_WITH_OUTER=1` lets it
split the thread budget with the outer workers; set it only after measuring
that this is faster for your workload.

## Persistent performance hints

Persistent hints let the adaptive inner scheduler record wall-clock
measurements and reuse them on later calls, which saves calibration on the
second and later runs. The state file defaults to
`output/parallel_policy_state/inner_policy_state_<profile>_<machine>_t<threads>.toml`
under the working directory; set a fixed path to keep the hints across Julia
sessions and working directories:

```bash
export SPACEAGORA_PARALLEL_POLICY_STATE_PATH=/path/to/hints.toml
export SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS=1
export SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST=1
export SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD=1
```

`SPACEAGORA_PARALLEL_POLICY_STATE_RESET=1` discards the accumulated hints and
restarts calibration.

## Auditing active controls

For reproducible studies, record the `SPACEAGORA_*` controls that define the
route and inner-policy behaviour:

```julia
for key in sort(collect(keys(ENV)))
    if startswith(key, "SPACEAGORA_")
        println("$key = $(ENV[key])")
    end
end
```

## Optional direct derivative assembly

For a performance experiment on a translational constellation in vacuum, you
can enable direct writes into the derivative buffer (`args` is your
configuration):

```julia
withenv("SPACEAGORA_RHS_FINAL_ASSEMBLY_DIRECT_LAYOUT" => "1") do
    run_simulation(args)
end
```

The option is off by default. It changes only the final assembly step, and only
when the run already uses the flat constellation-effector route; it does not
select that route. It supports the standard position, velocity, mass and
heat-load layout with contiguous `Float64` storage, a `NoAtmosphereModel`, and
no attitude propagation, control effectors or robot arm. Other configurations
and custom layouts take the usual assembly path. The option changes neither the
force calculation nor its accumulation order.

Keep it off for the baseline, then compare results, elapsed time and
allocations with it on, on the same warmed-up workload; the correctness tests
do not measure its effect on run time. Set the value to `"0"` or remove the
variable to restore the default. As with the other runtime controls, set it
before the simulation starts.
