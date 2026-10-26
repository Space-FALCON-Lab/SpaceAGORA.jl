# Parallel runtime internals

This page is for contributors and maintainers. It documents the controls
behind `SolverConfig(parallel=true)` that serve the paper benchmark harness
and performance diagnosis. Users should start from the published
[Parallel Execution](../src/user/parallel_execution.md) page; nothing here is
a supported interface.

## Internal and benchmark controls

Everything below is the machinery `SolverConfig(parallel=true)` drives, kept
for the paper benchmark harness (which needs fixed routes and older profiles
to reproduce its measurements) and for diagnosing performance. It is not the
supported way to parallelize and may change without notice. The profile
functions are not exported; reach them as `SpaceAGORA.ParallelProfiles.<name>`.

### Runtime Controls

SpaceAGORA's parallel runtime is controlled by route and policy settings. The
important axes are:

| Axis | Values | Purpose |
|---|---|---|
| Outer route | `none`, `threads`, `process`, `auto` | Chooses whether campaign or constellation work runs serially, on Julia threads, or on process workers |
| Callback and effector modes | `off`, `auto`, `on` | Controls density, control, thermal, multibody, and dynamic-effector inner parallelism |
| RHS execution mode | `auto`, `serial`, `satellite`, `per_satellite`, `flat` | Chooses the dominant RHS threading strategy for a simulation |
| Inner scheduler | `static`, `dynamic` | Chooses fixed worker assignment or dynamic work stealing |
| Adaptive policy | `0`, `1` | Allows measured runtime feedback to adjust worker counts and route choices |
| Persistent hints | `0`, `1` | Reuses measured policy choices across repeated runs when state persistence is enabled |

For a machine with many cores running a multi-satellite campaign, start with
`auto` outer routing, `auto` inner modes, and the default `auto` RHS execution
mode. Enable adaptive policy and persistent hints for repeated performance
studies after validating correctness against a forced serial run.

### Parallel profiles

A profile is a named bundle of the routing and policy settings above.
`SolverConfig(parallel=true)` runs under `R7`. Apply a profile explicitly with
`SpaceAGORA.ParallelProfiles.with_parallel_profile(f, "R6")` or with
`ParallelConfig(profile="R6")` in a `SimulationEngineConfig`; both expand
the whole bundle. Setting `SPACEAGORA_PARALLEL_PROFILE` on its own does
**not** apply a profile: only the two routes above expand it, and alone it
merely names the persisted policy state files (and, for `R0`, forces a serial
right-hand side).

| Profile | Outer route | Inner modes | What it adds |
|---|---|---|---|
| `R0` | serial | off | The true serial baseline; nothing is threaded, at any thread count |
| `R1_a` | threads | off | Outer parallelism only, on Julia threads |
| `R1_b` | process | off | Outer parallelism only, on worker processes |
| `R2` | serial | auto | Inner parallelism only |
| `R3` | auto | auto | Heuristic outer route with static inner widths |
| `R4` | auto, adaptive | auto | The outer-route bandit and pre-solve RHS plan calibration |
| `R5` | auto, adaptive | auto | R4 plus persistent hints and the measured-reward width chooser |
| `R6` | auto, adaptive | auto | R5 plus `SPACEAGORA_PARALLEL_POLICY_V2` |
| `R7` | auto, adaptive | auto | R6 plus `SPACEAGORA_CAMPAIGN_PLANNER=predictive` |

`R6` is the bandit-planner counterpart of `R7`. Everything it changes sits behind one switch,
`SPACEAGORA_PARALLEL_POLICY_V2`, so that with the switch off it is exactly
`R5`; the switch is what a paired comparison between the two measures.
Under the switch:

- A Monte Carlo campaign routes to whichever outer axis has more cores:
  process workers when the pool can match or exceed the thread count
  (`SPACEAGORA_PERF_PROCS` caps the pool), threads otherwise, on any machine
  size. The router does not spend campaigns trying the other route; it
  switches only on history it already holds.
- The process pool is sized by memory as well as cores. Each worker is priced
  at the coordinator's own resident set (never under 2 GB) plus 2 MB per
  spacecraft for native GRAM, against the machine's memory (or its cgroup
  limit) less a reserve and less what the coordinator already holds; the
  route is offered only when at least two workers fit, and the split ladder
  stops at the number that fit. `SPACEAGORA_MEMORY_BUDGET_GB`,
  `SPACEAGORA_PERF_WORKER_MEMORY_GB` and `SPACEAGORA_GRAM_SAT_MEMORY_MB`
  override the three terms. The per-spacecraft figure was measured on
  missions of at most 1800 s that did not return their solutions
  (`docs/architecture/gram_memory_footprint.md`); for longer missions, denser
  saving, or samples that return their solutions it may be too low, so
  measure a worker and set the overrides on a tight memory budget.
- The outer split width is learned by racing the candidate widths inside the
  first campaign of a new workload, with at least three samples per worker
  per width. A campaign too small to race takes the widest split. Widths are
  never explored one per campaign.
- A cached "the heuristic won" verdict from RHS plan calibration is honoured
  on long solves instead of re-running the sweep in every process.
- Callback and effector widths take the static `min(items, budget)`. The
  density callback additionally gets a short pre-solve width sweep, and a
  narrower width is pinned only if it beats the static one by the calibration
  margin.
- The coordinator collects garbage before dispatching a campaign, so a
  threaded campaign's heap does not stall the next campaign's process
  dispatch.

The [September 2026 session record](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/ab545b85c99cd795c61a29dc0ae928d8b6195770/docs/policy_v2_session_record_20260903.md)
reports historical paired measurements on a Ryzen 9 9900X workstation with
Julia 1.12.1. It records improvements for selected workloads and launch
configurations, as well as slower cases. Those results describe the recorded
source and measurement conditions; they do not establish that `R6` is fastest
for every workload or on the current revision. Compare routes for the target
workload and machine with `scripts/calibration/paired_profile_probe.jl` and
`scripts/calibration/paired_campaign_probe.jl --src-runner --profile=full`. Retain the source
revision, configuration, calibration state, cold/warm distinction and run records
with any reported timing; use paired measurements and their observed variability
rather than assuming a fixed benchmark resolution.

`R7` is `R6` with the campaign route chosen by a predictive planner instead of
the bandit. The planner scores the candidate routes for the campaign in front
of it -- its sample count, spacecraft count, force model and mission length --
against per-machine calibration constants, and dispatches the route it prices
as cheapest. The initial plan is chosen before the first sample. A mixed plan can then
be narrowed by the dispatch guard: it closes consumers before their next job,
without interrupting an active sample or starting a second dispatch. The
[planner and guard contract](../architecture/predictive_routing_r7.md#the-guard)
explains those conditions and the default-disabled route switch. Everything else is exactly `R6`, and the single switch
`SPACEAGORA_CAMPAIGN_PLANNER` (`bandit`, the default, or `predictive`) is the
only difference between them, so a paired run of the two measures the planner
and nothing else.

The calibration constants are created automatically on the first
`parallel=true` run, or explicitly with

```bash
julia --project=. --threads=<T> scripts/calibration/calibrate_machine.jl
```

which is run once per machine, at the thread count you intend to simulate at,
and writes `output/parallel_policy_state/cost_constants_<fingerprint>.toml`.
Every later run on that machine reads the file instead of measuring anything
itself. `R7` applied without the flag (through a profile) does not calibrate;
it then prices the routes with no model of contention, so run the script
before comparing `R7` against `R6` on a new machine. Historical R6/R7 comparisons
are recorded in the [predictive-routing study](https://github.com/Space-FALCON-Lab/SpaceAGORA.jl/blob/ab545b85c99cd795c61a29dc0ae928d8b6195770/docs/architecture/predictive_routing_r7.md).
Their source revisions, machines, workloads and calibration state define their
scope; they do not establish a relative performance guarantee for this revision.

### Environment variables

For CLI runs or scripted batch execution, set the controls before launching.
The examples run in the examples environment; prepare it once with
`julia --startup-file=no examples/setup.jl`.

```powershell
$env:SPACEAGORA_PERF_PARALLEL_BACKEND = "auto"
$env:SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE = "1"
$env:SPACEAGORA_DENSITY_CALLBACK_PARALLEL = "auto"
$env:SPACEAGORA_CONTROL_CALLBACK_PARALLEL = "auto"
$env:SPACEAGORA_THERMAL_CALLBACK_PARALLEL = "auto"
$env:SPACEAGORA_MULTIBODY_PARALLEL = "auto"
$env:SPACEAGORA_EFFECTOR_PARALLEL = "auto"
$env:SPACEAGORA_RHS_BATCH_PARALLEL = "auto"
$env:SPACEAGORA_RHS_EXECUTION_MODE = "auto"
$env:SPACEAGORA_RHS_CALIBRATE = "auto"
$env:SPACEAGORA_PARALLEL_POLICY_ADAPTIVE = "1"
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

`cmd.exe`:

```bat
set SPACEAGORA_PERF_PARALLEL_BACKEND=auto
set SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE=1
set SPACEAGORA_DENSITY_CALLBACK_PARALLEL=auto
set SPACEAGORA_CONTROL_CALLBACK_PARALLEL=auto
set SPACEAGORA_THERMAL_CALLBACK_PARALLEL=auto
set SPACEAGORA_MULTIBODY_PARALLEL=auto
set SPACEAGORA_EFFECTOR_PARALLEL=auto
set SPACEAGORA_RHS_BATCH_PARALLEL=auto
set SPACEAGORA_RHS_EXECUTION_MODE=auto
set SPACEAGORA_RHS_CALIBRATE=auto
set SPACEAGORA_PARALLEL_POLICY_ADAPTIVE=1
julia --project=examples examples/AGORA_Basic_Quickstart.jl
```

POSIX shells:

```bash
export SPACEAGORA_PERF_PARALLEL_BACKEND=auto
export SPACEAGORA_PERF_OUTER_ROUTE_ADAPTIVE=1
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

To force a serial baseline for validation:

```bash
export SPACEAGORA_PERF_PARALLEL_BACKEND=none
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

Automatic mode evaluates the effector loop serially when at most one satellite
is active, avoiding per-step task overhead. `SPACEAGORA_EFFECTOR_PARALLEL=on`
can force threading when the remaining eligibility checks allow it. Threaded
and serial evaluation give bit-identical results either way.

That guarantee holds within an RHS execution route. The flat constellation
route (`SPACEAGORA_RHS_EXECUTION_MODE=flat`, also one of the candidates the
start-up calibration sweep may select for two or more satellites) evaluates the
batched gravity terms through SIMD kernels and sums contributions in a
different order, so its trajectories agree with the serial and per-satellite
routes only to the last bit of the force, not bit for bit. For runs that must
be reproducible across machines and thread counts, pin the route and skip the
sweep:

```bash
export SPACEAGORA_RHS_CALIBRATE=off
export SPACEAGORA_RHS_EXECUTION_MODE=per_satellite   # or serial
```

The same environment variables can be scoped in Julia with `withenv` when you
want one process to run several scenarios with different settings.

Runtime parallelism and cache knobs are resolved once at `run_simulation`
start into a typed, run-scoped snapshot (together with any active
`SimulationEngineConfig` overrides), so the RHS and callback hot paths never
touch process-global `ENV` during integration. Wrap the whole
`run_simulation` call in `withenv` — changing a variable while a run is in
flight does not affect that run.

When nothing reads the trajectory (`return_solution=false`,
`results=false`, no solver metadata requested), the solver skips per-step
solution storage entirely (`save_on=false`, endpoints kept). Skipping this storage can reduce memory use. The effect on a campaign depends
on its workload and output requirements; it does not establish a scaling
guarantee. Explicitly
set `SPACEAGORA_SOLVER_SAVE_EVERYSTEP` / `SPACEAGORA_SOLVER_SAVE_ON`
values override this default in either direction. With storage on, a run
holds one copy of the state (plus the interpolant stages) per accepted step
and the before/after saves of continuous events; the per-step housekeeping
callbacks (planet frame, density and thermal samples, quaternion projection)
add no saves of their own.

### Adaptive campaign routing (`threads=:auto`)

`parallel=true` passes `threads=:auto` to the campaign runners under `R7`.
Passed directly, without the flag, `threads=:auto` uses whichever campaign
planner the environment selects (`SPACEAGORA_CAMPAIGN_PLANNER`, default the
outer-route bandit `SpaceAGORA.ParallelProfiles.select_outer_route!`, which
keeps empirical runtime statistics per workload signature):

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
campaign shape (sample count, per-sample satellite count, density-model
family, mission length); the `SimulationConfiguration` method derives those
fields for you. After every campaign the runner records per-sample success and
amortized wall-clock feedback via
`SpaceAGORA.ParallelProfiles.record_outer_route_feedback!`, so repeated
campaigns with the same shape can adapt their choices from observed outcomes.
The result depends on workload, host load and available history; it is not a
guarantee of convergence to the fastest route. History accumulates in the process-global
`SpaceAGORA.SimulationCampaigns.campaign_outer_route_state()`; inspect it with
`SpaceAGORA.ParallelProfiles.outer_route_stats_snapshot`, reset it with
`SpaceAGORA.ParallelProfiles.reset_outer_route_state!`, or pass an isolated
`SpaceAGORA.ParallelProfiles.OuterRouteState()` via `route_state` (useful for
tests and one-off studies). These names are internal: they are not exported
and may change without notice.

While the adaptive route runs threaded workers, the runner sets
`SPACEAGORA_OUTER_PARALLEL_ACTIVE=1` and — unless you exported one yourself —
an `SPACEAGORA_INNER_THREAD_BUDGET` of `nthreads() ÷ workers`, so per-sample
inner threading and the outer campaign split the thread pool instead of
oversubscribing it. Nested adaptive campaigns (an `:auto` campaign running
inside another campaign's worker, where `SPACEAGORA_OUTER_PARALLEL_ACTIVE` is
already set) yield to the enclosing split: they execute serially and record no
feedback, so contended timings never poison the shared route statistics.

The spherical-harmonics gravity SIMD batch route is a separate hot path that
fires on every RHS/ODE step rather than once per sample, so historically it
did not check `SPACEAGORA_OUTER_PARALLEL_ACTIVE` before spawning its own
nested worker batch — under `:threads` outer parallelism this produced severe
nested contention (anywhere from ~1x overhead to multi-minute hangs on
otherwise-identical runs). It now defaults to running serially whenever outer
parallelism is active. Only set `SPACEAGORA_HARMONICS_BATCH_ALLOW_WITH_OUTER=1`
if you have measured that splitting the thread budget between outer workers
and the harmonics batch is actually faster for your workload; the default is
the safe choice.

### Persistent Performance Hints

Persistent performance hints let the adaptive inner scheduler record wall-clock
measurements and reuse them across repeated calls. This reduces calibration
overhead on the second and subsequent runs.

Set a stable hint state path to make hints persist across Julia sessions:

```powershell
$env:SPACEAGORA_PARALLEL_POLICY_STATE_PATH = "C:\\path\\to\\hints.toml"
$env:SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS = "1"
$env:SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST = "1"
$env:SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD = "1"
```

`cmd.exe`:

```bat
set SPACEAGORA_PARALLEL_POLICY_STATE_PATH=C:\path\to\hints.toml
set SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS=1
set SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST=1
set SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD=1
```

POSIX shells:

```bash
export SPACEAGORA_PARALLEL_POLICY_STATE_PATH=/path/to/hints.toml
export SPACEAGORA_PARALLEL_POLICY_PERSISTENT_HINTS=1
export SPACEAGORA_PARALLEL_POLICY_STATE_PERSIST=1
export SPACEAGORA_PARALLEL_POLICY_MEASURED_REWARD=1
```

Use `SPACEAGORA_PARALLEL_POLICY_STATE_RESET=1` to discard accumulated hints
and restart calibration.

### Auditing Active Controls

For reproducible studies, record the `SPACEAGORA_*` controls that define the
route and inner-policy behavior:

```julia
for key in sort(collect(keys(ENV)))
    if startswith(key, "SPACEAGORA_")
        println("$key = $(ENV[key])")
    end
end
```

### Optional direct derivative assembly

For a performance experiment on a translational constellation in vacuum, you
can enable direct writes into the derivative buffer:

```julia
withenv("SPACEAGORA_RHS_FINAL_ASSEMBLY_DIRECT_LAYOUT" => "1") do
    run_simulation(args)
end
```

This option is off by default. It only changes the final assembly step when
the runtime already uses the flat constellation-effector route; it does not
select that route itself. It supports the standard position, velocity, mass
and heat-load layout with contiguous `Float64` storage, a `NoAtmosphereModel`,
and no attitude propagation, control effectors or robot arm. Unsupported
configurations and custom layouts continue through the existing assembly path.
The option does not change the force calculation or its accumulation order.

Keep it off for your baseline, then compare results, elapsed time and allocations
with it enabled on the same warmed-up workload. No whole-simulation speedup is
established by the correctness tests. Set the value to `"0"` or remove the
variable to restore the default. As with other runtime controls, set it before
starting the simulation.
