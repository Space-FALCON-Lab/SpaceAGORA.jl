# Wall-Clock Runtime Comparison: 50 and 100 Orbits

For the subsequent comparison including `8_SpaceAGORA.jl-main_v5`, see
[Runtime Comparison Including v5](RUNTIME_COMPARISON_WITH_V5.md).
That snapshot lacks inclination scheduling, so the follow-up reruns all
four codes using SMA scheduling rather than mixing policies in this table.

## Setup

The earlier quoted times (27.92 s prototype, 3.30 s v1, 11.27 s v9)
already corresponded to **50 orbits**, not a shorter simulation.
This experiment repeats 50 orbits and adds **100 orbits**, while separating
first-call effects from warmed-up execution.

Physical settings are the same as the
[inclination comparison](INCLINATION_COMPARISON_50_HELPERS_50_ORBITS.md):

- 50 evenly phased circular equatorial helpers at 1000 km altitude.
- One initially circular target at 1050 km altitude and 1 degree inclination.
- Initial target anomaly, RAAN, and argument of periapsis: zero.
- Inclination-raising scheduler, `:gve_inc`.
- 227 kg spacecraft, 10 kW input power, amplification 100, and 200 km range.
- Central gravity only, without J2 or drag.
- 50 initial target periods: 318,567.018025 simulation seconds (88.49 hours).
- 100 initial target periods: 637,134.036050 simulation seconds (176.98 hours).

The prototype retains its Vern9 solver and RHS-level scheduling. v1 retains
its lightweight central-gravity wrapper and actual scheduler source; v9
retains native `run_simulation`. v1 and v9 use Tsit5 with a 10 s maximum step.
The numerical tolerances and output behavior were not standardized beyond
the settings used in the existing comparison.

## Benchmark method

For each implementation and duration:

1. Launch a fresh Julia process and load its existing runner.
2. Execute one first simulation call.
3. Execute three more calls in the same process, each constructing fresh
   initial states and models.
4. Report the median and range of the three warmed-up simulation-call times.

All six groups were run **sequentially**, not concurrently, to avoid
competition between these benchmark processes. Each used one Julia thread
and one OpenBLAS thread. A full garbage collection was requested before each
trial, outside its timer; garbage collection occurring during a trial remains
included.

Environment: Julia 1.12.1, Linux, AMD Ryzen 7 9700X 8-Core Processor.
The machine is shared, so unrelated system load was not controlled.
Measurements were taken on 2026-10-18 using the current workspace sources.

## Results: simulation-call wall time

All numbers in this section are **wall-clock seconds**, not simulation time.

| Code | Orbits | First call | Warm median | Warm minimum | Warm maximum |
|---|---:|---:|---:|---:|---:|
| Prototype | 50 | 27.158 | **23.088** | 23.072 | 23.200 |
| v1 | 50 | 3.068 | **0.923** | 0.895 | 0.928 |
| v9 | 50 | 10.926 | **1.633** | 1.625 | 1.654 |
| Prototype | 100 | 49.747 | **45.506** | 45.386 | 45.517 |
| v1 | 100 | 4.011 | **1.719** | 1.706 | 1.904 |
| v9 | 100 | 12.129 | **2.894** | 2.469 | 3.290 |

### Scaling when duration doubles

| Code | 100-orbit / 50-orbit warm median |
|---|---:|
| Prototype | 1.971x |
| v1 | 1.863x |
| v9 | 1.772x |

Runtime grows roughly with simulated duration in these two measurements.
This is not proof of a general scaling law: only two durations were tested,
and the v9 100-orbit samples show noticeable variability.

### Interpretation

- **v1 is fastest in this particular runner comparison.**
- **v9 is much closer to v1 after warm-up** than the earlier single-call
  measurements suggested.
- First-call overhead, including compilation and initialization within the
  timed call, is substantial for v1 and especially v9. The difference
  between first and warm timings is not a measurement of compilation alone.
- The prototype remains substantially slower after warm-up. It uses a
  different solver/tolerances and performs scheduling during RHS evaluation.
  A profile would be needed to attribute its runtime to individual components.
- This does **not** establish that v1's scheduling algorithm alone is faster
  than v9's matching algorithm. The measured scopes include different
  simulation frameworks, state layouts, and diagnostics.
- No runtime optimization or change to scheduler mathematics was made here.

## What the timers include

The primary metric reuses each runner's existing timer:

| Code | Timed call |
|---|---|
| Prototype | Original `run_open_cavity_multi` |
| v1 | ODE `solve` |
| v9 | Native `run_simulation` |

These are useful practical simulation-call comparisons, but not identical
low-level scopes. For example, the prototype's original runner can perform
internal setup/diagnostics that are outside v1's `solve` call.

The primary timer excludes loading the runner beforehand, the subsequent
comparison-runner CSV output and analysis, and Julia process startup.
It includes any compilation or framework work inside the timed call.
First call means **first call in that process**, not a cold package-cache
or clean-machine measurement.

An additional outer timer measures the complete `run_case` wrapper, including
configuration, simulation, post-processing, CSV writing, and result printing:

| Code | 50-orbit warm wrapper median (s) | 100-orbit warm wrapper median (s) |
|---|---:|---:|
| Prototype | 23.104 | 45.543 |
| v1 | 0.935 | 1.743 |
| v9 | 1.648 | 2.921 |

Runner-loading times are saved separately in the raw summaries. They do not
include the whole Julia process startup or the benchmark harness's own imports.
Neither timer should be described as full command-start-to-process-exit time.

## Verification

**All 24 simulation calls returned `Success`** and reached their requested
duration. The harness checked the helper count, policy, initial inclination,
positive finite runtime, and consistency of inner/outer timers.

Within each implementation/duration group, all four calls produced exactly
the same final SMA, inclination, eccentricity, activation count, active-update
count, and final helper. Repeated calls therefore did not reuse an already
evolved simulation state.

| Code | Final inclination after 50 orbits (deg) | Final inclination after 100 orbits (deg) |
|---|---:|---:|
| Prototype | 1.002199742725 | 1.004336908250 |
| v1 | 1.002199107272 | 1.004335950649 |
| v9 | 1.002202016132 | 1.004342054134 |

The 50-orbit results reproduce the earlier inclination comparison.
Differences between implementations remain numerical/behavioral comparisons,
not certified error bounds. No timestep-convergence study was added here.

## Reproduce and inspect

Run from the workspace root, sequentially:

```sh
export GKSwstype=100 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for orbits in 50 100; do
  julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/benchmark_scheduler_runtime.jl prototype 50 "$orbits" gve_inc 1 || exit
  julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/benchmark_scheduler_runtime.jl 4_v1 50 "$orbits" gve_inc 1 || exit
  julia --project="9_SpaceAGORA.jl-main_full integration" SciTech_Paper/benchmark_scheduler_runtime.jl 9_full 50 "$orbits" gve_inc 1 || exit
done
```

The [benchmark harness](benchmark_scheduler_runtime.jl) loads the existing
runners without executing them automatically, then calls `run_case` four
times with an isolated output root. Running the original runner files
directly still performs one simulation as before.

| Duration | Prototype raw timings | v1 raw timings | v9 raw timings |
|---|---|---|---|
| 50 orbits | [prototype.toml](runtime_benchmarks/50_helpers_50_orbits_gve_inc_target_inc_1.0_deg/prototype.toml) | [4_v1.toml](runtime_benchmarks/50_helpers_50_orbits_gve_inc_target_inc_1.0_deg/4_v1.toml) | [9_full.toml](runtime_benchmarks/50_helpers_50_orbits_gve_inc_target_inc_1.0_deg/9_full.toml) |
| 100 orbits | [prototype.toml](runtime_benchmarks/50_helpers_100_orbits_gve_inc_target_inc_1.0_deg/prototype.toml) | [4_v1.toml](runtime_benchmarks/50_helpers_100_orbits_gve_inc_target_inc_1.0_deg/4_v1.toml) | [9_full.toml](runtime_benchmarks/50_helpers_100_orbits_gve_inc_target_inc_1.0_deg/9_full.toml) |

Every trial's summary and schedule are saved under
`runtime_benchmarks/trial_1/` through `runtime_benchmarks/trial_4/`, within
their case-specific `results/` subdirectories. Trial 1 is the first-call
measurement; trials 2-4 are warm measurements. Earlier comparison results
were not overwritten.
