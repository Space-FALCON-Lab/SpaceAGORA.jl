# Runtime Comparison Including v5: 50 Helpers, 50 and 100 Orbits

An additional **100-orbit, first-call-only** measurement for
`2_SpaceAGORA.jl` is included in the final section. The original four-code
warm benchmark below remains unchanged.

For a separate test of compilation reuse when changing constellation size,
see [v9: 10, 20, and 30 Helpers in One Process](V9_HELPER_COUNT_WARMUP_SEQUENCE.md).
Each helper count was run exactly once.

## Important distinction from the previous benchmark

The current `8_SpaceAGORA.jl-main_v5` snapshot does **not** implement an
inclination scheduler. Its [policy constructor and dispatcher](../8_SpaceAGORA.jl-main_v5/I_src/1.4_satellite_interlink/2_dynamic_programming_models.jl)
support `:gve_sma`, `:gve_eccentricity`, and `:gve_sma_first_index`, not `:gve_inc`.

Therefore, this is a **new four-code SMA benchmark**, not an extra row added
to the previous [inclination runtime comparison](RUNTIME_COMPARISON_50_100_ORBITS.md).
All four codes were rerun using `:gve_sma`. The physical geometry is unchanged:
50 helpers at 1000 km/0 degrees and one target at 1050 km/1 degree.
No inclination scheduler was added to v5, and no production scheduler or
simulation implementation was modified for this comparison.

## Experiment and measurement method

- Circular initial orbits; helpers evenly phased by `360 * (slot - 1) / 50`.
- Target initial anomaly, RAAN, and argument of periapsis: zero.
- 227 kg per spacecraft; laser input power 10 kW, amplification 100, range 200 km.
- Central gravity only, without J2 or drag.
- Earth radius 6,378,136.6 m; propagation gravitational parameter
  3.98600436233e14 m^3/s^2.
- Target-helper links only, with one terminal per spacecraft. At most one
  helper can serve the target at a time; no helper-helper links.
- Objective: maximize the target's instantaneous SMA increase.
- Initial target period: 6,371.340360499 s.
- 50 periods: 318,567.018024937 simulation seconds.
- 100 periods: 637,134.036049874 simulation seconds.

The prototype retains its original Vern9 solver, 1e-12 relative/absolute
tolerances, and scheduling inside the RHS. v1 retains its actual legacy
scheduler in the existing central-gravity wrapper. v5 and v9 use their
respective native `run_simulation` implementations. v1, v5, and v9 use
Tsit5, relative tolerance 1e-9, absolute tolerance 1e-11, and a 10 s maximum step.
These are practical runner comparisons, not identical numerical workloads.

For each code and duration, a fresh Julia process performs one first call
followed by three warmed calls. Each trial builds fresh initial states and
models. The warm result is the median of calls 2-4. A full garbage collection
is requested before each call, outside its timer. All benchmark groups ran
sequentially, with one Julia thread and one OpenBLAS thread.

Environment: Julia 1.12.1, Linux, AMD Ryzen 7 9700X 8-Core Processor.
Measurements were taken on 2026-10-18. This is a shared machine;
unrelated system load was not controlled.

## Simulation-call wall time

All times below are **wall-clock seconds**, not simulated seconds.

| Code | 50 orbits: first call | 50 orbits: warm median | 100 orbits: first call | 100 orbits: warm median |
|---|---:|---:|---:|---:|
| Prototype | 19.361 | 15.315 | 35.350 | 31.320 |
| v1 | 3.048 | **0.927** | 3.996 | **1.693** |
| v5 | 32.695 | **4.424** | 38.071 | **9.521** |
| v9 | 10.875 | **1.628** | 12.162 | **3.297** |

| Code | 50-orbit warm range (s) | 100-orbit warm range (s) | Warm 100 / 50 ratio |
|---|---:|---:|---:|
| Prototype | 15.190-15.405 | 31.295-31.484 | 2.045x |
| v1 | 0.905-0.928 | 1.648-1.702 | 1.827x |
| v5 | 4.418-4.850 | 9.517-9.570 | 2.152x |
| v9 | 1.625-1.648 | 2.932-3.344 | 2.025x |

### Main conclusions

- Warm runtime ordering in this experiment is **v1 < v9 < v5 < prototype**.
- v5 takes **2.72x** v9's warm time at 50 orbits and **2.89x** at 100 orbits.
- v5 takes **4.77x** v1's warm time at 50 orbits and **5.62x** at 100 orbits.
- The prototype takes **3.46x** v5's warm time at 50 orbits and **3.29x**
  at 100 orbits.
- v5's first call is slower than the prototype's first call, but v5 is
  substantially faster after warm-up. First-call overhead is significant.
- Doubling the duration roughly doubles runtime here. Two durations and
  three warm samples per group do not establish a general scaling law.

The earlier prototype inclination benchmark took about 23.09 s warm at
50 orbits; this SMA benchmark takes about 15.31 s. These use different
objectives and trajectories, so that difference is **not** evidence of a
code optimization.

## What these timers do and do not measure

| Code | Primary timed call |
|---|---|
| Prototype | Original `run_open_cavity_multi` |
| v1 | ODE `solve` |
| v5 | Native `run_simulation` |
| v9 | Native `run_simulation` |

The primary timer excludes prior package/runner loading, Julia startup,
and the comparison runner's subsequent CSV writing. It includes framework
work and compilation occurring inside the timed call. First call means the
first simulation call in that process, not a clean-machine or cold-cache run.
First-call minus warm time is not a measurement of compilation alone.

An outer timer also measures configuration, simulation, post-processing,
CSV writing, and result printing:

| Code | 50-orbit warm wrapper median (s) | 100-orbit warm wrapper median (s) |
|---|---:|---:|
| Prototype | 15.326 | 31.342 |
| v1 | 0.938 | 1.716 |
| v5 | 4.437 | 9.547 |
| v9 | 1.643 | 3.324 |

v5's runner-loading time was 21.872 s in the 50-orbit process and 7.020 s
in the 100-orbit process. Its first package load automatically provisioned
the declared Conda environment and precompiled the package. That work is
outside the simulation-call timer. These measurements must not be described
as command-start-to-process-exit times.

### Why the v5/v9 gap is not a scheduler-only speed comparison

The current v5 and v9 both use physical force/mass SMA scoring and
maximum-weight terminal matching. Both have a tree dynamic-programming
path and a NetworkX Blossom path for cyclic candidate components.
This target-helper star is acyclic, so it exercises the tree path, not
Blossom. It is not a dense all-to-all constellation benchmark.

v5 and v9 share the same comparison configuration through the
[native runner](native_interlink_comparison.jl), but load different native
simulation packages. Their measured calls still include the whole simulation
framework, RHS, callbacks, diagnostics, and allocations. The observed v5/v9
gap cannot be attributed to one component without profiling; identical
output does not imply identical execution cost.

## Physical results and verification

All **32 benchmark calls succeeded** and reached the requested duration.
The saved summaries and schedule CSVs were checked for consistent trial
metrics, valid columns/helper indices, monotonic times, correct sample
counts, and exact final times. Within each code/duration group, all four
saved schedules are byte-identical and the final orbital metrics agree exactly.

| Code | SMA increase after 50 orbits (m) | SMA increase after 100 orbits (m) |
|---|---:|---:|
| Prototype | 1010.765077 | 2042.506581 |
| v1 | 1008.015803 | 2036.733413 |
| v5 | 1007.532924 | 2035.498480 |
| v9 | 1007.532924 | 2035.498480 |

**v5 and v9 produce exactly equal saved orbital/counter metrics and
byte-identical schedule CSVs at both durations**, not merely equal rounded
values. This is equality for these cases, not proof that the two packages
are equivalent for every configuration.

The four-code SMA-gain spread is 3.232153 m at 50 orbits and 7.008101 m at
100 orbits, respectively 0.3205% and 0.3439% of the mean gain.

| Code | Final inclination after 50 orbits (deg) | Final inclination after 100 orbits (deg) |
|---|---:|---:|
| Prototype | 1.000066178648 | 0.999938574455 |
| v1 | 1.000066512530 | 0.999942330308 |
| v5 | 1.000069707517 | 0.999946218719 |
| v9 | 1.000069707517 | 0.999946218719 |

Inclination is not the optimization objective here, so it need not increase.
The prototype/v1 differences remain subject to their different solver,
scoring, and scheduling behavior. No timestep-convergence study was added.
The prototype CSV records rescoring at saved states, not every internal
RHS scheduling decision.

The native runner extraction was also regression-tested against the
previous v9 50-orbit **inclination** result: all checked orbital metrics
and scheduler counters reproduced exactly. Existing inclination results
were not overwritten.

## Code, data, and reproduction

- [v5 entry point](run_8_v5.jl) and [v9 entry point](run_9_full.jl)
  load the [shared native configuration](native_interlink_comparison.jl).
  The shared runner verifies the loaded package entry path to prevent
  accidentally benchmarking the wrong SpaceAGORA snapshot.
- [Benchmark harness](benchmark_scheduler_runtime.jl) supports all four codes.
- [50-orbit raw timing summaries](runtime_benchmarks/50_helpers_50_orbits_gve_sma_target_inc_1.0_deg/)
  and [100-orbit raw timing summaries](runtime_benchmarks/50_helpers_100_orbits_gve_sma_target_inc_1.0_deg/)
  each contain `prototype.toml`, `4_v1.toml`, `8_v5.toml`, and `9_full.toml`.
- Per-trial summaries and schedules are under
  `runtime_benchmarks/trial_1/` through `trial_4/`, then
  `results/50_helpers_<orbits>_orbits_gve_sma_target_inc_1.0_deg/<code>/`.

From the workspace root, run sequentially:

```sh
export GKSwstype=100 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
for orbits in 50 100; do
  julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/benchmark_scheduler_runtime.jl prototype 50 "$orbits" gve_sma 1 || exit
  julia --project="4_SpaceAGORA.jl-main_v1" SciTech_Paper/benchmark_scheduler_runtime.jl 4_v1 50 "$orbits" gve_sma 1 || exit
  julia --project="8_SpaceAGORA.jl-main_v5" SciTech_Paper/benchmark_scheduler_runtime.jl 8_v5 50 "$orbits" gve_sma 1 || exit
  julia --project="9_SpaceAGORA.jl-main_full integration" SciTech_Paper/benchmark_scheduler_runtime.jl 9_full 50 "$orbits" gve_sma 1 || exit
done
```

The recorded run executed v5's 50-orbit group first, then the remaining
groups in the order above. Reproduction overwrites these case-specific
benchmark files, but not the earlier inclination benchmark.

## Additional snapshot: `2_SpaceAGORA.jl`, 100 orbits, first call only

At the user's request, this snapshot was run **once in a fresh Julia
process**, without a warm-up simulation or subsequent repeat. The physical
case is the same SMA case above: 50 equatorial helpers at 1000 km and one
initially circular target at 1050 km/1 degree, propagated for 100 initial
target periods.

### First simulation-call comparison

| Code | 100-orbit first simulation call (wall seconds) |
|---|---:|
| Prototype | 35.350 |
| v1 | 3.996 |
| `2_SpaceAGORA.jl` | **29.096** |
| v5 | 38.071 |
| v9 | 12.162 |

The added snapshot took **2.39x v9's first-call time**, **7.28x v1's**,
about 18% less time than the prototype, and about 24% less than v5.
These compare first calls only; the other four entries reuse the first-call
measurements already recorded above. They were not rerun for this addition.
One sample provides no variability estimate and does not establish warm
performance.

For the added snapshot:

- Native `run_simulation` wall time: **29.095600568 s**.
- Outer configuration/simulation/export wrapper: **29.753448391 s**.
- Prior common-runner/package load: **46.284683945 s**, including package
  precompilation. This is excluded from the simulation-call time.
- No command-start-to-process-exit timer was collected.

### Preserved native behavior and matched settings

The [new first-call runner](run_2_first_call.jl) loads the actual package
from `2_SpaceAGORA.jl` and calls its native simulation engine. It does not
substitute v1's lightweight ODE wrapper or v9's scheduler.

This snapshot's [laser implementation](../2_SpaceAGORA.jl/src/dynamics/coupled/force_torque_models/laser_link_effectors.jl)
uses direction-based GVE scoring and applies **discrete velocity kicks**.
The runner preserves the native example's callback order:

1. Apply the laser impulse for the elapsed step using the previous selection.
2. Update the scheduler for the next step.
3. Record the selected helper without changing the state or scheduler.

The original native [ORACLE example](../2_SpaceAGORA.jl/examples/ORACLE_LaserLinks/oracle_runner.jl)
uses J2 and 1e-12 tolerances. For this comparison only, its configuration
was replaced with central gravity, Tsit5, relative tolerance 1e-9, absolute
tolerance 1e-11, and a 10 s maximum step, matching the requested benchmark.
The laser parameters and equal 227 kg spacecraft masses are also matched.
No production source was modified.

`isolate_state=false` preserves the native example's shared laser-model
identity between the simulation and callbacks. The runner checks that
identity after the solve. State isolation, callbacks, state layout, and
framework costs still differ from other implementations, so these timings
are not scheduler-only comparisons.

### Result and validation

| Metric | `2_SpaceAGORA.jl` |
|---|---:|
| Return code | Success |
| Final simulation time (s) | 637134.036049874 |
| SMA increase (m) | 2029.830701171 |
| Final inclination (deg) | 0.999944662600 |
| Scheduler records, including initialization | 63723 |
| Native activation count, including initialization | 53 |
| Active scheduler records | 11674 |
| Final helper, raw index | 50 |

Raw helper index 50 means physical helper slot 49, because the target is
spacecraft 1. The initial recorded selection is raw helper index 2.
Activation counts include initialization and should not be interpreted
as 53 sustained encounter windows without interval analysis.

For reference, the other 100-orbit SMA increases were 2042.506581 m
(prototype), 2036.733413 m (v1), and 2035.498480 m (both v5 and v9).
The new snapshot's SMA gain is about 5.668 m below v9's. Its impulse-based
propagation differs from v9's continuous-force model; this experiment
does not isolate the contribution of that difference or establish an
accuracy ranking.

Validation checked `Success`, exact requested final time, initial
inclination, finite final metrics, callback/model identity, impulse
completion, all 63,723 schedule rows, time ordering, helper indices,
normalized orbit times, and agreement between recorded and native
activation/active-update counters. No repeatability claim is made for
this single-call addition.

### Saved data and reproduction

- [First-call timing and metadata](runtime_benchmarks/50_helpers_100_orbits_gve_sma_target_inc_1.0_deg/2_spaceagora.toml)
- [Physical summary](runtime_benchmarks/first_call/results/50_helpers_100_orbits_gve_sma_target_inc_1.0_deg/2_spaceagora/summary.toml)
- [Recorded schedule](runtime_benchmarks/first_call/results/50_helpers_100_orbits_gve_sma_target_inc_1.0_deg/2_spaceagora/schedule.csv)

From the workspace root, launch a fresh process:

```sh
GKSwstype=100 JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  julia --project="2_SpaceAGORA.jl" SciTech_Paper/run_2_first_call.jl
```

This executes exactly one fixed-case simulation and overwrites only its
own first-call result files.
