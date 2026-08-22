# Paper Parallelization Benchmarks — Test Cases

This study (`benchmarks/studies/paper_parallelization_benchmarks/`) runs the
benchmark suite backing the parallelization paper. It wraps the lower-level
`parallelization_performance` study (`main.jl`, `cli.jl`, `cases.jl`,
`modes.jl` in `benchmarks/studies/parallelization_performance/`), grouping its
cases and routing profiles into phases (B1–B8). Each phase isolates a
specific axis of SpaceAGORA's parallelization stack so the paper can make a
narrow, defensible claim per figure/table.

Two measurement properties are load-bearing for every phase below:

- **Warm-up is mandatory.** Each `(case, mode, threads, repeat)` combination
  runs in a fresh Julia subprocess, so without a warm-up solve the timed run is
  the process's *first* run and its wall time is dominated by JIT compilation of
  the RHS/solver stack, not by the simulation (measured: 21.4 s first solve vs
  0.017 s steady-state solve for `gravity_16sat_l20_vacuum`). Phase `warmup`
  counts are forwarded to the worker via `--warmup`; process-pool provisioning
  and per-worker warm-up likewise happen before the clock starts.
- **Thread ladders are in physical cores.** The RHS hot path is SIMD/FP-bound,
  gains nothing from SMT, and regresses sharply when SMT siblings are
  oversubscribed, so a ladder topped out at `Sys.CPU_THREADS` puts its worst
  data point at the end of the curve. See `SPACEAGORA_PPC_PHYSICAL_CORES`.

Pre-solve RHS-plan auto-calibration (`SPACEAGORA_RHS_CALIBRATE`) is forced off
for every mode in this harness: it pins a plan via `rhs_plan_override`, which
short-circuits the routing policy the profile ladder exists to compare.

## LaTeX table for the paper

```latex
\begin{table}[htbp]
\centering
\caption{Summary of the paper parallelization benchmark phases (B1--B8) and the parallelization axis each isolates.}
\label{tab:paper_benchmark_phases}
\begin{tabular}{lp{6.0cm}p{4.4cm}}
\toprule
Phase & Workload & Axis exercised \\
\midrule
B1 & Gravity-only constellation scaling, 1--1024 sat, inverse-square gravity & Outer-loop thread parallelism with minimal RHS cost (scaling baseline) \\
B2 & Constellation scaling with L20 spherical harmonics gravity, 1--1024 sat & Outer-loop parallelism combined with a harmonics-batched inner RHS \\
B3 & Atmosphere and GRAMTrackCache surrogate, 16 and 64 sat, drag plus density & Density-callback parallelism and surrogate caching under scaling \\
B4 & Monte Carlo throughput, single-spacecraft aerobraking and high-accuracy trials, 1--1024 samples & Process-level throughput scaling across 1--32 workers \\
B5 & Routing profile comparison on three representative workloads at a fixed thread budget & Full routing-profile ladder as the independent variable \\
B6 & Expanded router evaluation across spacecraft count, atmosphere/GRAM, actuators, propagation mode, duration/cadence & Router regret across workload axes \\
B7 & Heavy constellation scaling, 256--4096 sat, multi-hour missions, $\sim$8\,s serial baselines & Thread ladder on workloads large enough for scaling to be observable \\
B8 & Heavy Monte Carlo throughput, seconds-long aerobraking arcs, 1--64 process workers & Process-level throughput free of per-worker startup bias \\
\bottomrule
\end{tabular}
\end{table}
```

## Phases (B1–B8)

### B1 — Constellation Scaling, Gravity-Only Baseline
Cases: `single_inverse_square_vacuum` through `gravity_1024sat_inverse_square_vacuum`
(1, 4, 16, 64, 256, 1024 spacecraft, inverse-square gravity, no atmosphere).
Modes: `serial`, `outer_threads`, `full_smart`. Thread ladder: full
machine-scaled ladder.

Exercises **outer-loop (across-spacecraft) parallelism** in isolation, using
the cheapest possible right-hand side (a single `1/r²` gravity term) so the
per-satellite work is minimal and any speedup is attributable almost entirely
to outer-loop scheduling overhead/efficiency rather than to inner-loop or
callback parallelism. This is the scaling baseline every other phase is
compared against.

### B2 — Constellation Scaling, L20 Harmonics
Cases: `single_harmonics_l20_vacuum` through `gravity_1024sat_l20_vacuum`
(same constellation sizes as B1, but gravity is a degree/order-20 spherical
harmonics model instead of inverse-square). Modes: `serial`, `outer_threads`,
`outer_inner_adaptive`, `full_smart`.

Exercises **outer-loop parallelism combined with a heavier, harmonics-batched
inner right-hand side**. L20 harmonics evaluation is the first case where
inner-loop (within-satellite RHS) batching/parallelism has real work to do,
so this phase is what shows whether `outer_inner_adaptive`/`full_smart`
routing beats pure outer-thread parallelism once the RHS itself is
non-trivial.

### B3 — Constellation Scaling, Atmosphere and GRAMTrackCache Surrogate
Cases: `multi_16_aero_surrogate_cached` (16 sat, aerodynamic drag with a
cached/surrogate exponential density model) and `multi_64_high_fidelity` (64
sat, harmonics + SRP + aero + density). Modes: `serial`, `outer_threads`,
`full_smart`.

Exercises **atmosphere/density-callback parallelism and the
`GRAMTrackCache` surrogate path** under constellation scaling — i.e. whether
caching/parallelizing the density callback keeps pace once drag effectors and
a non-trivial density model are added on top of gravity.

### B4 — Monte Carlo Throughput Scaling
Cases: `montecarlo_mars_aerobraking`, `montecarlo_high_accuracy` (single
spacecraft, MC-jittered initial conditions). Modes: `serial`,
`outer_process`, `full_smart`. MC sample ladder: 1, 4, 16, 64, 256, 1024.
Thread mode: `:single` (1 thread per Julia process — all parallelism comes
from process count). Worker ladder: 1, 2, 4, 8, 16, 32 process workers.

Exercises **process-level (multi-process) throughput scaling for
independent Monte Carlo trials**, separate from thread-based parallelism.
Each worker count is a distinct controller run, producing a
throughput-vs-workers curve; this is the phase that validates SpaceAGORA's
ability to farm out embarrassingly-parallel MC ensembles across processes
(relevant to aerobraking/uncertainty-quantification workloads) rather than
across threads within one process.

### B5 — Routing Profile Comparison, R0–R5
Cases: `gravity_16sat_l20_vacuum`, `multi_16_aero_surrogate_cached`,
`articulated_1sat_fullstack` (three representative workloads: pure gravity,
gravity+aero, and a fully articulated single spacecraft with harmonics, aero,
thermal, and attitude). Modes: all six routing profiles — `serial`,
`outer_threads`, `inner_only`, `outer_inner_static`, `outer_inner_adaptive`,
`full_smart` (= R0–R5). Thread mode: `:max_only` (fixed at the machine's
maximum thread count for every mode, so all profiles are compared at the
same total thread budget).

Exercises **the full routing-profile ladder itself** — this is the phase
that isolates the *policy* (which combination of outer/inner/callback
parallelism is active, and whether routing is static or adaptive) as the
independent variable, holding thread budget and workload fixed. It is the
basis for the paper's "which profile should you actually use" comparison.

### B6 — Expanded Router Evaluation
Cases: the `many_sat_high_fidelity` spacecraft-count ladder (`multi_4`,
`multi_16`, `multi_64`, `multi_256`), `gravity_16sat_l20_vacuum`,
`montecarlo_mars_gram_live`, `montecarlo_multi_sat`,
`articulated_1sat_fullstack`, `multi_8sat_magnetorquer_attitude`, and the
duration/cadence pair `gravity_16sat_l20_vacuum_longmission` /
`_sparse_output`. Modes: `serial`, `outer_threads`, `outer_process`,
`outer_inner_adaptive`, `full_smart`. Thread mode: `:low_high` (1 and the
machine maximum). MC samples: 1 and 8.

Sweeps the workload axes B5 holds fixed — spacecraft count, atmosphere/GRAM
usage, force and actuator model count, interacting versus independent
propagation, thread and process budget, mission duration and output cadence —
so router regret can be reported across a matrix rather than three points.
`multi_16_gram_live` is defined in `cases.jl` but deliberately excluded: it
triggers an unbounded per-call memory leak in the vendored GRAMSuite native
binding (see the comment on the phase definition in `cli.jl`).

Caveat on the output-cadence axis: `num_steps_to_save` is currently inert as
far as the solver is concerned — nothing outside
`analysis/verification/telemetry_verification/` reads it, and the solve runs
with `save_everystep` and no `saveat`, so `gravity_16sat_l20_vacuum_sparse_output`
and `gravity_16sat_l20_vacuum` produce an identical number of saved steps.
That axis needs a real `saveat` plumbing change before it measures anything.

### B7 — Heavy Constellation Scaling
Cases: `heavy_1024sat_l50_6hr`, `heavy_4096sat_l50_1hr`,
`heavy_1024sat_fullstack_1hr`, `heavy_256sat_coupled6dof_2hr`. Modes:
`serial`, `outer_threads`, `outer_inner_adaptive`, `full_smart`. Thread ladder:
full machine-scaled ladder (physical cores).

Exercises **the thread ladder on workloads that are actually large enough to
scale**. Every case B1–B6 draws on resolves in a fraction of a second of real
integration once JIT compilation is excluded — measured post-warm-up on a
12-physical-core box: 0.017 s for `gravity_16sat_l20_vacuum`, 0.18 s for
`multi_64_high_fidelity`, 2.9 s for `multi_256_high_fidelity`. At that size the
wall time is fixed per-solve setup plus the integrator's serial spine, so the
ladder is flat by construction and no routing profile can distinguish itself.
B7's cases have ~8 s serial baselines with the parallelisable RHS as the
majority of that.

The four cases are a contrast rather than a sweep: `heavy_1024sat_l50_6hr` is a
single heavy effector with no callbacks (cleanest read on outer-loop scaling);
`heavy_4096sat_l50_1hr` is the same physics with 4× the satellites per RHS call
(reference measurement on 12 physical cores: 1024 satellites saturate at ~4
workers / 3.9×, 4096 keeps scaling past 12 / 5.1×, which localises where
per-RHS fork/join overhead stops being amortised); `heavy_1024sat_fullstack_1hr`
adds a heterogeneous effector mix plus density and thermal callbacks at scale;
`heavy_256sat_coupled6dof_2hr` is the only many-satellite case in the catalog
that propagates attitude.

### B8 — Heavy Monte Carlo Process Throughput
Case: `montecarlo_heavy_aerobraking`. Modes: `serial`, `outer_process`,
`full_smart`. MC sample ladder: 16, 64. Thread mode: `:single`. Worker
ladder: 1, 2, 4, 8, 16, 32, 64.

Same axis as B4, but each sample is a 12× longer aerobraking arc — ~1.05 s of
integration per trial (measured), against which pmap dispatch is roughly a 1%
term. B4's samples are short enough that its throughput-vs-workers curve was
largely per-worker process startup and JIT rather than the trials themselves.

The sample ladder stops at 64 because every worker-ladder entry is a separate
controller run that re-executes the whole mode × sample × repeat grid, and
`serial`/`full_smart` produce identical numbers at every worker count (both run
the batch in-process at `thread_mode = :single`), so their cost is paid seven
times over. 64 samples across 64 workers is already exactly one dispatch round.

## Underlying case families (`parallelization_performance/cases.jl`)

The phases above draw from a shared case catalog, grouped into families:

- **`gravity_only`** — `single_inverse_square_vacuum`,
  `gravity_{4,16,64,256,1024,2048}sat_inverse_square_vacuum`,
  `single_harmonics_l20_vacuum`, `gravity_{n}sat_l20_vacuum`. Pure orbital
  dynamics (inverse-square or L20 harmonics), no atmosphere, no attitude.
  Isolates outer-loop and gravity-RHS parallelism from every other subsystem.
- **`few_sat_high_fidelity`** — `single_harmonics_l50_vacuum` (1 sat, L50
  harmonics), `srp_heavy_high_area` (1 articulated sat, stacked SRP,
  attitude), `articulated_1sat_fullstack` (1 articulated sat: harmonics,
  aero, thermal, attitude). Stresses per-satellite (inner-loop/callback)
  cost with a single spacecraft, so outer-loop parallelism has nothing to
  exploit — these cases isolate inner-loop/effector/callback parallelism.
- **`many_sat_high_fidelity`** — `multi_16_aero_surrogate_cached`,
  `multi_{64,128}_high_fidelity`, `callback_128_aero_thermal`. Constellations
  with harmonics, SRP, aero, and density models together, exercising outer-
  and inner-loop parallelism simultaneously plus callback-heavy aero/thermal
  paths at scale.
- **`monte_carlo`** — `montecarlo_high_accuracy` (1 sat, jittered orbit,
  L50 harmonics), `montecarlo_multi_sat` (4 sat, jittered orbits, L20
  harmonics), `montecarlo_mars_aerobraking` (Mars aerobraking, jittered
  periapsis/apoapsis). Independent-trial workloads that isolate
  process-level throughput scaling from thread-level RHS parallelism.
- **`heavy_scaling`** — `heavy_1024sat_l50_6hr`, `heavy_4096sat_l50_1hr`,
  `heavy_1024sat_fullstack_1hr`, `heavy_256sat_coupled6dof_2hr`,
  `montecarlo_heavy_aerobraking`. Sized so the serial baseline is several
  seconds of integration rather than a fraction of one, which is the
  precondition for any of the other families' routing questions to have a
  measurable answer. See B7/B8 above.

## Routing profiles / modes (`parallelization_performance/modes.jl`)

Each mode corresponds to a `ParallelProfile` (R0–R5) and toggles a specific
combination of outer-loop backend and inner-loop/callback parallelism:

| Mode | Profile | Outer backend | Inner/callback parallelism | Adaptive routing | Exercises |
|---|---|---|---|---|---|
| `serial` | R0 | none | off | no | True serial baseline, no parallelism anywhere |
| `outer_threads` | R1_a | threads | off | no | Outer-loop (across-spacecraft) thread parallelism only |
| `outer_process` | R1_b | process | off | no | Outer-loop (across-spacecraft/trial) process parallelism only |
| `inner_only` | R2 | none | auto (density/control/thermal/multibody/effector) | no | Inner-loop/callback parallelism only, no outer parallelism |
| `outer_inner_static` | R3 | threads | auto, statically scheduled | no | Combined outer+inner parallelism with a fixed (static) schedule |
| `outer_inner_adaptive` | R4 | threads | auto, dynamically scheduled | yes | Combined outer+inner parallelism with adaptive/dynamic routing |
| `full_smart` | R5 | threads | auto + thermal always on, persistent hints | yes | The full adaptive routing policy SpaceAGORA ships with, including persistent scheduling hints and control-tail guarding |

## Preview mode

`--preview` (or `SPACEAGORA_PPB_PREVIEW=1`) caps constellation size at 64
spacecraft, MC samples at 16, process workers at 4, and repeats at 2, so the
same phase structure can be smoke-tested on a laptop before committing to a
full run on the benchmark machine.
