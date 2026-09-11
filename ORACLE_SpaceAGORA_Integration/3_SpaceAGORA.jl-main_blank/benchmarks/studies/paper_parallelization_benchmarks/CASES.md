# Paper Parallelization Benchmarks — Test Cases

This study (`benchmarks/studies/paper_parallelization_benchmarks/`) runs the
benchmark suite backing the parallelization paper. It wraps the lower-level
`parallelization_performance` study (`main.jl`, `cli.jl`, `cases.jl`,
`modes.jl` in `benchmarks/studies/parallelization_performance/`), grouping its
cases and routing profiles into five phases (B1–B5). Each phase isolates a
specific axis of SpaceAGORA's parallelization stack so the paper can make a
narrow, defensible claim per figure/table.

## LaTeX table for the paper

```latex
\begin{table}[htbp]
\centering
\caption{Summary of the paper parallelization benchmark phases (B1--B5) and the parallelization axis each isolates.}
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
\bottomrule
\end{tabular}
\end{table}
```

## Phases (B1–B5)

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
