# Paper Parallelization Benchmarks — Test Cases

This study (`benchmarks/studies/paper_parallelization_benchmarks/`) runs the
benchmark suite backing the parallelization paper. It wraps the lower-level
`parallelization_performance` study (`main.jl`, `cli.jl`, `cases.jl`,
`modes.jl` in `benchmarks/studies/parallelization_performance/`), grouping its
cases and routing profiles into phases (B1–B8). Each phase isolates a
specific axis of SpaceAGORA's parallelization stack so the paper can make a
narrow, defensible claim per figure/table.

Two measurement properties are load-bearing for every phase below:

- **Warm-up is mandatory.** Each `(case, mode, threads, mc_samples)` point runs
  in a fresh Julia subprocess, so without a warm-up solve the timed run is the
  process's *first* run and its wall time is dominated by JIT compilation of the
  RHS/solver stack, not by the simulation (measured: 21.4 s first solve vs
  0.017 s steady-state solve for `gravity_16sat_l20_vacuum`). Phase `warmup`
  counts are forwarded to the worker via `--warmup`; process-pool provisioning
  and per-worker warm-up likewise happen before the clock starts.
- **Repeats run inside one worker process, not one process each.** A worker
  spends ~80 s on Julia startup plus JIT of the RHS/solver stack before it can
  time anything — measured at 79 s per row for a case whose solve is 0.016 s — so
  a subprocess per repeat paid that cost `repeats` times to collect `repeats`
  samples of the same point, which at three repeats was two thirds of the whole
  run. The controller now launches one subprocess per
  `(case, mode, threads, mc_samples)` point and the worker loops the repeats
  inside it (`PPCConfig.worker_repeats`), re-running warm-up before each so every
  repeat is measured identically. Two consequences: `--resume` granularity is now
  a whole point rather than a single repeat, and worker CSVs are named
  `perf_<case>_<mode>_t<N>_mc<M>.csv` with no `_r<n>` suffix.

  This also *improves* measurement validity, which was not the motivation.
  A single warm-up solve does not fully warm the stack: in a worker running three
  repeats of `stack32_e5_6dof` after one warm-up, repeat 1 measured 0.139 s
  against 0.019 s and 0.016 s for repeats 2 and 3. Under per-repeat subprocesses
  *every* row was a "repeat 1" and carried that residue, and the median of three
  cold subprocesses inherited it wholesale; a median over three in-process
  repeats does not. B9–B14 use `warmup = 2` so even repeat 1 is clean.
- **The router ladder stops at 32 threads** (`PPB_ROUTER_LADDER_MAX_THREADS`),
  even on a 64-core host. B9–B14 are about which route the router picks, not how
  far the machine scales — that is B1/B2/B7's job, and they run the full ladder.
  A top rung that spends every core mostly measures oversubscription: §9 of the
  findings doc records every workload regressing past the physical-core count,
  and the benchmark host is shared with other users' jobs, so the widest rung is
  the one most likely to be contending with something else rather than measuring
  routing. B13's budget grid is rebased to the same 32-core total.
- **Thread ladders are in physical cores.** The RHS hot path is SIMD/FP-bound,
  gains nothing from SMT, and regresses sharply when SMT siblings are
  oversubscribed, so a ladder topped out at `Sys.CPU_THREADS` puts its worst
  data point at the end of the curve. See `SPACEAGORA_PPC_PHYSICAL_CORES`.

Pre-solve RHS-plan auto-calibration (`SPACEAGORA_RHS_CALIBRATE`) follows each
mode's own adaptive flag: `off` for the static profiles R0–R3, `auto` for the
adaptive ones R4/R5. It was previously forced off for *every* mode, on the
reasoning that pinning a plan via `rhs_plan_override` short-circuits the routing
policy the profile ladder exists to compare. That holds for the static profiles,
which are meant to measure one fixed route, but applying it to R4/R5 disabled the
mechanism those profiles are defined by — calibration is part of the adaptive
policy SpaceAGORA ships with. The correction is worth a factor of five: on
`heavy_1024sat_fullstack_1hr` at 12 threads, 29.90 s with calibration off vs
5.71 s with it on. Every "multi-effector constellations do not scale" number
taken under the blanket-off setting understates the shipped configuration by
that factor. See `modes.jl` and §8b of
`docs/architecture/parallelization_measurement_and_allocation_findings.md`.

**Two things that do not measure what their names suggest**, both relevant when
reading older results:

- `num_steps_to_save` is inert. Nothing outside
  `analysis/verification/telemetry_verification/` reads it, and the solve runs
  `save_everystep` with no `saveat`, so varying it changes no work at all. The
  case built on it (`gravity_16sat_l20_vacuum_sparse_output`) has been removed
  rather than repaired. Output cadence is really controlled by
  `simulation_settings.results`, which gates `get_data_saving_callback`'s
  `SavingCallback(...; saveat=data_rate)` — and this harness hardcoded
  `results=false`, so no benchmark case before B14 produced any trajectory
  output. `ppc_build_config` now takes `results` and `data_rate`.
- `rhs_execution_mode` in the raw CSV is the mode *requested*, always `"auto"`,
  not the route selected. The selected route now appears separately in the
  `rhs_plan_*` and `policy_*` columns (see "Router-selection telemetry" below).

## LaTeX table for the paper

```latex
\begin{table}[htbp]
\centering
\caption{Summary of the paper parallelization benchmark phases (B1--B14) and the parallelization axis each isolates.}
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
B6 & Small-workload control: sub-second constellation and MC cases & Demonstrates that below $\sim$3\,s serial no routing profile is distinguishable \\
B7 & Heavy constellation scaling, 256--4096 sat, multi-hour missions, $\sim$8\,s serial baselines & Thread ladder on workloads large enough for scaling to be observable \\
B8 & Heavy Monte Carlo throughput, seconds-long aerobraking arcs, 1--64 process workers & Process-level throughput free of per-worker startup bias \\
B9 & Spacecraft count 256--4096, L50 harmonics, vacuum, 1\,h & Router regret vs.\ satellites per RHS call \\
B10 & Atmosphere fidelity at $N{=}256$: none, exponential, GRAM surrogate, live GRAM, live GRAM + third-body & Router regret vs.\ native-library serialisation \\
B11 & Force/actuator model count, models added one at a time; separate $N{=}32$ sub-ladder for the actuator rung & Router regret vs.\ effector-tuple heterogeneity \\
B12 & Interacting $N$-satellite propagation vs.\ $N$ independent trials at matched work & Router regret vs.\ workload coupling \\
B13 & One fixed 32-core budget split six ways between processes and threads & Router regret vs.\ the shape of the parallel budget \\
B14 & Mission duration 15\,min--6\,h and trajectory output off/10\,s/1\,s & Router regret vs.\ duration and the serial output spine \\
\bottomrule
\end{tabular}
\end{table}
```

## Phases (B1–B14)

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

### B6 — Small-Workload Control
Cases: the `many_sat_high_fidelity` spacecraft-count ladder (`multi_4`,
`multi_16`, `multi_64`, `multi_256`), `gravity_16sat_l20_vacuum`,
`montecarlo_mars_gram_live`, `montecarlo_multi_sat`,
`articulated_1sat_fullstack`, `multi_8sat_magnetorquer_attitude`, and
`gravity_16sat_l20_vacuum_longmission`. Modes: `serial`, `outer_threads`,
`outer_process`, `outer_inner_adaptive`, `full_smart`. Thread mode:
`:low_high` (1 and the machine maximum). MC samples: 1 and 8.

**This was the first attempt at the expanded router evaluation, and its numbers
are not reportable as router performance.** Every case it draws on resolves in
0.017–2.9 s post-warm-up, so fixed per-solve setup dominates and no route can
distinguish itself: the run of 2026-08-18 measured 1.0× speedup on 9 of 11 cases
and regret of 0.0–3.3%, all of it inside the ~16% run-to-run variance §7.4 of
the findings doc records for measurements that size. B9–B14 are the actual
expanded evaluation.

It is retained because "below roughly a 3 s serial baseline no routing profile
is distinguishable" is a real boundary on the paper's claim and this phase is
the evidence for it. Reporting flags its rows `below_noise_floor` and excludes
them from the per-axis regret statistics.

`multi_16_gram_live` is defined in `cases.jl` but excluded here for a historical
unbounded native memory leak. §8.1 of the findings doc reports that neither that
leak nor `multi_4_gram_live`'s step-size collapse reproduces on the current
branch after the ephemeris-sharing fixes; B10's GRAM rungs are the replacement
either way, and are sized to be measurable where these were not.

### Running order of B9–B14

`--phases` is order-respecting: the phases run in the order named on the command
line, and fall back to catalog order only when the flag is omitted. The expanded
set is roughly 20 h of wall clock on a 12-core box against an overnight budget,
so which phases run first decides which axes exist if the run is stopped early
(`--resume` picks up the rest). The default order for a short session is

```
--phases=B9,B10,B12,B11,B13,B14
```

which is by evidence-per-hour rather than by phase number:

| Order | Phase | Est. wall clock | Why here |
|---|---|---:|---|
| 1 | B9 | ~1.3 h | Cheapest phase, and it exercises the whole pipeline (six modes, three budgets, parity, regret, telemetry) on real data before anything expensive starts. A harness fault shows up in the first 20 min instead of two hours into a GRAM rung. |
| 2 | B10 | ~3.5 h | The load-bearing axis. Review point 7's "R5 is ~25 % slower than the best static route on GRAM" is the strongest argument against the contribution, and this is the only phase that isolates density-model fidelity. |
| 3 | B12 | ~2.3 h | The manuscript's third contribution (coupled multi-spacecraft vs. process-isolable campaigns) at matched satellite-hours. |
| 4 | B11 | ~3.9 h | Model-count axis. Valuable but the most expensive per rung, and §7.1 already characterises the effector-tuple cost qualitatively. |
| 5 | B13 | ~2 h | Budget-split axis. Self-contained (one case), so it can be run on its own later at no loss. |
| 6 | B14 | ~7 h | Last because calibration already answered half of it: cadence is flat (`_10s` vs `_1s` differ by 3 %, both 62x `_none`), so the two most expensive rungs in the whole set are confirming a known-null result. The duration half is the part worth the time. |

Within each phase the case list is already ordered cheapest-first, so a phase cut
off part-way still yields its low rungs.

### B9 — Router Evaluation: Spacecraft Count
Cases: `gravity_{64,256,1024,4096}sat_l50_vacuum_1hr`. All six modes.
Thread mode: `:router_ladder` (`[1, 8, 32]`).

Varies satellites per RHS call at fixed physics (L50 harmonics, vacuum, no
callbacks) and fixed 1 h duration, so outer-loop routing is the only thing that
can move the number. The ladder straddles the measurability floor on purpose —
64 satellites solve in ~0.15 s serial and 4096 in ~8 s — because *where* routing
stops being decidable is part of the answer.

### B10 — Router Evaluation: Atmosphere and GRAM Usage
Cases: `atmo256_vacuum_10min`, `atmo256_exponential_10min`,
`atmo256_gram_surrogate_10min`, `atmo256_gram_live_10min`,
`atmo256_gram_live_nbody_10min`. All six modes. Thread mode: `:router_ladder` (`[1, 8, 32]`).

Five rungs of increasing native-library serialisation at fixed N=256, fixed
600 s duration and fixed L50 harmonics: no atmosphere → analytic exponential →
precomputed GRAM surrogate table → live native GRAM (one shared instance,
serialised by its own `instance_lock`) → live GRAM plus Sun/Moon third-body
gravity, which adds CSPICE under `SPICE_LOCK` as a second serialisation point.

This is the load-bearing phase. Review point 7 reports R5 running ~25% slower
than the best tested route on the GRAM workload and slower than serial; nothing
in the harness previously swept GRAM fidelity as an isolated variable, so that
finding could be measured but not localised. 600 s rather than 1 h because that
is the duration the `heavy_<N>sat_gram_nbody_l50` family is validated at up to
N=256.

### B11 — Router Evaluation: Force and Actuator Model Count
Cases: the `stack<N>_e1..e6` ladder. All six modes. Thread mode:
`:router_ladder`.

Six rungs at fixed N and fixed 600 s duration, each adding exactly one model:
L20 harmonics → +SRP → +aero and exponential density → +Sun/Moon third body →
attitude propagation on → +LVLH controller and per-satellite magnetorquer
actuators. §7.1 of the findings doc identifies the heterogeneous effector tuple
as a first-order scheduling cost, but no previous case swept model count in
isolation. `e5` and `e6` are split so that enabling attitude propagation and
adding the actuator path are measured separately rather than confounded.

**The ladder runs at a smaller N than the other axes because the actuator rung
is quadratic in satellite count.** `get_control_callbacks` builds one
`PeriodicCallback` per control effector and each loops over *all* satellites, so
N per-satellite `MagneticMomentumManagerModel`s cost N² `calcControlEffect!`
calls per tick — N intended, N(N−1) against a `sat_idx` the model does not own.
Measured post-warm-up: `stack256_e6_actuated` takes 33.4 s for a *10 s* mission
against 0.13 s for the otherwise identical `stack256_e5_6dof`. That is a real
defect in the control path (the shared-model-with-per-spacecraft-vectors shape
`BaseThrusterModel` enforces avoids it), not a benchmark artifact, but fixing it
is out of scope; the ladder is sized so every rung is measurable and the
quadratic is reported instead.

### B12 — Router Evaluation: Interacting vs. Independent Propagation
Cases: `interact_{16,64,256}sat_1hr` and `independent_1sat_1hr`. Modes: the
static ladder plus `outer_process`, minus `inner_only`. MC samples: 16/64/256.
Thread mode: `:router_ladder` (`[1, 8, 32]`).

A matched pair at equal satellite-hours: `interact_<N>sat_1hr` propagates N
satellites as one coupled simulation, `independent_1sat_1hr` at `mc_samples=N`
propagates the same work as N independent solves. Same planet, orbit, effectors,
density model, tolerances and duration, so coupling is the only variable. B6
gestured at this axis by putting constellation and Monte Carlo cases in one
phase, but those differ in physics, duration and size as well, so nothing there
isolated coupling. This pair is the direct test of the manuscript's third
contribution.

`independent_1sat_1hr` is deliberately *not* jittered the way the `montecarlo_*`
cases are: an initial-condition spread gives each sample a different step count,
which appears as load imbalance and would confound the comparison.

### B13 — Router Evaluation: Thread vs. Process Budget Split
Case: `montecarlo_heavy_aerobraking` at 64 samples. Modes: `outer_process`,
`outer_threads`, `outer_inner_adaptive`, `full_smart`. Budget grid:
`(workers, threads)` ∈ {(1,32), (2,16), (4,8), (8,4), (16,2), (32,1)}.

Every entry multiplies out to 32 cores, so the sweep varies the *shape* of the
parallelism and not its amount. Nothing else in the harness tests hybrid splits:
B4 and B8 both pin `thread_mode = :single` and vary workers alone, and every
constellation phase runs one process. The nested case is what
`SPACEAGORA_OUTER_PARALLEL_ACTIVE` / `SPACEAGORA_INNER_THREAD_BUDGET` exist to
arbitrate — and §3 of the findings doc records asserting `outer_active`
unconditionally pinning the RHS to one worker — which makes it the configuration
the router is most likely to get wrong and the least covered by evidence. Uses
the `budget_grid` phase field rather than `worker_ladder`, whose cross product
with the thread ladder would also vary the total budget.

### B14 — Router Evaluation: Mission Duration and Output Cadence
Cases: `gravity_1024sat_l50_vacuum_15min`, `gravity_1024sat_l50_vacuum_1hr`,
`heavy_1024sat_l50_6hr` (duration), and `cadence_1024sat_{none,60s,10s,1s}`
(cadence). All six modes. Thread mode: `:router_ladder` (`[1, 8, 32]`).

Two sub-ladders at N=1024, L50, vacuum, sharing the 1 h rung. The cadence half
is new work rather than a re-run of B6's: that case varied `num_steps_to_save`,
which reaches nothing, and this harness had trajectory output switched off
entirely. The real path is `simulation_settings.results` gating
`get_data_saving_callback`'s `SavingCallback(...; saveat=data_rate)`, both now
exposed through `ppc_build_config`. At 1024 satellites the callback's
per-satellite snapshot is a genuinely serial spine, so this ladder measures an
Amdahl term the router has no visibility into — which is the result either way
it comes out.

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
- **`atmosphere_ladder`** — `atmo256_{vacuum,exponential,gram_surrogate,
  gram_live,gram_live_nbody}_10min`. N=256, 600 s, L50 harmonics held fixed;
  only the density path varies. B10.
- **`effector_ladder`** — `stack<N>_e1..e6` for N ∈ {32, 64, 256}. Fixed N and
  duration; one model added per rung. Registered at several N because the
  actuator rung is quadratic in N (see B11); a phase picks one per sub-ladder.
- **`coupling_mode`** — `interact_{16,64,256}sat_1hr` and
  `independent_1sat_1hr`. Matched satellite-hours, arranged as one coupled
  simulation or as N independent trials. B12.
- **`duration_cadence`** — `gravity_16sat_l20_vacuum_longmission`,
  `gravity_1024sat_l50_vacuum_15min`, `cadence_1024sat_{none,10s,1s}`. Mission
  length and trajectory-output volume. B14.

### Measured case costs (space-falcon-1, serial, post-warm-up)

Seconds of solve per 10 s of simulated mission, which is what the phase sizing
above was derived from. Multiply by `mission_seconds / 10` for the full-profile
serial baseline.

| Case | s / 10 s mission | Case | s / 10 s mission |
|---|---:|---|---:|
| `atmo256_vacuum_10min` | 0.010 | `stack32_e1_harm` | 0.0009 |
| `atmo256_exponential_10min` | 0.063 | `stack32_e3_aero` | 0.0045 |
| `atmo256_gram_surrogate_10min` | 0.239 | `stack32_e5_6dof` | 0.016 |
| `atmo256_gram_live_10min` | 0.623 | `stack32_e6_actuated` | 0.583 |
| `interact_16sat_1hr` | 0.0035 | `stack64_e5_6dof` | 0.031 |
| `interact_64sat_1hr` | 0.016 | `stack64_e6_actuated` | 2.728 |
| `interact_256sat_1hr` | 0.063 | `stack256_e5_6dof` | 0.130 |
| `independent_1sat_1hr` | 0.0006 | `stack256_e6_actuated` | 33.4 |
| `gravity_1024sat_l50_vacuum_15min` | 0.077 | `cadence_1024sat_none` | 0.051 |
| `cadence_1024sat_10s` | 3.09 | `cadence_1024sat_1s` | 3.17 |

Two things this table is load-bearing for:

- **The actuator rung is O(N^2.2).** `stack{32,64,256}_e6_actuated` at
  0.583 / 2.728 / 33.4 s is a 4.7× step for each doubling of N, against a
  near-linear `e5`. That is why B11 runs two sub-ladders instead of one.
- **Output cadence is flat.** `cadence_1024sat_10s` and `_1s` differ by 3%
  while both are 62× `_none`. Turning output on is the variable; how often it
  saves is not.

## Router-selection telemetry

The raw CSV records what the router *did*, not just what it was asked for. Until
B9–B14 the only routing column was `rhs_execution_mode`, hardcoded to `"auto"`,
which meant a finding like review point 8's "the router selects poorly for the
GRAM case" could be measured but never attributed to a decision.

| Column | Meaning |
|---|---|
| `rhs_plan_source` | `none` (calibration never ran), `cache`, or `sweep` |
| `rhs_plan_mode` | the RHS plan installed, e.g. `satellite_batch` or `flat` |
| `rhs_plan_allotment` | that plan's worker allotment |
| `policy_last_mode` / `policy_last_allotment` | last inner-policy decision |
| `policy_last_outer_active` | whether an enclosing outer split was claimed |
| `policy_decisions_total` / `policy_adaptive_decisions_total` | decision counts |
| `policy_threads_enabled_total` | decisions that actually threaded |
| `policy_discarded_by_route_total` | policy decisions the route then overrode |

The selected plan lives in per-solve `shared_buffers.rhs_plan_override` and is
gone once the solve returns, so `record_rhs_plan_selection!`
(`src/parallel/policy/policy_telemetry.jl`) copies it into the process-level
policy telemetry, which `ppc_policy_snapshot` reads back after the timed batch.
Telemetry is reset immediately before the timed region so warm-up decisions are
not mixed in.

**Process-route caveat.** These counters are per-process. Under `outer_process`
the solves run on Distributed workers, so the controller's snapshot is empty and
the columns come back as `none`/0. That is accurate rather than missing — the
controller genuinely made no routing decision — but it means the selected-route
columns are only informative for the in-process (serial/threads) routes.

## Outer-backend resolution for the adaptive profiles

`ppc_run_sample_batch` dispatches the outer split itself, branching on the
mode's `backend` string. R4 (`outer_inner_adaptive`) and R5 (`full_smart`) were
declared with `backend="threads"`, which meant the two profiles whose purpose is
to *choose* an outer route never reached `select_outer_route!` — the harness
answered the question it existed to ask. The shipped profile definitions say
otherwise: `profile_definitions.jl` declares `outer_backend=:auto` for R3, R4 and
R5.

R4/R5 now declare `backend="auto"`, and `ppc_resolve_outer_backend` builds
`OuterRouteFeatures` from the case's own `SimulationConfiguration` (via the
shipped `campaign_route_features`), calls `select_outer_route!`, and substitutes
the answer before warm-up so warm-up, the recorded env and the timed dispatch all
agree. R3 keeps its explicit `threads` pin on purpose: it is the fixed hybrid
baseline the adaptive profiles are scored against, so it must not move with the
router.

Single-sample batches resolve to `"threads"` unconditionally — with one
simulation there is no outer split to route, so the resolution is moot and the
literal string stays comparable with the constellation phases measured before the
change.

**What this invalidated.** Only phases with `mc_samples > 1` have an outer route
to select, so B12 and B13 were affected and B9/B10/B11/B14 were not. The
2026-08-23 B12 results were discarded (see
`output/performance/paper_benchmarks/20260823_130939_invalidated/README.md`): on
`independent_1sat_1hr` at 256 samples the router selects `:process` (8.12x) while
the pin forced `:threads` (2.03x), so the "+301% router regret" read off those
rows described the mode table rather than the router.

## Effective core budget in regret

Regret is only comparable between points given the same machine, and
`thread_count` is the wrong proxy for that on the process route:
`ppc_run_sample_batch` spreads it across `min(mc_samples, process_workers)`
single-threaded worker processes regardless of `thread_count`, so an
`outer_process` row recorded at `thread_count=1` actually ran on 12 cores.

`_ppb_add_effective_cores!` therefore adds `effective_cores` — `1` for `serial`,
`min(process_workers, mc_samples)` for a point whose `outer_backend_actual` is
`process`, and `thread_count` otherwise — and both regret figures key on it
instead of on `thread_count`. Thread-backed points are unchanged, because there
`effective_cores == thread_count`; B9's and B10's numbers are identical before
and after.

Without it the metric charged the router for losing races it entered with a
twelfth of the hardware, which is where the discarded "+1428%" on
`independent_1sat_1hr` came from.

## Reported regret

Two definitions, both in the aggregated CSV and the report:

- `regret_vs_best_static` — best static route at the **same** thread/worker/MC
  point. This is review point 8's formula literally.
- `regret_vs_best_static_oracle` — best static route at **any** budget tested
  for that case. R4/R5 pick the parallel width as well as the route, so the
  matched-budget figure credits them for a budget someone else chose; the oracle
  charges an adaptive route for running wide on a workload that saturates
  narrow.

`best_static_mode` names which static route set the bar. `below_noise_floor`
marks points whose serial baseline is under 3 s; those are excluded from the
per-axis summary and never quoted as router performance. The summary itself
(worst, median, and fraction of points within 10%) is written to
`router_regret_by_axis_<stamp>.csv`.

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
