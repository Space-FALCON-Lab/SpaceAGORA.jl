# Parallelization Performance Cases

This study defines its case catalog in `cases.jl` and selects default run sets in
`cli.jl`. Cases are grouped by the kind of parallel work they are meant to
exercise: simple gravity, higher-fidelity single/few-spacecraft dynamics,
larger constellations, and Monte Carlo sweeps.

## Run Commands

Run the minimal integration test with:

```bash
julia --project=. benchmarks/studies/parallelization_performance.jl test
```

Run the short smoke check with:

```bash
julia --project=. benchmarks/studies/parallelization_performance.jl smoke
```

Run the full benchmark suite with the machine-scaled default thread ladder:

```bash
julia --project=. benchmarks/studies/parallelization_performance.jl full --outdir=output/performance/parallelization_performance_manual
```

Run the standard benchmark-machine protocol wrapper with the same
machine-scaled defaults:

```bash
bash benchmarks/studies/parallelization_performance/protocol.sh
```

The protocol wrapper changes into the repo root, defaults to the `full` profile,
creates a timestamped output directory under `output/performance/`, sets
`OPENBLAS_NUM_THREADS=1` and `GKSwstype=100`, and then launches the Julia study
entrypoint.

To override the thread ladder explicitly, pass `--threads`:

```bash
julia --project=. benchmarks/studies/parallelization_performance.jl full --threads=1,2,4,8,12,16,24
```

Or set the protocol wrapper's thread ladder through the environment:

```bash
SPACEAGORA_PPC_THREADS=1,2,4,8,12,16,24 bash benchmarks/studies/parallelization_performance/protocol.sh
```

## Default Profiles

### Test

The test profile is a fast execution-path check, not a performance-valid study.
It runs `montecarlo_multi_sat` for two 10-second samples with one repeat, no
warmup, 16 trajectory parity samples, one process worker, and at most two Julia
threads. It covers `serial`, `outer_threads`, `outer_process`, `inner_only`, and
`full_smart`.

Its pass criteria are successful routing, solves, parity, and artifact
generation. Speedup values from this profile should not be used in reports.

### Smoke

The smoke profile is a short correctness and wiring check. It runs:

| Case | Purpose |
| --- | --- |
| `single_inverse_square_vacuum` | Baseline 1-spacecraft inverse-square gravity case. |
| `gravity_4sat_inverse_square_vacuum` | Small constellation gravity-only case. |
| `montecarlo_high_accuracy` | One sampled high-accuracy gravity case with randomized orbital altitudes. |

Smoke trajectory parity checks run:

| Case | Purpose |
| --- | --- |
| `single_harmonics_l20_vacuum` | Single-spacecraft L20 harmonics parity. |
| `gravity_16sat_l20_vacuum` | Multi-spacecraft L20 harmonics parity. |
| `montecarlo_high_accuracy` | Monte Carlo high-accuracy parity. |

### Full

The full profile runs all cataloged performance cases by default and runs parity
on a representative subset:

| Parity Case | Purpose |
| --- | --- |
| `single_harmonics_l20_vacuum` | Single-spacecraft harmonics correctness. |
| `gravity_16sat_l20_vacuum` | Multi-spacecraft harmonics correctness. |
| `articulated_1sat_fullstack` | Attitude, aero, thermal, and harmonics correctness. |
| `multi_16_aero_surrogate_cached` | Multi-spacecraft aero correctness. |
| `montecarlo_high_accuracy` | High-accuracy Monte Carlo correctness. |
| `montecarlo_multi_sat` | Multi-spacecraft Monte Carlo correctness. |
| `montecarlo_mars_aerobraking` | Mars aerobraking Monte Carlo correctness. |

## Case Catalog

### Gravity-Only Cases

These cases use Earth-centered trajectories with no atmosphere. They isolate
gravity/RHS costs and constellation scaling without aero or attitude coupling.

| Case | Spacecraft | Effectors | Notes |
| --- | ---: | --- | --- |
| `single_inverse_square_vacuum` | 1 | Inverse-square gravity | Minimal baseline for serial overhead and solver cost. |
| `gravity_4sat_inverse_square_vacuum` | 4 | Inverse-square gravity | Small multi-spacecraft scaling case. |
| `gravity_16sat_inverse_square_vacuum` | 16 | Inverse-square gravity | Medium constellation scaling case. |
| `gravity_64sat_inverse_square_vacuum` | 64 | Inverse-square gravity | Large constellation scaling case. |
| `gravity_256sat_inverse_square_vacuum` | 256 | Inverse-square gravity | Very large constellation scaling case. |
| `gravity_1024sat_inverse_square_vacuum` | 1024 | Inverse-square gravity | Extreme constellation overhead and scheduling case. |
| `single_harmonics_l20_vacuum` | 1 | L20 gravity harmonics | Single high-cost gravity model. |
| `gravity_4sat_l20_vacuum` | 4 | L20 gravity harmonics | Small constellation with expensive gravity. |
| `gravity_16sat_l20_vacuum` | 16 | L20 gravity harmonics | Medium constellation with expensive gravity. |
| `gravity_64sat_l20_vacuum` | 64 | L20 gravity harmonics | Large constellation with expensive gravity. |
| `gravity_256sat_l20_vacuum` | 256 | L20 gravity harmonics | Very large constellation with expensive gravity. |
| `gravity_1024sat_l20_vacuum` | 1024 | L20 gravity harmonics | Extreme expensive-gravity scheduling case. |

The harmonics cases use `data/Gravity_harmonics_data/EarthGGM05C.csv` when it is
available. If that file is missing, the builder falls back to J2 gravity.

### Few-Spacecraft High-Fidelity Cases

These cases keep spacecraft count low and increase per-spacecraft physics cost.

| Case | Spacecraft | Effectors | Orientation | Notes |
| --- | ---: | --- | --- | --- |
| `single_harmonics_l50_vacuum` | 1 | L50 gravity harmonics | Off | Tight-tolerance, expensive gravity baseline. |
| `srp_heavy_high_area` | 1 | Inverse-square gravity plus two SRP models | On | High-area articulated spacecraft with stacked SRP work. |
| `articulated_1sat_fullstack` | 1 | L20 harmonics plus aero | On | Includes analytic atmosphere, attitude, thermal path, and many panels. |

### Many-Spacecraft High-Fidelity Cases

These cases combine constellation size with more expensive environmental models.

| Case | Spacecraft | Effectors | Atmosphere | Notes |
| --- | ---: | --- | --- | --- |
| `multi_16_aero_surrogate_cached` | 16 | Inverse-square gravity plus aero | Exponential Earth atmosphere | Exercises multi-spacecraft aero and density-cache behavior. |
| `multi_64_high_fidelity` | 64 | L20 harmonics, SRP, and aero | Exponential Earth atmosphere | Stresses mixed expensive effectors across a larger constellation. |

### Monte Carlo Cases

Monte Carlo cases run multiple seeded samples when `--mc-samples` is greater
than one. Each sample perturbs the initial conditions while keeping the case
structure fixed.

| Case | Spacecraft | Effectors | Randomized Inputs | Notes |
| --- | ---: | --- | --- | --- |
| `montecarlo_high_accuracy` | 1 | L50 gravity harmonics | Apoapsis/periapsis altitude jitter | Tight-tolerance high-accuracy gravity sweep. |
| `montecarlo_multi_sat` | 4 | L20 gravity harmonics | Per-spacecraft altitude and true-anomaly jitter | Multi-spacecraft Monte Carlo gravity sweep. |
| `montecarlo_mars_aerobraking` | 1 | Mars inverse-square gravity plus aero | Apoapsis, periapsis, and true-anomaly jitter | Mars aerobraking-style Monte Carlo case with exponential atmosphere. |

## Shared Setup

Test-profile missions are 10 s. Most Earth cases use a 120 s smoke mission and
a 1800 s full mission. The
full-stack, aero-surrogate, and high-fidelity constellation cases use shorter
90 s smoke and 1200 s full missions. Mars aerobraking uses 120 s smoke and
1800 s full missions.

The default solver mode is `auto_stiff`, unless overridden through
`--solver-mode` or `SPACEAGORA_PPC_SOLVER_MODE`.

The full-profile thread ladder is scaled to the machine's `Sys.CPU_THREADS`.
On a 24-thread machine, the default ladder is `1,2,4,8,12,16,24`. Machines with
64 or more threads keep the original `1,2,4,8,16,32,64` ladder. Smaller
machines use fewer unique thread points rather than oversubscribing by default.

## Mode Coverage

Test runs `serial`, `outer_threads`, `outer_process`, `inner_only`, and
`full_smart`. Smoke runs `serial` and `full_smart`.

Full runs:

| Mode | Intent |
| --- | --- |
| `serial` | Baseline with inner and outer parallelism disabled. |
| `outer_threads` | Outer threaded parallel route with inner parallelism disabled. |
| `outer_process` | Outer process-based route with inner parallelism disabled. |
| `inner_only` | Inner RHS/callback/effector parallelism without outer parallelism. |
| `outer_inner_static` | Outer threaded plus inner parallelism with static policy. |
| `outer_inner_adaptive` | Outer threaded plus inner parallelism with adaptive policy. |
| `full_smart` | Adaptive threaded route with persistent policy hints and thermal parallelism enabled. |
