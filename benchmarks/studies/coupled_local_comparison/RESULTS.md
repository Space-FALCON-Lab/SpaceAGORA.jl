# Coupled-case comparison: results

Measured with this study's `run.sh`, `compare.jl` and `profile_allocations.jl`
(see `README.md`). Every ratio below is between runs on one host. Outputs stay
local, under `output/coupled_local_comparison/` in the run checkout.

## Stage A: 32 spacecraft, before and after PR #222

Host: workstation, AMD Ryzen 9 9900X (12 cores, 24 threads), run 2026-10-06.
Case `stack32_e6_actuated_saved`, 120 s mission, 3 warm-up and 6 timed solves
per point; times are medians of the timed solves. `pre` is `ea42d2946`, `post`
is `402018b98`.

| | pre | post |
|---|---|---|
| Serial wall time | 29.40 s | 5.32 s (5.53x faster) |
| Allocation per solve | 140,981 MiB | 4,577 MiB (-96.8%) |
| GC time per serial solve | 8.20 s | 1.23 s |

Speedup over the same code's serial run, range over the five parallel routes:

| Threads | pre | post |
|---|---|---|
| 2 | 1.21-1.27x | 1.63-1.72x |
| 4 | 1.36-1.43x | 2.24-2.38x |
| 8 | 1.20-1.25x | 2.89-2.98x |
| 12 | 1.08-1.15x | 2.86-3.06x |

- PR #222's speedup at a fixed route and thread count grows with threads (about
  7x at 2, 9x at 4, 13x at 8, 15x at 12) because the pre-#222 code slows down
  past 4 threads; it is not a scaling result.
- The adaptive route (`predictive`, R7) was 1.03-1.05x the time of the fastest
  static route at the same thread count on `post` and 1.02-1.06x on `pre`. It
  was never faster. It allocates exactly what `inner_only` allocates at every
  thread count, so it ends on that route.
- The routes differ by 3-8% at one thread count, about the spread of the
  repeats, so six repeats do not resolve which static route is best.
- Step times and final states agree byte for byte across all 252 timed solves:
  within each point, across every route and thread count, and between `pre` and
  `post`.

## Allocation profile after PR #222

Host and case as Stage A, code `a5d3fccc9` (main), serial, one solve after two
warm-ups, `Profile.Allocs` at sample rate 0.001 (24,729 samples).

- One solve allocates 4,577.5 MiB with trajectory output and 4,561.7 MiB
  without it, so saving is not where the allocation comes from.
- Innermost frame on each sampled stack, in any package:

  | Share | Site |
  |---|---|
  | 47.2% | `src/simulation/engine/dynamics_rhs.jl:293`, the call into `_accumulate_control_effectors_from_tuple!` |
  | 14.9% | `src/simulation/engine/dynamics_rhs.jl:162`, `any(_wrench_method_available(effector) for effector in dynamic_effectors)` |
  | 4.4% | `src/simulation/engine/setup.jl:547`, `_dynamic_effectors_parallel_supported` |

- By allocated type: `ODEParams{...}` 41.6% (about 2.4 KB each), the
  dynamic-effector tuple 17.7%, `NTuple{32, MagneticMomentumManagerModel}`
  4.2%. The profiler's own buffer is 12.4% and is charged to the `solve` call
  site, which therefore reads higher than its real share.

The allocations sit at the call lines themselves, not in the functions they
call. The likely mechanism, not yet confirmed by a change: the control-effector
tuple's type is erased in `ControlModel`, so the PR #222 function barrier is
still a dynamic call, and passing the inline immutable `ODEParams` through a
dynamic call copies it to the heap on every call, once per spacecraft per RHS
evaluation. Line 162 probably boxes the heterogeneous effector tuple through
the generator. GC is 23% of serial time and about 30% at 12 threads after
PR #222, so these sites are the first candidates for the remaining thread
ceiling.

## Stage B: spacecraft-count sweep

Pending: TRX50 (AMD Threadripper PRO 9985WX, 64 cores), N = 32, 64, 128, 256,
serial and five routes at 12 and 32 threads, code `402018b98`.
