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

Host: TRX50, AMD Threadripper PRO 9985WX (64 cores), run 2026-10-07, remote
job `20261007-101748-2161710`. Code `402018b98` (after PR #222, before the
allocation fixes below). Same case family and protocol as Stage A, at
N = 32, 64, 128. Ratios below are within this host only.

| N | Serial | Allocation | GC (serial) | Best at 12 threads | Best at 32 threads |
|---|---|---|---|---|---|
| 32 | 5.70 s | 4,577 MiB | 1.36 s (24%) | 1.75 s (3.26x) | 1.83 s (3.12x) |
| 64 | 13.32 s | 11,642 MiB | 3.44 s (26%) | 3.97 s (3.35x) | 4.00 s (3.33x) |
| 128 | 90.74 s | 32,682 MiB | 66.26 s (73%) | 20.45 s (4.44x) | 18.38 s (4.94x) |

- Serial time grows faster than N: 2.3x from 32 to 64 and 6.8x from 64 to
  128. Allocation grows 2.5x and 2.8x over the same steps, so most of the jump
  at 128 is collection, not work: GC is 73% of the serial run and 68% of the
  fastest 32-thread run (12.46 of 18.38 s).
- Threading pays off more at larger N, 3.1-3.3x at 32 and 64 against 4.4-4.9x
  at 128, and 32 threads beats 12 only at 128.
- The adaptive route was 1.00-1.08x the fastest static route's time at each
  (N, threads) point, never faster.
- Every N is byte-identical across its routes and thread counts.

N = 256 did not run under the 16 GB per-point cap. Julia's type inference
overflows its stack compiling `initialize_callbacks!` over the callback set,
which holds one periodic callback per control effector (256 of them), and the
process exceeded the cap about four minutes in, before any solve. A rerun with
`MEMORY_MAX=96G` (job `20261007-131000-2637476`) continued past the same
overflow warning and completed:

| N | Cap | Serial | Allocation | GC (serial) | Best at 12 threads | Best at 32 threads |
|---|---|---|---|---|---|---|
| 256 | 96G | 92.52 s | 101,105 MiB | 31.89 s (34%) | 25.22 s (3.67x) | 23.80 s (3.89x) |

Adaptive/best static was 1.01 at 12 threads and 1.00 at 32; all 11 points are
byte-identical.

**The 16 GB cap distorts the N = 128 row.** N = 256 under a 96 GB cap runs
serially in about the same time as N = 128 under 16 GB (92.5 s against
90.7 s), allocates 3.1x as much, and spends half as long collecting (31.9 s
against 66.3 s). Julia sizes its heap target to the cgroup memory limit, so a
solve whose live heap approaches the cap collects far more often. The
workstation sweep below, under a 36 GB cap, confirms it: the two hosts agree
within about 13% at N = 32 and 64, where the cap does not bind, but N = 128
takes 33.8 s serially on the workstation against 90.7 s on TRX50. Read the
TRX50 N = 128 row as a capped measurement, not as scaling. The N = 256 code was
compiled after an inference failure, which may also change how it runs.

## Stage B on the workstation, 36 GB cap

Host and code as in the allocation-fix comparison (main `a5d3fccc9`), run
2026-10-08, every point under a 36 GB cap and started only with at least 40 GB
free. Serial and five routes at 12 threads; N = 128 is the main side of the
fix pair below (serial and four routes at 12 threads).

| N | Serial | Allocation | GC (serial) | Best at 12 threads |
|---|---|---|---|---|
| 32 | 5.03 s | 4,577 MiB | 1.12 s (22%) | 1.62 s (3.10x) |
| 64 | 12.94 s | 11,642 MiB | 2.82 s (22%) | 3.57 s (3.63x) |
| 128 | 33.78 s | 32,682 MiB | 10.06 s (30%) | 9.23 s (3.66x) |
| 256 | 129.61 s | 101,104 MiB | 71.46 s (55%) | 33.54 s (3.87x) |

- Serial time grows 2.6x per doubling from 32 to 128 and 3.8x from 128 to 256;
  allocation grows 2.5-3.1x per doubling. The per-spacecraft control scan is
  quadratic in N, which fits growth above 2x.
- GC rises from 22% to 55% of serial time at N = 256 even with a 36 GB cap; on
  TRX50 under 96 GB it was 34%. At this size the cap is still a variable, so
  N = 256 numbers hold for their stated cap only.
- The best speedup at 12 threads rises slowly with N, from 3.1x to 3.9x.
- Adaptive/best static was 1.06, 1.03, 1.00 and 1.00 at N = 32, 64, 128 and 256.
- Each N is byte-identical across its routes.

## Allocation fixes (branch `perf/coupled-alloc-fixes`)

The profile above led to a separate change: `ODEParams` became a mutable
struct with `const` fields, so dynamic calls pass it by reference instead of
copying it, and three smaller tuple boxes were removed. Against main
`a5d3fccc9` on the Stage A host and case (serial plus four routes at 4, 8 and
12 threads, same protocol):

| | main | fix |
|---|---|---|
| Serial | 5.155 s | 2.589 s (1.99x faster) |
| Allocation per solve | 4,577 MiB | 1,383 MiB (-69.8%) |
| GC (serial) | 1.119 s | 0.098 s |
| Fastest parallel point | 1.683 s | 0.939 s |

Every route and thread count measured was 1.67-1.85x faster, step times and
final states are byte-identical to main, and the full test suite passes. The
fix speeds the serial path more than the parallel one, so speedup over serial
falls (about 3.0x to 2.7x at 8-12 threads) while every absolute time drops.
On TRX50 at N = 128 (remote job `20261007-142758-2901849`, 16 GB cap like its
baseline, the Stage B N = 128 points from `402018b98`; the fix branch's base
also includes #224, which removed two unused aerodynamic methods):

| | baseline | fix |
|---|---|---|
| Serial | 90.74 s | 32.14 s (2.82x faster) |
| Allocation per solve | 32,682 MiB | 13,198 MiB (-59.6%) |
| GC (serial) | 66.26 s | 17.44 s |
| Best at 12 threads | 20.45 s | 9.04 s (2.26x faster) |
| Best at 32 threads | 18.38 s | 8.17 s (2.25x faster) |

Every route was 2.22-2.35x faster at 12 and 32 threads, and the fix's states
are byte-identical to the baseline's. Both sides ran under the 16 GB cap
described above, so this gain includes the cap's effect on collection.

The same pair on the workstation under a 36 GB cap (main `a5d3fccc9` against
the fix, serial and four routes at 4, 8 and 12 threads, 2026-10-08):

| | main | fix |
|---|---|---|
| Serial | 33.78 s | 18.07 s (1.87x faster) |
| Allocation per solve | 32,682 MiB | 13,198 MiB (-59.6%) |
| GC (serial) | 10.06 s | 3.36 s |
| Best at 12 threads | 9.23 s | 5.33 s (1.73x faster) |

Every route was 1.72-1.81x faster at 4, 8 and 12 threads, and states are
byte-identical. So the fix's gain at N = 128 without cap pressure is about
1.8x, close to the 1.7-2.0x at N = 32; the larger factor under 16 GB was
partly relief from the cap.
