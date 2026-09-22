# The native-GRAM memory footprint, measured

## The question

`src/parallel/routing/machine_topology.jl` charges every process worker a
per-spacecraft memory term whenever the workload is native GRAM point density,
and `effective_process_workers` / `mixed_local_slots` turn that charge into a
hard routing decision: fewer than two affordable workers withdraws the process
route entirely.

Until this study the term was a single constant, 90 MB per spacecraft,
introduced in `bcda828cd` (2026-09-04) with the justification "from the 22 GB
single-process footprint measured at 256 spacecraft on an 18 GB laptop". At 256
spacecraft that is 22 GB charged to each worker; at 4096 it is 360 GB, more than
any machine in the lab, so the process route was foreclosed for every large GRAM
constellation.

Meanwhile the committed `benchmarks/studies/paper_scenarios` results record peak
resident set for the same kind of workload and show something two orders of
magnitude smaller. This document reports what the footprint actually is.

## Where the 90 MB came from

It cannot be reconstructed from anything in the repository.

- The constant and its justification appear in exactly two places: the commit
  message of `bcda828cd` and §10 of `docs/policy_v2_session_record_20260903.md`.
  Both give the same sentence and no data.
- The case named as its source, `atmo256_gram_live_10min`, is run by
  `benchmarks/studies/parallelization_performance`. That harness has **no**
  resident-set instrumentation at all -- no `Sys.maxrss`, no `maxrss` column, no
  RSS field in any of its CSVs. There is therefore no recorded measurement
  behind the number, only a reported observation.
- The machine is described as an 18 GB laptop. A single process cannot hold a
  22 GB resident set on an 18 GB machine, so whatever was observed was either
  virtual size, a swap-inclusive figure, or the sum of several processes.

So the honest answer to "which density path was the 90 MB for" is: **none of
them**. It is not a per-path value that this study failed to reproduce; it is
not a measurement of a per-spacecraft term at all. The predicate it feeds,
`_is_native_gram_point_density`, matches the direct native path, the
freeze-per-step path and the look-ahead cache path alike, and all three measure
about eighty times cheaper than the constant claims.

## Method

`benchmarks/studies/gram_memory_footprint` runs **one subprocess per (density
path, constellation size) point**. A shared process would carry the previous
point's heap into the next point's baseline, which is precisely the quantity
being measured.

Every point runs `julia --threads=1`, because that is how a process-pool worker
is started and the worker is the thing the routing estimate prices. The worker
first solves a one-spacecraft mission to pay the solver's compilation and the
native GRAM image's first-call allocations -- neither is a per-spacecraft cost
-- then reads its resident set, builds the N-spacecraft configuration, solves,
and reports `Sys.maxrss` along with live RSS before and after a full collection
and the `VmRSS` line from `/proc/self/status`.

The constellation geometry, effector stack, tolerances and save settings are the
same as `benchmarks/studies/paper_scenarios/scenario_worker.jl`, so a row here
is directly comparable with an S1 or S2 row there.

Paths measured, selected exactly as `paper_scenarios/common.jl` selects them:

| path | what it is | switches |
|---|---|---|
| `none` | no atmosphere: the control that isolates what GRAM itself costs | -- |
| `gram_point` | direct native point density, the path the predicate is named for | `FREEZE_PER_STEP=0`, `VACUUM_GRAM_CACHE=0` |
| `gram_freeze_step` | density frozen per accepted step (S2's `serial_standard`) | `FREEZE_PER_STEP=1` |
| `gram_lookahead` | vacuum-predicted look-ahead cache (S2's `threads_lookahead`) | `VACUUM_GRAM_CACHE=1`, horizon and deviation loosened past the mission |
| `gram_surrogate` | offline surrogate table | `GRAMAtmosphereModelSurrogate` |

`isolate_state` is left at its production default, so each point pays
`run_simulation`'s own deep copy of the configuration.

Two arcs were run: a short one (60 s, degree-20 harmonics) across
N = 16, 64, 256, 1024, 4096, and an anchor arc matching the S1/S2 convention
(1800 s, degree-50 harmonics) across N = 256, 1024, 4096.

`rss_base_mb` is recorded per point but is **not** subtracted from the peak: it
is read right after a collection whose page release the kernel may or may not
have completed, and it was seen to move by 260 MB between two points of the same
path. The peak is the stable reading, and the intercept of its fit is the honest
baseline.

## Cross-machine check

Before trusting any of this, the method was checked against data that already
existed in the repository. `paper_scenarios/results/space-falcon-lab-TRX50-AERO-D/s1_constellation_scaling_l50_none.csv`
records 1712.3 MB for 256 spacecraft, no atmosphere, degree-50 harmonics, an
1800 s arc, serial. This study's matching point on this workstation measures
**1713.4 MB** -- a 1.1 MB difference across two different machines and two
independently written harnesses.

## What the footprint is

Fitted peak resident set against spacecraft count, ordinary least squares, one
solve per point (`fit_footprint.jl` reproduces this table from the committed
CSVs):

| path | arc | points | intercept | slope | R² |
|---|---|---|---|---|---|
| `none` | 60 s, L20 | 5 | 1699 MB | 0.0063 MB/sat | 0.81 |
| `gram_point` | 60 s, L20 | 5 | 1909 MB | 0.0099 MB/sat | 0.56 |
| `gram_freeze_step` | 60 s, L20 | 5 | 1882 MB | 0.0133 MB/sat | 0.57 |
| `gram_lookahead` | 60 s, L20 | 5 | 1903 MB | 0.0112 MB/sat | 0.36 |
| `gram_surrogate` | 60 s, L20 | 5 | 2210 MB | 0.0052 MB/sat | 0.23 |
| `none` | 1800 s, L50 | 3 | 1715 MB | 0.0092 MB/sat | 0.94 |
| `gram_point` | 1800 s, L50 | 3 | 1907 MB | -0.0023 MB/sat | 0.02 |
| `gram_freeze_step` | 1800 s, L50 | 3 | 1914 MB | 0.0067 MB/sat | 0.78 |
| `gram_lookahead` | 1800 s, L50 | 3 | 1913 MB | 0.0067 MB/sat | 0.45 |
| `gram_surrogate` | 1800 s, L50 | 3 | 2185 MB | 0.0212 MB/sat | 0.85 |

`gram_surrogate` is not one of the three native calling methods the routing
predicate matches (see "Not covered" below for why it is charged nothing
regardless); it is carried in this table because it was measured, not because
it changes the fit the constants come from.

Three things follow, about the three native paths.

**The three native paths are indistinguishable.** Their intercepts agree to
about 30 MB and their slopes are all within the noise of each other and of zero.
Nothing in the measurement separates direct point density from freeze-per-step
from the look-ahead cache, at any size from 16 to 4096 spacecraft.

**What GRAM costs is a fixed image, not a per-spacecraft term.** Loading the
native GRAM library and its data adds about 200 MB to the process -- the gap
between the `none` intercept and the GRAM intercepts -- and that gap does not
grow with the constellation. This is what the code says it should be: the
extension builds one `GRAMAtmosphereModel` wrapping one native core behind one
`ReentrantLock` per process, not one per spacecraft.

**The per-spacecraft term that does exist is the solver's, not GRAM's.** Over a
256-fold change in constellation size the peak moves by tens of megabytes, and
it moves by the same amount with and without an atmosphere model. The quantity
that scales with spacecraft count is retained trajectory state, and a run that
solves repeatedly in one process accumulates more of it: the S1 rows, which run
four solves per point, reach 5410 MB at 4096 spacecraft where a single solve
here measures 1752 MB.

## The slope the routing constant is set from

A pool worker runs many samples in its lifetime, so it is the repeated-solve
figure, not the first-solve peak, that a worker's memory estimate has to cover.
Measured slopes for four solves of the same size in one process:

| source | workload | slope |
|---|---|---|
| `s1_constellation_scaling_l50_none.csv` (TRX50) | no atmosphere, L50, 1800 s | 0.963 MB/sat (256 → 4096) |
| `s2_gram_atmosphere_modes.csv` (TRX50), `serial_standard` | native GRAM, freeze-per-step, 600 s | 1.140 MB/sat (256 → 1024) |
| `s2_gram_atmosphere_modes.csv` (TRX50), `threads_standard` | native GRAM, freeze-per-step, 64 threads | 1.097 MB/sat (256 → 1024) |
| `s2_gram_atmosphere_modes.csv` (TRX50), `threads_lookahead` | native GRAM, look-ahead cache | 1.015 MB/sat (256 → 1024) |

The largest slope measured anywhere, across two machines, two harnesses and
four workloads, is **1.14 MB per spacecraft**.

## The constants

| constant | value | status |
|---|---|---|
| `_WORKER_MEMORY_FLOOR_BYTES` | 2 GB | DERIVED: 1903.8 MB (1.86 GiB) measured at `gram_point`, N = 16, rounded up |
| `_GRAM_SAT_MEMORY_BYTES_BY_PATH.point` | 2 MB/sat | DERIVED: 1.14 MB/sat largest measured, carried up with a safety factor |
| `.freeze_per_step` | 2 MB/sat | same; the measurement does not separate it from `point` |
| `.lookahead` | 2 MB/sat | same |
| `.surrogate` | 0 | SOURCED: the offline table holds no per-spacecraft native state, and the router already excludes it via `gram_surrogate_enabled` |

The per-path table is kept even though the three native entries are equal: it
records that they were measured separately and found equal, and a future
re-measurement that does separate them changes one line rather than the shape of
the model. `SPACEAGORA_GRAM_SAT_MEMORY_MB` still overrides every path at once,
as before, so an operator who had pinned it keeps that behavior.

`OuterRouteFeatures` was **not** given a field to say which native calling
method a workload is on. `density_family` alone cannot distinguish them --
`gram_point`, `gram_freeze_step` and `gram_lookahead` are all selected through
`SPACEAGORA_*` environment switches on top of the same `gram_point` density
family, so `_is_native_gram_point_density` already matches all three -- but
since the measured per-path constants are equal, `native_gram_worker_extra_bytes`
keeps its default `path = :point` at its one call site in
`effective_process_workers` and nothing downstream needs to know which method is
running. If a future re-measurement does separate the paths, that call site is
where a features field (and a way to fill it from the env switches in
`campaign_route_features`) would need to be added.

Raising the worker floor from 1.5 GB to 2 GB applies to **every** workload, not
only GRAM ones. It makes the estimate more conservative, and it is the direction
the measurement points: a bare `--threads=1` worker's measured peak ranges
1697-1752 MB (1.66-1.71 GiB) with no atmosphere and 1866-1946 MB (1.82-1.90 GiB)
with native GRAM, across every size and both arcs, so the old 1.5 GB floor
started workers that did not fit.

## What changes in routing

Machine facts are supplied through the existing override hooks, so the table
answers for the profile rather than for the host that computed it:
`SPACEAGORA_CORE_BUDGET`, `SPACEAGORA_MEMORY_BUDGET_GB`,
`SPACEAGORA_MEMORY_AVAILABLE_GB` and `SPACEAGORA_PERF_WORKER_MEMORY_GB`. The
budget is the machine's memory less the 10 % reserve `machine_topology` applies;
the headroom is that budget less the coordinator's own 1.5 GB image. "Before" is
the shipped-before 90 MB reproduced through `SPACEAGORA_GRAM_SAT_MEMORY_MB=90`;
"after" is the measured defaults. The coordinator's thread count is held at 8 in
both columns, because the local-slot term `min(nthreads - 1, usable_cores - W)`
is unchanged by this work.

Reproduce with:

```
julia --project=. --threads=8 benchmarks/studies/gram_memory_footprint/route_fit_table.jl
```

| profile | N | charge | extra/worker | workers | route fits | mixed slots | capacity |
|---|---|---|---|---|---|---|---|
| workstation, 60 GB / 12 cores | 256 | before (90 MB/sat) | 22.50 GB | 2 | yes | 0 | 2 |
| workstation, 60 GB / 12 cores | 256 | after (measured) | 0.50 GB | 12 | yes | 0 | 12 |
| workstation, 60 GB / 12 cores | 1024 | before (90 MB/sat) | 90.00 GB | 0 | **NO** | 0 | 0 |
| workstation, 60 GB / 12 cores | 1024 | after (measured) | 2.00 GB | 12 | yes | 0 | 12 |
| workstation, 60 GB / 12 cores | 4096 | before (90 MB/sat) | 360.00 GB | 0 | **NO** | 0 | 0 |
| workstation, 60 GB / 12 cores | 4096 | after (measured) | 8.00 GB | 5 | yes | 0 | 5 |
| TRX50-class, 250 GB / 64 cores | 256 | before (90 MB/sat) | 22.50 GB | 9 | yes | 0 | 9 |
| TRX50-class, 250 GB / 64 cores | 256 | after (measured) | 0.50 GB | 64 | yes | 0 | 64 |
| TRX50-class, 250 GB / 64 cores | 1024 | before (90 MB/sat) | 90.00 GB | 2 | yes | 0 | 2 |
| TRX50-class, 250 GB / 64 cores | 1024 | after (measured) | 2.00 GB | 63 | yes | 1 | 64 |
| TRX50-class, 250 GB / 64 cores | 4096 | before (90 MB/sat) | 360.00 GB | 0 | **NO** | 0 | 0 |
| TRX50-class, 250 GB / 64 cores | 4096 | after (measured) | 8.00 GB | 23 | yes | 0 | 23 |

`extra/worker` is `native_gram_worker_extra_bytes(N)`, the per-spacecraft term
alone (not the full worker estimate, which also carries the pinned 1.5 GB
package image). "route fits" is `effective_process_workers >= 2`.

The old constant withdraws the process route outright at 1024 spacecraft on
this workstation and at 4096 on both profiles; the measured constants leave it
affordable everywhere in this table, including a wider pool than the old
constant left even where the old one still technically fit (256 spacecraft on
the workstation: 2 workers before, 12 after).

## Re-running this

```
julia --project=. benchmarks/studies/gram_memory_footprint/run_footprint.jl
julia --project=. benchmarks/studies/gram_memory_footprint/fit_footprint.jl \
    benchmarks/studies/gram_memory_footprint/results/gram_memory_footprint_<host>.csv
```

`GMF_PATHS`, `GMF_SIZES`, `GMF_MISSION_S`, `GMF_GRAVITY`, `GMF_ISOLATE`,
`GMF_SOLVE_REPEATS`, `GMF_TIMEOUT_S`, `GMF_MIN_FREE_GB` and `GMF_SUFFIX` are
documented at the top of the controller. Results land under `results/` named by
host, and the committed CSVs from this workstation are what the constants above
are read off.

## Not covered

- **The offline surrogate path was measured, but is not what the constants above
  are read from.** `GRAMAtmosphereModelSurrogate` carries a fixed table and
  `with_density_model_epoch` refuses to re-epoch it, so the worker has to build
  the surrogate model at the same epoch the run configuration uses rather than
  letting the engine realign it; `gmf_initial_time()` does that for both. With
  that in place the path solved successfully at every size on both arcs (the
  `gram_surrogate` rows above). Its base image is the largest measured, roughly
  250-340 MB above the native paths' (comparing `rss_base_mb` at matching
  points), because it holds both the offline table and a native fallback; its
  fitted slope (0.005-0.021 MB/sat) is in the same small range as every other
  path's. Its routing charge is zero regardless of
  what this measured, because `_is_native_gram_point_density` already excludes
  a surrogate workload through `gram_surrogate_enabled` -- the surrogate makes
  no native call per spacecraft, so there is nothing to price per spacecraft --
  and this work does not change that exclusion.
- **The GRAM density-service pool and the isolated GRAM pool are not priced.**
  Both spawn their own processes with their own native GRAM images. Those are
  per-process costs, not per-spacecraft ones, and neither is enabled by any
  shipped profile.
- **Inner thread width is not priced.** Every point here ran `--threads=1`.
- **Wall times in the CSVs are not benchmark results.** The measurements were
  taken while another agent was profiling on the same workstation, by the
  user's explicit instruction; the `solve_s` column is recorded so a re-run can
  be recognized as comparable, not as a performance claim.
