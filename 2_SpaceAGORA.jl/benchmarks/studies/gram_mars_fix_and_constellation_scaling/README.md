# GRAM Mars Fix and Constellation Size Scaling

Two related pieces of work from the same investigation, kept together because
the scaling study depends on the Mars fix to be able to run Mars scenarios at
all.

## The Mars GRAM fix

`libGRAM.dylib` (the vendored native atmosphere library, `data/GRAMSuite.jl`)
statically links its own private copy of CSPICE, completely isolated from the
CSPICE instance SpaceAGORA's own `SPICE.jl` uses. GRAM's internal
`SpiceLoader` populates that isolated instance from its own default per-body
SPK kernels, which are Sun-inclusive only for Earth/Venus (a full planetary
ephemeris); Mars's default (`mar097_GRAM.bsp`) is Mars-system-only. As a
result, every Mars density query used to fail immediately with:

```
GRAM update failed (code=1): Error: A Spice error occurred.
       Error in longitude of the sun calculation.
```

This is not a missing-kernel-data problem -- the same `lspcn_c`/`ltime_c`
calls succeed when replayed against the identical kernel set via `SPICE.jl`.
It is that GRAM's own isolated CSPICE instance never gets fed the data.

The fix (`data/GRAMSuite.jl/src/GRAMSuite.jl`'s `_GRAM_EPHEMERIS_STATE_FN`
hook, populated by `ext/SpaceAGORAGRAMSuiteExt.jl`'s
`_gram_spice_ephemeris_state`) computes GRAM's required ephemeris state
(solar longitude, one-way light time, subsolar point, local solar time, solar
zenith angle, seconds-per-sol) using SpaceAGORA's own working `SPICE.jl`
bindings and feeds it to GRAM directly via `set_ephemeris_state!`, bypassing
GRAM's own broken internal SPICE calls entirely. The hook defaults to
`nothing` (GRAM's native behavior, unaffected) unless the extension is
loaded, and is registered for every body GRAM supports, not just Mars -- the
same isolated-CSPICE limitation plausibly affects Jupiter/Saturn/Uranus/
Neptune/Titan too (same satellites-only default kernel pattern), though only
Earth and Mars have been exercised end-to-end.

`mars_150km_gram_scaling.jl` is both a working example and a regression
check: if it starts failing again with the "longitude of the sun" error
above, the ephemeris bypass has broken.

```bash
julia --project=. --threads=4 benchmarks/studies/gram_mars_fix_and_constellation_scaling/mars_150km_gram_scaling.jl parallel
```

## Constellation size scaling

`leo_constellation_size_scaling.jl` sweeps LEO constellation size (1, 2, 4,
..., 1024 satellites -- powers of 2) **and** thread count (default 1, 2, 4,
8, 16, 32, 64 -- capped at 64) at a fixed 600s mission, comparing:

- **standard**: real/native GRAM, no vacuum-predicted lookahead cache
  (`SPACEAGORA_DENSITY_FREEZE_PER_STEP=1`, density calls serialized)
- **surrogate**: the GRAM offline surrogate (`GRAMAtmosphereModelSurrogate`)
  -- precomputed interpolation grid, no native calls, no lock
- **no_gram**: `NoAtmosphereModel()` with the aerodynamic effector dropped
  entirely -- baseline gravity-only dynamics cost, isolating
  constellation-size scaling with zero atmosphere-model overhead

All 33 `(N_SATS, mode)` points for one thread count run back-to-back in ONE
process (`leo_constellation_size_scaling_point.jl`'s `run_scaling_point`) --
`ODEParams`/`SharedBuffers`'s `N_sats` is a runtime field
(`src/core/types/runtime_types.jl`), so every `N_sats` value for a given mode
shares one compiled specialization and only the 3 modes still differ in
type. Thread count, however, is fixed at Julia process startup, so each
thread count in the ladder still runs as its own subprocess (this script
re-invoking itself with `--threads=<N>`); a crashed thread-count subprocess
only drops that thread count from the merged results rather than aborting
the rest of the ladder.

Each point also records its own resource footprint -- peak RSS and mean/peak
CPU%, sampled from that process only via `ps -p <pid>` polling
(`resource_monitor.jl`) so other load on the machine never pollutes the
numbers.

```bash
julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_constellation_size_scaling.jl
```

Override the thread ladder (comma-separated, values must be <= 64):

```bash
SPACEAGORA_SCALING_THREADS=1,8,64 julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_constellation_size_scaling.jl
```

Produces `leo_constellation_size_scaling_with_gram.png` (standard vs.
surrogate, both GRAM-derived, one small-multiple panel per thread count) and
`leo_constellation_size_scaling_without_gram.png` (the no-GRAM baseline, one
line per thread count), each with a dashed O(N) linear-scaling reference line
anchored at N=1 for visually judging super-/sub-linear scaling, plus
`leo_constellation_size_scaling_resource_ram.png` and
`leo_constellation_size_scaling_resource_cpu.png` (peak RSS / mean CPU%, one
panel per mode, one line per thread count).

**Key finding**: the lookahead/vacuum-predicted-cache GRAM mode -- used
elsewhere in this repo's benchmarks as the default "fast path" -- is actually
*slower* than the "standard" (no cache) mode at every constellation size
tested here (1 to 1024 satellites), because a 600s mission is too short to
amortize the cache's upfront drag-free-trajectory + spline-fit cost. The
lookahead cache's benefit is framed around much longer missions elsewhere in
this repo (e.g. a 1-hour mission) -- check actual mission duration before
assuming the lookahead cache is the fast option.

Companion single-point scripts used to establish this before building the
full sweep:

- `earth_standard_gram_scaling.jl` -- standard (no-lookahead) real GRAM, one
  constellation size at a time (`LEO_SCALING_N_SATS` env override)
- `earth_surrogate_gram_scaling.jl` -- GRAM surrogate, same pattern

## Thread allocation, GRAM concurrency, and the `:process` outer-parallel route

Full detail in `THREAD_ALLOCATION_AND_GRAM_CONCURRENCY_HANDOFF.md` (Findings
1-7) and `PARALLELIZATION_CURRENT_STATE.md` (current-state snapshot); summary
here.

A later investigation found that concurrent `GRAMAtmosphereModel`
*construction* (not sampling -- sampling was already correctly locked) is not
thread-safe, crashing any outer-parallel campaign
(`run_monte_carlo`/`run_constellation_ensemble`) that builds more than one
native-GRAM model at a time on different threads. Fixed by extending
`_with_gram_lock` to also guard construction. Verified: `mc_multisat_thread_allocation.jl`'s
`mode="standard"` (native GRAM) points, previously excluded from that study
entirely because they crashed, now complete 8/8 successfully across the full
outer x inner thread grid.

That fix unblocked a second piece of work: `SpaceAGORA.ParallelProcess`, a
shared `Distributed`-backed worker pool wired into `run_monte_carlo`/
`run_constellation_ensemble`'s adaptive `threads=:auto` routing as a real
`:process` outer-parallel route (separate OS processes, no shared native
GRAM/CSPICE state to lock in the first place -- see
`project_gram_cspice_symbol_collision` memory for why GRAM's *own* vendored
CSPICE isn't actually isolated from SpaceAGORA's `SPICE.jl` at the
threaded-in-process level, which is what made the lock necessary at all).

`leo_ensemble_outer_backend_scaling.jl`/`_worker.jl` compares `:threads` vs.
`:process` at a fixed N=64-satellite native-GRAM ensemble:

```bash
julia --project=. benchmarks/studies/gram_mars_fix_and_constellation_scaling/leo_ensemble_outer_backend_scaling.jl
```

| outer_workers | threads | process |
|---|---|---|
| 1 | 1.00x | 1.00x |
| 2 | 0.98x | **1.20x** |
| 4 | 0.98x | **1.44x** |
| 8 | 0.93x | 0.75x (regresses -- oversubscription past this 11-core machine's headroom) |

`:threads` stays pinned near 1x regardless of worker count -- GRAM's
process-wide lock caps it, exactly as intended (safe, not fast). `:process`
breaks past that ceiling up to 4 workers on this machine. Results:
`leo_ensemble_outer_backend_scaling_summary.csv` /
`leo_ensemble_outer_backend_scaling.png`.

`mc_multisat_thread_allocation.jl` was extended the same way at the
Monte-Carlo-campaign level (8 samples, 4 sats/sample), adding `mode="standard"`
and an `outer_backend` axis. Its own full-sweep CSV's `process`-backend rows
at `outer_workers>=2` are **not reliable** -- a benchmark-harness precompile-cache
issue hit mid-session made most of those points record as `TIMEOUT` rather than
real measurements. That surfaced a bigger, now-fixed finding: this workload's
`:process` backend originally showed no scaling with `outer_workers` at all,
not because of per-sample model reconstruction (that hypothesis was tested
directly and refuted -- GRAM construction is sub-millisecond, no cross-process
contention), but because a worker's first-ever full `run_simulation` call pays
a large one-time JIT/specialization cost (measured directly at ~58-70s cold vs.
a few seconds warm, same workload) that `SpaceAGORA.ParallelProcess`'s generic
bootstrap warmup didn't cover. **Fixed**: `ensure_process_workers!` gained an
optional `warmup_fn` that campaign call sites populate with their own real
per-sample closure, so a newly-added worker's JIT cost is paid once during
pool growth instead of inside the caller's first real, measured dispatch --
see handoff doc Finding 7 for the fix and its (surprisingly involved)
validation story.

The `no_gram`/`surrogate`/`standard` nested outer+inner `threads`-backend
points in this same sweep, previously wildly non-deterministic (some ranging
2.7s-300s+ across identical repeated runs) or hanging for 16+ minutes, are now
fixed and bounded -- see the next section.

## Nested outer+inner thread contention: root-caused and fixed

The nested-parallelism pathology found while building the study above (severe,
non-deterministic slowdowns whenever a campaign's outer split, e.g.
`run_monte_carlo`'s `Threads.@spawn` workers, was combined with per-sample
inner parallelism) was originally characterized as a fundamental Julia
scheduler limitation, out of scope to fix. It wasn't: `outer_parallel_active()`
is meant to be a uniform kill-switch against exactly this kind of nesting, and
correctly gates most call sites in this codebase (density/control/thermal
callbacks, dynamic-effector threading) -- but three multi-worker routing
branches inside `_rhs_execution_plan` (`src/simulation/engine/setup.jl`),
including the harmonics-batch/flat-constellation route this study's own env
config enables by default, never checked it. Since that route fires on every
RHS/ODE step (not once per sample), an outer worker already blocked waiting on
its own children would repeatedly spawn a further nested batch throughout the
whole integration -- the actual mechanism behind the observed livelock-like
behavior. Fixed by adding the same guard already used elsewhere in the file to
all three branches (new `SPACEAGORA_HARMONICS_BATCH_ALLOW_WITH_OUTER` env var,
default off). The two worst previously-observed repro points now run
consistently at ~0.04-0.07s across repeated trials, matching the outer-only
baseline exactly -- full detail and the before/after numbers in the handoff
doc's Finding 6.

Separately, `SPACEAGORA_DENSITY_CALLBACK_PARALLEL=auto`'s 16-thread floor
(Finding 1) applied uniformly regardless of density model, needlessly gating
the lock-free GRAM surrogate (native GRAM's process-wide lock is the actual
reason the floor exists; the surrogate has none). Fixed with a separate,
lower-floor gate category for lock-free density models; native GRAM's
existing threshold is unchanged.
