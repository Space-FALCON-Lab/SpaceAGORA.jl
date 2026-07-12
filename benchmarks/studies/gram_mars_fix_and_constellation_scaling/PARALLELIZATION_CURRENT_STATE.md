# Parallelization: Current State of the Repository

Snapshot of what's implemented, proven, and still open on the
`parllelization_paper` branch, focused specifically on the parallelization
work (GRAM/CSPICE thread-safety and the new `:process` outer-parallel route).
Written as a single up-to-date reference instead of reconstructing the story
from commit history and the two other docs in this folder
(`README.md`'s "Thread allocation..." section and
`THREAD_ALLOCATION_AND_GRAM_CONCURRENCY_HANDOFF.md`'s Findings 1-7), which
remain the detailed, dated record this document summarizes.

**Nothing described below as "uncommitted" has been pushed.** Only the GRAM
constructor lock and CSPICE/SPICE lock unification (first bullet below) are
committed and pushed.

## 1. What's fixed

### 1a. GRAM/CSPICE lock unification (committed: submodule `b66e1b8`, parent `444a23c5`)

`libGRAM.dylib` statically links its own copy of CSPICE. It isn't actually
isolated from SpaceAGORA's own `SPICE.jl` CSPICE instance at the OS
symbol-export level (`nm -gU` shows shared symbol names) -- two separate
locks (`GRAM_LOCK`, `SPICE_LOCK`) were silently protecting against two
threads racing on what is, underneath, the *same* global CSPICE state.
`src/simulation/runtime_services.jl` now defines `const GRAM_LOCK = SPICE_LOCK`
so every GRAM density call and every SPICE ephemeris call serialize against
each other, not just against their own kind.

### 1b. GRAM model construction thread-safety (uncommitted, in the `data/GRAMSuite.jl` submodule)

Density *sampling* (`getDensity`/`density_state`/`point_density_state`) was
already correctly locked. Model *construction*
(`GRAMAtmosphereModel`'s keyword constructor, which calls
`gram.initialize!`/`gram.create_atmosphere`/`gram.set_start_time!`/etc.) was
not -- two threads each building a fresh GRAM model concurrently crashed with
`SPICE(BADSUBSCRIPT)` or a `libc++abi` abort, 100% reproducible. This is what
caused every native-GRAM Monte Carlo/ensemble campaign to crash outright at
`outer_workers >= 2` before this work. Fixed by wrapping the constructor body
in the same `_with_gram_lock`/`GRAM_LOCK` used for sampling.

### 1c. World-age / `invokelatest` fixes (uncommitted, same submodule)

A second, unrelated bug surfaced only by the new `:process` route (section 2):
GRAM's native wrapper module is dynamically `Base.include`d on first-ever
model construction in a process, defining new types/methods at a *later*
Julia world age. A Distributed campaign closure dispatched to a fresh worker
-- whose call frame was compiled before that dynamic load -- threw
`MethodError(..., world_age=...)` on unwrapped calls into the newly-loaded
module, even ones textually after the load. Two fixes:

- `_gram_apply_user_ephemeris_state!`'s `EphemerisStateC(state...)` call now
  goes through `Base.invokelatest` (fetching the type via an
  `invokelatest`'d `getfield` was already correct; *calling* the fetched
  type object needed its own wrap too).
- `_with_gram_lock`'s construction-only branch (`lock_obj === nothing`, the
  branch the constructor uses) now runs the whole locked body via
  `Base.invokelatest(f)` instead of `f()`, so every `gram.xxx(...)` call
  inside the constructor -- there are about a dozen -- sees the current
  world age without each one needing its own wrap. Deliberately **not**
  applied to the sampling branch (hot per-density-query path); construction
  always precedes sampling in a given process, so the dynamic load has
  already settled by the time sampling runs, and `invokelatest` has a real
  per-call dispatch cost not worth paying there.

Custom `Serialization.serialize`/`deserialize` methods were also added for
both `GRAMAtmosphereModel` and `GRAMAtmosphereModelSurrogate` (serialize the
plain-data constructor arguments, reconstruct via the now-locked constructor
on the receiving side) -- needed for a `SimulationConfiguration` carrying a
GRAM density model to cross a `Distributed` process boundary at all.

### 1d. Density-callback 16-thread auto-gate made density-model-aware (uncommitted)

Unrelated to the GRAM/CSPICE work above, but the same "parallelization
correctness/tuning" family. `SPACEAGORA_DENSITY_CALLBACK_PARALLEL=auto`'s
16-thread floor (`auto_thread_min_budget(:density_callback)`,
`src/parallel/policy/env_config.jl`) applied uniformly regardless of density
model, gating the lock-free GRAM surrogate for no reason (native GRAM's
process-wide lock is the actual reason the floor exists; the surrogate has no
such lock). Fixed with a new `:density_callback_lockfree` source category
(falls through to the general 4-thread default) that
`_density_callback_thread_decision` (`density_callbacks/config.jl`) now
selects for non-`GRAMAtmosphereModel` density models, leaving native GRAM's
existing 16-thread floor completely unchanged. Directly verified: a surrogate
model at `--threads=4` under `auto` mode now gives `use_threads=true,
allotment=4` (previously would have been `false`, since 4 < the old
unconditional 16-thread floor).

## 2. What's new: `SpaceAGORA.ParallelProcess`

A shared, auto-bootstrapping `Distributed` worker pool
(`src/parallel/process/`, uncommitted), wired into `run_monte_carlo` and
`run_constellation_ensemble`'s adaptive `threads=:auto` routing as a real
`:process` outer-parallel route -- separate OS processes, so there's no
shared native GRAM/CSPICE global state to contend over in the first place.
Previously `:process` was a route name that silently downgraded to
`:threads`/`:none`; it now actually dispatches.

- `campaign_process_pool()` / `ensure_process_workers!(pool, n)`: grow a
  process-global pool to at least `n` workers, bootstrapping each new one
  with `using SpaceAGORA`, best-effort `import GRAMSuite`, default SPICE
  kernels (a fresh `Earth("")` construction), and a GRAM warmup (construction
  + one density query, to force the world-age settling from section 1c
  before any real campaign work is dispatched).
- Each worker runs `--threads=1` -- process workers don't share the
  coordinator's thread pool, so inner per-worker thread budget isn't
  currently configurable (flagged as a known limitation, not a bug).
- `OuterRouteTuning` gained a `process_max_workers` field (default
  `Sys.CPU_THREADS`) so a caller can pin the process-route worker count
  directly instead of relying on the adaptive bandit's exploration.
- **Done**: `export`/`@doc` lines for
  `ProcessPool`/`campaign_process_pool`/`ensure_process_workers!`/
  `shutdown_process_pool!` added to `src/SpaceAGORA.jl` (mirroring the
  existing `ParallelProfiles` pattern), `SpaceAGORA.ParallelProcess` added to
  `docs/make.jl`'s `checkdocs_ignored_modules`, and its 4 symbols added to
  `docs/public_api_symbols.jl`'s generated API page (a second, separate gap
  from `checkdocs_ignored_modules` -- the auto-generated `public_api.md` page
  has its own manually-maintained symbol list that also needed the new
  entries, or `makedocs` fails with "N docstrings not included in the
  manual" even with `checkdocs_ignored_modules` set correctly). `julia
  --project=. docs/make.jl` now completes clean end to end, including
  linkcheck.

## 3. What's proven, with numbers

### Constellation ensembles (`leo_ensemble_outer_backend_scaling.jl`, new benchmark)

Fixed 64-satellite native-GRAM ensemble, sweeping `outer_workers ∈ {1,2,4,8}`
× `{threads, process}`:

| outer_workers | threads | process |
|---|---|---|
| 1 | 1.00x | 1.00x |
| 2 | 0.98x | **1.20x** |
| 4 | 0.98x | **1.44x** |
| 8 | 0.93x | 0.75x (regresses -- likely this 11-core machine's own ceiling: 8 process workers + 1 coordinator = 9 concurrent OS processes) |

`:threads` stays pinned near 1x regardless of worker count -- the lock from
section 1a/1b is doing its job (safe, not fast). `:process` clearly breaks
past that ceiling up to 4 workers.

### Monte Carlo campaigns (`mc_multisat_thread_allocation.jl`, extended)

Added a `mode="standard"` (native GRAM) option -- previously excluded from
this study entirely because it crashed -- and an `outer_backend` axis.

- **Regression check passed**: every `threads`-backend combo for
  `mode="standard"`, across the full outer×inner grid up to 8×1, now
  completes 8/8 successfully. No crash, at any combo. This is the direct,
  concrete proof that sections 1a-1c fixed the crash.
- **Process-backend cold-start JIT cost: root-caused and fixed.** Direct,
  uncapped single-point measurements originally showed every `outer_workers`
  value clustering around the serial baseline (~2.2-2.45s), no scaling with
  worker count. The original hypothesis (per-sample `GRAMAtmosphereModel`/MERRA2
  reconstruction serializing across OS processes) was tested directly and
  **refuted**: construction alone is sub-millisecond whether run serially or
  concurrently across 8 processes, no contention at all. The real cause: a
  worker's **first-ever** full `run_simulation` call pays a large one-time
  JIT/specialization cost -- directly measured at ~58-70s for a cold worker vs.
  a few seconds for an already-warm one, on the exact same workload.
  `SpaceAGORA.ParallelProcess`'s per-worker bootstrap warmup only exercised
  GRAM's native wrapper (construction + one density query), not the much
  larger JIT surface of a full `run_simulation` call. **Fixed** with a new
  `warmup_fn` parameter on `ensure_process_workers!` (`worker_pool.jl`):
  campaign call sites now pass their own real per-sample closure as the
  warmup, so a newly-added worker's JIT cost is paid once during pool growth
  instead of landing inside the caller's first real, measured dispatch.
  Directly verified: all 4 newly-added workers' warmup calls succeeded
  (~57s each, matching the known cold cost) during pool growth, after which
  all 4 real, timed dispatch samples came back uniformly fast (~3.8s each, no
  asymmetry). Full test suite re-run clean. See handoff doc Finding 7 for the
  fix and, notably, the debugging story behind validating it (three stacked
  test-harness mistakes -- a benchmark-local closure dependency, a bare-name
  resolution gap, and a wrong `density_family` string -- each of which made
  early verification attempts silently bypass the process pool entirely and
  run everything serially on the coordinator instead, producing a misleading
  "one slow sample, rest fast" pattern that had nothing to do with the actual
  fix).
  
  Separately (not fixed, and not fixable in code): the *originally reported*
  flat process-backend numbers are, independently, also explained by CPU
  oversubscription on this 11-core dev machine once several workers are doing
  real concurrent compute -- see handoff doc Finding 6's pass-1/pass-2
  instrumented data. That's a hardware-capacity fact, distinct from the JIT
  cost this fix addresses.

### Phase 3 (nested outer+inner thread contention) -- root-caused and fixed

Was characterized as a "genuine scheduler-level, out of scope" livelock in the
prior write-up of this document. A follow-up session found it is not that --
`outer_parallel_active()` is supposed to be a uniform kill-switch against
nested outer+inner parallelism, but three multi-worker routing branches inside
`_rhs_execution_plan` (`src/simulation/engine/setup.jl`) -- including the
harmonics-batch/flat-constellation route this study's own env config enables
by default -- never checked it, unlike every other call site in the codebase.
Fixed by adding the same guard already used elsewhere in the same file to all
three branches. Directly re-verified: the two worst repro points (previously
2.7s-300s+ non-deterministic, and once left running 16+ minutes before being
force-killed) now run consistently at ~0.04-0.07s across repeated trials --
matching the outer-only baseline exactly, no residual slowdown. Full test
suite re-run clean (4059 pass, 1 expected-broken, matching the pre-fix
baseline). See handoff doc Finding 6 for the full writeup. The existing
"never split the thread budget between outer and inner parallelism" guidance
is now the *automatic default behavior* rather than something a caller has to
manually ensure -- a nested request degrades gracefully to serial instead of
pathologically hanging.

The full 42-point sweep was re-run after this fix landed and confirms it
comprehensively, not just at the two directly-repro'd points -- every one of
the 30 `threads`-backend rows across all three modes (`no_gram`/`surrogate`/
`standard`) now falls in a normal, bounded range with no timeouts:
`no_gram` outer=4x2 (previously timed out on every run, once left running
16+ min) is now 0.050s; `standard` outer=4x2 (previously 146.3s) is now
2.286s, matching that mode's flat ~2.2-2.4s baseline. Updated CSVs in this
folder; full before/after table in the handoff doc's Finding 6.

## 4. Known rough edges in the benchmark tooling itself (not library bugs)

- `mc_multisat_thread_allocation.jl`'s full-sweep CSV has `TIMEOUT` entries
  for every `outer_backend=process, outer_workers>=2` row. This is a
  benchmark-harness artifact: a Julia precompile-cache invalidation
  triggered mid-session (from a `kill -9` on a hung process during
  debugging) made every fresh process-backend point pay a ~190-380s
  one-time recompile cost inside the harness's 300s per-point timeout. The
  direct, uncapped measurements in section 3 are the real numbers for that
  combo; the CSV itself should not be read as-is for those rows.
- The sweep controller (`run_worker` in `mc_multisat_thread_allocation.jl`)
  now polls `process_running` + a wall-clock deadline instead of
  `wait(proc)`, and reaps a timed-out subprocess's own children via `pkill
  -9 -P <pid>` before killing it -- `wait(proc)` was observed to hang
  indefinitely past a killed subprocess's own exit when that subprocess left
  orphaned `SpaceAGORA.ParallelProcess` pool workers holding its
  stdout/stderr pipe open. The worker scripts also now call
  `shutdown_process_pool!` on normal exit for the same reason. Both fixes
  are in the benchmark scripts, not the library.

## 5. Commit / push status

| Change | Location | Status |
|---|---|---|
| GRAM/SPICE lock unification | parent repo + submodule | **committed, pushed** |
| GRAM constructor lock | submodule | uncommitted |
| World-age/`invokelatest` fixes | submodule | uncommitted |
| Serialization for GRAM models | submodule | uncommitted |
| `SpaceAGORA.ParallelProcess` | parent repo (`src/parallel/process/`) | uncommitted |
| Adaptive routing `:process` wiring | parent repo | uncommitted |
| New/extended benchmarks (this folder) | parent repo | uncommitted |
| `export`/`@doc`/`checkdocs`/generated-API-page entries for `ParallelProcess` | parent repo (`src/SpaceAGORA.jl`, `docs/make.jl`, `docs/public_api_symbols.jl`) | uncommitted, **done and verified** (`docs/make.jl` builds clean) |
| Density-callback lock-free auto-gate | parent repo (`src/parallel/policy/env_config.jl`, `src/core/types/runtime_types.jl`, `src/simulation/callbacks/density_callbacks/config.jl`) | uncommitted, **done and verified** |
| Phase 3 fix: `outer_parallel_active` gating for `_rhs_execution_plan`'s 3 multi-worker branches | parent repo (`src/simulation/engine/setup.jl`, `src/core/types/runtime_types.jl`) | uncommitted, **done and verified** (direct repro + full test suite) |

`test/suites/01_contract_and_api_tests.jl` has the user's own unrelated
pre-existing WIP mixed into the same file -- committing the GRAM_LOCK
contract-test fix needs patch-level staging (`git add -p`), not a whole-file
add.

## 6. Where to look next

- `THREAD_ALLOCATION_AND_GRAM_CONCURRENCY_HANDOFF.md` (this folder) --
  full dated findings log, Findings 1-7, with the "Suggested next checks"
  list this document's section 3/4 draws from.
- `README.md` (this folder) -- per-script usage/invocation reference.
- Memory: `project_parallelization_tier1_status`,
  `feedback_julia_distributed_gotchas`,
  `project_gram_concurrent_construction_crash`,
  `project_gram_cspice_symbol_collision`.
