# One-heap contention above 16 threads: what could and couldn't be measured here

A TRX50 cold run showed the pinned-threads route on the P3/P4 Monte Carlo
resource-ladder phases (`benchmarks/studies/paper_parallelization_benchmarks/`)
getting *slower* as thread count rises past 8 (1.55 s -> 2.06 s -> 2.72 s at
8/16/32 threads on P3), while the pinned pool route keeps getting faster
(1.28 -> 0.76 -> 0.43 s). That is the signature of single-heap GC contention,
and it needs 16+ threads to show up.

**This workstation has 12 physical cores / 24 threads and cannot reproduce the
16+-thread regime.** Everything in sections 1-4 is either (a) allocation
attribution and reduction, which is thread-count-independent and fully
measurable here, or (b) an 8-thread collector-settings grid, which is the
closest this box gets to the question and is reported as exactly that -- a
data point below the regime that actually matters, not a stand-in for it.
Sections 5-6 cover a related but separate diagnostic -- process-pool worker
heap growth, unbounded across many samples rather than across threads -- and
the same 16+-thread/8+-worker reach limit applies there too, for the same
reason. Section "What to run on TRX50" at the bottom says what closes the
gap for both.

Every Julia process launched by the scripts in this study runs under a hard
memory cap:

```bash
systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
    julia --project=. ...
```

`collector_grid.sh` wraps it automatically; `alloc_profile.jl`'s and
`dump_state.jl`'s usage headers show the wrapped invocation directly. This is
non-optional for `alloc_profile.jl` in particular: `Profile.Allocs` captures a
full stack per *sampled* allocation, so its own bookkeeping grows with the
sampled allocation count, not the run's byte total, and an unrestrained
sampling rate on a large enough case can consume far more memory than the run
being profiled.

## 1. Allocation attribution (P3, P4)

Tool: `benchmarks/studies/heap_contention/alloc_profile.jl`. Warms up one
solve (JIT), times a second with `@timed` for the headline bytes/accepted-step
number, then profiles a third with `Profile.Allocs` and buckets every sampled
allocation by the `src/` domain that owns its innermost SpaceAGORA.jl stack
frame -- split finely enough to separate `engine/execution.jl` and
`engine/persistence.jl` (the two files this study modifies) from the rest of
`src/simulation/engine/` (`solver_policy.jl`, `dynamics_rhs.jl`,
`rhs_calibration.jl`, ...).

P4 (`montecarlo_heavy_aerobraking`) is a 6 h mission at a 1 s max step, ~365
MiB and several million allocations per sample -- far too large to sample
exhaustively. It runs at `--sample-rate=0.005`, well below Julia's own 0.01
practical ceiling for this kind of profiling; the script warns if a caller
passes a higher rate on anything other than the small P3 case. P3
(`independent_1sat_1hr`, ~15 MiB/sample, ~180k allocations) is small enough to
sample exhaustively at `--sample-rate=1.0`. A uniform sub-sample still gives a
valid *fraction* of allocation per bucket (all this study needs); it does not
give exact byte/count totals, and the script's CSV output says so in its
column names (`estimated_bytes_at_1x` is explicitly an estimate, not a
measurement).

Measured single-threaded (`--threads=1`), space-falcon-1, under the memory
cap, `SPACEAGORA_PARALLEL_POLICY` unset (serial baseline):

| Case | naccept | total alloc | bytes / accepted step |
|---|---|---|---|
| P3 (`independent_1sat_1hr`) | 250 | 14.99 MiB | 62,860 |
| P4 (`montecarlo_heavy_aerobraking`) | 22,487 | 365.79 MiB | 17,057 |

(P4's 365.79 MiB matches the ~364 MB per-sample figure measured on TRX50
closely enough to confirm this is the same case/shape.)

Allocation attribution by bucket (see `results/alloc_p3_independent_1sat_1hr.csv`,
`results/alloc_p4_montecarlo_heavy_aerobraking.csv`):

| Bucket | P3 (sample_rate=1.0) | P4 (sample_rate=0.005, estimate) |
|---|---|---|
| `engine/other` (solver_policy, dynamics_rhs, rhs_calibration) | 51.5% | 56.1% |
| `callbacks` (SavedValues/SaveData/SavingCallback) | 28.2% | 12.1% |
| `external` (OrdinaryDiffEq/DiffEqCallbacks/Base) | 9.0% | 0.3% |
| `parallel` (routing/policy) | 5.6% | 13.7% |
| `environment` (gravity/atmosphere/ephemerides) | 3.4% | -- |
| `dynamics` (RHS terms) | 2.0% | -- |
| `core` | 0.0% | 17.7% |
| **`engine/execution.jl`** | **0.4%** | **0.0%** |
| `engine/persistence.jl` | 0.0% (no samples) | 0.0% (no samples) |

**`execution.jl` and `persistence.jl` are not where P3/P4's per-sample
allocation lives.** `execution.jl` accounts for 0.4% of P3's sampled
allocation and effectively none of P4's; `persistence.jl` (46 lines of thin
re-export wrappers, see its header) accounts for none at either scale. The
dominant sources -- `solver_policy.jl`/`dynamics_rhs.jl` (the solve driver
and RHS dispatch) and the `SavedValues`/`SaveData` callback storage -- belong
to the solver/RHS-routing and callback code respectively (see the
domain-ownership table in `CLAUDE.md`). "What fraction is the saved solution"
is therefore mostly the `callbacks` bucket (SavingCallback building and
growing `saved_values.t`/`saved_values.saveval`), not anything inside
`execution.jl`'s own accumulation loops -- those loops
(`_append_backbone_saved_segment!`, `_append_checkpoint_saved_segment!`) are
not exercised by P3/P4 at all, since neither case checkpoints or uses
`gravity_backbone_split`.

## 2. Reductions made, with dump-and-cmp proof

Given (1), the one change inside `execution.jl`/`persistence.jl` that is both
real and correctness-neutral is: **stop building the debug `ODEProblem` when
nothing reads it.**

`run_simulation` used to construct `prob_debug = ODEProblem(spacecraft_dynamics!,
...)` unconditionally, on every single solve, even though its only reader is
the `SPACEAGORA_DEBUG_INITIAL_DERIVATIVE` NaN-probe branch that almost no run
takes. The construction is now inside the `if` that gates that branch, so the
default (debug flag off) path never builds it. Nothing outside the two debug
branches ever read `prob_debug`, so this changes no observable behavior --
confirmed below by bit-for-bit trajectory comparison, not just by inspection.

No other reduction survived scrutiny at this workload scale: pre-sizing the
checkpoint/backbone accumulation vectors (`checkpoint_saved_times`, etc.) was
considered and **skipped** -- none of the three reference cases below exercise
checkpointing or `gravity_backbone_split`, so there was no way to dump-and-cmp
a change to that code path without adding a fourth case beyond the three
already covered here, and a change that cannot be proven bit-identical is not
shipped.

### Dump-and-cmp evidence

Tool: `benchmarks/studies/heap_contention/dump_state.jl` (same pattern as
`benchmarks/studies/third_body_cost/variants.jl --dump`: raw `Float64.(sol.t)`
followed by every `getdata(u)`, no tolerance). Base commit `adb343566`
(this branch's start point), tip = this branch with the debug-`ODEProblem`
skip applied. Both runs: `--threads=1`, space-falcon-1, `mode=serial`, warmed
once before the timed/dumped solve.

Cases (see `results/dump_cmp_debug_problem_skip.csv`):

| Case | Why this one | Dump size | `cmp` |
|---|---|---|---|
| `independent_1sat_1hr` (P3) | the case backing the P3 Monte Carlo resource-ladder phase in `paper_parallelization_benchmarks/CASES.md` -- one spacecraft per sample, 1 h mission | 18,072 bytes | **identical** |
| `montecarlo_heavy_aerobraking` (P4) | the case backing the P4 phase -- the same resource ladder on a compute-bound 6 h Mars aerobraking sample | 1,619,136 bytes | **identical** |
| `atmo256_vacuum_10min` | a 256-spacecraft, L50-harmonics, no-atmosphere constellation from `parallelization_performance/cases.jl`'s atmosphere-fidelity ladder -- same physics family as `gravity_256sat_l50_vacuum_124000s` but a 10 min mission instead of 124,000 s, chosen so it dumps twice per iteration instead of once per multi-minute run | 2,032,608 bytes | **identical** |

All three: `cmp base_<case>.bin tip_<case>.bin` exits 0. Total allocation moved
by less than the MiB-display's rounding (P3: 15.00 -> 14.99 MiB; P4: 365.81 ->
365.79 MiB; 256-vacuum: 34.76 -> 34.75 MiB) -- consistent with section 1's
finding that this file was never a large contributor, and with the change
itself being "skip one small unconditional allocation on every solve."

## 3. Opt-in state thinning

Two new env-driven settings, read inside
`_save_simulation_results_if_enabled!` in `execution.jl`, apply ONLY to the
*written results output* (CSV / results bundle) -- never to the `ODESolution`
`run_simulation` returns under `return_solution=true`, which is unaffected by
either setting:

- `SPACEAGORA_RESULTS_THIN_STRIDE=k` (default `1`, i.e. off): keep every
  k-th saved state, and always the final saved state even if it does not
  land on a multiple of k (a results file that silently dropped the
  mission-end value would be a correctness hazard for anything reading it).
- `SPACEAGORA_RESULTS_FINAL_ONLY=1` (default off, takes priority over the
  stride setting when both are set): keep only the final saved state.

Both default off, and off is the identity: `_thin_results_segment` returns
the *same* `times`/`data` vector objects it was given, unchanged, so the
default write path is byte-identical to before this feature existed (proven
in `test/unit/simulation/solution_storage_tests.jl`'s
`"_thin_results_segment: pure-function contract"` testset -- `out_t === times`,
not just `==`).

With thinning on, every retained state is an *element of* the untouched
default set -- never recomputed, copied field-by-field, or interpolated -- so
a retained state is always byte-identical to what the default (unthinned) run
would have written for that same saved time. `solution_storage_tests.jl`'s
second testset proves this end to end: it runs a small simulation once with
`SPACEAGORA_RESULTS_THIN_STRIDE` unset (baseline CSV) and once with it set to
`3` (thinned CSV), and asserts every thinned row is a member of the baseline
row set, the thinned CSV is strictly shorter, and the mission-end time is
present in both. A third run with `SPACEAGORA_RESULTS_FINAL_ONLY=1` asserts
exactly one row, and that it matches the baseline's final row.

This is a real but *partial* reduction, stated plainly rather than
overstated: thinning shrinks what gets written and what the caller holds
onto after the solve, but section 1 showed the solver's own per-step
accumulation (which the solve has already paid for by the time this code
runs) is not inside `execution.jl` or `persistence.jl` to begin with, so this
setting does not reduce the solve's own peak allocation -- only the
results-writing tail past it.

## 4. Collector settings (`--gcthreads`, `--heap-size-hint`)

### What was measured, and what wasn't

**Only 8 threads was measured here; 12 threads was not run on this machine.**
Neither was CPU pinning applied (this box has no `--cpu-list` configured for
this study the way TRX50 does), so this is the `outer_threads` route
unpinned, not the pinned-threads configuration TRX50 runs -- the closest
available approximation on this machine.

Tooling: `benchmarks/studies/heap_contention/collector_grid.sh` launches one
Julia subprocess per grid point (under the memory cap above), each running
`collector_grid.jl`'s 5 back-to-back `outer_threads`-mode solves and
reporting the median wall time. Grid: `--gcthreads` in {unset (Julia
default), `2`}; `--heap-size-hint` in {unset, `4G`}. The machine was **not**
quiet during this run (`uptime` load average 6.6-11.2, ~17 other
`julia`-named processes running concurrently on this shared machine
throughout) -- ratios below are reported as measured, not cleaned up, and
should be read with that noise floor in mind.

| Case | gcthreads | heap-size-hint | median wall (s) | ratio vs. both-default |
|---|---|---|---|---|
| P3 | (default) | (default) | 0.042347 | 1.000 |
| P3 | (default) | 4G | 0.042385 | 1.001 |
| P3 | 2 | (default) | 0.038160 | 0.901 |
| P3 | 2 | 4G | 0.038002 | 0.897 |
| P4 | (default) | (default) | 0.813219 | 1.000 |
| P4 | (default) | 4G | 0.844239 | 1.038 |
| P4 | 2 | (default) | 0.903386 | 1.111 |
| P4 | 2 | 4G | 0.831965 | 1.023 |

`--gcthreads=2` moved P3 about 10% faster and P4 about 11% *slower* in the
same run. That is not a consistent winner -- it is two workloads disagreeing
under a noisy shared machine, at a thread count (8) well below the 16+ regime
where this behavior actually appears on TRX50. **No default was changed.**
The mechanism was wired as an opt-in only (`ppc_worker_cmd` in
`benchmarks/studies/parallelization_performance/execution.jl`, gated by
`SPACEAGORA_PPC_WORKER_GCTHREADS` / `SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT`,
both unset by default -- unset means no flag is added to the worker's `julia`
invocation at all, not "pass Julia's own default explicitly"), so a future
run on a quiet machine can exercise the grid without a source change:

```bash
SPACEAGORA_PPC_WORKER_GCTHREADS=2 SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT=4G \
  julia --project=. benchmarks/studies/parallelization_performance.jl \
  --phases=P3,P4 --threads=8,16,32 --process-workers=32
```

Raw data: `results/collector_grid_p3_t8.csv`, `results/collector_grid_p4_t8.csv`.

## 5. Process-pool workers: an unbounded heap, and a default heap-size-hint

A separate diagnostic (P6p, a 4096-sample native-GRAM one-satellite campaign
on TRX50, process route, job `20260923-111926-55080`, memory sampled every
5 s under a 200 GB cgroup cap) found a second, larger instance of the same
kind of problem: at the 8-worker rung, during three consecutive 10 s
samples, system memory went from 107 GB to 181 GB to 211 GB, with thirteen
Julia processes totaling 103/174/203 GB and the largest single process
reaching 13/24/30 GB, before the cap killed the job. Earlier, at the
1-worker rung, the harness's own per-point subprocess reached 46 GB by
itself in serial mode. The sample function returns a small named tuple
(retcode, wall time, allocation, terminal metrics, metadata), not the
solution, so this is not result retention -- it is Julia's collector
letting each process's heap grow toward the machine's total with nothing
giving it a target to aim for.

Three separate worker-launch sites needed the same fix, because there are
three separate places a process worker gets its `julia` command line built:

1. **`SpaceAGORA.ParallelProcess`'s own pool** (`src/parallel/process/worker_pool.jl`,
   `_process_worker_exeflags`/`ensure_process_workers!`) -- the general
   campaign-dispatch process pool.
2. **The parallelization_performance harness's per-point controller
   subprocess** (`ppc_worker_cmd` in
   `benchmarks/studies/parallelization_performance/execution.jl`) -- the
   worker that runs one `(case, mode, threads, mc)` point end to end (already
   carried the opt-in `SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT` from section 4,
   which was OFF by default; it is now ON by default).
3. **The same harness's own Distributed sub-pool for `outer_process`/GRAM-live
   batches** (`ppc_ensure_process_workers!`/`_ppc_pool_worker_exeflags` in the
   same file) -- a second, independent `addprocs` call with its own hardcoded
   exeflags, separate from (1) and not reached by fixing (1) or (2). This is
   very likely the actual worker-spawn path P6p exercised, since it is the
   one built specifically for GRAM-live Monte Carlo batches.

All three now default to a derived `--heap-size-hint`, computed once in (1)
and reused by (2) and (3) rather than copied three times:

```julia
const _POOL_WORKER_HEAP_HINT_FLOOR_BYTES = 2 * 1024^3    # 2 GiB
const _POOL_WORKER_HEAP_HINT_CEIL_BYTES = 64 * 1024^3    # 64 GiB

const _POOL_WORKER_HEAP_HINT_SHARE_DEFAULT = 0.5           # SPACEAGORA_POOL_WORKER_HEAP_HINT_SHARE

_default_pool_worker_heap_hint_bytes(pool_size) =
    clamp(floor(Int, SHARE * Sys.total_memory() / (pool_size + 1)), FLOOR, CEIL)
```

DERIVED: a share of total memory split `pool_size + 1` ways -- the pool's
workers plus the coordinator process sharing the same machine. The share is
ASSUMED at 0.5: a heap-size hint is a soft target a process can overshoot,
and hinting every process at its equal split of the whole machine would put
the pool plus the coordinator at 100 percent of physical memory with nothing
held back for the native GRAM images, MERRA2 buffers, the OS or another job.
Half the machine goes to Julia heaps in aggregate; the rest is headroom.
On TRX50's roughly 250 GB, the formula gives about 13.9 GB per worker at the
diagnostic's 8-worker rung and about 7.4 GB at 16 workers, in both cases well
below the 30 GB single-process peak the diagnostic observed right before the
kill, so the collector is made to collect long before that point. The floor
and cap are ASSUMED, not measured: 2 GiB so the
hint never fights a single sample's own working set (GRAM tables, harmonics
buffers) rather than bounding growth across many samples, and 64 GiB because
it is comfortably above every per-process peak the diagnostic recorded, so
it only engages for a small pool on a very large machine, where it adds no
protection anyway.

Each of the three sites keeps its own override knob, none of them new
except in default value: `SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT` for (1) and
(3) (the same knob, since (3)'s workers are the same kind of thing under the
same per-process memory-share reasoning), `SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT`
for (2) (the per-point controller's own heap, a different process). `off`
(case-insensitive) removes the flag entirely; any other value overrides the
derived default and is passed through verbatim; unset uses the derived
default.

Byte-identity of the dynamics is untouched by construction: `--heap-size-hint`
is a GC collector setting read by the Julia runtime before any simulation
code runs, and moves no floating-point operation in the RHS, the solver, or
anywhere else the trajectory is computed. Nothing in this workstream's
existing dump-and-cmp evidence (section 2) exercises a process-pool worker,
so there is no new pair to re-prove here -- the claim rests on what
`--heap-size-hint` is (a runtime flag consumed by the garbage collector, with
no code path into the dynamics), not on a repeated trajectory comparison.

Tests: `test/unit/parallel/pool_worker_heap_hint_tests.jl` covers all three
sites' exeflags/argv builders directly (no worker actually spawned) --
default carries the hint, `off` removes it, an explicit value passes
through verbatim, for `_process_worker_exeflags` (1), `_ppc_worker_gc_flags`
(2), and `_ppc_pool_worker_exeflags` (3) -- plus the shared formula's
floor/cap clamping and monotonicity.

## 6. Verification: peak RSS with and without the pool-worker hint

Tooling: `benchmarks/studies/heap_contention/pool_rss_probe.jl` runs a
P6p-shaped campaign (`montecarlo_mars_gram_live`, native-GRAM one-satellite
Monte Carlo samples, `outer_process` route) through the real dispatch path
(`ppc_run_sample_batch` -> `ppc_ensure_process_workers!`, site 3 above) and
prints markers; `pool_rss_probe.sh` runs it twice back to back (once with
`SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT=off`, once at the default), always
under the `systemd-run` memory cap, sampling the RSS of every process
belonging to that run (the coordinator plus its `--worker` children,
distinguished from unrelated `julia`-named processes on the shared machine
by parentage -- see `summarize_pool_rss.py`) every 2 s via `ps`.

This uses the `outer_process` static mode, not the diagnostic's "predictive"
adaptive profile -- both dispatch through the same worker-spawn path this
fix touches, and `outer_process` is simpler to reproduce deterministically.
GRAM is available on this workstation via the `data/GRAMSuite.jl` symlink
into the main checkout, so no substitution was needed.

Three runs, all at 4 workers on space-falcon-1, `--threads=1`:

| Run | Profile | Samples | Arm | Peak total RSS | Max per-worker peak | Wall (s) | Ratio hint/no-hint |
|---|---|---|---|---|---|---|---|
| n1500_smoke | smoke (120 s/sample) | 1500 | no-hint | 8.00 GiB | 1.76 GiB | 55.056 | -- |
| n1500_smoke | smoke | 1500 | hint | 8.22 GiB | 1.81 GiB | 55.426 | 1.0067 |
| n400_full | full (1800 s/sample) | 400 | no-hint | 8.13 GiB | 1.80 GiB | 53.336 | -- |
| n400_full | full | 400 | hint | 8.15 GiB | 1.79 GiB | 53.717 | 1.0071 |
| n4000_full | full | 4000 | no-hint | 8.21 GiB | 1.77 GiB | 95.575 | -- |
| n4000_full | full | 4000 | hint | 8.35 GiB | 1.84 GiB | 91.432 | 0.9567 |

(Full data: `results/pool_rss_verification_summary.csv`. Wall-time ratios
scatter within about ±4% across runs, consistent with noise on a shared
machine rather than a real cost from the extra GC pressure the hint can
introduce -- no run showed a directional, repeatable slowdown.)

**No difference was measurable, at any of the three scales tried, and that
has a specific, checkable explanation rather than being a shrug.** The
computed default hint under this workstation's own `systemd-run
-p MemoryMax=16G` cap was, with the formula as it stood when these runs were
made (no share factor), `Sys.total_memory() ÷ (4 + 1)` = 3.2 GiB per worker
(confirmed directly: the hint arm's worker argv carried
`--heap-size-hint=3.2G`, exactly matching that formula). With the share
factor the same scope gives 8 GiB / 5 = 1.6 GiB, clamped up to the 2 GiB
floor; the runs were not repeated at that value. Every worker's
actual peak RSS across all three runs stayed at 1.7-1.84 GiB -- well under
that 3.2 GiB target -- so the hint never became a binding constraint at any
sample count tried; there was nothing for it to bound. This is the expected
outcome of testing at 1/12.5 the memory budget (16 GB vs. TRX50's ~250 GB
visible total) and 1/2 the worker count (4 vs. 8) the diagnostic ran at: the
formula scales the hint down with the smaller visible memory and the
smaller pool, and growth here never approached even that smaller target.

**This does not verify that the fix prevents the TRX50-scale failure --
only that the mechanism computes and applies the intended flag correctly and
costs no measurable time when it does not bind.** Confirming the actual
mitigation needs one of: (a) running this same probe on TRX50 at the
diagnostic's own scale (8 workers, thousands of samples, predictive mode),
where the formula's ~13.9 GB/worker hint sits well below the 30 GB peak the
diagnostic recorded right before the kill; or (b) forcing an artificially
tight override (e.g. `SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT=500M`) on this
box to at least confirm Julia's collector visibly respects a hint set below
where RSS would otherwise land -- neither was run here.

**Correction, later measurement.** The P6p growth this section attributes to
Julia's collector was not on the Julia heap: every trace-6 sample rebuilt a
native GRAM atmosphere that only a finalizer frees, about 106 MB resident per
sample, and a worker with a 7.4 GB hint ran 147 samples without a single
collection. The hint is passed but cannot bound that memory. See
`gram_thread_scaling.md`, "Pool-worker memory growth", for the measurements and
the fix.

## What to run on TRX50 (16+ threads)

This box cannot reproduce the 16+-thread regime where the crossover in the
opening paragraph appears. On TRX50, quiet of other users' jobs
(`pgrep -u <user> julia` empty):

1. Re-run `collector_grid.sh` (or the `SPACEAGORA_PPC_WORKER_GCTHREADS`/
   `_HEAP_SIZE_HINT` env pair above through the real harness) at 8, 16, and
   32 threads on P3 and P4, pinned-threads (`outer_threads`) mode, with
   `--cpu-list` set the way the archived TRX50 cold run was. If GC threads
   scale with the widening ladder the way OS-thread contention does, the
   16/32-thread rows should show the same crossover the archive's raw
   wall-clock numbers show for `outer_threads` vs. `outer_process`.
2. Repeat the allocation attribution (`alloc_profile.jl`, remembering the
   P4 `--sample-rate` ceiling and the memory-cap wrapper above) at the same
   thread counts if `Profile.Allocs` proves safe to run concurrently with
   the solve's own thread pool -- this document's P3/P4 numbers are all
   single-threaded and say nothing about whether the *bucket* split shifts
   under contention, only about where the bytes originate.
3. Run `pool_rss_probe.sh montecarlo_mars_gram_live 8 4096` (matching the
   P6p diagnostic's own worker count and sample count) under whatever memory
   cap TRX50's own job used, with and without
   `SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT=off`. At that scale the derived
   hint (~13.9 GB/worker on a ~250 GB machine) sits well below the 30 GB
   single-process peak the diagnostic recorded right before its cap killed
   the job, which is where section 6's local runs (hint always well above
   the observed peak, so never binding) could not reach.
