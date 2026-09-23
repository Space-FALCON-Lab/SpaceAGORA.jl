# One-heap contention above 16 threads: what could and couldn't be measured here

This is WS11c's record. The problem this workstream exists for is a TRX50
observation: the pinned-threads route on P3/P4 gets *slower* as thread count
rises past 8 (1.55 s -> 2.06 s -> 2.72 s at 8/16/32 threads on P3), while the
pinned pool route keeps getting faster (1.28 -> 0.76 -> 0.43 s). That is the
signature of single-heap GC contention, and it needs 16+ threads to show up.

**This workstation has 12 physical cores / 24 threads and cannot reproduce the
16+-thread regime.** Everything below is either (a) allocation attribution and
reduction, which is thread-count-independent and fully measurable here, or (b)
an 8-thread collector-settings grid, which is the closest this box gets to the
question and is reported as exactly that -- a data point below the regime that
actually matters, not a stand-in for it. Section "What to run on TRX50" at the
bottom says what closes the gap.

## Incident: `Profile.Allocs` at `sample_rate=1.0` on P4 took the whole
## session down

While attributing P4's allocation (see below), an early version of
`alloc_profile.jl` ran `Profile.Allocs` at `sample_rate=1.0` on
`montecarlo_heavy_aerobraking`. That case is a 6 h mission at a 1 s max step,
~365 MiB and several million allocations per sample; `Profile.Allocs`
captures a full stack per *sampled* allocation, so its own bookkeeping grows
with the sampled allocation count, not the run's byte total. At full sampling
that bookkeeping reached 45.5 GB resident on a 60 GB box and the kernel OOM
killer took out the process's cgroup -- which also killed the orchestrator,
every other WS11 agent's running job, and the remote-job drivers sharing the
machine.

Two changes came out of that, and both are load-bearing, not optional
polish:

1. **Every Julia process this study launches is wrapped in a hard memory
   cap**: `systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q
   -- julia ...`. A run that exceeds 16 GB is killed alone instead of taking
   the machine down. Every measurement in this document was run this way;
   `collector_grid.sh` wraps it automatically, and `alloc_profile.jl`'s and
   `dump_state.jl`'s usage headers show the wrapped invocation.
2. **`alloc_profile.jl` refuses full sampling on anything P4-scale.** P3
   (~15 MiB/sample) is small enough to sample exhaustively
   (`--sample-rate=1.0`, ~180k allocations). P4 runs at `--sample-rate=0.005`,
   two orders of magnitude below the coordinator's 0.01 ceiling. A uniform
   sub-sample still gives a valid *fraction* of allocation per bucket (that is
   all this study needs); it does not give exact byte/count totals, and the
   script's CSV output says so in its column names
   (`estimated_bytes_at_1x` is explicitly an estimate, not a measurement).

## 1. Allocation attribution (P3, P4)

Tool: `benchmarks/studies/heap_contention/alloc_profile.jl`. Warms up one
solve (JIT), times a second with `@timed` for the headline bytes/accepted-step
number, then profiles a third with `Profile.Allocs` and buckets every sampled
allocation by the `src/` domain that owns its innermost SpaceAGORA.jl stack
frame -- split finely enough to separate WS11c's own two owned files
(`engine/execution.jl`, `engine/persistence.jl`) from the rest of
`src/simulation/engine/` (`solver_policy.jl`, `dynamics_rhs.jl`,
`rhs_calibration.jl`, ...), which other workstreams own.

Measured single-threaded (`--threads=1`), space-falcon-1, under the memory
cap, `SPACEAGORA_PARALLEL_POLICY` unset (serial baseline):

| Case | naccept | total alloc | bytes / accepted step |
|---|---|---|---|
| P3 (`independent_1sat_1hr`) | 250 | 14.99 MiB | 62,860 |
| P4 (`montecarlo_heavy_aerobraking`) | 22,487 | 365.79 MiB | 17,057 |

(P4's 365.79 MiB matches the contract's cited 364 MB closely enough to confirm
this is the same case/shape the TRX50 archive measured.)

Allocation attribution by bucket (see `results/alloc_p3_independent_1sat_1hr.csv`,
`results/alloc_p4_montecarlo_heavy_aerobraking.csv`):

| Bucket | P3 (sample_rate=1.0) | P4 (sample_rate=0.005, estimate) |
|---|---|---|
| `engine/other` (solver_policy, dynamics_rhs, rhs_calibration -- **not WS11c**) | 51.5% | 56.1% |
| `callbacks` (SavedValues/SaveData/SavingCallback -- **not WS11c**) | 28.2% | 12.1% |
| `external` (OrdinaryDiffEq/DiffEqCallbacks/Base) | 9.0% | 0.3% |
| `parallel` (routing/policy) | 5.6% | 13.7% |
| `environment` (gravity/atmosphere/ephemerides) | 3.4% | -- |
| `dynamics` (RHS terms) | 2.0% | -- |
| `core` | 0.0% | 17.7% |
| **`engine/execution.jl` (WS11c-owned)** | **0.4%** | **0.0%** |
| `engine/persistence.jl` (WS11c-owned) | 0.0% (no samples) | 0.0% (no samples) |

**The honest finding: the two files WS11c owns are not where P3/P4's
per-sample allocation lives.** `execution.jl` accounts for 0.4% of P3's
sampled allocation and effectively none of P4's; `persistence.jl` (46 lines of
thin re-export wrappers, see its header) accounts for none at either scale.
The dominant sources -- `solver_policy.jl`/`dynamics_rhs.jl` (the solve driver
and RHS dispatch) and the `SavedValues`/`SaveData` callback storage -- are
owned by other workstreams (solver/RHS routing, and callbacks, respectively;
see the domain-ownership table in `CLAUDE.md`). "What fraction is the saved
solution" is therefore mostly the `callbacks` bucket (SavingCallback building
and growing `saved_values.t`/`saved_values.saveval`), not anything inside
`execution.jl`'s own accumulation loops -- those loops (`_append_backbone_saved_segment!`,
`_append_checkpoint_saved_segment!`) are not exercised by P3/P4 at all, since
neither case checkpoints or uses `gravity_backbone_split`.

## 2. Reductions made, with dump-and-cmp proof

Given (1), the one change inside WS11c's own files that is both real and
correctness-neutral is: **stop building the debug `ODEProblem` when nothing
reads it.**

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
a change to that code path without inventing a fourth case outside this
contract's list, and per the standing rule ("a change that cannot be proven
bit-identical is not shipped"), it stays out.

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
| `independent_1sat_1hr` (P3) | the contract's P3 sample | 18,072 bytes | **identical** |
| `montecarlo_heavy_aerobraking` (P4) | the contract's P4 sample | 1,619,136 bytes | **identical** |
| `atmo256_vacuum_10min` | the contract's "256-spacecraft vacuum case" -- 256 sat, L50 harmonics, no atmosphere, same physics family as `gravity_256sat_l50_vacuum_124000s` but a 10 min mission instead of 124,000 s, chosen so the case is fast enough to dump twice per iteration instead of once per multi-minute run | 2,032,608 bytes | **identical** |

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

This is a real but *partial* reduction, and the doc says so plainly rather
than overstating it: thinning shrinks what gets written and what the caller
holds onto after the solve, but section 1 showed the solver's own per-step
accumulation (which the solve has already paid for by the time this code
runs) is not inside WS11c's files to begin with, so this setting does not
reduce the solve's own peak allocation -- only the results-writing tail past
it.

## 4. Collector settings (`--gcthreads`, `--heap-size-hint`)

### What was measured, and what wasn't

The contract asked for 8 and 12 threads. Mid-task, after the OOM incident
above, the coordinator's recovery message reiterated a hard "at most 8
threads" rule for every Julia process launched from this workstream, without
qualification. That supersedes the contract's 12-thread ask for this specific
measurement, so **only 8 threads was measured here**; 12 was not run. Neither
was CPU pinning applied (this box has no `--cpu-list` configured for this
study the way TRX50 does), so this is the `outer_threads` route unpinned, not
literally the "pinned-threads" route the contract describes -- the closest
available approximation on this machine.

Tooling: `benchmarks/studies/heap_contention/collector_grid.sh` launches one
Julia subprocess per grid point (under the same memory cap as everywhere
else in this study), each running `collector_grid.jl`'s 5 back-to-back
`outer_threads`-mode solves and reporting the median wall time. Grid:
`--gcthreads` in {unset (Julia default), `2`}; `--heap-size-hint` in {unset,
`4G`}. The machine was **not** quiet during this run (`uptime` load average
6.6-11.2, ~17 other `julia`-named processes from the other concurrently
running WS11 agents throughout) -- ratios below are reported as measured,
not cleaned up, and should be read with that noise floor in mind.

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
under a noisy shared machine, at a thread count (8) below the regime this
workstream is actually about. **No default was changed.** The mechanism was
wired as an opt-in only (`ppc_worker_cmd` in
`benchmarks/studies/parallelization_performance/execution.jl`, gated by
`SPACEAGORA_PPC_WORKER_GCTHREADS` / `SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT`,
both unset by default -- unset means no flag is added to the worker's `julia`
invocation at all, not "pass Julia's own default explicitly"), so a future
run on a quiet machine (TRX50, or this box between other agents' jobs) can
exercise the grid without a source change:

```bash
SPACEAGORA_PPC_WORKER_GCTHREADS=2 SPACEAGORA_PPC_WORKER_HEAP_SIZE_HINT=4G \
  julia --project=. benchmarks/studies/parallelization_performance.jl \
  --phases=P3,P4 --threads=8,16,32 --process-workers=32
```

Raw data: `results/collector_grid_p3_t8.csv`, `results/collector_grid_p4_t8.csv`.

## What to run on TRX50 (16+ threads, this workstream's actual question)

This box cannot produce the signal this workstream exists to explain. On
TRX50, quiet (`pgrep -u <user> julia` empty, per the standing rule):

1. Re-run `collector_grid.sh` (or the `SPACEAGORA_PPC_WORKER_GCTHREADS`/
   `_HEAP_SIZE_HINT` env pair above through the real harness) at 8, 16, and
   32 threads on P3 and P4, pinned-threads (`outer_threads`) mode, with
   `--cpu-list` set the way the archived TRX50 cold run was. If GC threads
   scale with the widening ladder the way OS-thread contention does, the
   16/32-thread rows should show the same crossover the archive's raw
   wall-clock numbers show for `outer_threads` vs. `outer_process`.
2. Repeat the allocation attribution (`alloc_profile.jl`, remembering the
   P4 `--sample-rate` ceiling and the memory-cap wrapper) at the same thread
   counts if `Profile.Allocs` proves safe to run concurrently with the
   solve's own thread pool -- this document's P3/P4 numbers are all
   single-threaded and say nothing about whether the *bucket* split shifts
   under contention, only about where the bytes originate.
