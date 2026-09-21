# paper_scenarios — JAIS Parallelization Paper Benchmark Suite (S1–S5)

Self-contained benchmark harness backing the JAIS paper's parallelization claims.
Five scenarios tell one story: **on a fixed, resource-constrained machine, SpaceAGORA
extracts near-full utilization from whatever cores and RAM exist**, across two
orthogonal workload axes — constellation size and Monte Carlo sample count —
including both at once. All scenarios are Earth-based (Earth GRAM, GGM05C L20/L50
harmonics) so no cross-planet model mismatch confounds comparisons.

| Scenario | Script | Claim it isolates |
|---|---|---|
| S1 | `s1_constellation_scaling.jl` | Outer thread routing scales a coupled 1→256-sat propagation near-linearly; marginal memory per satellite is a state vector, not a process image |
| S2 | `s2_gram_atmosphere_modes.jl` | With a non-thread-safe native atmosphere model, the look-ahead density cache recovers most of the process-route speedup with **one** GRAM instance — the speed-per-GB centerpiece |
| S3 | `s3_montecarlo_process_scaling.jl` | Process-route Monte Carlo throughput scales near-linearly with workers; per-worker RSS quantifies the RAM price of process parallelism |
| S4 | `s4_hybrid_mc_constellation.jl` | A fixed thread budget can be partitioned between outer (sample) and inner (constellation) parallelism, safely, under one resource governor |
| S5 | `s5_routing_profile_ladder.jl` | Adaptive profiles R4/R5 match or approach the best hand-picked profile at fixed budget — no expert tuning, portable across machines |

Story arc: S1 establishes the clean scaling mechanism → S2 shows it surviving the
real obstacle (native GRAM) without buying RAM → S3 covers the orthogonal MC axis →
S4 composes the two axes → S5 shows the whole thing is automatic and portable.

## Running

Each point is a separate `julia --threads=N` subprocess (thread count is fixed at
Julia startup), spawned by the controller and measured by `scenario_worker.jl`.

```bash
# Individual scenarios
julia --project=. benchmarks/studies/paper_scenarios/s1_constellation_scaling.jl
julia --project=. benchmarks/studies/paper_scenarios/s2_gram_atmosphere_modes.jl
julia --project=. benchmarks/studies/paper_scenarios/s3_montecarlo_process_scaling.jl
julia --project=. benchmarks/studies/paper_scenarios/s4_hybrid_mc_constellation.jl
julia --project=. benchmarks/studies/paper_scenarios/s5_routing_profile_ladder.jl

# Everything (paper protocol order; S2/S3 are the long ones)
bash benchmarks/studies/paper_scenarios/run_all.sh
```

Results land in `results/<hostname>/sN_*.csv` plus one `.log` per point. Run the
suite on each paper machine (M3 Pro MacBook Pro, Zen 5 Ryzen 9 / 64 GB,
Threadripper 64-thread / 256 GB); the hostname keys the per-machine tables.

**Status:** `space-falcon-1` (Zen 5 Ryzen 9, 12c/24t) and
`space-falcon-lab-TRX50-AERO-D` (Threadripper, 64c/128t, run at
`PS_THREADS=64`) both have full S1–S5 data. The M3 Pro MacBook Pro is still
outstanding.

### Running on a remote box

`scripts/remote/spaceagora-remote` mirrors this repo to a remote host and runs
a job there detached (config: `scripts/remote/remotes.conf`, gitignored — copy
`scripts/remote/remotes.conf.example`). It syncs the full local working tree
(gitignore-aware, so uncommitted edits and untracked scripts ship too — not
just `git ls-files`) and runs `Pkg.instantiate()` before every job, so a stale
shared package depot fails fast with one clear log instead of every scenario
point silently erroring on a missing package. GRAM/SPICE data (~20G) is opt-in
via `spaceagora-remote seed-data --remote <name>` (or `push --sync-data`) — run
it once per remote before the first job that needs it. Its auto-pull only
captures a release's top-level `output/`/`results/` dirs, not the nested
`benchmarks/studies/paper_scenarios/results/<hostname>/` path this suite
writes to, so pull that manually once the job finishes:

```bash
rsync -az <ssh-alias>:<remote-base>/releases/<job-id>/benchmarks/studies/paper_scenarios/results/<hostname>/ \
  benchmarks/studies/paper_scenarios/results/<hostname>/
```

## Knobs

| Env var | Meaning | Default |
|---|---|---|
| `PS_THREADS` | Parallel thread budget per point — **set to physical core count for paper numbers**. S1 also accepts a comma-separated ladder (e.g. `1,2,4,8,12`): serial runs once per size regardless of ladder length (thread-count-independent), and each ladder value adds one more `parallel` row for that size, distinguished by `julia_threads` — a full thread-count sweep is one invocation/CSV write, not one per thread value | `Sys.CPU_THREADS` |
| `PS_REPEATS` / `PS_WARMUP` | Timed repeats / warmup solves per point | 3 / 1 |
| `PS_MISSION_S` | Mission length override | 1800 s (S1), 600 s (S2–S5) |
| `PS_SIZES` | Constellation size ladder (S1, S2) | `1,4,16,64,256,1024,2048,4096` / `4,16,64` |
| `PS_WORKERS` | S3 process-worker ladder | `1,2,4,8` (extend to `16,32` on the Threadripper) |
| `PS_MC_SAMPLES` | Samples per campaign (S3, S4) | 32 / 8 |
| `PS_PROC_WORKERS` | S2 process-mode worker count | `min(8, PS_THREADS)` |
| `PS_TIMEOUT_S` | Per-point kill deadline | 3600 s |
| `PS_GRAVITY` / `PS_DENSITY` | S1-only: force-model override (`invsq\|l20\|l50` / `none\|gram_standard\|gram_lookahead\|gram_surrogate`), isolates per-satellite cost as a variable | `l20` / `none` |
| `PS_RESULTS_SUFFIX` | Appended to `results/<hostname>` — use for exploratory runs (e.g. a thread-count sweep, or a one-off thread-budget override) so they land in their own pseudo-host directory instead of overwriting the canonical per-host baseline CSV | (empty) |

Mission lengths were chosen so each point is long enough that per-step work
dominates dispatch overhead (hundreds of accepted steps) but a full scenario stays
in the minutes-to-an-hour range on the workstation-class machines. `--preview`-style
smoke: set `PS_SIZES=1,4 PS_MC_SAMPLES=4 PS_WORKERS=1,2 PS_MISSION_S=60 PS_REPEATS=1`.

## Interpreting the CSVs

- `median_s`/`min_s`/`max_s`/`times_s` — wall time over `PS_REPEATS` (median is the
  paper metric); warmups are excluded, so JIT and GRAM data loads never pollute timing.
- `maxrss_mb` — measured process peak RSS (`Sys.maxrss`), the memory column.
- `workers_rss_mb` — summed peak RSS of Distributed pool workers (process modes
  only): the RAM cost of process-level parallelism that S2's shared-instance
  thread route avoids.
- Speedup for S1/S2 = serial median / mode median at the same size. S3's headline
  is per-sample time (`median_s / samples`) vs. workers, which should stay ~flat.
- `PS_GRAVITY`/`PS_DENSITY` overrides in S1 write a suffixed CSV
  (`s1_constellation_scaling_<gravity>_<density>.csv`) instead of overwriting the
  (`l20`, vacuum) baseline; `plot_results.jl` picks up every such file automatically
  and tags rows with a `variant` column. A logged finding from an `l50` comparison
  run: at low N the baseline's cheap per-satellite cost makes fixed dispatch
  overhead dominate (including a specific dip at N=16, where
  `SPACEAGORA_RHS_BATCH_THREAD_THRESHOLD` first turns batched RHS threading on);
  a heavier force model shrinks that overhead fraction and produces *superlinear*
  speedup at N >= 1024 (up to 21.6x on a 12-thread budget) — reproducible but not
  confirmed against hardware cache-miss counters. See "S1 Force-Model Sensitivity
  Finding" in `docs/architecture/parallelization_paper_notes.md` for the full data
  and the caveat to state alongside it in the paper.

## S1 finding: TRX50 thread-count sensitivity at N=2048/4096

Both `l20` and `l50` full ladders (N=1…4096) are now collected on both paper
machines at their native thread budget (12 on `space-falcon-1`, 64 on
`space-falcon-lab-TRX50-AERO-D`). At those native budgets, TRX50 badly
underutilizes its larger thread count at N=2048/4096 relative to
`min(N, threads)` ideal — e.g. `l20`/N=2048 is 3.86x speedup on TRX50 (64
threads, 6% of ideal) vs. 8.96x on `space-falcon-1` (12 threads, 75% of
ideal).

Two follow-up runs isolate why:

- **Matched thread count.** Re-running TRX50 capped to `PS_THREADS=12`
  (matching `space-falcon-1`) shows TRX50 *matches or exceeds* local's
  efficiency at every (N, gravity) point tested — e.g. `l50`/N=4096: 27.6x vs.
  21.6x. This rules out a raw hardware/per-core-speed deficit: serial timing
  is consistent between the 64-thread and 12-thread TRX50 runs (as expected,
  since `PS_THREADS` only affects parallel mode).
- **Full thread-count sweep.** `PS_THREADS=1,2,4,8,16,32,64` (TRX50) /
  `1,2,4,8,12` (`space-falcon-1`) at N=2048/4096, both gravities, shows TRX50's
  speedup curve peaking around **T=16–32** and *degrading* out to T=64 in
  every case (e.g. `l20`/N=4096: 13.90x at T=16 → 8.54x at T=64), while
  `space-falcon-1`'s curve is still rising or flat at its 12-thread ceiling.

Conclusion: N=2048/4096 doesn't carry enough work to keep more than ~16–32
threads fed on TRX50 at this per-satellite cost; past that point,
dispatch/synchronization overhead across the extra threads exceeds the
marginal parallel benefit. **Any TRX50/64-thread S1 number used in the paper
should carry this caveat** — the headline speedup depends on matching thread
budget to problem size, not simply "more threads is better."

Raw data:

- `results/<hostname>/s1_constellation_scaling{,_l50_none}.csv` — full ladder
  at each machine's native/paper thread budget (the baseline table).
- `results/space-falcon-lab-TRX50-AERO-D_PS_THREADS12/` — TRX50 re-run at
  `PS_THREADS=12`, N=2048/4096 only, for the matched-thread-count comparison.
- `results/<hostname>_thread_ladder/` — the full `PS_THREADS` sweep,
  N=2048/4096 only, both gravities, on each machine.

The `_PS_THREADS12`/`_thread_ladder` directories are exploratory pseudo-hosts
(via `PS_RESULTS_SUFFIX`), not part of the paper-baseline table, and
`plot_results.jl`'s auto-generated per-host S1 plot isn't built for multiple
`parallel` rows sharing one `n_sats` (it draws a connect-the-dots line across
thread values at that N) — treat its output for these directories as
diagnostic only, not a paper figure.

## GRAM safety conventions (inherited from prior studies)

- `gram_standard` (direct native GRAM): density frozen per accepted step
  (`SPACEAGORA_DENSITY_FREEZE_PER_STEP=1`; per-RK-stage perturbation noise
  otherwise collapses adaptive dt) and density callback kept serial.
- `gram_lookahead`: vacuum-predicted cache with horizon/deviation loosened past
  mission reach so only the proven-safe initial build runs (native rebuild with
  2+ satellites hangs — see `[[project_gram_multisat_rebuild_bug]]`).
- Process modes: each pool worker is `--threads=1` with its own native GRAM; the
  coordinator only dispatches.

`FABLE_FINDINGS.md` in this folder documents the performance fixes implemented
alongside this suite and the deferred optimization opportunities found during the
same code review.
