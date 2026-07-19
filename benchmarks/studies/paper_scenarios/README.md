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

## Knobs

| Env var | Meaning | Default |
|---|---|---|
| `PS_THREADS` | Parallel thread budget per point — **set to physical core count for paper numbers** | `Sys.CPU_THREADS` |
| `PS_REPEATS` / `PS_WARMUP` | Timed repeats / warmup solves per point | 3 / 1 |
| `PS_MISSION_S` | Mission length override | 1800 s (S1), 600 s (S2–S5) |
| `PS_SIZES` | Constellation size ladder (S1, S2) | `1,4,16,64,256,1024,2048,4096` / `4,16,64` |
| `PS_WORKERS` | S3 process-worker ladder | `1,2,4,8` (extend to `16,32` on the Threadripper) |
| `PS_MC_SAMPLES` | Samples per campaign (S3, S4) | 32 / 8 |
| `PS_PROC_WORKERS` | S2 process-mode worker count | `min(8, PS_THREADS)` |
| `PS_TIMEOUT_S` | Per-point kill deadline | 3600 s |
| `PS_GRAVITY` / `PS_DENSITY` | S1-only: force-model override (`invsq\|l20\|l50` / `none\|gram_standard\|gram_lookahead\|gram_surrogate`), isolates per-satellite cost as a variable | `l20` / `none` |

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
