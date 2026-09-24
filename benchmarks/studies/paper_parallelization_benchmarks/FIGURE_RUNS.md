# The paper's figure runs

Which benchmark goes on which machine, in what order, from what calibration-store
state, and — for every number that had to be chosen rather than measured — where
that number came from.

`CASES.md` says what each phase measures. `protocol.sh` runs the suite.
`paper_figure_runs.sh` in this directory runs *the figures*: it prints the
ordered sequence per machine and, with `--execute`, performs it, refusing to
start a timed run on a busy box.

Every number below is marked **SOURCED** (read out of a named file), **DERIVED**
(computed from sourced numbers by a stated rule) or **ASSUMED** (a judgement,
with the reason). Nothing here is an estimate with no provenance.

---

## Figure F2 — P6 and P6p: thread scaling across force models and density paths

Six traces at one constellation size over the thread ladder `1, 2, 4, 8, 16, 32`.
Trace 2 is P1/P2's own 4096 rung, reused rather than duplicated, so the figure
and those tables share a serial baseline. See `CASES.md` for the trace table, for
why the figure needs two phases, and for how the GRAM density path is selected.

### Choosing the constellation size

The user asked for 4096 spacecraft. The routing layer's own memory model says
that cannot work for the native-GRAM traces:
`native_gram_worker_extra_bytes` in `src/parallel/routing/machine_topology.jl`
charges `_GRAM_SAT_MEMORY_BYTES = 90 MB` per spacecraft ("~22 GB / 256"), which
puts 4096 spacecraft at about 360 GB against the benchmark box's 250 GB.

The measurements say otherwise, and the measurements are what this decision is
made on.

**Installed memory, benchmark box: 250 GB.**
SOURCED — `docs/adaptive_policy_r6_nested_campaign_record_20260906.md` line 656,
"The 250 GB TRX50 never hit the cap". The budget this decision is held to is the
contract's 200 GB, i.e. 80% of that.

#### Resident memory of a live-GRAM constellation, in one process

SOURCED — `benchmarks/studies/paper_scenarios/results/space-falcon-lab-TRX50-AERO-D/s2_gram_atmosphere_modes.csv`,
column `maxrss_mb` (measured `Sys.maxrss`, decimal MB), scenario S2 on the
benchmark box, 600 s mission, degree-20 gravity, live native GRAM:

| `n_sats` | `serial_standard` | `threads_standard` | `threads_lookahead` |
|---:|---:|---:|---:|
| 1 | 1920.9 | 1734.9 | 1737.1 |
| 4 | 1759.3 | 1799.0 | 1797.1 |
| 16 | 1831.6 | 1751.8 | 1713.2 |
| 64 | 1851.9 | 1709.4 | 1818.8 |
| 256 | 1966.2 | 1938.8 | 2039.8 |
| 1024 | 2842.0 | 2781.0 | 2819.6 |

DERIVED — ordinary least squares of `maxrss_mb` on `n_sats` over those six
points, per mode:

| Mode | Fit (MB) | R² | Predicted at 4096 |
|---|---|---:|---:|
| `serial_standard` | 1803.0 + 0.9917·n | 0.969 | 5865 MB |
| `threads_standard` | 1720.1 + 1.0215·n | 0.984 | 5904 MB |
| `threads_lookahead` | 1748.5 + 1.0510·n | 0.994 | 6053 MB |

So a live-GRAM constellation costs **~1.7 GB fixed plus ~1.0 MB per
spacecraft**, not 90 MB per spacecraft. At 4096 that is ~5.9 GB, not ~360 GB.

That fit is confirmed independently, and at the exact size in question, by the
S1 scenario, which ran the same constellation ladder at degree 50 in vacuum all
the way to 4096 on the same box. SOURCED —
`.../space-falcon-lab-TRX50-AERO-D/s1_constellation_scaling_l50_none.csv`:
`maxrss_mb` = **5410.5 MB** at `n_sats = 4096, mode = serial` and **5469.2 MB**
at `mode = parallel` (64 threads). The degree-20 baseline CSV
(`s1_constellation_scaling.csv`) gives 5609.6 / 5343.9 MB at the same size, i.e.
harmonic degree does not move the memory; per-spacecraft state does. Subtracting
S1's vacuum figure from S2's GRAM figure at matching sizes puts the whole native
GRAM + aero addition at 0.05–0.3 GB, flat in `n_sats`.

**Verdict for traces 1–5 (one process): 4096 fits, at ~6 GB, with a 30×
margin.** The `_GRAM_SAT_MEMORY_BYTES` constant over-predicts this workload by
roughly 46× at 256 spacecraft. That discrepancy is recorded, not resolved: the
constant may have been fitted to a configuration that holds per-spacecraft GRAM
state, which none of these cases do. It is not this workstream's to change, and
nothing here depends on it — the P6p worker count is set by the phase, not by the
router.

#### Resident memory of the process route

SOURCED — same S2 CSV, `mode = process_members`: one coordinator plus 8
`--threads=1` pool workers, each with its own native GRAM image, `n_sats`
one-spacecraft samples spread across them. `workers_rss_mb` is the *sum* over
the pool.

| `n_sats` | samples/worker | `workers_rss_mb` | per worker | coordinator `maxrss_mb` |
|---:|---:|---:|---:|---:|
| 1 | 0.125 | 16353.6 | 2044.2 | 1216.6 |
| 4 | 0.5 | 17532.4 | 2191.6 | 1217.2 |
| 16 | 2 | 21582.7 | 2697.8 | 1225.8 |
| 64 | 8 | 26727.8 | 3341.0 | 1217.2 |
| 256 | 32 | 30057.9 | 3757.2 | 1222.8 |
| 1024 | 128 | 43225.7 | **5403.2** | 1224.8 |

Per-worker RSS grows with the number of samples that worker runs, sublinearly —
about 823 MB per doubling over the top two points. The binding rung for P6p is
the widest one: **4096 samples over 32 workers is 128 samples per worker, which
is exactly the measured `n_sats = 1024` row**, so the per-worker figure needs no
extrapolation in the direction that is uncertain. The only extrapolation is in
worker *count*, and each worker is a separate process image, so that is a sum of
independent measured quantities.

DERIVED — total resident memory per rung of P6p at 4096 samples, using the
measured 5403.2 MB where `samples/worker ≤ 128` and the 823 MB-per-doubling
slope above it:

| Workers | samples/worker | per worker | total (with coordinator) |
|---:|---:|---:|---:|
| 1 | 4096 | 9518 MB | 10.7 GB |
| 2 | 2048 | 8695 MB | 18.6 GB |
| 4 | 1024 | 7872 MB | 32.7 GB |
| 8 | 512 | 7049 MB | 57.6 GB |
| 16 | 256 | 6226 MB | 100.8 GB |
| 32 | 128 | 5403 MB | **174.1 GB** |

**Verdict: 4096 fits the process trace too, at 174.1 GB against the 200 GB
budget and 250 GB installed.** It is the one rung in this figure with a margin
under 1.5×, so P6p is the phase to watch in `dmesg` for the OOM killer, and the
`--process-workers 32` cap must not be raised on this case.

**Decision: N = 4096 for all six traces. No second phase at a smaller size is
needed**, so the `P6g` fallback the contract allowed for is not defined. If a
future host has less memory, the rung table above says where it has to stop.

### Choosing the mission lengths

A single mission length across force models this different is not usable: it puts
the light traces under the harness's 3 s measurability floor and the GRAM traces
into the hours. So each trace is sized against **trace 2's measured serial
baseline**, which is the same rule `PPC_L50_ISO_MISSION_S` applies across
spacecraft counts, and the mission length becomes a column of the figure's table.

**Trace 2's serial baseline: 11.856630086898804 s.**
SOURCED — `.claude/worktrees/paper-routing/output/performance/paper_benchmarks_trx50_cold11/20260918_162845/paper_benchmarks_aggregated_20260918_162845.csv`,
row `phase_id = P2, case = gravity_4096sat_l50_vacuum_5800s, mode = serial,
thread_count = 1`, `wall_time_median_s`, `n_repeats = 11`, machine
`space-falcon-lab-TRX50-AERO-D` (from the sibling raw CSV's `machine` column).

#### Trace 1 — degree 20, vacuum: 19 700 s

SOURCED — `s1_constellation_scaling_l50_none.csv` and
`s1_constellation_scaling.csv`, `n_sats = 4096, mode = serial`, 1800 s mission,
same box: **210.4728 s** at degree 50 and **62.0167 s** at degree 20.
DERIVED — ratio 3.3938; mission = 5800 × 3.3938 = 19 684 s, registered as
**19 700 s**, predicted serial baseline 11.9 s.

The caveat, and why the rounding is not the uncertainty that matters: S1's
`serial` mode runs with `SPACEAGORA_HARMONICS_BATCH_ENABLED=0`
(`paper_scenarios/common.jl`) while the ppb `serial` mode always runs with it set
to 1 (`ppc_mode_env_pairs`), so 3.3938 is the *unbatched* ratio. The batched
kernel's cost is closer to the coefficient count, whose ratio is
(51·52/2)/(21·22/2) = 5.7403. The true ratio is bracketed by [3.39, 5.74], and
19 700 s puts the predicted serial baseline at **7.0–11.9 s across that whole
bracket** — above the 3 s floor either way. Trace 2's own 5800 s would put it at
2.1–3.5 s, on or under the floor, which is why this trace is iso-work against
trace 2 rather than iso-duration.

No in-harness degree-20/degree-50 pair could be used instead: the only ppb run
containing `gravity_1024sat_l20_vacuum` is `20260701_190405`, from before the
batched harmonics kernel existed, and its serial median (24.698 s at 1800 s)
is 5.5× *more* than the degree-50 case measures today. It is not a usable ratio.

#### Trace 3 — degree 50 + SRP + Sun/Moon third body, vacuum: 5800 s

**Trace 2's mission, unchanged.** Adding effectors can only raise the serial
baseline, so this trace is above the floor by construction and needs no
extrapolation to size — and trace 3 against trace 2 becomes an exact iso-mission,
iso-size, single-variable comparison.

DERIVED, for the cost estimate only — SOURCED from
`output/performance/paper_benchmarks/trx50_quiet_20260829/paper_benchmarks_aggregated_20260829_195452.csv`:
`atmo256_gram_live_nbody_10min` serial 46.34575700759888 s against
`atmo256_gram_live_10min` serial 42.05648183822632 s, i.e. the Sun/Moon third
body costs +10.2% at 256 spacecraft. Predicted baseline **13–15 s**.

#### Traces 4, 5, 6 — the density ladder: 100 s

These three share one mission length on purpose: they are the `atmo256_*`
ladder's design (fixed spacecraft count, fixed duration, fixed harmonic degree,
only the density model varying) moved to the figure's constellation size, and
that design is the only thing that makes a difference between them attributable
to the density path. They also share `dt_max_orbit = 5 s` with that ladder,
rather than the 20 s the vacuum traces use.

SOURCED — same `trx50_quiet_20260829` aggregated CSV, `mode = serial`,
`thread_count = 1`, 256 spacecraft, 600 s mission, degree-50 harmonics:

| Case | serial median |
|---|---:|
| `atmo256_exponential_10min` | 4.853820085525513 s |
| `atmo256_gram_surrogate_10min` | 14.83997893333435 s |
| `atmo256_gram_live_10min` | 42.05648183822632 s |
| `atmo256_gram_live_nbody_10min` | 46.34575700759888 s |

(The second archived benchmark-box run, `trx50_20260828_193820`, gives 4.771 /
13.065 / 42.380 / 46.084 s for the same four — within 12% on the surrogate and
within 2% on the rest, which is the run-to-run spread these predictions carry.)

DERIVED, assuming cost is linear in spacecraft count and in mission length:

- exponential: 4.8538 / (256 × 600) = 3.160×10⁻⁵ s per spacecraft-second.
  A 11.86 s baseline at 4096 spacecraft wants **91.6 s** of mission.
- live GRAM: 42.0565 / (256 × 600) = 2.738×10⁻⁴ s per spacecraft-second.
  A 11.86 s baseline at 4096 spacecraft wants **10.6 s** of mission.

ASSUMED — **100 s for all three**, rounded up from the exponential trace's
91.6 s. Predicted serial baselines: trace 4 **12.9 s** (9% above trace 2's
measured value, the closest match the shared duration allows), trace 5
**112 s**, trace 6 **~357 s** over 4096 samples.

Trace 5 is deliberately not sized to 11.86 s. Doing so would demand a 10.6 s
mission at 4096 spacecraft, at which point the measurement is dominated by
building 4096 spacecraft and assembling a solver rather than by propagating
anything, and the look-ahead cache — whose entire subject is a spline over a
look-ahead horizon — would never be exercised. 100 s is the shortest mission at
which all three density traces are still simulations, and trace 5's 9.4× cost
over trace 2 is the single largest term in P6's wall clock. That is a real
property of native GRAM at constellation scale, and it is the thing the figure
reports.

Trace 6's per-sample cost is SOURCED rather than scaled: S2's `process_members`
mode ran 128 one-spacecraft samples per worker in 11.1525 s median at a 600 s
mission, i.e. **87.1 ms per sample** including per-solve fixed cost, so 4096
samples is ~357 s of serial-equivalent work. That is held flat against the
shorter mission rather than scaled down by 6, because at 100 s the per-solve
fixed cost is the dominant term and does not shrink — a conservative choice for
a duration estimate.

**These are per-machine predictions, not constants of the physics**, exactly as
`PPC_L50_ISO_MISSION_S` is. Before a paper run on a host they were not derived
for:

```bash
bash benchmarks/studies/paper_parallelization_benchmarks/paper_figure_runs.sh calibrate-p6 --execute
```

which times each of the six in serial at the full profile and prints the mission
length that would hit the target, the way `calibrate_iso_ladder.sh` does for the
constellation-count ladder. Move a suggestion into
`PPB_P6_{L20,VACUUM,AERO}_MISSION_S` in `cli.jl` and into the case-catalog names
in `parallelization_performance/cases.jl` if it moves a rung by more than ~20%.

### Expected duration

DERIVED, from a per-point cost model calibrated against a measured P-phase run.

SOURCED — the P2 phase of `20260918_162845` took **42m56s** for 25 points at 11
repeats plus 1 warm-up. Summing that run's own 25 medians gives 126.63 s per
repeat-set, so 12 × 126.63 = 1519.6 s is solving and the remaining 1056 s is
fixed per-point cost: **42 s per point** of Julia start-up and JIT on the
benchmark box. (The 12-core reference box's figure is ~80 s; see `CASES.md`.)

Applying `cost = 42 s + 4 × solve` per point, with the predicted baselines above
and trace 2's own measured ladder shape:

| Phase | Points | Estimate |
|---|---:|---:|
| P6 trace 1 (`l20_vacuum`) | 13 | ~14 min |
| P6 trace 2 (`l50_vacuum`) | 13 | ~14 min |
| P6 trace 3 (`srp_nbody`) | 13 | ~14 min |
| P6 trace 4 (`expatm`) | 13 | ~14 min |
| P6 trace 5 (`gram_lookahead`) | 13 | ~72 min |
| **P6 total** | **65** | **~2 h 10 m** |
| P6p `serial` (6 grid entries) | 6 | ~2 h 27 m |
| P6p `outer_process` | 6 | ~55 min |
| P6p `predictive` | 6 | ~55 min |
| **P6p total** | **18** | **~4 h 20 m** |
| **Figure F2 total** | **83** | **~6 h 30 m** |

P6p's `serial` mode is over half of the figure's cost: it is budget-independent,
so those six points are the same measurement six times, but the aggregation joins
the serial baseline on `(phase_id, case, mc_samples, process_workers)` and a
point with no serial row at its own worker count gets no `speedup` column at all.
That is the same trade every P phase makes (see `CASES.md`, "Serial is in every P
phase's mode list").

### Dry runs

Both pasted verbatim. Run on `space-falcon-1` on 2026-09-21 from the
`pv2c-ws9c-paper-runs` worktree.

**P6 and P6p as the benchmark box will run them**
(`SPACEAGORA_PPB_PAPER_BUDGET=32` stands in for the box's core count, which is
what sets P6p's worker ladder):

```
$ SPACEAGORA_PPB_PAPER_BUDGET=32 julia --project=. \
    benchmarks/studies/paper_parallelization_benchmarks.jl \
    --phases=P6,P6p --threads=1,2,4,8,16,32 --process-workers=32 --dry-run

[paper-benchmarks] phases           = P6, P6p
[paper-benchmarks] thread_ladder    = 1,2,4,8,16,32
[paper-benchmarks] process_workers  = 32
[paper-benchmarks] mc_samples_max   = unset
[paper-benchmarks] solver_mode      = auto_stiff
[paper-benchmarks] seed             = 20260615
[paper-benchmarks] dry_run          = true

[phase P6] Paper — Thread Scaling at 4096 Spacecraft, Force-Model and Atmosphere Variants
[dry-run] phase=P6 — single run
[dry-run]   cases           = gravity_4096sat_l20_vacuum_19700s, gravity_4096sat_l50_vacuum_5800s, gravity_4096sat_l50_srp_nbody_vacuum_5800s, aero_4096sat_l50_expatm_100s, aero_4096sat_l50_gram_lookahead_100s
[dry-run]   parity_cases    =
[dry-run]   modes           = serial, inner_only, predictive
[dry-run]   threads         = [1, 2, 4, 8, 16, 32]
[dry-run]   mc_samples      = [1]
[dry-run]   process_workers = 32
[dry-run]   repeats         = 3, warmup = 1

[phase P6p] Paper — Process Route at 4096 Spacecraft, Native GRAM
[dry-run] phase=P6p — workers=1 x threads=1 (per-route budget 1)
[dry-run]   cases           = aero_4096sat_l50_gram_process_100s
[dry-run]   modes           = serial, outer_process, predictive
[dry-run]   threads         = [1]
[dry-run]   mc_samples      = [4096]
[dry-run]   process_workers = 1
[dry-run]   repeats         = 3, warmup = 1
[dry-run] phase=P6p — workers=2 x threads=1 (per-route budget 2)     ... process_workers = 2
[dry-run] phase=P6p — workers=4 x threads=1 (per-route budget 4)     ... process_workers = 4
[dry-run] phase=P6p — workers=8 x threads=1 (per-route budget 8)     ... process_workers = 8
[dry-run] phase=P6p — workers=16 x threads=1 (per-route budget 16)   ... process_workers = 16
[dry-run] phase=P6p — workers=32 x threads=1 (per-route budget 32)   ... process_workers = 32

Phase Status   Runs  Elapsed     Label
P6    ok       1     0s          Paper — Thread Scaling at 4096 Spacecraft, Force-Model and Atmosphere Variants
P6p   ok       6     0s          Paper — Process Route at 4096 Spacecraft, Native GRAM
  2/2 phases completed successfully.
```

(The six P6p entries differ only in `process_workers`; the repeated lines are
elided above and are in full in the run's own output.)

**P6 and P6p on this 12-core workstation**, which is what a smoke of the phase
structure here will show — the worker ladder rescales to the host and stops at
12, as every P phase's host-sized grid does:

```
$ julia --project=. benchmarks/studies/paper_parallelization_benchmarks.jl \
    --phases=P6,P6p --dry-run

[paper-benchmarks] thread_ladder    = auto (24 CPU threads)
[paper-benchmarks] process_workers  = 12
[dry-run] phase=P6 — single run
[dry-run]   threads         = [1, 2, 4, 6, 9, 12]
[dry-run]   modes           = serial, inner_only, predictive
[dry-run]   repeats         = 3, warmup = 1
[dry-run] phase=P6p — workers=1 x threads=1 (per-route budget 1)
[dry-run] phase=P6p — workers=2 x threads=1 (per-route budget 2)
[dry-run] phase=P6p — workers=4 x threads=1 (per-route budget 4)
[dry-run] phase=P6p — workers=8 x threads=1 (per-route budget 8)
[dry-run] phase=P6p — workers=12 x threads=1 (per-route budget 12)
[dry-run]   mc_samples      = [4096]
```

### Smoke

Every new case run once at N = 16 under `--profile=test` on `space-falcon-1`
(which pins the mission at 10 s regardless of the duration in the name), one
Julia process at a time, `timeout 600`:

| Case | Mode | `success` | `wall_time_s` | `retcode` |
|---|---|---|---:|---|
| `gravity_16sat_l20_vacuum_19700s` | serial | true | 0.206 | Success |
| `gravity_16sat_l50_srp_nbody_vacuum_5800s` | serial | true | 0.211 | Success |
| `aero_16sat_l50_expatm_100s` | serial | true | 0.123 | Success |
| `aero_16sat_l50_gram_lookahead_100s` | serial | true | 0.195 | Success |
| `aero_16sat_l50_gram_process_100s` | serial, 16 samples | true | 0.732 | Success |
| `aero_16sat_l50_gram_process_100s` | outer_process, 16 samples, 2 workers | true | 1.446 | Success |

The process-route row is the one that had to be checked rather than assumed: its
CSV records `outer_backend_actual = process` and `outer_tasks = 16`, i.e. the
pool really was provisioned and the sixteen one-spacecraft samples really were
dispatched across it, with each `--threads=1` worker loading its own native GRAM
image. The `serial` row of the same case records `outer_backend_actual = serial`
and the same sixteen `outer_tasks`, which is the baseline that trace pairs with.

The density-path selector was checked separately against the four argument
shapes it has to distinguish:

```
--case=aero_4096sat_l50_gram_lookahead_100s -> FREEZE_PER_STEP=0 VACUUM_GRAM_CACHE=1 HORIZON_S=600.0 DEVIATION_M=1e8
--case=aero_4096sat_l50_gram_process_100s   -> FREEZE_PER_STEP=1 VACUUM_GRAM_CACHE=0 (horizon, deviation unset)
--case=gravity_4096sat_l50_vacuum_5800s     -> all unset
--phases=P6 --dry-run                       -> all unset
```

i.e. the look-ahead horizon is the mission length plus 500 s so no rebuild can
occur mid-run, the process trace gets the freeze instead, and neither a
non-GRAM case nor the controller process (which has no `--case=` in `ARGS`) is
touched.

---

## Figure F4 — the workstation arm

F4 is cross-machine portability: the same two Monte Carlo phases, P3 and P5, on
the benchmark box and on this 12-core workstation, every mode including
`predictive` (R7) and `policy_v2` (R6), 11 repeats, cold calibration store.

```bash
bash benchmarks/studies/paper_parallelization_benchmarks/paper_figure_runs.sh workstation            # describe
bash benchmarks/studies/paper_parallelization_benchmarks/paper_figure_runs.sh workstation --execute  # run it
```

What that does, in order: refuse if the machine is busy; run
`scripts/calibrate_machine.jl` at 12 threads so the predictive planner has this
machine's cost constants; move `output/parallel_policy_state/` aside to
`output/parallel_policy_state_backup_<UTC stamp>` (moved, never deleted — a
converged store is hours of measurement and is the input to the converged arm);
re-run the calibration into the now-empty store, because emptying it takes the
constants with it; then

```bash
SPACEAGORA_PPB_MIN_REPEATS=11 OPENBLAS_NUM_THREADS=1 GKSwstype=100 \
  julia --project=. benchmarks/studies/paper_parallelization_benchmarks.jl \
    --phases=P3,P5 --threads=1,2,4,8,12 --process-workers=12
```

and finally the archive call. **Eleven repeats** because five put the 90%
interval of an adaptive point's median at ~7%, the same width as the band used to
decide whether two runs differ at all; eleven takes it to ~2.2% and it plateaus
(`_ppb_min_repeats`, `main.jl`). A cross-machine figure is exactly the comparison
that needs the tighter interval.

**The busy guard.** It mirrors the harness's own `ppc_assert_machine_quiet!`
(50% load headroom) but runs in the launcher, so the refusal happens before the
store is rotated and the calibration spent rather than forty minutes into a
phase. It also refuses if any Julia process other than an editor language server
is running: two benchmark jobs on one box do not merely add noise to each other,
they make both sets of numbers unusable and the machine unusable for its owner.

### Dry run

Verbatim, `space-falcon-1`, 2026-09-21. The rungs are what the ppb budget logic
produces on 12 physical cores — P3's `(b, b)` resource ladder and P5's six splits
of 12 — and the repeats floor is visible rewriting both phases from 5 to 11:

```
$ SPACEAGORA_PPB_MIN_REPEATS=11 julia --project=. \
    benchmarks/studies/paper_parallelization_benchmarks.jl \
    --phases=P3,P5 --threads=1,2,4,8,12 --process-workers=12 --dry-run

[paper-benchmarks] phase P3: repeats 5 -> 11
[paper-benchmarks] phase P5: repeats 5 -> 11
[paper-benchmarks] phases           = P3, P5
[paper-benchmarks] thread_ladder    = 1,2,4,8,12
[paper-benchmarks] process_workers  = 12

[phase P3] Paper — Monte Carlo Resource Ladder, One Spacecraft per Sample
[dry-run] phase=P3 — workers=1 x threads=1 (per-route budget 1)
[dry-run]   cases           = independent_1sat_1hr
[dry-run]   modes           = serial, outer_threads, outer_process, policy_v2, predictive
[dry-run]   threads         = [1]
[dry-run]   mc_samples      = [256]
[dry-run]   process_workers = 1
[dry-run]   repeats         = 11, warmup = 1
[dry-run] phase=P3 — workers=2 x threads=2   (per-route budget 2)
[dry-run] phase=P3 — workers=4 x threads=4   (per-route budget 4)
[dry-run] phase=P3 — workers=8 x threads=8   (per-route budget 8)
[dry-run] phase=P3 — workers=12 x threads=12 (per-route budget 12)

[phase P5] Paper — Monte Carlo over Constellations, Worker/Thread Split at a Fixed Budget
[dry-run] phase=P5 — workers=1 x threads=12  (per-route budget 1)
[dry-run]   cases           = mcgrid_16sat_8mc, mcgrid_8sat_16mc
[dry-run]   modes           = serial, outer_threads, outer_process, outer_inner_static, policy_v2, predictive
[dry-run]   mc_samples      = [1]
[dry-run]   repeats         = 11, warmup = 1
[dry-run] phase=P5 — workers=2 x threads=6   (per-route budget 2)
[dry-run] phase=P5 — workers=3 x threads=4   (per-route budget 3)
[dry-run] phase=P5 — workers=4 x threads=3   (per-route budget 4)
[dry-run] phase=P5 — workers=6 x threads=2   (per-route budget 6)
[dry-run] phase=P5 — workers=12 x threads=1  (per-route budget 12)

Phase Status   Runs  Elapsed     Label
P3    ok       5     0s          Paper — Monte Carlo Resource Ladder, One Spacecraft per Sample
P5    ok       6     0s          Paper — Monte Carlo over Constellations, Worker/Thread Split at a Fixed Budget
  2/2 phases completed successfully.
```

### Expected duration

SOURCED — P5 on this machine took **1h53m29s** at 5 repeats
(`.claude/worktrees/paper-routing/output/performance/paper_benchmarks/20260915_181642/paper_benchmarks_report_20260915_181642.md`,
`Machine: space-falcon-1`). DERIVED — at 11 repeats the timed portion roughly
doubles (12 solves per point against 6) while the per-point start-up does not, so
**3–4 h**. P3 has never been run at the full profile on this machine — only under
`--preview` (9m14s, same report) — so its duration here is UNMEASURED; the
benchmark box did it in 1h3m43s at 11 repeats with 32 cores, which bounds this
box's figure from below. Budget **5–7 h for the workstation arm**, i.e. one
overnight slot, and do not schedule anything else on the machine during it.

---

## Precompile workload

Every point is a fresh Julia process, and a campaign point's pool workers are
too; each compiles its solver specializations before it can time anything.
`workload/` builds a package image holding them for every P1-P6p point except
the two native-GRAM traces, and the harness loads it into every worker when
asked. It is off by default: with the image loaded the timed repeats measured 1 to 4
percent slower on four of five points (table below), so it is off by default, and the default for
benchmark timing is the configuration the archived runs used. `workload/README.md` has the details; in short:

- **Build**: `bash benchmarks/studies/paper_parallelization_benchmarks/workload/build_workload.sh`
  (3 m 40 s and 9.3 GB peak on the workstation, MEASURED). Rebuild after any
  change to `src/` or to the harness files; a stale image is never used.
- **Use**: opt-in. `SPACEAGORA_PPB_WORKLOAD=auto` uses a current image and
  falls back without one; `SPACEAGORA_PPB_WORKLOAD=1` makes a missing or stale
  image an error. The controller prints `precompile_workload=<env>` or
  `precompile_workload=off` at the start of each run. Remote jobs build it once
  per job only when asked (`spaceagora-remote push --workload auto|on|off`,
  default `off`).
- **Default**: off (`SPACEAGORA_PPB_WORKLOAD=0`). Rows measured with and without
  the image are not comparable and must not share a figure.

Measured on the workstation (MEASURED, medians over three alternating launches
per variant, `workload/results/space-falcon-1_20260924/workload_validation.csv`):

| Point (threads x workers) | Worker wall, stock / workload | Timed median shift |
|---|---|---|
| P1 1 spacecraft, serial, 8 threads | 76.5 / 47.1 s (1.6x) | +0.7% |
| P2 4096 spacecraft, inner_only, 8 threads | 49.2 / 25.3 s (1.9x) | +4.6% (launch spread -9.5% to +14.1%) |
| P3 outer_process, 256 samples, 4 x 4 | 151.8 / 42.9 s (3.5x) | +3.6% |
| P3 policy_v2, 256 samples, 4 x 4 | 150.1 / 41.0 s (3.7x) | +3.9% |
| P5 mcgrid_16sat_8mc, outer_inner_static, 4 x 2 | 78.3 / 28.8 s (2.7x) | +1.9% |

Final states and step times are byte-identical with and without the image. The
timed repeats run 1-4% slower with it (same allocations, same GC; the time is in
the samples' compute), so **a run's rows are comparable only with rows measured
the same way**: do not compare a run that used the workload against one that did
not, including the archived runs cited below, which did not. Set
`SPACEAGORA_PPB_WORKLOAD=0` for a run that must line up with them. The
durations quoted below were measured without it.

---

## The benchmark-box sequence

```bash
bash benchmarks/studies/paper_parallelization_benchmarks/paper_figure_runs.sh trx50            # print the sequence
bash benchmarks/studies/paper_parallelization_benchmarks/paper_figure_runs.sh trx50 --execute  # offer each step in turn
```

Six steps, ordered, non-overlapping. The ordering is load-bearing: the
calibration is itself a timed measurement, the targeted points gate the
eleven-hour run, and the converged arm consumes the store the cold arm produced.

Three facts about the interface that are easy to get wrong, all confirmed by
reading `scripts/remote/spaceagora-remote`:

- **Flags take a space, not an `=`.** `--threads 1,2,4,8,16,32` parses;
  `--threads=1,2,4,8,16,32` hits the parser's catch-all and aborts with
  `error: unexpected push argument`.
- **The remainder after `--` is exec'd directly, not through a shell.** It is
  shell-quoted once and pasted into the generated `run.sh`, so `&&`, `;`, `$VAR`
  and redirects do not work there — use `env K=V ... cmd` for environment, and
  `bash -c '...'` when a shell is genuinely needed.
- **The calibration store is not shipped and not pulled.** `sync_code` excludes
  `/output/` outright and filters by `.gitignore` besides. On the remote the
  store is a persistent per-remote directory, `~/spaceagora_remote/policy_state`,
  symlinked into each release's `output/parallel_policy_state`. Its state is
  therefore set over `ssh`, out of band, before the job that depends on it —
  which is what the `ssh ... && spaceagora-remote push ...` shape below is doing.

### Step 1 — machine calibration

```bash
scripts/remote/spaceagora-remote push --remote trx50 --threads 32 \
  -- julia --project=. --threads=32 scripts/calibrate_machine.jl
```

At the thread count the jobs will use: the dispatch and atomic constants are
properties of the inner dispatch primitive at that width
(`scripts/calibrate_machine.jl` header). Duration UNMEASURED — no archived run on
disk records it. It writes
`output/parallel_policy_state/cost_constants_<fingerprint>.toml`, which
`targeted_points.sh` warns about the absence of; confirm the file exists before
step 2.

### Step 2 — targeted points, one job per POINT

`targeted_points.sh` documents the store state each point needs. `finding8` is
the converged-store point; the other three are cold.

```bash
# finding8 — converged store
ssh trx50 'rm -rf ~/spaceagora_remote/policy_state && cp -a ~/spaceagora_remote/policy_state_converged_20260918 ~/spaceagora_remote/policy_state' \
  && scripts/remote/spaceagora-remote push --remote trx50 --threads 1,2,4,8,16,32 --process-workers 32 \
       -- env POINT=finding8 bash benchmarks/studies/paper_parallelization_benchmarks/targeted_points.sh

# p5_16sat, finding9, defectA — cold store (repeat for each POINT)
ssh trx50 'mv ~/spaceagora_remote/policy_state ~/spaceagora_remote/policy_state_bak_$(date -u +%Y%m%d_%H%M%S); mkdir -p ~/spaceagora_remote/policy_state' \
  && scripts/remote/spaceagora-remote push --remote trx50 --threads 1,2,4,8,16,32 --process-workers 32 \
       -- env POINT=p5_16sat bash benchmarks/studies/paper_parallelization_benchmarks/targeted_points.sh
```

Duration UNMEASURED — `targeted_points.sh` has never been run on the box, so no
archived timing exists. Its own bound is three runs per POINT under `timeout
7200`, i.e. at most 6 h per POINT; the three points it measures are single
`(case, mode)` points at 11 repeats, so the realistic figure is far below that
cap and should be observed on the first job rather than predicted here.

**`policy_state_converged_20260918` does not exist on this machine.** A search of
the whole repository (including `remote_runs/`, `output/` and every worktree)
finds no file or directory by that name and no tarball that could contain one;
the only stores present are plain `parallel_policy_state` directories. It may
exist on the remote — this was not checked, because the remote was not touched in
this workstream. **Confirm `ssh trx50 'ls -d ~/spaceagora_remote/policy_state_converged_20260918'`
before running step 2's `finding8` or step 4 at all**; both commands above
*delete* the live store before copying, and if the snapshot is absent they
destroy the store and leave nothing in its place.

### Step 3 — the cold 11-repeat P1–P5 run

```bash
ssh trx50 'mv ~/spaceagora_remote/policy_state ~/spaceagora_remote/policy_state_bak_$(date -u +%Y%m%d_%H%M%S); mkdir -p ~/spaceagora_remote/policy_state' \
  && scripts/remote/spaceagora-remote push --remote trx50 --threads 1,2,4,8,16,32 --process-workers 32 \
       -- env SPACEAGORA_PPB_MIN_REPEATS=11 julia --project=. \
          benchmarks/studies/paper_parallelization_benchmarks.jl \
          --phases=P1,P2,P3,P4,P5 --threads=1,2,4,8,16,32 --process-workers=32

```

**Duration: 11h49m.** SOURCED — the per-phase elapsed table of the last such run,
`.claude/worktrees/paper-routing/output/performance/paper_benchmarks_trx50_cold11/20260918_162845/paper_benchmarks_report_20260918_162845.md`,
`Machine: space-falcon-lab-TRX50-AERO-D`:

| Phase | Elapsed |
|---|---|
| P1 Constellation Size Scaling at a Fixed Budget | 1h27m3s |
| P2 Thread Scaling at 4096 Spacecraft | 42m56s |
| P3 Monte Carlo Resource Ladder, One Spacecraft per Sample | 1h3m43s |
| P4 Monte Carlo Resource Ladder, Compute-Bound Samples | 1h30m54s |
| P5 Monte Carlo over Constellations, Worker/Thread Split | 6h44m20s |
| **sum** | **11h48m56s** |

Corroborated by that run's raw CSV `timestamp_utc` span, 2026-09-18T16:29:02 to
2026-09-19T03:55:07 (11h26m) — the difference is the final aggregation and
plotting, which happens after the last phase reports.

This run's `predictive` arm is the one addition over the archived one, which
carried `serial, outer_threads, inner_only, outer_inner_static, policy_v2`. R7
adds one mode to each phase's ladder, so budget roughly **+20–25%, i.e. 14–15 h**.

### Step 4 — the converged P1/P5 arm

```bash
ssh trx50 'rm -rf ~/spaceagora_remote/policy_state && cp -a ~/spaceagora_remote/policy_state_converged_20260918 ~/spaceagora_remote/policy_state' \
  && scripts/remote/spaceagora-remote push --remote trx50 --threads 1,2,4,8,16,32 --process-workers 32 \
       -- env SPACEAGORA_PPB_MIN_REPEATS=11 julia --project=. \
          benchmarks/studies/paper_parallelization_benchmarks.jl \
          --phases=P1,P5 --threads=1,2,4,8,16,32 --process-workers=32
```

**Duration: 8h11m.** DERIVED — P1 (1h27m3s) plus P5 (6h44m20s) from the same
measured table; the converged arm re-runs the same points, and only the adaptive
modes' own timings move. See the warning in step 2 about the snapshot's existence
before running the `rm -rf`.

### Step 5 — P6 and P6p (figure F2)

If the box has not run P6 before, recalibrate its mission lengths first, as its
own job:

```bash
scripts/remote/spaceagora-remote push --remote trx50 --threads 1 \
  -- bash benchmarks/studies/paper_parallelization_benchmarks/paper_figure_runs.sh calibrate-p6 --execute
```

Then:

```bash
ssh trx50 'mv ~/spaceagora_remote/policy_state ~/spaceagora_remote/policy_state_bak_$(date -u +%Y%m%d_%H%M%S); mkdir -p ~/spaceagora_remote/policy_state' \
  && scripts/remote/spaceagora-remote push --remote trx50 --threads 1,2,4,8,16,32 --process-workers 32 \
       -- julia --project=. benchmarks/studies/paper_parallelization_benchmarks.jl \
          --phases=P6,P6p --threads=1,2,4,8,16,32 --process-workers=32
```

**Duration: ~6h30m**, DERIVED — see "Expected duration" under F2 above, from the
42 s-per-point fixed cost measured on this box in run `20260918_162845` and the
predicted per-trace serial baselines. The calibration job is ~20 min on the same
model (six serial solves plus start-up, trace 5 dominating).

Watch P6p's widest rung: 32 workers × 5.4 GB is 174 GB of the box's 250 GB.
Do not raise `--process-workers` above 32 for this phase.

### Step 6 — pull, then archive

`spaceagora-remote`'s auto-pull takes a release's top-level `output/` and
`results/` only, and the watcher runs on the *launching* machine — if it slept or
rebooted while the job ran, nothing was pulled. Name the path either way:

```bash
rsync -az trx50:~/spaceagora_remote/releases/<job-id>/output/performance/paper_benchmarks/ \
  output/performance/paper_benchmarks/
```

Then archive each run into the benchmark-data directory, which is a plain local
directory with no git remote — nothing in it goes into the paper repository:

```bash
python3 scripts/archive_paper_run.py output/performance/paper_benchmarks/<stamp> \
  --archive "${SPACEAGORA_PAPER_ARCHIVE:-../SpaceAGORA-paper-data/data/raw}" \
  --machine trx50 --store cold      --notes 'P1-P5, 11 repeats, cold store'

python3 scripts/archive_paper_run.py output/performance/paper_benchmarks/<stamp> \
  --archive "${SPACEAGORA_PAPER_ARCHIVE:-../SpaceAGORA-paper-data/data/raw}" \
  --machine trx50 --store converged --notes 'P1+P5, 11 repeats, converged store'

python3 scripts/archive_paper_run.py output/performance/paper_benchmarks/<stamp> \
  --archive "${SPACEAGORA_PAPER_ARCHIVE:-../SpaceAGORA-paper-data/data/raw}" \
  --machine trx50 --store cold      --notes 'P6+P6p, 3 repeats, figure F2'

python3 scripts/archive_paper_run.py --verify \
  --archive "${SPACEAGORA_PAPER_ARCHIVE:-../SpaceAGORA-paper-data/data/raw}"
```

`--store` is required for ppb runs and is the one field the CSV cannot supply:
only the person who set the store knows which state the job ran under. The
workstation arm passes `--machine workstation --store cold`.

A benchmark run lands in a gitignored `output/` directory inside a worktree, so
a reboot, a `git clean` or a deleted worktree takes it with it. Archive each run
as soon as it is pulled, not at the end of the campaign.

---

## Provenance ledger

Every number in this document, and what it rests on.

| Quantity | Value | Status | Source |
|---|---|---|---|
| Trace 2 serial baseline | 11.8566 s | SOURCED | `paper_benchmarks_trx50_cold11/20260918_162845/paper_benchmarks_aggregated_20260918_162845.csv`, P2 / serial / t1 |
| Trace 2 thread ladder | 11.80 / 6.77 / 3.63 / 2.19 / 1.99 / 2.24 s | SOURCED | same CSV, `inner_only` rows |
| P1–P5 phase elapsed | 1h27m / 42m56s / 1h3m43s / 1h30m54s / 6h44m20s | SOURCED | `.../20260918_162845/paper_benchmarks_report_20260918_162845.md` |
| Per-point fixed cost, benchmark box | 42 s | DERIVED | P2's 42m56s minus 12 × the sum of its own 25 medians, ÷ 25 |
| `atmo256_exponential_10min` serial | 4.8538 s | SOURCED | `trx50_quiet_20260829/paper_benchmarks_aggregated_20260829_195452.csv` |
| `atmo256_gram_live_10min` serial | 42.0565 s | SOURCED | same |
| `atmo256_gram_live_nbody_10min` serial | 46.3458 s | SOURCED | same |
| Third-body cost | +10.2% | DERIVED | ratio of the two above |
| Degree-50 / degree-20 serial ratio at 4096 | 3.3938 | SOURCED | `s1_constellation_scaling_l50_none.csv` (210.4728 s) / `s1_constellation_scaling.csv` (62.0167 s) |
| …batched-kernel bracket | 5.7403 | DERIVED | coefficient counts (51·52/2)/(21·22/2) |
| Trace 1 mission | 19 700 s | DERIVED | 5800 × 3.3938, rounded |
| Trace 4/5/6 mission | 100 s | ASSUMED | rounded up from the exponential trace's derived 91.6 s; see the argument above |
| Single-process GRAM RSS | 1.7 GB + 1.0 MB/spacecraft | DERIVED | OLS over `s2_gram_atmosphere_modes.csv`, `maxrss_mb` |
| …confirmed at 4096 | 5410.5 / 5469.2 MB | SOURCED | `s1_constellation_scaling_l50_none.csv`, `n_sats = 4096` |
| Per-worker GRAM RSS at 128 samples/worker | 5403.2 MB | SOURCED | `s2_gram_atmosphere_modes.csv`, `process_members`, `n_sats = 1024`, `workers_rss_mb` ÷ 8 |
| P6p worst-rung memory | 174.1 GB | DERIVED | 32 × 5403.2 MB + 1224.8 MB coordinator |
| Benchmark box installed memory | 250 GB | SOURCED | `docs/adaptive_policy_r6_nested_campaign_record_20260906.md:656` |
| Per-sample process-route cost | 87.1 ms | SOURCED | `s2_gram_atmosphere_modes.csv`, `process_members`, `n_sats = 1024`: 11.1525 s ÷ 128 |
| P6 / P6p duration | ~2h10m / ~4h20m | DERIVED | 42 s per point + 4 × the predicted solves |
| P5 on the workstation | 1h53m29s at 5 repeats | SOURCED | `paper_benchmarks/20260915_181642/paper_benchmarks_report_20260915_181642.md` |
| Workstation arm duration | 5–7 h | DERIVED | P5 doubled for 11 repeats, plus an unmeasured P3 bounded below by the box's 1h3m43s |

## What is not settled

- **`policy_state_converged_20260918` was not found** anywhere in this
  repository. Steps 2 and 4 depend on it and both delete the live store before
  copying it. Verify it exists on the remote first.
- **Step 1 and step 2 durations are unmeasured.** Neither
  `scripts/calibrate_machine.jl` nor `targeted_points.sh` has an archived run on
  disk. Observe them on their first job rather than planning around a guess.
- **P6's mission lengths are predictions for one machine.** Traces 1 and 3 in
  particular rest on ratios taken from a different harness configuration (S1's
  unbatched serial kernel) or from a different spacecraft count (the `atmo256_*`
  ladder). `calibrate-p6` settles both, and should be run on the box before the
  figure is.
- **`_GRAM_SAT_MEMORY_BYTES` (90 MB per spacecraft) disagrees with the measured
  ~1.0 MB per spacecraft by about 46×** at 256 spacecraft. Recorded, not
  resolved; it does not affect this figure, because P6p's worker count is set by
  the phase rather than by the router's memory model. It would affect any run
  that lets the router size a native-GRAM process pool.
- **Trace 5 is expected to be slower than serial at every thread count.** S2
  measured `threads_lookahead` at 80.19 s against `serial_standard` at 29.91 s on
  1024 spacecraft, and the archived `atmo256_gram_live_10min` rows top out at
  2.07× on 12 threads. That is the figure's result, not a defect in the setup,
  but it means trace 5's wall time does not fall with the thread count the way
  traces 1–4's do and the phase's cost estimate assumes it will not.
