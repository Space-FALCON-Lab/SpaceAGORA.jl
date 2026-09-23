# The default RHS execution plan for a constellation

> Status: measurement record and the derivation of the shipped defaults.
> Measurement: `benchmarks/studies/rhs_heuristic_defaults/`.
> Code: `_rhs_execution_plan_uncached` in `src/simulation/engine/setup.jl`.

## What was wrong

A solve picks its RHS execution plan one of two ways. With calibration enabled
(`SPACEAGORA_RHS_CALIBRATE=auto`, which the adaptive profiles set) a short
pre-solve sweep times a candidate set and pins the fastest plan for the solve.
With calibration off -- which is what every pinned static profile runs, and what
a user who has not opted into an adaptive profile gets -- the plan comes from
the heuristic in `_rhs_execution_plan_uncached`, re-derived on every RHS call.

For a constellation whose only effector is a gravitational-harmonics model, the
heuristic took the flat constellation route and sized its worker team as
`min(thread_budget, fld(active_sats, 4))`. The 4 there is
`SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER`, whose own docstring derives it
from SIMD register width: "one AVX2 / half AVX-512 register width -- enough for
`@turbo` to have a full SIMD iteration". That is a correct floor for *sizing a
pre-pass slice* and it is not a floor for *deciding whether to open a worker
team at all*, which is the question the routing was using it to answer.

## The diagnosis

Phase P1 of the archived TRX50 cold run
`trx50_ppb_cold_20260922_002429` (11 repeats, 32 Julia threads, the iso-work
L50 vacuum constellation ladder) carries `rhs_plan_source`, `rhs_plan_mode`,
`rhs_plan_allotment` and `rhs_plan_scheduler` per row. Tabulated per size and
per mode:

| N | mode (profile) | calibrate | plan source | plan the solve ran | median wall s | reps |
|---|---|---|---|---|---|---|
| 1 | serial (R0) | off | none | heuristic (no telemetry) | 8.92 | 11 |
| 1 | inner_only (R2) | off | none | heuristic (no telemetry) | 8.64 | 11 |
| 1 | outer_threads (R1a) | off | none | heuristic (no telemetry) | 8.62 | 11 |
| 1 | outer_inner_static (R3) | off | none | heuristic (no telemetry) | 8.67 | 11 |
| 1 | policy_v2 (R6) | auto | none | heuristic (no telemetry) | 8.68 | 11 |
| 1 | predictive (R7) | auto | none | heuristic (no telemetry) | 8.57 | 11 |
| 16 | serial (R0) | off | none | heuristic (no telemetry) | 11.36 | 11 |
| 16 | inner_only (R2) | off | none | heuristic (no telemetry) | 2.13 | 11 |
| 16 | outer_threads (R1a) | off | none | heuristic (no telemetry) | 2.17 | 11 |
| 16 | outer_inner_static (R3) | off | none | heuristic (no telemetry) | 2.14 | 11 |
| 16 | policy_v2 (R6) | auto | cache | heuristic retained | 2.18 | 11 |
| 16 | predictive (R7) | auto | cache,sweep | heuristic retained | 2.09 | 11 |
| 64 | serial (R0) | off | none | heuristic (no telemetry) | 11.22 | 11 |
| 64 | inner_only (R2) | off | none | heuristic (no telemetry) | 35.13 | 11 |
| 64 | outer_threads (R1a) | off | none | heuristic (no telemetry) | 34.74 | 11 |
| 64 | outer_inner_static (R3) | off | none | heuristic (no telemetry) | 34.42 | 11 |
| 64 | policy_v2 (R6) | auto | cache,sweep | satellite_batch@1/static | 3.31 | 11 |
| 64 | predictive (R7) | auto | cache,sweep | heuristic retained; satellite_batch@1/static; satellite_batch@16/static | 3.30 | 11 |
| 256 | serial (R0) | off | none | heuristic (no telemetry) | 11.96 | 11 |
| 256 | inner_only (R2) | off | none | heuristic (no telemetry) | 21.16 | 11 |
| 256 | outer_threads (R1a) | off | none | heuristic (no telemetry) | 21.10 | 11 |
| 256 | outer_inner_static (R3) | off | none | heuristic (no telemetry) | 21.08 | 11 |
| 256 | policy_v2 (R6) | auto | cache,sweep | satellite_batch@1/static | 2.73 | 11 |
| 256 | predictive (R7) | auto | cache,sweep | satellite_batch@1/static | 2.73 | 11 |
| 1024 | serial (R0) | off | none | heuristic (no telemetry) | 12.15 | 11 |
| 1024 | inner_only (R2) | off | none | heuristic (no telemetry) | 5.00 | 11 |
| 1024 | outer_threads (R1a) | off | none | heuristic (no telemetry) | 4.94 | 11 |
| 1024 | outer_inner_static (R3) | off | none | heuristic (no telemetry) | 4.97 | 11 |
| 1024 | policy_v2 (R6) | auto | cache,sweep | satellite_batch@1/static | 1.97 | 11 |
| 1024 | predictive (R7) | auto | cache,sweep | heuristic retained; satellite_batch@1/static | 2.04 | 11 |
| 4096 | serial (R0) | off | none | heuristic (no telemetry) | 11.66 | 11 |
| 4096 | inner_only (R2) | off | none | heuristic (no telemetry) | 2.78 | 11 |
| 4096 | outer_threads (R1a) | off | none | heuristic (no telemetry) | 2.15 | 11 |
| 4096 | outer_inner_static (R3) | off | none | heuristic (no telemetry) | 2.21 | 11 |
| 4096 | policy_v2 (R6) | auto | cache | heuristic retained | 2.17 | 11 |
| 4096 | predictive (R7) | auto | cache,sweep | heuristic retained | 2.17 | 11 |

The pinned routes (`R1a`, `R2`, `R3`) and the serial baseline (`R0`) all run
`SPACEAGORA_RHS_CALIBRATE=off`, so they record no plan telemetry and ran the
heuristic's plan. The adaptive routes (`R6`, `R7`) ran
`SPACEAGORA_RHS_CALIBRATE=auto` and recorded what the sweep chose. Reading the
two against each other:

- At 1 and 16 spacecraft the heuristic does not take the flat route at all
  (`active_sats < SPACEAGORA_EFFECTOR_FLAT_MIN_SATS`, default 24), and the two
  families agree to within the repeat spread.
- At 64, 256 and 1024 spacecraft the heuristic took
  `flat_constellation_effector_queue` at 16, 32 and 32 workers, and the sweep
  pinned `satellite_batch`. The pinned routes ran 10.6x, 7.7x and 2.5x slower
  than the adaptive ones. At 64 spacecraft the default was also 3.1x slower than
  the *serial* baseline.
- At 4096 spacecraft the sweep retained the heuristic and the two families agree
  again, so the flat route does win at the top of the ladder.

The same shape, at that machine's own crossover, is in the archived workstation
runs at 12 threads (`workstation_ppb_cold_20260915_181642`, 3 repeats): only the
64-spacecraft point is off, at 14.1 s pinned against 5.73 s adaptive, and the
sweep retains the heuristic from 256 spacecraft up.

## The mechanism

The flat route pays a fixed per-RHS-call cost: it dispatches a worker team for
the harmonics pre-pass (`_accumulate_harmonics_flat_batch!`), and the team is
woken through a persistent channel pool whose own dispatch is documented at
~1-5 microseconds per call (`threaded_foreach_worker_persistent`,
`src/parallel/policy/thread_execution.jl`). The satellite batch route pays no
such cost: its Polyester `@batch` over the spacecraft is the only split there
is.

The P1 ladder is iso-work, so the number of derivative evaluations falls
roughly as 1/N while the work inside each one stays constant. The flat route's
fixed per-call cost is therefore charged to a shrinking number of calls as N
grows, which is exactly the shape the archive shows: the default's penalty is
largest at the small end of the constellation ladder and disappears at the
large end.

## The measurement

`benchmarks/studies/rhs_heuristic_defaults/plan_ladder.jl`, one process, 8 Julia
threads, 3 repeats, every plan the sweep would consider at that budget plus the
heuristic's own, all pinned through the production calibration cache. Case: the
P1 iso-work L50 vacuum ladder. Medians, in seconds
(`results/plan_ladder_vacuum_space-falcon-1_t8.csv`):

| plan | 16 | 32 | 64 | 128 | 256 | 512 | 1024 | 4096 |
|---|---|---|---|---|---|---|---|---|
| heuristic (= flat@8/dynamic from 32 up) | 3.13 | 11.62 | 11.03 | 9.04 | 6.31 | 4.12 | 3.24 | 3.19 |
| satellite_batch@1 | **2.82** | **4.20** | **7.09** | **7.28** | 7.44 | 6.37 | 5.91 | 8.49 |
| satellite_batch@4 | 4.01 | 6.23 | 10.84 | 11.79 | 12.58 | 10.88 | 10.77 | 15.91 |
| satellite_batch@2 | 6.62 | 10.99 | 20.14 | 22.23 | 22.54 | 20.68 | 21.62 | 30.46 |
| flat@1/static | 5.08 | 7.73 | 12.06 | 11.96 | 11.73 | 10.50 | 11.79 | 11.92 |
| flat@2/static | 10.09 | 10.87 | 13.03 | 12.43 | 10.55 | 7.69 | 6.66 | 7.54 |
| flat@2/dynamic | 9.71 | 10.83 | 12.96 | 12.22 | 10.26 | 7.66 | 6.88 | 7.13 |
| flat@4/static | 11.37 | 9.57 | 10.75 | 9.62 | 7.13 | 5.32 | 4.42 | 4.70 |
| flat@4/dynamic | 11.22 | 10.20 | 10.93 | 9.53 | 7.08 | 5.51 | 4.22 | 4.39 |
| flat@8/static | - | 11.53 | 11.00 | 9.01 | **6.36** | **4.02** | **3.20** | **2.84** |
| flat@8/dynamic | - | 11.93 | 11.18 | 9.12 | 6.36 | 4.08 | 3.22 | 2.87 |

Three things the table says.

The satellite batch is the fastest plan through 128 spacecraft and the widest
flat plan is the fastest from 256 up, so there is a crossover and it sits
between 16 and 32 satellites per worker at this budget. The default is on the
wrong side of it from 32 through 128 spacecraft, by 2.77x, 1.56x and 1.24x.

Narrowing the flat route does not help. At 128 spacecraft flat@8 is 1.24x
behind the best plan, flat@4 is 1.32x, flat@2 is 1.68x and flat@1 is 1.64x --
the route is worse the *narrower* it is, everywhere on the ladder. Its cost is
therefore a fixed charge per RHS call rather than a per-worker charge, and the
lever that matters is the route, not its width. A width floor that admitted
flat@2 at a size where flat@8 loses would make the default worse, not better.

The heuristic costs what its plan costs. Where the heuristic takes flat@8 its
median tracks the pinned flat@8 row to within the repeat spread (9.04 against
9.01 at 128, 6.31 against 6.36 at 256), so re-deriving the plan on every RHS
call -- which a pinned plan skips entirely, returning at the top of
`_rhs_execution_plan_uncached` -- is not a measurable part of the gap at these
sizes.


### The two-effector stack

The same ladder over L50 harmonics plus `AerodynamicCoefficientfM` in an
analytic exponential atmosphere (`--case=aero`, 3 repeats, missions 30000 /
10000 / 4000 s at 64 / 256 / 1024 spacecraft, DERIVED to keep the rungs near
the same wall; `results/plan_ladder_aero_space-falcon-1_t8_n*.csv`) reaches the
generic flat branch rather than the single-harmonics one, and shows the same
defect at the same place:

| N | default (flat@8) | satellite_batch@1 | flat@8/static | best |
|---|---|---|---|---|
| 64 | 27.3 | **9.43** | 26.7 | batch, default 2.89x behind |
| 256 | 12.11 | **9.11** | 12.06 | batch, default 1.33x behind |
| 1024 | **9.23** | 9.83 | 9.59 | the default |

## The decision

A new routing-only floor, `SPACEAGORA_HARMONICS_FLAT_MIN_SATS_PER_WORKER`
(default **64**), in `src/simulation/engine/setup.jl`: the default opens a
multi-worker flat team only when `active_sats >= floor * thread_budget`, i.e.
when every worker of a full-budget team gets at least 64 spacecraft. Below
that it returns `satellite_batch` at the full budget. Above it nothing changes:
the flat route at the width it always had.

It gates the *route*, not the width, because the ladder shows narrowing the
flat route only makes it slower. It is a new knob rather than a new default for
`SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER`, because that one also sizes
the pre-pass slice and caps the kernel's worker count for every plan, including
calibrated ones; raising it would have changed what the adaptive routes run.

Where 64 comes from (SOURCED). Crossovers in satellites per worker at the full
budget:

| data | flat loses at | flat wins at |
|---|---|---|
| this ladder, vacuum, 8 threads | 16 (N=128) | 32 (N=256) |
| this ladder, + aero, 8 threads | 32 (N=256) | 128 (N=1024) |
| TRX50 P1, 32 threads | 32 (N=1024, 2.5x) | 128 (N=4096) |

Every value in (32, 64] separates every point in all three rows except the
8-thread vacuum N=256 point, which no value above 32 can keep; 64 is the
smallest satellites-per-worker at which the flat route was measured to be the
best plan. The cost of that choice is the one point: at 256 vacuum spacecraft on
8 threads the default is now 1.18x behind the best plan (7.44 s against 6.31 s).
That is the largest residual anywhere in the data.

Applied to both flat branches, with three exemptions that keep unmeasured paths
exactly as they were: a thread budget of 1 (no tasks are spawned there; the
flat route is taken for the batched coefficient sweep, WS10c's admission);
an enclosing outer split (already clamped to one worker; only single
simulations were measured); and, on the generic branch, any stack containing a
batched pre-pass kernel (n-body, SRP, inverse-square), where the flat route
shares ephemeris samples the batch cannot and nothing here measured it above
one thread. The opt-in spin-barrier dispatch bypasses the floor too (ASSUMED:
that dispatch exists to remove exactly this per-call charge, and the ladder did
not run it).

### What the default now does, and what it costs

8 threads, this machine, from the committed CSVs (seconds):

| N (vacuum) | old default | new default | change |
|---|---|---|---|
| 16 | 3.13 | same plan | none |
| 32 | 11.62 | 4.20 | 2.77x faster |
| 64 | 11.03 | 7.09 | 1.56x faster |
| 128 | 9.04 | 7.28 | 1.24x faster |
| 256 | 6.31 | 7.44 | 1.18x slower |
| 512, 1024, 4096 | | same plan | none |

Aero at 8 threads: 2.89x faster at 64, 1.33x faster at 256, unchanged at 1024.

For the archived runs, the new default picks the plan the sweep pinned, so the
expected pinned-route P1 timings are the adaptive routes' own measurements in
the same run (TRX50, 32 threads: new default flat only from 2048 spacecraft):

| N | pinned routes, old | expected pinned, new (R6 median, same run) |
|---|---|---|
| 16 | 2.13-2.17 | unchanged |
| 64 | 34.4-35.1 | ~3.3 |
| 256 | 21.1-21.2 | ~2.7 |
| 1024 | 4.94-5.00 | ~2.0 |
| 4096 | 2.15-2.78 | unchanged |

At 12 threads on the workstation archive (flat from 768 spacecraft) the 64
point moves from 14.1 s toward the adaptive 5.7 s; the 256 point now takes the
batch, where the sweep had retained flat within its 10% override margin, so
that point may be up to about 1.1x slower -- not measured here.

The serial baseline (`R0`) is untouched: it forces the budget to 1.


## Bit-identity

The plans differ only in execution route; each satellite's accumulation runs in
the same order on every route. Proven, not assumed, with full state histories
(every saved time, every component of every spacecraft, raw `Float64`) and
`cmp`, 8 threads for every dump:

- At base commit `adb343566`, all 11 plans on each reference case are byte
  identical to one another: L50 vacuum, 256 spacecraft, 30000 s (24,702,744
  bytes each); L50 + aero over the exponential atmosphere, 64 spacecraft,
  30000 s (17,552,808 bytes each).
- The new default's trajectory against the old default's: `cmp` reports no
  difference on both cases (MD5 `6da2ea8e6bd7e6223c4cbcbb9df69985` vacuum,
  `9069c41ff7ff7f4a84b280529b0e1b77` aero, same at base and tip), and every
  plan at the tip is byte identical to the base old-default file.
- `test/unit/simulation/harmonics_batch_parity_tests.jl` (22 pass) and
  `test/unit/simulation/third_body_route_parity_tests.jl` (84 pass) pass at the
  tip.

## Open

- The single inverse-square branch sizes a flat team as `min(budget,
  active_sats)` from 8 spacecraft with no work floor. It is the same shape of
  decision and was not measured here.
- The TRX50 expectations above are the sweep's pinned-batch medians from the
  same run, not a rerun of the new default there.
- The 12-thread workstation point at 256 spacecraft is the one place the new
  default may be slower by more than the 8-thread residual; a 12-thread ladder
  would settle it.


