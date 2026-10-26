# RHS heuristic defaults

What the default RHS execution plan costs against every plan the calibration
sweep would consider, size by size, on one machine at one thread budget.

## Why

The pinned parallel routes (`R1a`/`R2`/`R3` in
the paper study (SpaceAGORAParallelizationStudies)) run whatever
`_rhs_execution_plan_uncached` (`src/orchestration/simulation/engine/setup.jl`) decides; the
adaptive routes (`R6`/`R7`) run whatever the pre-solve calibration sweep
measured. Compare the heuristic and calibrated plans at the same resource
budget before changing either. `docs/architecture/rhs_heuristic_defaults.md`
describes the shipped rule and its scope; this directory holds the measurement
procedure.

## Running it

```bash
julia --project=. --threads=8 benchmarks/studies/rhs_heuristic_defaults/plan_ladder.jl \
    --sizes=16,32,64,128,256,512,1024,4096 --repeats=3
```

Useful flags:

- `--case=vacuum|aero` — L50 harmonics in vacuum (the P1 shape), or L50
  harmonics plus an aerodynamic coefficient model over an analytic exponential
  atmosphere.
- `--dump-dir=<dir>` — also write each timed solve's full state history (every
  saved time, every component of every spacecraft, raw `Float64`) so two plans,
  or two commits, can be compared with `cmp`.
- `--mission=<s>` — functional smoke only; overrides every rung's mission.

Results land in `results/` and belong in the private research record. Compare
plans within one invocation and record host load; running back to back reduces
some variability but does not eliminate contention or make a ratio portable.

## How a plan is pinned

Through the production calibration cache, not a private hook: each plan gets a
one-row TOML store, handed to the engine with
`SPACEAGORA_RHS_CALIBRATION_PATH`, with
`SPACEAGORA_RHS_CALIBRATE_MIN_SOLVE_S` raised past any solve length here so the
cached verdict is honored rather than re-swept. The `heuristic` row runs
`SPACEAGORA_RHS_CALIBRATE=off` instead, which is exactly what the pinned routes
in the paper harness run. Every row records the plan the engine reported having
applied (`applied_source`/`applied_mode`/`applied_allotment`/`applied_scheduler`),
so a pin that silently failed is visible in the CSV rather than being averaged
into a verdict.
