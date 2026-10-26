# RHS heuristic defaults

This study compares the default right-hand-side (RHS) execution plan with every
plan the calibration sweep would consider, size by size, on one machine at one
thread budget.

## Why

The pinned parallel routes (`R1a`, `R2` and `R3` in the paper study, which lives
in the SpaceAGORAParallelizationStudies repository) run whatever
`_rhs_execution_plan_uncached` in `src/orchestration/simulation/engine/setup.jl`
decides. The adaptive routes (`R6`, `R7`) run whatever the calibration sweep
before the solve measured. Compare the heuristic and calibrated plans at the
same resource budget before you change either one.
`docs/architecture/rhs_heuristic_defaults.md` describes the shipped rule and
where it applies; this folder holds the measurement procedure.

## Running it

```bash
julia --project=. --threads=8 benchmarks/studies/rhs_heuristic_defaults/plan_ladder.jl \
    --sizes=16,32,64,128,256,512,1024,4096 --repeats=3
```

Other flags:

- `--case=vacuum|aero`: L50 harmonics in vacuum (the paper study's P1 case),
  or L50 harmonics plus an aerodynamic coefficient model over an analytic
  exponential atmosphere.
- `--dump-dir=<dir>`: also write each timed solve's full state history (every
  saved time, every component of every spacecraft, raw `Float64`), so two plans
  or two commits can be compared with `cmp`.
- `--mission=<s>`: overrides every rung's mission length; for a quick check
  that the script works, not for measurements.

Results go to `results/` and belong in the private research record. Compare
plans within one invocation and note the host load; running them back to back
reduces some variability, but it does not remove contention or make a ratio
hold on another machine.

## How a plan is pinned

Through the production calibration cache, not a private hook. Each plan gets a
one-row TOML store, handed to the engine with `SPACEAGORA_RHS_CALIBRATION_PATH`,
and `SPACEAGORA_RHS_CALIBRATE_MIN_SOLVE_S` is raised above any solve length
here so the cached choice is used instead of being swept again. The `heuristic`
row runs with `SPACEAGORA_RHS_CALIBRATE=off` instead, exactly as the pinned
routes in the paper harness do. Every row records the plan the engine reports
it applied (`applied_source`, `applied_mode`, `applied_allotment`,
`applied_scheduler`), so a pin that silently failed shows up in the CSV instead
of being averaged into a result.
