# RHS heuristic defaults (WS11a)

What the default RHS execution plan costs against every plan the calibration
sweep would consider, size by size, on one machine at one thread budget.

## Why

The pinned parallel routes (`R1a`/`R2`/`R3` in
`benchmarks/studies/paper_parallelization_benchmarks`) run whatever
`_rhs_execution_plan_uncached` (`src/simulation/engine/setup.jl`) decides; the
adaptive routes (`R6`/`R7`) run whatever the pre-solve calibration sweep
measured. On the archived TRX50 P1 ladder the two answers differ by up to ten
times in the middle of the constellation-size range, so what needs fixing is
the default, not the sweep. `docs/architecture/rhs_heuristic_defaults.md` is
the write-up; this directory holds the measurement that backs it.

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

Results land in `results/` and are committed. Only ratios within one invocation
mean anything: the machine is shared, and the whole point of running every plan
back to back in one process is that the comparison does not depend on the
machine being quiet.

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
