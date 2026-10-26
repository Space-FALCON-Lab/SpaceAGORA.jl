# Default RHS execution plans for constellations

`_rhs_execution_plan_uncached` in
`src/orchestration/simulation/engine/setup.jl` supplies the heuristic plan.
With `SPACEAGORA_RHS_CALIBRATE=auto`, the pre-solve calibration can instead
pin a plan through the production calibration cache. Static profiles and the
heuristic comparison in this study disable calibration so that the two choices
can be examined separately.

Historical plan ladders and the derivation of shipped defaults are retained
under the [research-record policy](../src/user/studies_benchmarks.md#research-records-and-generated-results).
This page keeps the rule, its scope and the procedure for checking it.

## Route admission and worker width

The flat route has dispatch and reduction costs as well as the work performed
inside its pre-passes. `SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER` sizes a
pre-pass slice; it is not the same decision as whether to open a flat worker
team. `SPACEAGORA_HARMONICS_FLAT_MIN_SATS_PER_WORKER` supplies the latter heuristic
floor, with a default of 64 spacecraft per worker.

The floor gates route admission rather than narrowing the flat team's width.
The supported exceptions include a budget of one, an enclosing outer split,
the generic branch's batchable pre-pass kernels and opt-in spin-barrier
dispatch. These exceptions retain their own support checks. They do not make
the heuristic a guarantee of the fastest route on every machine or workload.

The generic flat branch also has
`SPACEAGORA_EFFECTOR_FLAT_MIN_THREAD_BUDGET` (default 4, from
`SPACEAGORA_AUTO_THREAD_MIN_BUDGET`) to amortize queue overhead. A stack whose
effectors are all served by pre-passes does not build that queue and is exempt
from this budget floor. It still passes the generic branch's work and
satellites-per-worker admission checks. A stack containing queue-only
effectors retains the floor.

At allotment one, supported pre-pass-only stacks can reuse batched kernels
without spawning tasks. See [third-body and SRP paths](third_body_cost.md) and
[aerodynamic batching](../../benchmarks/studies/aero_batch/FINDINGS.md) for the associated ownership constraints.

## Comparing plans

The [plan-ladder study](../../benchmarks/studies/rhs_heuristic_defaults/README.md)
pins candidate plans through one-row production calibration stores. It records
the applied source, mode, allotment and scheduler so a failed pin is visible.

```bash
julia --project=. --threads=8 benchmarks/studies/rhs_heuristic_defaults/plan_ladder.jl --sizes=16,32,64,128,256,512,1024,4096 --repeats=3
```

Use `--case=vacuum|aero` for the L50 vacuum or analytic-exponential-aero
configuration and `--dump-dir=<dir>` to save full state histories. Compare
repeated plans within the same invocation and record host load, revision and
thread budget. Re-measure on the target host before changing a default.

Route changes must retain each spacecraft's accumulation order. Compare every
saved time and state component, not just the final state or a rounded table.
`test/unit/simulation/rhs_heuristic_defaults_tests.jl`,
`harmonics_batch_parity_tests.jl` and `third_body_route_parity_tests.jl` retain
regression coverage; archived state identity only qualifies its recorded
revision and cases.
