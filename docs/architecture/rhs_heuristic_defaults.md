# Default RHS execution plans for constellations

`_rhs_execution_plan_uncached` in
`src/orchestration/simulation/engine/setup.jl` supplies the heuristic plan for a
constellation's right-hand side. With `SPACEAGORA_RHS_CALIBRATE=auto` (the
default), the pre-solve calibration sweep can instead pin a plan through the
production calibration store. Static profiles, and the heuristic arm of the
study below, run with calibration off so that the two choices can be examined
separately.

This page covers the rule, its scope and how to check it. The plan ladders and
the derivation of the shipped defaults are kept in the lab's research records
(see [Record a result](../src/user/studies_benchmarks.md#record-a-result)).

## Route admission and worker width

The flat route has dispatch and reduction costs on top of the work done in its
pre-passes. Two settings answer different questions:

- `SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER` (default 4) sizes a
  harmonics pre-pass slice for SIMD and caps the pre-pass's worker count, for
  any plan.
- `SPACEAGORA_HARMONICS_FLAT_MIN_SATS_PER_WORKER` (default 64) is the heuristic
  floor for opening a flat worker team at all: the default plan takes the flat
  route at the full budget only when every worker gets at least this many
  spacecraft.

The second floor gates route admission; it does not narrow the flat team's
width. It does not apply at a budget of one, under an enclosing outer split, to
stacks with a batchable pre-pass kernel (n-body, SRP, inverse-square) on the
generic branch, or when the opt-in spin-barrier dispatch is active
(`SPACEAGORA_HARMONICS_BATCH_SPIN_BARRIER`). Those cases keep their own support
checks. The floor is a heuristic, not a guarantee that the chosen route is the
fastest on every machine or workload.

The generic flat branch also requires `SPACEAGORA_EFFECTOR_FLAT_MIN_THREAD_BUDGET`
threads (default `SPACEAGORA_AUTO_THREAD_MIN_BUDGET`, itself 4 by default), to
amortize the overhead of its work queue. A stack whose effectors are all served
by pre-passes (`_rhs_all_prepass_effectors`) never builds that queue and is
exempt from this floor, though it still has to pass the generic branch's work
and spacecraft-per-worker checks. A stack with any queue-only effector keeps
the floor.

At a budget of one, a supported pre-pass-only stack with at least
`SPACEAGORA_EFFECTOR_FLAT_MIN_SATS` spacecraft (default 24) takes the flat route
and runs the batched kernels inline, without spawning tasks. See
[third-body and SRP paths](third_body_cost.md) and
[aerodynamic batching](../../benchmarks/studies/aero_batch/FINDINGS.md) for the
related ownership constraints.

## Comparing plans

The [plan-ladder study](../../benchmarks/studies/rhs_heuristic_defaults/README.md)
pins each candidate plan through a one-row production calibration store and
records the plan the engine reports having applied (source, mode, allotment and
scheduler), so a pin that failed is visible.

```bash
julia --project=. --threads=8 benchmarks/studies/rhs_heuristic_defaults/plan_ladder.jl --sizes=16,32,64,128,256,512,1024,4096 --repeats=3
```

`--case=vacuum` (the default) flies degree-50 harmonics in vacuum, and
`--case=aero` adds an aerodynamic coefficient model over an analytic
exponential atmosphere. `--dump-dir=<dir>` also saves full state histories.
Compare plans within one invocation, and record the host load, revision and
thread budget. Measure again on the target host before changing a default.

A route change must keep each spacecraft's accumulation order. Compare every
saved time and state component, not only the final state or a rounded table.
`test/unit/simulation/rhs_heuristic_defaults_tests.jl`,
`test/unit/simulation/harmonics_batch_parity_tests.jl` and
`test/unit/simulation/third_body_route_parity_tests.jl` hold the regression
coverage. A recorded state-identity result qualifies only its own revision and
cases.
