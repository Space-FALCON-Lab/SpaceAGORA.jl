# Third-body and solar-radiation-pressure execution paths

What third-body gravity and solar radiation pressure (SRP) cost depends on the
solver, the RHS route and who owns the environment samples, as well as on the
force calculation itself. This page describes the current implementation choices
and the studies that exercise them. Profiles and timing tables are kept in the
lab's research records (see
[Record a result](../src/user/studies_benchmarks.md#record-a-result)).

## Shared samples and name resolution

`_spice_query_name` in
`src/models/dynamics/effectors/perturbations/gravity/models_and_workspaces.jl`
interns canonical SPICE body names. The table is published copy-on-write behind
a lock, so readers use a dictionary that is never mutated afterwards. The cached
result must equal the uncached string transformation
(`_spice_query_name_uncached`).

`sample_third_body_ephemerides` in
`src/orchestration/simulation/engine/effector_sampling.jl` maps over the typed
tuple of body names the effector declares
(`environment_requirements(model).third_body_names`) rather than building a
tuple from a runtime length. Body order and values are unchanged, and the
compiler can infer the tuple length. The flat
route's shared prefill and its n-body and SRP pre-passes reuse the
time-dependent body samples across spacecraft.

## Flat-route admission

`_rhs_all_prepass_effectors` in `src/orchestration/simulation/engine/setup.jl`
identifies stacks whose effectors are all served by the flat driver's
pre-passes. Such a stack with at least `SPACEAGORA_EFFECTOR_FLAT_MIN_SATS`
spacecraft (default 24) takes the flat route at an allotment of one. That route
performs its batches and reduction inline and builds no work queue when no
queue-only effector remains.

This admission is separate from the worker-team size floor: a serial pre-pass
can reuse coefficients and shared samples without spawning tasks. Stacks with
other effectors keep the route's support and budget checks; see
[RHS heuristic defaults](rhs_heuristic_defaults.md).

## Solver scope

`_auto_stiff_smooth_gravity_effector` in
`src/orchestration/simulation/engine/solver_policy.jl` admits
`SolarRadiationPressureModel` to the smooth-gravity Tsit5 fast path alongside
the gravity models, and `_auto_stiff_smooth_gravity_reject_reason` allows the
solar sample SRP needs. The other eligibility checks still apply: no attitude
simulation, no active control, guidance or navigation effectors, and no
effector that needs atmosphere samples. `SolverConfig.auto_stiff_gravity_tsit5`
(default `true`) set to `false`, or `SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5=0`,
turns the fast path off and keeps automatic stiff-solver selection.

An explicit solver and an auto-switching solver need not take the same steps or
reach the same states. Compare trajectory error as well as runtime, especially
when the step cap changes across eclipse transitions. Agreement on one case does
not establish solver equivalence for every mission.

## Validation and measurement

`test/unit/dynamics/third_body_cost_tests.jl` checks name-resolution identity,
allocation behavior and planet-frame reuse. The route comparisons are in
`test/unit/simulation/third_body_route_parity_tests.jl`.

Run each comparison with an explicit thread count, and keep the case, revision,
solver, route, timings and complete state dumps together:

```bash
julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl --n=256 --mission=5800 --repeats=3
julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl --n=256 --mission=5800 --variants=vacuum,srp_nbody --solver=tsit5
julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl --n=256 --mission=5800 --variants=vacuum,nbody --rhs=satellite
julia --project=. --threads=8 benchmarks/studies/third_body_cost/variants.jl --n=256 --mission=5800 --mode=inner_only --repeats=2
julia --project=. --threads=1 benchmarks/studies/third_body_cost/attribution.jl --n=256 --mission=5800 --profile
SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5=0 julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl --n=256 --mission=5800 --variants=vacuum,srp_nbody
```

`variants.jl --dump=<prefix>` writes each variant's complete `Float64` state
history to `<prefix>_<variant>.bin`. Compare full histories when checking an
execution change that should not change results. Keep solver-change error
comparisons separate from exact route-parity checks, and keep generated
measurements out of the source tree.
