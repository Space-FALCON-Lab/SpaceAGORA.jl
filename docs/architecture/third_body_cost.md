# Third-body and solar-radiation-pressure execution paths

Third-body and SRP cost depends on the solver, RHS route and environment-sample
ownership as well as the force calculation itself. This page describes the
current implementation choices and the studies that exercise them. Historical
profiles and timing tables are retained under the
[research-record policy](../src/user/studies_benchmarks.md#research-records-and-generated-results).

## Shared samples and name resolution

`_spice_query_name` in `src/models/dynamics/effectors/perturbations.jl` interns
canonical body names. Its table is published copy-on-write behind a lock;
readers use a dictionary that is not subsequently mutated. The cached result
must equal the uncached string transformation.

`sample_third_body_ephemerides` in
`src/orchestration/simulation/engine/effector_sampling.jl` maps over the typed
`model.body_names` tuple. This retains body order and values while allowing the
compiler to infer the tuple length. The flat route's shared prefill and n-body
and SRP pre-passes reuse time-dependent body samples across spacecraft.

## Flat-route admission

`_rhs_all_prepass_effectors` in
`src/orchestration/simulation/engine/setup.jl` identifies stacks served by the
flat driver's pre-passes. A supported stack above `flat_min_sats` can take the
flat route at allotment one. That route performs its batches and reduction
inline and does not build the work queue when no queue-only effectors remain.

This admission is different from the worker-team size threshold: a serial
pre-pass can reuse coefficients and shared samples without spawning tasks.
Stacks with other effectors retain the route's support and budget checks.
See [RHS heuristic defaults](rhs_heuristic_defaults.md) for those constraints.

## Solver scope

`_auto_stiff_smooth_gravity_effector` and
`_auto_stiff_smooth_gravity_reject_reason` in
`src/orchestration/simulation/engine/solver_policy.jl` admit
`SolarRadiationPressureModel` to the smooth-gravity Tsit5 fast path and allow
the solar sample it requires. Other eligibility checks still apply.
`SolverConfig.auto_stiff_gravity_tsit5=false`, or
`SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5=0`, disables that fast path.

An explicit solver and an autoswitching solver need not have identical step
sequences or states. Compare trajectory error as well as runtime, especially
when changing the step cap across eclipse transitions. Historical agreement on
one case does not establish solver equivalence for every mission.

## Validation and measurement

`test/unit/dynamics/third_body_cost_tests.jl` checks name-resolution identity,
allocation behavior and planet-frame reuse. The route comparisons live in
`test/unit/simulation/third_body_route_parity_tests.jl`.

Run each comparison with an explicit thread count and preserve the case,
revision, solver, route, timings and complete state dumps together:

```bash
julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl --n=256 --mission=5800 --repeats=3
julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl --n=256 --mission=5800 --variants=vacuum,srp_nbody --solver=tsit5
julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl --n=256 --mission=5800 --variants=vacuum,nbody --rhs=satellite
julia --project=. --threads=8 benchmarks/studies/third_body_cost/variants.jl --n=256 --mission=5800 --mode=inner_only --repeats=2
julia --project=. --threads=1 benchmarks/studies/third_body_cost/attribution.jl --n=256 --mission=5800 --profile
SPACEAGORA_AUTO_STIFF_GRAVITY_TSIT5=0 julia --project=. --threads=1 benchmarks/studies/third_body_cost/variants.jl --n=256 --mission=5800 --variants=vacuum,srp_nbody
```

`--dump` writes complete `Float64` state histories. Compare full histories when
checking a behavior-preserving execution change. Keep solver-change error
comparisons distinct from exact route-parity checks, and keep generated
measurements outside the public source tree.
