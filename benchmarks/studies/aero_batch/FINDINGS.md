# Aerodynamics and density at constellation scale

The flat constellation route can evaluate supported aerodynamic effectors in a
pre-pass before reducing the force and torque slots. This page describes the
ownership and numerical constraints of that path. Historical profiles and
before/after measurements are retained under the
[research-record policy](../../../docs/src/user/studies_benchmarks.md#research-records-and-generated-results).

## Aerodynamic pre-pass

`_accumulate_aero_flat_batch!` uses the existing `wrench_caching!` call and
`_aero_pure_wrench` arithmetic for each spacecraft. It writes that spacecraft's
force and torque slot and its own drag/lift/cross cache entries. The flat queue
excludes effectors already handled by the pre-pass, so they are evaluated once.
The later slot reduction still visits effectors in their established order.

The pre-pass supports `AerodynamicCoefficientfM` without per-link atmosphere.
Per-link atmosphere, the constant/no-ballistic-flight/mesh-surrogate models and
IMEX-partitioned calls retain their existing paths. These exclusions preserve
model-specific sampling, warning and partition behavior.

Keep the arithmetic in the same compiled `_aero_pure_wrench` method. Moving it
into another caller can change compiler treatment of library `muladd` operations
even when the written expressions look identical. The pre-pass adds no
reassociation, FMA or SIMD annotation, and preserves link accumulation order.
Worker slices may change which spacecraft is evaluated first; they must not
change any accumulation within one spacecraft or share a writable cache slot.

## Shared analytic atmosphere sampling

`_uniform_light_density_model` and `_fill_uniform_light_atmosphere!` in
`src/orchestration/simulation/engine/effector_sampling.jl` allow one
`getDensityBatch!` call over the constellation's altitude/latitude/longitude
buffers. The per-spacecraft path remains the fallback. Admission requires:

- all spacecraft active and resolving to the same supported density object;
- GRAM runtime statistics off and no GRAM grid snapshot;
- a light density model, with freeze-per-step, vacuum look-ahead cache and
  GRAM track cache inactive;
- no selected density-service route, and sufficiently sized output/time buffers.

The specialized exponential and piecewise-exponential batch methods must retain
the scalar methods' arithmetic. This path does not change the atmosphere model
or native-GRAM locking.

## Route scope

Both pre-passes run on the flat constellation route. At a budget of one, the
auto route can admit a supported pre-pass-only stack through
`_rhs_all_prepass_effectors` once `flat_min_sats` is met. It spawns no worker
team. An explicit `SPACEAGORA_RHS_EXECUTION_MODE=flat` request at that budget
still follows its own fallback rules. A per-link-atmosphere fM model does not
qualify as a pre-pass effector.

## Validation and measurement

`test/unit/dynamics/aero_batch_parity_tests.jl` compares the pre-pass with the
per-spacecraft call, varies worker allotments and inactive spacecraft, checks
vacuum/zero-airspeed returns, compares scalar and batch density values, and
compares short solves across routes. Its exact comparisons are regression
checks for the supported cases, not a proof for excluded models.

Use `benchmarks/studies/aero_batch/attribution.jl` to profile the wrench,
atmosphere sampling, callbacks, queue and harmonics work separately. Use
`benchmarks/studies/aero_batch/variants.jl` with `--dump` for complete state
histories and `--repeats` for timings. Compare histories before interpreting
performance differences, and record the revision, route, thread count and case.
Serial and threaded routes need separate attribution. Preserve raw outputs
outside the public source tree; no archived ratio qualifies a new revision.
