# Aerodynamics and density at constellation scale

On the flat constellation route, supported aerodynamic effectors can be
evaluated for every spacecraft in a pre-pass before the force and torque slots
are summed. This page describes who owns that path and the numerical rules it
must keep. Earlier profiles and before-and-after measurements are kept outside
this repository, as [Record a result](../../../docs/src/user/studies_benchmarks.md#record-a-result)
describes.

## Aerodynamic pre-pass

`_accumulate_aero_flat_batch!` calls `_aero_prepass_satellite!` for each
spacecraft (`src/orchestration/simulation/engine/dynamics_rhs.jl`). That
function does, step by step, what `wrench_caching!` would do: it builds the
environment sample the effector asks for (planet frame and atmosphere only),
calls `_aero_pure_wrench(:fm, ...)` with the per-link hook off, and stores the
drag, lift and cross vectors in the spacecraft's own cache entries with the same
`_store_vector_cache!` helper. It then writes that spacecraft's force and torque
slot. The flat queue skips effectors the pre-pass has already handled, so each
is evaluated once, and the later slot sum still visits effectors in their usual
order.

The pre-pass supports `AerodynamicCoefficientfM` without per-link atmosphere.
Per-link atmosphere, the constant, no-ballistic-flight and mesh-surrogate
models, and IMEX-partitioned calls keep their existing paths, so their own
sampling, warnings and partition behavior do not change.

`_aero_pure_wrench` is called, not inlined, so the pre-pass runs the same
compiled method as the per-satellite route. Inlining it elsewhere can change
how the compiler treats the `muladd` operations inside StaticArrays even when
the written expressions are identical, and the last bit of the result moves.
The pre-pass adds no reassociation, FMA or SIMD annotation and keeps the order
in which link contributions are summed. Worker slices may change which
spacecraft is evaluated first; they must not change any sum within one
spacecraft or share a writable cache slot.

## Shared analytic atmosphere sampling

`_uniform_light_density_model` and `_fill_uniform_light_atmosphere!` in
`src/orchestration/simulation/engine/effector_sampling.jl` let one
`getDensityBatch!` call fill the constellation's altitude, latitude and
longitude buffers. The per-spacecraft path stays as the fallback. The batch path
is used only when:

- every spacecraft is active and resolves to the same supported density object;
- GRAM run-time statistics are off and there is no GRAM grid snapshot;
- the density model is a light one, with the freeze-per-step, vacuum look-ahead
  and GRAM track caches inactive;
- no density-service route is selected, and the output and time buffers are
  large enough.

The batch methods of the exponential and piecewise-exponential models must keep
the arithmetic of their scalar methods. This path changes neither the
atmosphere model nor the locking around native GRAM.

## Route scope

Both pre-passes run on the flat constellation route. With a thread budget of
one, the automatic route can still choose a stack made only of pre-pass
effectors, through `_rhs_all_prepass_effectors`, once there are at least
`flat_min_sats` spacecraft; it starts no worker team. An explicit
`SPACEAGORA_RHS_EXECUTION_MODE=flat` request at that budget follows its own
fallback rules. An fM model with per-link atmosphere does not count as a
pre-pass effector.

## Validation and measurement

`test/unit/dynamics/aero_batch_parity_tests.jl` compares single derivative
evaluations, bit for bit: the flat route against the per-satellite route,
different worker counts, a spacecraft outside the atmosphere (which must
contribute exactly zero), `getDensityBatch!` against `getDensity`, and the flat
route's atmosphere buffers against the per-satellite samples. These exact
comparisons cover the supported cases; they prove nothing for the excluded
models.

`benchmarks/studies/aero_batch/attribution.jl` profiles the wrench, the
atmosphere sampling, callbacks, the queue and the harmonics separately.
`benchmarks/studies/aero_batch/variants.jl` writes complete state histories with
`--dump` and times runs with `--repeats`. Compare the histories before you read
anything into a timing difference, and record the revision, route, thread count
and case. Serial and threaded routes need separate attribution. Keep raw
outputs outside the source tree; no archived ratio qualifies a new revision.
