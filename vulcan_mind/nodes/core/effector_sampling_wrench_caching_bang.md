---
id: core.effector_sampling_wrench_caching_bang
label: wrench_caching!
kind: function
source:
  file: src/core/types/effector_sampling.jl
  symbol: wrench_caching!
  lines:
  - 172
  - 172
inputs:
- id: model
  type: Any
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: env
  type: Any
  units: n/a
  required: true
  description: Positional argument `env`.
- id: t
  type: Any
  units: n/a
  required: true
  description: Positional argument `t`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Any
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: result
  type: Any
  units: n/a
  description: Return value of `wrench_caching!`; mutates `model` in place. Returns
    `wrench(model, x, env, t)`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# wrench_caching!

## Purpose
Variant of the `wrench` hook that also lets an effector record per-component diagnostics (for example separate drag, lift and cross-wind force vectors for aerodynamic models) into the integrator's `p.save_cache` for a given satellite, while still returning the total force and torque.

## Design & Implementation
Declared with `function wrench_caching! end` and given one generic `@inline` fallback `wrench_caching!(model, x, env, t, p, sat_idx) = wrench(model, x, env, t)` that ignores `p` and `sat_idx`. The expected signature is `(model, x::StateSample, env::EnvironmentSample, t::Float64, p, sat_idx::Int)` returning `(force_ii, torque_body)` in newtons and newton-metres. The engine calls this rather than `wrench` when a save cache is present, so models that override it write to `p.save_cache` and models that do not pay no cost beyond the extra call layer.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | Any | n/a | yes | Positional argument `model`. |
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `env` | Any | n/a | yes | Positional argument `env`. |
| in | `t` | Any | n/a | yes | Positional argument `t`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Any | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `wrench_caching!`; mutates `model` in place. Returns `wrench(model, x, env, t)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.effector_sampling_environment_requirements|environment_requirements]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:165-165`
- [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:17-17`

**Downstream**

- `callees` → [[core.effector_sampling_solver_partition|solver_partition]] · `callers` · call · `src/core/types/effector_sampling.jl:175-175`
- `callees` → [[dynamics.aerodynamic_wrench_models_solver_partition|solver_partition]] · `callers` · call · `src/core/types/effector_sampling.jl:175-175`
- `callees` → [[dynamics.aerodynamic_wrench_models_wrench|wrench]] · `callers` · call · `src/core/types/effector_sampling.jl:172-172`
- `callees` → [[dynamics.perturbations_wrench|wrench]] · `callers` · call · `src/core/types/effector_sampling.jl:172-172`
- `callees` → [[environment.gravity_models_wrench|wrench]] · `callers` · call · `src/core/types/effector_sampling.jl:172-172`
<!-- vulcan:connections:end -->

## Limitations
Writing to `p.save_cache` from inside an ODE stage means diagnostics reflect the last evaluated stage, not the accepted step, unless the engine copies them at save time. The fallback is untyped on every argument, so it participates in dispatch for any six-argument call and can mask a mistyped override. Because it is a `!` function, an override is expected to mutate `p.save_cache`, but nothing guards against concurrent writes when satellites are evaluated in parallel threads.

## Provenance
Mapped from `src/core/types/effector_sampling.jl` line 172.
