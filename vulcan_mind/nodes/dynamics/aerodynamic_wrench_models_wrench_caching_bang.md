---
id: dynamics.aerodynamic_wrench_models_wrench_caching_bang
label: wrench_caching!
kind: function
source:
  file: src/dynamics/coupled/aerodynamic_wrench_models.jl
  symbol: wrench_caching!
  lines:
  - 542
  - 542
inputs:
- id: model
  type: AerodynamicCoefficientConstant
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: StateSample
  units: n/a
  required: true
  description: Positional argument `x`.
- id: env
  type: EnvironmentSample
  units: n/a
  required: true
  description: Positional argument `env`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
- id: p
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
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
  type: Tuple{SVector{3,
  units: n/a
  description: Return value of `wrench_caching!`; mutates `model` in place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynamics
origin: agent
---

# wrench_caching!

## Purpose
ODE-aware variant of `wrench` that optionally samples the atmosphere per link and stores drag/lift/cross components into the save cache for output logging.

## Design & Implementation
Signature `(model, x, env, t, p::ODEParams, sat_idx::Int)`. Builds `link_atmosphere_fn = pos -> _aero_link_atmosphere_query(p, sat_idx, t, pos, env.planet)` when `_per_link_enabled(model)`, else `nothing`. Calls `_aero_pure_wrench` with the model's coefficient mode (and `fixed_attitude_incidence` for fM), then `_store_aero_caches!(p, sat_idx, drag_ii, lift_ii, cross_ii)` and returns `(force, torque)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | AerodynamicCoefficientConstant | n/a | yes | Positional argument `model`. |
| in | `x` | StateSample | n/a | yes | Positional argument `x`. |
| in | `env` | EnvironmentSample | n/a | yes | Positional argument `env`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `p` | ODEParams | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{SVector{3, | n/a | — | Return value of `wrench_caching!`; mutates `model` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.effector_sampling_environment_requirements|environment_requirements]] · `callees` → `callers` · call · `src/core/types/effector_sampling.jl:165-165`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl`
- [[simulation.dynamics_rhs__evaluate_dynamic_effector|_evaluate_dynamic_effector]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:17-17`

**Downstream**

- `callees` → [[dynamics.aerodynamic_wrench_models__aero_link_atmosphere_query|_aero_link_atmosphere_query]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:551-551`
- `callees` → [[dynamics.aerodynamic_wrench_models__store_aero_caches_bang|_store_aero_caches!]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:553-553`
- `callees` → [[dynx.coupled_aerodynamic_wrench_models_aero_pure_wrench|_aero_pure_wrench]] · `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:552-552`
<!-- vulcan:connections:end -->

## Limitations
The closure allocates per call when per-link sampling is enabled. Cache writes happen on every RHS evaluation, including rejected adaptive steps, so logged components may lag the accepted state. Per-link sampling applies only when `orientation_sim` is on and only for non-root links.

## Provenance
Mapped from `src/dynamics/coupled/aerodynamic_wrench_models.jl` line 542.
