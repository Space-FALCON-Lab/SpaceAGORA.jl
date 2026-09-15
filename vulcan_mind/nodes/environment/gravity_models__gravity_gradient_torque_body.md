---
id: environment.gravity_models__gravity_gradient_torque_body
label: _gravity_gradient_torque_body
kind: function
source:
  file: src/environment/gravity/gravity_models.jl
  symbol: _gravity_gradient_torque_body
  lines:
  - 63
  - 63
inputs:
- id: model
  type: Any
  units: n/a
  required: true
  description: Positional argument `model`.
- id: pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_ii`.
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: param
  type: ODEParams
  units: n/a
  required: true
  description: Positional argument `param`.
- id: i
  type: Int64
  units: n/a
  required: true
  description: Positional argument `i`.
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
  type: SVector{3,
  units: n/a
  description: Return value of `_gravity_gradient_torque_body`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- environment
origin: agent
---

# _gravity_gradient_torque_body

## Purpose
Legacy-path (line 63) computation of the body-frame gravity-gradient torque for satellite `i` from the `ComponentVector` state, returning zero unless the model has `gravity_gradient = true`, orientation simulation is enabled, and the state carries a quaternion. A second method on `(model, x::StateSample, planet)` serves the sampled `wrench` path.

## Design & Implementation
Marked `@inline`. Early-outs return `SVector(0,0,0)` when `!model.gravity_gradient`, `!param.args.mission_configuration.orientation_sim`, `i` is outside `1:length(param.args.dynamics_model.spacecraft)`, `x` has no `:q` property, or `r = norm(pos_ii)` is non-finite or non-positive. Otherwise it reads `inertia_tensor` from `spacecraft[i]`, builds `q_body::SVector{4,Float64}` from `x.q`, rotates the inertial position into the body frame with `rot(q_body) * pos_ii`, and returns `gravity_gradient(inertia_tensor, r_body, planet.μ)` in N m. The `StateSample` overload does the same with `x.q_ib`, `x.spacecraft.inertia_tensor` and `env.planet`, returning zero when either is `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | Any | n/a | yes | Positional argument `model`. |
| in | `pos_ii` | SVector{3, Float64} | n/a | yes | Positional argument `pos_ii`. |
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `param` | ODEParams | n/a | yes | Positional argument `param`. |
| in | `i` | Int64 | n/a | yes | Positional argument `i`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector{3, | n/a | — | Return value of `_gravity_gradient_torque_body`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[environment.gravity_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:189-189`
- [[environment.gravity_models_environment_requirements|environment_requirements]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:270-270`
- [[environment.gravity_models_gravity_backbone_acceleration_ii|gravity_backbone_acceleration_ii]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:221-221`
- [[environment.gravity_models_wrench|wrench]] · `callees` → `callers` · call · `src/environment/gravity/gravity_models.jl:201-201`
- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/gravity/gravity_models.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/environment/gravity/gravity_models.jl:87-87`
- `callees` → [[core.quaternion_utils_rot|rot]] · `callers` · call · `src/environment/gravity/gravity_models.jl:88-88`
- `callees` → [[environment.gravity_models_gravity_gradient|gravity_gradient]] · `callers` · call · `src/environment/gravity/gravity_models.jl:89-89`
- `callees` → [[environment.gravity_models_j2_secular_rates|j2_secular_rates]] · `callers` · call · `src/environment/gravity/gravity_models.jl:112-112`
<!-- vulcan:connections:end -->

## Limitations
The quaternion convention must match `rot`: the code treats `x.q` as inertial-to-body, and a body-to-inertial quaternion would yield a torque of the wrong sign with no error. `hasproperty(x, :q)` on a `ComponentVector` is a dynamic check performed every derivative evaluation. Silently returning zero for an out-of-range `i` hides indexing bugs. Only a central field is used for the gradient even when the force model includes J2.

## Provenance
Mapped from `src/environment/gravity/gravity_models.jl` line 63.
