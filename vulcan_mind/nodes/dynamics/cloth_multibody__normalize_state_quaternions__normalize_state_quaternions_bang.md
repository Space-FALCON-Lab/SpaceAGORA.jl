---
id: dynamics.cloth_multibody__normalize_state_quaternions__normalize_state_quaternions_bang
label: _normalize_state_quaternions!
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: _normalize_state_quaternions!
  lines:
  - 469
  - 469
inputs:
- id: x
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `x`.
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
  description: Return value of `_normalize_state_quaternions!`; mutates `x` in place.
    Returns `x`.
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

# _normalize_state_quaternions!

## Purpose
Renormalises every body quaternion in a state vector in place, correcting the drift that integrating the quaternion rate accumulates.

## Design & Implementation
Divides the length by thirteen for the body count and overwrites slots four to seven of each block with `_unit_quat` of themselves. Returns the same vector. Called after every RK4 step and every implicit-midpoint trial.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | AbstractVector | n/a | yes | Positional argument `x`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_normalize_state_quaternions!`; mutates `x` in place. Returns `x`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.cloth_multibody_residual|residual]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:667-667`
- [[dynamics.cloth_multibody_simulate_compliant_multibody|simulate_compliant_multibody]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:718-718`
- [[dynamics.cloth_multibody_step_compliant_multibody_rk4|step_compliant_multibody_rk4]] · `callees` → `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:642-642`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:473-473`
- `callees` → [[dynamics.cloth_robot_arm_dynamics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:473-473`
- `callees` → [[vehicle.robotics__unit_quat|_unit_quat]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:473-473`
<!-- vulcan:connections:end -->

## Limitations
Projection rather than a constraint-preserving integrator; the projected state is no longer exactly the integrator's output, which the implicit stepper accounts for by re-evaluating the residual after projection.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 469.
