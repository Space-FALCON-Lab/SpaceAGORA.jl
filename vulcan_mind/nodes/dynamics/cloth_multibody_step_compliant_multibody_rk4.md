---
id: dynamics.cloth_multibody_step_compliant_multibody_rk4
label: step_compliant_multibody_rk4
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: step_compliant_multibody_rk4
  lines:
  - 627
  - 627
inputs:
- id: model
  type: CompliantMultibodyModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `x`.
- id: t
  type: Real
  units: n/a
  required: true
  description: Positional argument `t`.
- id: dt
  type: Real
  units: n/a
  required: true
  description: Positional argument `dt`.
- id: dynamics_kwargs
  type: Vararg{Any}
  units: n/a
  required: false
  description: Keyword argument `dynamics_kwargs` (variadic).
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
  description: Return value of `step_compliant_multibody_rk4`. Returns `_normalize_state_quaternions!(xn)`.
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

# step_compliant_multibody_rk4

## Purpose
Advances the compliant state one step with classical fourth-order Runge-Kutta, for smooth or lightly stiff configurations.

## Design & Implementation
Copies the state, evaluates the derivative at the four standard stages with any `dynamics_kwargs` forwarded, combines with weights 1, 2, 2, 1 over six, and renormalises the quaternions. Returns the new state vector.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | CompliantMultibodyModel | n/a | yes | Positional argument `model`. |
| in | `x` | AbstractVector | n/a | yes | Positional argument `x`. |
| in | `t` | Real | n/a | yes | Positional argument `t`. |
| in | `dt` | Real | n/a | yes | Positional argument `dt`. |
| in | `dynamics_kwargs` | Vararg{Any} | n/a | no | Keyword argument `dynamics_kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `step_compliant_multibody_rk4`. Returns `_normalize_state_quaternions!(xn)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:634-634`
- `callees` → [[dynamics.cloth_multibody__normalize_state_quaternions__normalize_state_quaternions_bang|_normalize_state_quaternions!]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:642-642`
- `callees` → [[dynx.multibody_cloth_compliant_multibody_dynamics|compliant_multibody_dynamics]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:636-636`
- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:637-637`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:637-637`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:637-637`
<!-- vulcan:connections:end -->

## Limitations
Explicit, so the stiff joint springs impose a stability bound roughly `dt < 2 sqrt(m/k)`; at the default 5 kN/m and 1 kg that is a few milliseconds, and exceeding it diverges without warning.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 627.
