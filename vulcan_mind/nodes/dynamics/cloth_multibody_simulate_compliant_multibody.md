---
id: dynamics.cloth_multibody_simulate_compliant_multibody
label: simulate_compliant_multibody
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: simulate_compliant_multibody
  lines:
  - 701
  - 701
inputs:
- id: model
  type: CompliantMultibodyModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: x0
  type: AbstractVector
  units: n/a
  required: true
  description: Positional argument `x0`.
- id: dt_s
  type: Real
  units: n/a
  required: true
  description: Keyword argument `dt_s`.
- id: duration_s
  type: Real
  units: n/a
  required: true
  description: Keyword argument `duration_s`.
- id: integrator
  type: Symbol
  units: n/a
  required: false
  description: Keyword argument `integrator` (default `:implicit_midpoint`).
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
  type: CompliantMultibodyTrajectory
  units: n/a
  description: Return value of `simulate_compliant_multibody`.
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

# simulate_compliant_multibody

## Purpose
Runs a fixed-step simulation of the compliant model over a duration and returns the trajectory.

## Design & Implementation
Validates positive `dt_s` and non-negative `duration_s`, builds the time grid with `0:dt:duration` and appends the end time if the grid falls short, normalises the initial state, and steps with either the implicit-midpoint or RK4 stepper selected by `integrator`, raising for any other symbol. Forwards `dynamics_kwargs` to every step.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | CompliantMultibodyModel | n/a | yes | Positional argument `model`. |
| in | `x0` | AbstractVector | n/a | yes | Positional argument `x0`. |
| in | `dt_s` | Real | n/a | yes | Keyword argument `dt_s`. |
| in | `duration_s` | Real | n/a | yes | Keyword argument `duration_s`. |
| in | `integrator` | Symbol | n/a | no | Keyword argument `integrator` (default `:implicit_midpoint`). |
| in | `dynamics_kwargs` | Vararg{Any} | n/a | no | Keyword argument `dynamics_kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CompliantMultibodyTrajectory | n/a | — | Return value of `simulate_compliant_multibody`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:709-709`
- `callees` → [[dynamics.cloth_multibody__normalize_state_quaternions__normalize_state_quaternions_bang|_normalize_state_quaternions!]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:718-718`
- `callees` → [[dynamics.cloth_multibody_compliantmultibodytrajectory|CompliantMultibodyTrajectory]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:729-729`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:715-715`
<!-- vulcan:connections:end -->

## Limitations
The integrator symbol is re-dispatched inside the loop rather than resolved once, and there is no event handling, output thinning or adaptive stepping.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 701.
