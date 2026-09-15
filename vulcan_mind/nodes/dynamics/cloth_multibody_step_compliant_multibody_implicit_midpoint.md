---
id: dynamics.cloth_multibody_step_compliant_multibody_implicit_midpoint
label: step_compliant_multibody_implicit_midpoint
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: step_compliant_multibody_implicit_midpoint
  lines:
  - 646
  - 646
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
- id: tol
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `tol` (default `1.0e-9`).
- id: max_iters
  type: Int
  units: n/a
  required: false
  description: Keyword argument `max_iters` (default `12`).
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
  description: Return value of `step_compliant_multibody_implicit_midpoint`. Returns
    `z .- x0 .- h .* compliant_multibody_dynamics(model, mid, t + 0.5h; dynamics_kwar`
    or `_normalize_state_quaternions!(xn)`.
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

# step_compliant_multibody_implicit_midpoint

## Purpose
Advances the compliant state one step with the implicit midpoint rule solved by damped Newton, allowing much larger steps than RK4 on stiff joints.

## Design & Implementation
Starts from a forward-Euler predictor. Up to `max_iters` (12) times it checks the residual against `tol` (1e-9), builds a dense Jacobian by forward differences with step `sqrt(eps) max(1, |x_j|)`, solves for the Newton direction, and backtracks `α` from one by halving until the residual maximum decreases, renormalising quaternions in each trial. If no `α` down to 1e-3 improves, iteration stops. It accepts the result if the residual is within `100 tol`, otherwise raises `ErrorException` with the residual.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | CompliantMultibodyModel | n/a | yes | Positional argument `model`. |
| in | `x` | AbstractVector | n/a | yes | Positional argument `x`. |
| in | `t` | Real | n/a | yes | Positional argument `t`. |
| in | `dt` | Real | n/a | yes | Positional argument `dt`. |
| in | `tol` | Float64 | n/a | no | Keyword argument `tol` (default `1.0e-9`). |
| in | `max_iters` | Int | n/a | no | Keyword argument `max_iters` (default `12`). |
| in | `dynamics_kwargs` | Vararg{Any} | n/a | no | Keyword argument `dynamics_kwargs` (variadic). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `step_compliant_multibody_implicit_midpoint`. Returns `z .- x0 .- h .* compliant_multibody_dynamics(model, mid, t + 0.5h; dynamics_kwar` or `_normalize_state_quaternions!(xn)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:655-655`
- `callees` → [[dynx.multibody_cloth_compliant_multibody_dynamics|compliant_multibody_dynamics]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:657-657`
<!-- vulcan:connections:end -->

## Limitations
The Jacobian is `13n` by `13n` dense and rebuilt every iteration by `13n` dynamics evaluations, so cost grows as `n³`; the acceptance tolerance being a hundred times the target means a stalled solve can return a state that only loosely satisfies the midpoint equations.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 646.
