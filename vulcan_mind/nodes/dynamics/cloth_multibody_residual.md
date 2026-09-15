---
id: dynamics.cloth_multibody_residual
label: residual
kind: function
source:
  file: src/dynamics/multibody_cloth/cloth_multibody.jl
  symbol: residual
  lines:
  - 660
  - 660
inputs:
- id: z
  type: Any
  units: n/a
  required: true
  description: Positional argument `z`.
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
  description: Return value of `residual`. Returns `z .- x0 .- h .* compliant_multibody_dynamics(model,
    mid, t + 0.5h; dynamics_kwar`.
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

# residual

## Purpose
The nonlinear residual of the implicit-midpoint step, defined as a closure so Newton iteration can evaluate it at trial states.

## Theory & Math
$$
r(z) = z - x_0 - h\, f\!\left(\tfrac{x_0 + z}{2},\; t + \tfrac{h}{2}\right)
$$

## Design & Implementation
Given a candidate next state `z`, forms the midpoint `(x0 + z)/2`, evaluates the dynamics there at `t + h/2`, and returns `z - x0 - h f(mid)`. Captures `x0`, `h`, `t`, `model` and the dynamics keywords from the enclosing stepper.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `z` | Any | n/a | yes | Positional argument `z`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `residual`. Returns `z .- x0 .- h .* compliant_multibody_dynamics(model, mid, t + 0.5h; dynamics_kwar`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncx.heat_load_control__edg_solve_heat_load_switches|_edg_solve_heat_load_switches]] · `callees` → `callers` · call · `src/gnc/control/heat_load_control.jl:648-648`
- [[grp.src_gnc_control|gnc/control/]] · `members_out` → `callers` · call · `src/gnc/control/heat_load_control.jl:648-648`

**Downstream**

- `callees` → [[dynamics.cloth_multibody__normalize_state_quaternions__normalize_state_quaternions_bang|_normalize_state_quaternions!]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:667-667`
- `callees` → [[dynx.multibody_cloth_compliant_multibody_dynamics|compliant_multibody_dynamics]] · `callers` · call · `src/dynamics/multibody_cloth/cloth_multibody.jl:662-662`
<!-- vulcan:connections:end -->

## Limitations
Each evaluation is a full dynamics call; the finite-difference Jacobian calls it once per state dimension per Newton iteration.

## Provenance
Mapped from `src/dynamics/multibody_cloth/cloth_multibody.jl` line 660.
