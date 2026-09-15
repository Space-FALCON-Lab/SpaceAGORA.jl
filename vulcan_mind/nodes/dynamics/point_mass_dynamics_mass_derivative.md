---
id: dynamics.point_mass_dynamics_mass_derivative
label: mass_derivative
kind: function
source:
  file: src/dynamics/translational/point_mass_dynamics.jl
  symbol: mass_derivative
  lines:
  - 18
  - 18
inputs:
- id: mass_rate
  type: Real
  units: n/a
  required: true
  description: Positional argument `mass_rate`.
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
  type: Float64
  units: n/a
  description: Return value of `mass_derivative`.
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

# mass_derivative

## Purpose

`mass_derivative(mass_rate::Real)::Float64` sanitises the spacecraft mass rate of change before it is written into a state derivative. It converts the incoming value to `Float64` and substitutes `0.0` whenever the value is not finite.

## Design & Implementation

The body is the single expression `isfinite(mass_rate) ? Float64(mass_rate) : 0.0`, and the function is marked `@inline` so it disappears into the calling right-hand side. The guard means a `NaN` or `Inf` produced by a misbehaving propulsion or mass-flow model cannot poison the integrator's mass state; the mass simply stops changing for that step. The declared return type `::Float64` also pins inference for the `du_view.mass` assignment in the callers `assign_full_translational_rhs!` and `assign_control_only_translational_rhs!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `mass_rate` | Real | n/a | yes | Positional argument `mass_rate`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `mass_derivative`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.point_mass_dynamics_assign_control_only_translational_rhs_bang|assign_control_only_translational_rhs!]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:53-53`
- [[dynamics.point_mass_dynamics_assign_full_translational_rhs_bang|assign_full_translational_rhs!]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:30-30`
- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/translational/point_mass_dynamics.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:19-19`
<!-- vulcan:connections:end -->

## Limitations

Silently replacing a non-finite rate with zero hides the upstream fault: no warning is emitted and no counter is incremented, so a persistently `NaN` thruster mass flow looks like a coasting vehicle. The sign convention is not enforced either, so a positive `mass_rate` will grow the spacecraft mass without complaint. Units are whatever the caller supplies, nominally kilograms per second.

## Provenance
Mapped from `src/dynamics/translational/point_mass_dynamics.jl` line 18.
