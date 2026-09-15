---
id: dynamics.point_mass_dynamics_assign_force_only_translational_rhs_bang
label: assign_force_only_translational_rhs!
kind: function
source:
  file: src/dynamics/translational/point_mass_dynamics.jl
  symbol: assign_force_only_translational_rhs!
  lines:
  - 57
  - 57
inputs:
- id: du_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `du_view`.
- id: sc_view
  type: Any
  units: n/a
  required: true
  description: Positional argument `sc_view`.
- id: net_force
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `net_force`.
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
  type: Nothing
  units: n/a
  description: Return value of `assign_force_only_translational_rhs!`; mutates `du_view`
    in place. Returns `nothing`.
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

# assign_force_only_translational_rhs!

## Purpose

`assign_force_only_translational_rhs!` is the most restricted of the four translational assignment variants: it computes only the acceleration due to `net_force`, zeroing both the position derivative and the mass derivative.

## Design & Implementation

`du_view.pos .= zero_position_derivative()` and `du_view.mass = 0.0` blank the two channels this variant does not own, while `du_view.vel .= acceleration_from_force(net_force, sc_view.mass)` divides the summed force by the current spacecraft mass to give an acceleration in metres per second squared. There is no `mass_rate` parameter in the signature. The function is `@inline`, mutates `du_view`, and returns `nothing`. It pairs with `assign_control_only_translational_rhs!` as the non-propulsive force channel of the split.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `du_view` | Any | n/a | yes | Positional argument `du_view`. |
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `net_force` | AbstractVector{<:Real} | n/a | yes | Positional argument `net_force`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `assign_force_only_translational_rhs!`; mutates `du_view` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/translational/point_mass_dynamics.jl`
- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1352-1352`
- [[simulation.dynamics_rhs_spacecraft_dynamics_implicit_atmosphere_bang|spacecraft_dynamics_implicit_atmosphere!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2017-2017`

**Downstream**

- `callees` → [[dynamics.position_kinematics_zero_position_derivative|zero_position_derivative]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:62-62`
- `callees` → [[dynx.translational_point_mass_dynamics_acceleration_from_force|acceleration_from_force]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:63-63`
<!-- vulcan:connections:end -->

## Limitations

Only meaningful within an operator-splitting scheme; used alone it neither moves the vehicle nor consumes propellant. Because it shares `acceleration_from_force`, a `sc_view.mass` that is non-finite or below `eps(Float64)` produces a zero acceleration with no diagnostic. The untyped `du_view` and `sc_view` arguments mean the required `pos`, `vel` and `mass` fields are a convention rather than a checked contract.

## Provenance
Mapped from `src/dynamics/translational/point_mass_dynamics.jl` line 57.
