---
id: dynamics.point_mass_dynamics_assign_full_translational_rhs_bang
label: assign_full_translational_rhs!
kind: function
source:
  file: src/dynamics/translational/point_mass_dynamics.jl
  symbol: assign_full_translational_rhs!
  lines:
  - 22
  - 22
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
  type: Nothing
  units: n/a
  description: Return value of `assign_full_translational_rhs!`; mutates `du_view`
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

# assign_full_translational_rhs!

## Purpose

`assign_full_translational_rhs!` writes the complete translational derivative of a spacecraft into `du_view`: position rate, velocity rate and mass rate. It is the unrestricted variant used when position, velocity and propellant consumption all evolve on the same integration step.

## Design & Implementation

Three assignments make up the body. `du_view.pos .= position_derivative(sc_view.vel)` broadcasts the current velocity into the position slot; `du_view.vel .= acceleration_from_force(net_force, sc_view.mass)` divides the summed force by the current mass; `du_view.mass = mass_derivative(mass_rate)` stores the sanitised mass flow. The function is `@inline`, mutates `du_view` in place, and returns `nothing`. `acceleration_from_force` guards against a non-finite or effectively zero mass by returning a zero acceleration when `!isfinite(mass_f64) || abs(mass_f64) <= eps(Float64)`.

## Theory & Math

The assignments implement the point-mass translational state derivative

$$\dot{\mathbf{r}} = \mathbf{v}, \qquad \dot{\mathbf{v}} = \frac{\mathbf{F}_{\text{net}}}{m}, \qquad \dot{m} = \dot{m}_{\text{rate}}$$

where $\mathbf{r}$ is `sc_view` position, $\mathbf{v}$ is `sc_view.vel`, $\mathbf{F}_{\text{net}}$ is `net_force` in newtons, and $m$ is `sc_view.mass` in kilograms, so the acceleration is in metres per second squared.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `du_view` | Any | n/a | yes | Positional argument `du_view`. |
| in | `sc_view` | Any | n/a | yes | Positional argument `sc_view`. |
| in | `net_force` | AbstractVector{<:Real} | n/a | yes | Positional argument `net_force`. |
| in | `mass_rate` | Real | n/a | yes | Positional argument `mass_rate`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `assign_full_translational_rhs!`; mutates `du_view` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/translational/point_mass_dynamics.jl`
- [[simulation.dynamics_rhs__spacecraft_dynamics_flat_constellation_effector_queue__spacecraft_dynamics_flat_constellation_effector_queue_bang|_spacecraft_dynamics_flat_constellation_effector_queue!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1410-1410`
- [[simulation.dynamics_rhs_spacecraft_dynamics_explicit_remainder_bang|spacecraft_dynamics_explicit_remainder!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:2126-2126`
- [[simx.engine_dynamics_rhs_spacecraft_dynamics_bang|spacecraft_dynamics!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1757-1757`

**Downstream**

- `callees` → [[dynamics.point_mass_dynamics_mass_derivative|mass_derivative]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:30-30`
- `callees` → [[dynx.translational_point_mass_dynamics_acceleration_from_force|acceleration_from_force]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:29-29`
- `callees` → [[dynx.translational_position_kinematics_position_derivative|position_derivative]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:28-28`
<!-- vulcan:connections:end -->

## Limitations

`du_view` and `sc_view` are untyped arguments, so the function relies entirely on the caller supplying views with `pos`, `vel` and `mass` fields of matching length; a mismatch surfaces as a broadcast error at runtime rather than a method error. When mass falls to zero or goes non-finite the acceleration is silently zeroed instead of raising, which masks propellant-accounting bugs. The formulation is Newtonian in a single inertial frame and ignores the $\dot{m}\mathbf{v}$ momentum term, which is assumed already folded into `net_force`.

## Provenance
Mapped from `src/dynamics/translational/point_mass_dynamics.jl` line 22.
