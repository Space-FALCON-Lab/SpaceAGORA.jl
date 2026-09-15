---
id: dynx.translational_point_mass_dynamics_acceleration_from_force
label: acceleration_from_force
kind: function
source:
  file: src/dynamics/translational/point_mass_dynamics.jl
  symbol: acceleration_from_force
  lines:
  - 1
  - 16
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace through which the translational point-mass
    kernel is reached.
- id: net_force
  type: AbstractVector
  units: N
  required: true
  description: Sum of all effector forces on the spacecraft, resolved in inertial
    axes.
- id: mass
  type: Real
  units: kg
  required: true
  description: Current spacecraft wet mass.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: accel_ii
  type: SVector{3,Float64}
  units: m/s^2
  description: Inertial translational acceleration, or a zero vector when the mass
    is non-finite or effectively zero.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- dynamics
charts:
- dynx
origin: agent
---

# acceleration_from_force

## Purpose
`acceleration_from_force` converts the accumulated effector force into a translational acceleration. Every translational right-hand-side assignment in the file routes through it, so the guard against a degenerate mass exists in exactly one place rather than being repeated at four call sites.

## Theory & Math
The kernel is Newton's second law for a variable-mass point in the form used by the integrator,

$$\vec{a}_{ii} = \frac{\vec{F}_{net}}{m}$$

with $\vec{F}_{net}$ in N, $m$ in kg and $\vec{a}_{ii}$ in m/s^2. The net force aggregates gravity $-\mu m \vec{r}/r^{3}$, aerodynamic drag $-\tfrac{1}{2}\rho \lVert\vec{v}_{rel}\rVert^{2} S C_D \hat{v}_{rel}$, solar radiation pressure, third-body attraction and thrust. Note that this is the constant-mass form: for a rocket the momentum balance is

$$m\frac{d\vec{v}}{dt} = \vec{F}_{ext} + \dot{m}\,\vec{v}_{e}$$

and the exhaust term is folded into $\vec{F}_{net}$ by the thruster effector rather than appearing here, which is why mass depletion is handled separately by `mass_derivative`. Companion state propagation is $\dot{\vec{r}} = \vec{v}$ from `position_derivative`, closing the six-state translational system.

## Model & Assumptions
The kernel assumes the force is already resolved in the same inertial frame as the velocity state and that the mass is the instantaneous total. When the mass is non-finite or its magnitude is at or below `eps(Float64)` the function returns a zero acceleration instead of raising, which keeps a propellant-exhausted or misconfigured spacecraft from poisoning the whole ensemble state with NaNs mid-integration.

## Design & Implementation
It is `@inline`, accepts any `AbstractVector{<:Real}` force so state views and static arrays both work, and computes a single reciprocal `inv(mass_f64)` that is then applied component-wise rather than performing three divisions. The return type is a concrete `SVector{3,Float64}`, which keeps the caller's broadcast assignment `du_view.vel .= ...` allocation-free.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace through which the translational point-mass kernel is reached. |
| in | `net_force` | AbstractVector | N | yes | Sum of all effector forces on the spacecraft, resolved in inertial axes. |
| in | `mass` | Real | kg | yes | Current spacecraft wet mass. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `accel_ii` | SVector{3,Float64} | m/s^2 | — | Inertial translational acceleration, or a zero vector when the mass is non-finite or effectively zero. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.point_mass_dynamics_assign_control_only_translational_rhs_bang|assign_control_only_translational_rhs!]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:52-52`
- [[dynamics.point_mass_dynamics_assign_force_only_translational_rhs_bang|assign_force_only_translational_rhs!]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:63-63`
- [[dynamics.point_mass_dynamics_assign_full_translational_rhs_bang|assign_full_translational_rhs!]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:29-29`
- [[dynamics.point_mass_dynamics_assign_slow_translational_rhs_bang|assign_slow_translational_rhs!]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:40-40`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:5-5`
<!-- vulcan:connections:end -->

## Limitations
Silently returning zero on a degenerate mass hides configuration errors: a scenario with mass mistakenly set to zero produces a spacecraft that simply coasts rather than failing loudly. Only the first three components of the force vector are read, so a longer vector is truncated without warning. Relativistic and variable-mass thrust terms are not represented here.

## Provenance
Mapped from `src/dynamics/translational/point_mass_dynamics.jl:1-16`, used by the four right-hand-side assignments at lines 22, 34, 45 and 57 of the same file.
