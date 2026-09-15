---
id: dynx.translational_position_kinematics_position_derivative
label: position_derivative
kind: function
source:
  file: src/dynamics/translational/position_kinematics.jl
  symbol: position_derivative
  lines:
  - 1
  - 5
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace through which the translational kinematics
    kernel is reached.
- id: velocity
  type: AbstractVector
  units: m/s
  required: true
  description: Inertial velocity of the spacecraft mass centre.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: pos_dot
  type: SVector{3,Float64}
  units: m/s
  description: Time derivative of the inertial position vector.
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

# position_derivative

## Purpose
`position_derivative` supplies the kinematic half of the translational state propagation: the time derivative of inertial position is the inertial velocity. Its companion in the same file, `zero_position_derivative`, returns a frozen-position rate for the control-only and force-only right-hand-side variants used in reduced-order studies.

## Theory & Math
The kinematic relation is the definition of velocity in an inertial frame,

$$\dot{\vec{r}}_{ii} = \vec{v}_{ii}$$

with $\vec{r}_{ii}$ in m and $\vec{v}_{ii}$ in m/s. Paired with $\dot{\vec{v}}_{ii} = \vec{F}_{net}/m$ this closes the six-state translational system

$$\frac{d}{dt}\begin{bmatrix}\vec{r}\\ \vec{v}\end{bmatrix} = \begin{bmatrix}\vec{v}\\ \vec{F}_{net}/m\end{bmatrix}$$

which, with $\vec{F}_{net} = -\mu m \vec{r}/r^{3}$, reduces to the two-body problem whose conserved quantities are specific energy $\varepsilon = v^{2}/2 - \mu/r$ in J/kg and specific angular momentum $\vec{h} = \vec{r}\times\vec{v}$ in m^2/s. Because the relation holds only in an inertial frame, a rotating-frame formulation would instead require $\dot{\vec{r}} = \vec{v} + \vec{\omega}\times\vec{r}$, with $\vec{\omega}$ the frame rate in rad/s.

## Model & Assumptions
The kernel assumes both position and velocity states are expressed in the same non-rotating, planet-centred inertial frame, so no transport term appears. It reads exactly three components and converts them to `Float64`, which fixes the state layout as a contiguous three-element position block followed by a three-element velocity block.

## Design & Implementation
The function is `@inline` and returns a concrete `SVector{3,Float64}` built from explicit `Float64` conversions, so it accepts dual numbers or state views without allocating a heap array and keeps the caller's broadcast assignment type-stable. Keeping this trivial relation behind a named export means the frame contract has one documented location and the reduced-order variants can substitute `zero_position_derivative` at the same call site.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace through which the translational kinematics kernel is reached. |
| in | `velocity` | AbstractVector | m/s | yes | Inertial velocity of the spacecraft mass centre. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `pos_dot` | SVector{3,Float64} | m/s | — | Time derivative of the inertial position vector. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.point_mass_dynamics_assign_full_translational_rhs_bang|assign_full_translational_rhs!]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:28-28`
- [[dynamics.point_mass_dynamics_assign_slow_translational_rhs_bang|assign_slow_translational_rhs!]] · `callees` → `callers` · call · `src/dynamics/translational/point_mass_dynamics.jl:39-39`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/dynamics/translational/position_kinematics.jl:4-4`
<!-- vulcan:connections:end -->

## Limitations
There is no frame validation: supplying a planet-fixed velocity with an inertial position silently omits the $\vec{\omega}\times\vec{r}$ transport term and produces a wrong trajectory. Components beyond the third are ignored, and non-finite velocity values propagate unchecked into the position rate.

## Provenance
Mapped from `src/dynamics/translational/position_kinematics.jl:1-5`, used by the right-hand-side assignments at lines 28 and 39 of `src/dynamics/translational/point_mass_dynamics.jl`.
