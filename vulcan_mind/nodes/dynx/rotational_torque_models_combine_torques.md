---
id: dynx.rotational_torque_models_combine_torques
label: combine_torques
kind: function
source:
  file: src/dynamics/rotational/torque_models.jl
  symbol: combine_torques
  lines:
  - 9
  - 14
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors namespace through which torque marshalling is reached.
- id: dynamic_torque
  type: SVector{3,Float64}
  units: N*m
  required: true
  description: Sum of environmental and disturbance torques in body axes.
- id: control_torque
  type: SVector{3,Float64}
  units: N*m
  required: true
  description: Commanded actuator torque in body axes.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: net_torque
  type: SVector{3,Float64}
  units: N*m
  description: Net body-frame torque passed to Euler's rotational equations.
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

# combine_torques

## Purpose
`combine_torques` forms the net body-frame torque that drives Euler's equations by summing the environmental contribution and the actuator command. It sits alongside `body_torque` and `body_angular_velocity`, the two conversion helpers in the same file that coerce loosely typed state slices into static three-vectors.

## Theory & Math
Torque superposition follows directly from the linearity of the angular momentum balance about the mass centre:

$$\vec{\tau}_{net} = \vec{\tau}_{dyn} + \vec{\tau}_{ctrl} = \sum_k \vec{\tau}_k$$

with every term in N*m and expressed in the same body axes. Substituted into Euler's equation this gives

$$\mathbf{J}\dot{\vec{\omega}}_b = \vec{\tau}_{dyn} + \vec{\tau}_{ctrl} - \vec{\omega}_b \times \mathbf{J}\vec{\omega}_b$$

The environmental term typically aggregates aerodynamic torque $\vec{r}_{cp/cm}\times\vec{F}_{aero}$, magnetic torque $\vec{m}\times\vec{B}$, gravity-gradient torque $\tfrac{3\mu}{r^{3}}\,\hat{n}\times\mathbf{J}\hat{n}$ with $\hat{n}$ the body-frame nadir unit vector, and solar pressure torque. In low Earth orbit these are of order $10^{-6}$ to $10^{-4}$ N*m, while reaction wheels deliver $10^{-3}$ to $10^{-1}$ N*m, so the control term normally dominates by two to three orders of magnitude.

## Model & Assumptions
The helper assumes both arguments are already resolved in body axes about the same reference point, the mass centre. It performs no saturation, so an actuator command exceeding hardware limits passes through unchanged and must be clipped upstream by the actuator model. Superposition is exact for rigid bodies; it is an approximation once structural flexibility couples the load paths.

## Design & Implementation
The function is `@inline` and typed on `SVector{3,Float64}` for both arguments, so the sum compiles to three floating-point additions with no allocation and no bounds checking inside the integrator. Keeping the addition behind a named function rather than inlining a `+` at each call site gives the graph a single, greppable point where the net torque contract is defined.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors namespace through which torque marshalling is reached. |
| in | `dynamic_torque` | SVector{3,Float64} | N*m | yes | Sum of environmental and disturbance torques in body axes. |
| in | `control_torque` | SVector{3,Float64} | N*m | yes | Commanded actuator torque in body axes. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `net_torque` | SVector{3,Float64} | N*m | — | Net body-frame torque passed to Euler's rotational equations. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/rotational/torque_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No saturation, deadband, quantisation or actuator lag is applied, and no check confirms that the two torques share a frame or reference point. Only two contributions can be combined per call, so aggregating many disturbance sources requires the caller to pre-sum them. Non-finite inputs propagate silently into the angular acceleration.

## Provenance
Mapped from `src/dynamics/rotational/torque_models.jl:9-14`, exported by `DynamicsRotational` at line 14 of `src/dynamics/rotational/rotational_models.jl`.
