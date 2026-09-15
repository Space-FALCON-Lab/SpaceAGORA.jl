---
id: dynx.ftm_thruster_models_thrustermodels
label: ThrusterModels
kind: struct
source:
  file: src/dynamics/coupled/force_torque_models/thruster_models.jl
  symbol: ThrusterModels
  lines:
  - 1
  - 5
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: DynamicEffectors parent namespace under which the thruster type is
    made visible to effector lists.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: thruster_types
  type: Module
  units: n/a
  description: Namespace re-exporting BaseThrusterModel so propulsive effectors can
    be attached alongside gravity and aerodynamic models.
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

# ThrusterModels

## Purpose
`ThrusterModels` is the effector-side re-export of `BaseThrusterModel`. It exists so scenario code assembling a coupled effector list can name the propulsion type from the same namespace as the gravity, aerodynamic and perturbation effectors, while the thruster implementation itself stays in the vehicle propulsion package.

## Theory & Math
A thruster contributes both a body-frame force and the torque produced by its offset from the mass centre:

$$\vec{F}_b = T \hat{u}_b, \qquad \vec{\tau}_b = \vec{r}_{t/cm} \times \vec{F}_b, \qquad \dot{m} = -\frac{T}{I_{sp} g_0}$$

with $T$ the thrust magnitude in N, $\hat{u}_b$ the unit thrust direction in body axes, $\vec{r}_{t/cm}$ the thruster location relative to the mass centre in m, $I_{sp}$ the specific impulse in s and $g_0 = 9.80665$ m/s^2. For a pressure-fed or regulated engine the delivered thrust varies with ambient back-pressure as

$$T = \dot{m}\, v_e + (p_e - p_a) A_e$$

with $v_e$ the exhaust velocity in m/s, $p_e$ the nozzle exit pressure and $p_a$ the ambient pressure in Pa, and $A_e$ the exit area in m^2. In vacuum $p_a \to 0$ and thrust reaches its vacuum value.

## Model & Assumptions
The shim assumes the upstream thruster model owns throttle logic, mass-flow bookkeeping and duty-cycle state. Thrust is treated as instantaneously commanded, so ignition and shutdown transients shorter than the integrator step are not resolved. Mass depletion must be integrated through the translational mass state for the rocket equation to hold over long burns.

## Design & Implementation
The module body is a single `using ...ThrusterModels: BaseThrusterModel` plus a matching export. The submodule name intentionally shadows the upstream package name inside the effector namespace so effector code reads consistently, while the type identity is unchanged and dispatch continues to resolve to the upstream definitions.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | DynamicEffectors parent namespace under which the thruster type is made visible to effector lists. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `thruster_types` | Module | n/a | — | Namespace re-exporting BaseThrusterModel so propulsive effectors can be attached alongside gravity and aerodynamic models. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.dynamics|DynamicEffectors]] · `api` → `module_api` · call · `src/dynamics/coupled/force_torque_models/thruster_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No `calcForceTorque` method is defined here, so the propulsive wrench is produced upstream; this file only affects name resolution. Plume impingement on deployed structures, thrust-vector misalignment and pressurant blowdown are outside the scope of the re-exported base type.

## Provenance
Mapped from `src/dynamics/coupled/force_torque_models/thruster_models.jl:1-5`.
