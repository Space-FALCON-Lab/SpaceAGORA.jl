---
id: core.abstract_types_abstractforcetorquemodel
label: AbstractForceTorqueModel
kind: struct
source:
  file: src/core/types/abstract_types.jl
  symbol: AbstractForceTorqueModel
  lines:
  - 11
  - 11
inputs:
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
  type: AbstractForceTorqueModel
  units: n/a
  description: Abstract supertype `AbstractForceTorqueModel`; no fields.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# AbstractForceTorqueModel

## Purpose
Root abstract type for every dynamic effector that contributes force and/or torque to a spacecraft's translational or rotational equations of motion. It is the supertype declared for gravity models, aerodynamic coefficient models, magnetic actuators and structural reaction effectors, and it is the type the dynamics RHS dispatches on when it accumulates wrenches.

## Design & Implementation
Declared as `abstract type AbstractForceTorqueModel end` inside `module AbstractTypes` with no fields, no parameters and no default methods; it is re-exported from `SpaceAGORA` via `@doc`-forwarding. Concrete subtypes in the repo include `InverseSquaredGravityModel`, `InverseSquaredJ2GravityModel`, `GravitationalHarmonicsModel{P<:AbstractPlanet}`, `NBodyGravityModel`, `ConstantGravityModel`, `AerodynamicCoefficientConstant`, `AerodynamicCoefficientfM`, `AerodynamicCoefficientNoBallisticFlight`, `EddyCurrentDampingModel`, `MagneticTorqueRodModel`, `LVLHCascadeAttitudeControlModel` and `RobotArmReactionEffector`. The extension contract is documented in `src/SpaceAGORA.jl`: implement either the legacy `calcForceTorque(model, x, p, i) -> (force_n, torque_n_m)` or the preferred pure-function `wrench(model, x::StateSample, env::EnvironmentSample, t::Float64) -> (force_ii, torque_body)`, optionally with `environment_requirements(model)` and `solver_partition(model)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractForceTorqueModel | n/a | — | Abstract supertype `AbstractForceTorqueModel`; no fields. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/abstract_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Nothing is enforced at the type level: a subtype that implements neither `calcForceTorque` nor `wrench` compiles fine and fails only with a `MethodError` at the first RHS evaluation. Units (newtons, newton-metres, inertial force vs body torque) are a documentation convention, not checked. Effectors are stored in heterogeneous containers, so dispatch on this abstract type happens dynamically per effector per RHS call.

## Provenance
Mapped from `src/core/types/abstract_types.jl` line 11.
