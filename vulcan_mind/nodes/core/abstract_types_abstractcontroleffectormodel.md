---
id: core.abstract_types_abstractcontroleffectormodel
label: AbstractControlEffectorModel
kind: struct
source:
  file: src/core/types/abstract_types.jl
  symbol: AbstractControlEffectorModel
  lines:
  - 19
  - 19
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
  type: AbstractControlEffectorModel
  units: n/a
  description: Abstract supertype `AbstractControlEffectorModel`; no fields.
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

# AbstractControlEffectorModel

## Purpose
Supertype for control-layer effectors that turn navigation or guidance outputs into runtime actuation commands, as opposed to `AbstractForceTorqueModel` effectors that directly produce physical wrenches. Concrete examples are `RPOMPCControlModel`, `AerobrakingEnergyDepletionControlModel`, `SolarPanelAngleOfAttackControlModel`, `MagneticMomentumManagerModel` and `RobotArmControlEffector`.

## Design & Implementation
An empty `abstract type AbstractControlEffectorModel end` with no fields or parameters. The control pipeline in `src/gnc/control/propulsive_maneuvers.jl` uses it as a dispatch anchor for fallback hooks: `calcControlMassFlowRate(::AbstractControlEffectorModel, u::AbstractVector, p::ODEParams, i::Int64, t::Float64)::Float64` returns `0.0` and `calcReactionWheelTorque(...)` returns `nothing` unless a concrete subtype overrides them. Instances are collected in `control_model.control_effectors`, which `_validate_ensemble_uncoupled` inspects when deciding whether a constellation may be split into independent ensemble members.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractControlEffectorModel | n/a | — | Abstract supertype `AbstractControlEffectorModel`; no fields. |
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
Both fallback hooks also exist for an untyped `controlModel` argument, so subtyping this abstract type is optional in practice; a control effector that omits it still works, which erodes the value of the type as a contract. No method declares what a control effector must implement, so the required entry points are discoverable only by reading the control dispatch code. There is no way to declare whether an effector couples multiple spacecraft, which is why ensemble validation must reject all control effectors conservatively.

## Provenance
Mapped from `src/core/types/abstract_types.jl` line 19.
