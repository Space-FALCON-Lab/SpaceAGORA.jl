---
id: core.abstract_types_abstractthrustermodel
label: AbstractThrusterModel
kind: struct
source:
  file: src/core/types/abstract_types.jl
  symbol: AbstractThrusterModel
  lines:
  - 58
  - 58
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
  type: AbstractThrusterModel
  units: n/a
  description: Abstract supertype `AbstractThrusterModel`; no fields.
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

# AbstractThrusterModel

## Purpose
Supertype for thruster hardware descriptions used by actuator and control layers to convert a commanded impulse or force into mass flow, thrust and torque. Concrete subtypes are `BaseThrusterModel` (a single thrust axis with specific impulse and thrust level) and `SixAxisThrusterModel` (independent force and torque authority on all six axes).

## Design & Implementation
Declared as `abstract type AbstractThrusterModel end` with no fields. The concrete types live under `src/vehicle/actuators/thruster/` in the `ThrusterModels` module, are re-exported through `src/dynamics/coupled/force_torque_models.jl` and are consumed by `src/gnc/control/control_hooks.jl`. Nothing in `src` currently dispatches on the abstract type itself; propulsive maneuver code and scenario builders (for example `scenario_builders.jl:410`) construct and pass the concrete `BaseThrusterModel` directly.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractThrusterModel | n/a | — | Abstract supertype `AbstractThrusterModel`; no fields. |
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
The abstract type carries no method contract, so a custom thruster subtype gains nothing from subtyping except documentation intent; control hooks dispatch on `BaseThrusterModel` and `SixAxisThrusterModel` by name. Thrust and Isp units (N, s) and whether thrust is expressed in body or inertial frame are conventions fixed in the concrete types. Mass flow rate integration relies on the control effector calling `calcControlMassFlowRate`, not on anything the thruster type guarantees.

## Provenance
Mapped from `src/core/types/abstract_types.jl` line 58.
