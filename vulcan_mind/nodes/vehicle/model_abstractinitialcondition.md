---
id: vehicle.model_abstractinitialcondition
label: AbstractInitialCondition
kind: struct
source:
  file: src/vehicle/spacecraft/model.jl
  symbol: AbstractInitialCondition
  lines:
  - 15
  - 15
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
  type: AbstractInitialCondition
  units: n/a
  description: Abstract supertype `AbstractInitialCondition`; no fields.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
charts:
- vehicle
origin: agent
---

# AbstractInitialCondition

## Purpose
Abstract supertype under which every initial-state description for a spacecraft is grouped, so `SpacecraftModel.initial_condition` can hold either Keplerian or Cartesian inputs and the simulation setup can dispatch on the concrete type.

## Design & Implementation
Declared as `abstract type AbstractInitialCondition end` with two concrete subtypes in this file: `InitialCondition` (six classical orbital elements `a, e, i, ω, Ω, ν` plus attitude quaternion `q` and angular velocity `ang_vel`) and `CartesianInitialCondition` (inertial `pos`, `vel`, `q`, `ang_vel`). `SpacecraftModel` stores the field with the abstract type and defaults it to `InitialCondition()`. Downstream state-vector construction dispatches on the concrete subtype to produce `pos`/`vel`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractInitialCondition | n/a | — | Abstract supertype `AbstractInitialCondition`; no fields. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/model.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because `SpacecraftModel.initial_condition` is typed abstractly, accessing its fields inside hot loops incurs dynamic dispatch; this is acceptable since it is only read at setup. No interface functions are declared on the abstract type, so a new subtype must be wired into every consumer by hand.

## Provenance
Mapped from `src/vehicle/spacecraft/model.jl` line 15.
