---
id: vehicle.components_magnet
label: Magnet
kind: struct
source:
  file: src/vehicle/spacecraft/components.jl
  symbol: Magnet
  lines:
  - 33
  - 33
inputs:
- id: m
  type: MVector{3, Float64}
  units: n/a
  required: false
  description: Field `m` (default `MVector{3, Float64}(zeros(3))`).
- id: location
  type: MVector{3, Float64}
  units: n/a
  required: false
  description: Field `location` (default `MVector{3, Float64}(zeros(3))`).
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
  type: Magnet
  units: n/a
  description: Constructed `Magnet` (keyword constructor via @kwdef).
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

# Magnet

## Purpose
A fixed magnetic dipole mounted on the vehicle, the source term for magnetic torque against the local field.

## Design & Implementation
The smallest component in the module: a `@kwdef mutable struct` with `m`, the dipole moment in the body frame in ampere square metres, and `location`, the mounting point in the containing link's frame relative to that link's centre of mass in metres. Both default to zero three-vectors, so a default-constructed `Magnet` contributes nothing. Being mutable and `MVector`-backed lets a magnetorquer command be written into `m` each control tick without reallocating the component.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `m` | MVector{3, Float64} | n/a | no | Field `m` (default `MVector{3, Float64}(zeros(3))`). |
| in | `location` | MVector{3, Float64} | n/a | no | Field `location` (default `MVector{3, Float64}(zeros(3))`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Magnet | n/a | — | Constructed `Magnet` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/components.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only the dipole term is represented, so higher-order multipoles and any residual magnetism of the bus structure are outside this model; `location` is stored but does not enter a pure dipole torque, which depends on the moment and field alone.

## Provenance
Mapped from `src/vehicle/spacecraft/components.jl` line 33.
