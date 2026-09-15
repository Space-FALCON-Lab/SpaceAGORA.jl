---
id: vehx.actuators_thruster_models_module_thrustermodels
label: ThrusterModels
kind: struct
source:
  file: src/vehicle/actuators/thruster/thruster_models_module.jl
  symbol: ThrusterModels
  lines:
  - 1
  - 11
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Vehicle package namespace that loads this file and exports the symbol.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: exports
  type: Module
  units: n/a
  description: Namespace exporting BaseThrusterModel and SixAxisThrusterModel to the
    vehicle package.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- actuators
- thruster
charts:
- vehx
origin: agent
---

# ThrusterModels

## Purpose
`ThrusterModels` is the namespace wrapper that turns the raw thruster type definitions into a loadable unit of the vehicle package. It exists so that the concrete thruster descriptions can be included once, with a single controlled import surface, and referenced from actuator hooks, dynamic effectors and scenario builders without any of them reaching into a file path.

## Model & Assumptions
The module owns no state and performs no computation. Its contract is purely one of visibility: it brings `AbstractThrusterModel` into scope from the shared `AbstractTypes` module, adds `StaticArrays` and `LinearAlgebra` for the static matrix fields of the six-axis layout, and re-exports the two concrete types. Because the include happens inside the module body, the types are defined in this namespace and their subtype relationship to the shared abstract type is established at load time.

## Design & Implementation
Lines 1 through 11 form the whole file. Line 3 selectively imports only `AbstractThrusterModel` rather than the entire abstract type namespace, which keeps the symbol surface small and prevents accidental method capture. Line 7 declares the export list, and line 9 includes `thruster_models.jl` using `@__DIR__` so the path resolves correctly regardless of the working directory in which the package is loaded. Keeping the module shell separate from the type bodies also lets the type file be read or tested in isolation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `exports` | Module | n/a | — | Namespace exporting BaseThrusterModel and SixAxisThrusterModel to the vehicle package. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/actuators/thruster/thruster_models_module.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The module is a load-order dependency: it must be included after `AbstractTypes` or the subtype annotation fails, and this ordering is enforced only by the parent package file. It offers no constructor helpers, no registry of available thruster models, and no validation that the two exported types remain mutually consistent.

## Provenance
Mapped from `src/vehicle/actuators/thruster/thruster_models_module.jl:1-11`.
