---
id: vehicle.components_components
label: Components
kind: module
source:
  file: src/vehicle/spacecraft/components.jl
  symbol: Components
  lines:
  - 1
  - 1
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
  type: Any
  units: n/a
  description: Value produced by this symbol.
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

# Components

## Purpose
The module defining the physical hardware pieces a spacecraft model is assembled from — surface facets, thrusters, magnetic dipoles and reaction wheel assemblies.

## Design & Implementation
Imports StaticArrays and LinearAlgebra, defines the constant `I3` as a three-by-three static identity, and exports `Facet`, `Thruster`, `Magnet`, `ReactionWheelAssembly` and the `create_facet_list` builder. Every component is a `@kwdef mutable struct` with a complete set of defaults, so a caller can construct one naming only the fields that differ from a neutral part. Fixed-size `MVector` and `SMatrix` fields keep the components usable inside the integration hot path without heap traffic per evaluation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Value produced by this symbol. |
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
Components are mutable and hold live simulation state alongside fixed parameters — a `Thruster` carries both its `Isp` and its current `thrust` — so a component instance cannot be shared between two concurrently integrating vehicles without aliasing their state.

## Provenance
Mapped from `src/vehicle/spacecraft/components.jl` line 1.
