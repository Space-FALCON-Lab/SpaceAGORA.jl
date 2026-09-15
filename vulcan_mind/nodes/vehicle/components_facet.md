---
id: vehicle.components_facet
label: Facet
kind: struct
source:
  file: src/vehicle/spacecraft/components.jl
  symbol: Facet
  lines:
  - 9
  - 9
inputs:
- id: area
  type: Float64
  units: n/a
  required: false
  description: Field `area` (default `0.0`).
- id: attitude
  type: MVector{4, Float64}
  units: n/a
  required: false
  description: Field `attitude` (default `@MVector [0.0, 0.0, 0.0, 1.0]`).
- id: normal_vector
  type: MVector{3, Float64}
  units: n/a
  required: false
  description: Field `normal_vector` (default `@MVector [1.0, 0.0, 0.0]`).
- id: cp
  type: MVector{3, Float64}
  units: n/a
  required: false
  description: Field `cp` (default `@MVector [0.0, 0.0, 0.0]`).
- id: rho
  type: Float64
  units: n/a
  required: false
  description: Field `ρ` (default `0.0`).
- id: delta
  type: Float64
  units: n/a
  required: false
  description: Field `δ` (default `0.0`).
- id: name
  type: String
  units: n/a
  required: false
  description: Field `name` (default `""`).
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
  type: Facet
  units: n/a
  description: Constructed `Facet` (keyword constructor via @kwdef).
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

# Facet

## Purpose
One flat surface patch of the vehicle, carrying the geometry and optical coefficients that aerodynamic and radiation wrench models integrate over.

## Design & Implementation
A `@kwdef mutable struct` of seven fields. `area` is in square metres, `attitude` is a four-element quaternion relative to the body frame defaulting to identity, and `normal_vector` is the patch normal expressed in the facet's own frame, defaulting to the x direction because a flat plate is modelled with its normal along x. `cp` is the centre of pressure relative to the centre of mass of the containing `Link`, expressed in that link's frame, so torque contributions compose correctly through the multibody tree. `ρ` and `δ` are the diffuse and specular reflection coefficients, and `name` identifies the patch in diagnostics.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `area` | Float64 | n/a | no | Field `area` (default `0.0`). |
| in | `attitude` | MVector{4, Float64} | n/a | no | Field `attitude` (default `@MVector [0.0, 0.0, 0.0, 1.0]`). |
| in | `normal_vector` | MVector{3, Float64} | n/a | no | Field `normal_vector` (default `@MVector [1.0, 0.0, 0.0]`). |
| in | `cp` | MVector{3, Float64} | n/a | no | Field `cp` (default `@MVector [0.0, 0.0, 0.0]`). |
| in | `rho` | Float64 | n/a | no | Field `ρ` (default `0.0`). |
| in | `delta` | Float64 | n/a | no | Field `δ` (default `0.0`). |
| in | `name` | String | n/a | no | Field `name` (default `""`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Facet | n/a | — | Constructed `Facet` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/components.jl`
- [[vehicle.components_create_facet_list|create_facet_list]] · `callees` → `callers` · call · `src/vehicle/spacecraft/components.jl:82-82`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`ρ` and `δ` default to zero, giving a fully absorbing surface, and nothing checks that their sum stays at or below one, so a miscopied material table can produce a facet that reflects more than it receives.

## Provenance
Mapped from `src/vehicle/spacecraft/components.jl` line 9.
