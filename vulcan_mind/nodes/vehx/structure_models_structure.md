---
id: vehx.structure_models_structure
label: Structure
kind: struct
source:
  file: src/vehicle/structure/structure_models.jl
  symbol: Structure
  lines:
  - 1
  - 27
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
  description: Namespace exporting the assembly graph, mass property and geometry
    property helpers.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- structure
charts:
- vehx
origin: agent
---

# Structure

## Purpose
`Structure` is the namespace that gathers every routine which answers a question about the physical configuration of the vehicle: connectivity, centre of mass, inertia, mass budget, reference and projected areas, characteristic length, and surface normal and tangent directions. Dynamics, environment and analysis code import from this single module rather than from three separate files.

## Model & Assumptions
The module defines no types of its own. It closes over three source files that share a common contract: each takes a `SpacecraftModel` and optionally a `Link`, resolves the relevant assembly through the traversal helper, and reduces a property over that component. Because the assembly graph is included first, the mass and geometry files can call `traverse_bodies` without any forward declaration.

## Design & Implementation
Lines 6 and 7 use selective imports, pulling only `SpacecraftModel`, `Link` and `Joint` from the model module and only `rotate_to_inertial` from the kinematics module. That narrow surface is what keeps the structure layer free of a circular dependency on assembly, which is why `add_body!` cannot call the inertia update directly. The export list at lines 9 through 21 names twelve functions explicitly, which documents the public contract. The three includes at lines 23 to 25 are ordered so the graph traversal is defined before its consumers, and each uses `@__DIR__` for path safety.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `exports` | Module | n/a | — | Namespace exporting the assembly graph, mass property and geometry property helpers. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/structure/structure_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Because the helpers are split across included files rather than submodules, method names such as `get_spacecraft_reference_area` accumulate several methods whose applicability depends on model variant; the export list gives no indication of which overload is current. There is no runtime check that the includes were loaded in the intended order, and the module offers no caching layer, so repeated property queries repeat the full traversal.

## Provenance
Mapped from `src/vehicle/structure/structure_models.jl:1-27` and the three files it includes under `src/vehicle/structure/`.
