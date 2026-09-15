---
id: vehicle.assembly_assembly
label: Assembly
kind: module
source:
  file: src/vehicle/spacecraft/assembly.jl
  symbol: Assembly
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

# Assembly

## Purpose

`Assembly` is the module that mutates a `SpacecraftModel` and its `Link`s into a complete vehicle definition. It exports the five in-place builders `add_body!`, `add_joint!`, `add_facet!`, `add_magnet!` and `add_thruster!`, which append rigid bodies, articulation joints, SRP facets, magnetic dipoles and thrusters onto the growing model.

## Design & Implementation

The module pulls in `LinearAlgebra` and `StaticArrays` plus the sibling modules `SpacecraftModels` and `Components`, and every exported entry point works by `push!`ing into a vector field of the model or link, so the container types stay stable. A source comment records that `add_dynamic_effector!` was deliberately removed because appending to a tuple would change the parent struct's type; dynamic effectors must instead be passed as a tuple to the `SpacecraftModel` constructor. Inertia-tensor recomputation is likewise left out to avoid a circular dependency on `Analysis`, and is expected to be triggered by a higher-level assembly script.

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

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/assembly.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations

Assembly is purely additive: there is no remove or replace operation, and no validation that the resulting topology is a connected tree. Because `add_body!` only seeds the root inertia tensor with the identity, the model's mass properties are wrong until an external routine recomputes them. All functions mutate shared vectors without locking, so concurrent assembly of the same model is unsafe.

## Provenance
Mapped from `src/vehicle/spacecraft/assembly.jl` line 1.
