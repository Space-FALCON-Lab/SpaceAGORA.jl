---
id: vehicle.mass_properties_set_inertia_tensor_bang
label: set_inertia_tensor!
kind: function
source:
  file: src/vehicle/structure/mass_properties.jl
  symbol: set_inertia_tensor!
  lines:
  - 90
  - 90
inputs:
- id: model
  type: SpacecraftModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: body
  type: Link
  units: n/a
  required: true
  description: Positional argument `body`.
- id: inertia_tensor
  type: SMatrix{3, 3, Float64}
  units: n/a
  required: true
  description: Positional argument `inertia_tensor`.
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
  description: Return value of `set_inertia_tensor!`; mutates `model` in place.
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

# set_inertia_tensor!

## Purpose

`set_inertia_tensor!` lets a caller override the stored inertia tensor for the assembly that contains `body`, bypassing the parallel-axis computation in `update_inertia_tensor`. It is the escape hatch for models whose inertia comes from a CAD export or a measured value rather than from the link geometry.

## Design & Implementation

The function calls `traverse_bodies(model, body)` to resolve which assembly `body` belongs to, discards the returned `bodies` vector and keeps only `root_index`, then writes the supplied `inertia_tensor::SMatrix{3,3,Float64}` into `model.inertia_tensors[root_index]`. It mutates `model` and returns the assignment result rather than an explicit value.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SpacecraftModel | n/a | yes | Positional argument `model`. |
| in | `body` | Link | n/a | yes | Positional argument `body`. |
| in | `inertia_tensor` | SMatrix{3, 3, Float64} | n/a | yes | Positional argument `inertia_tensor`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `set_inertia_tensor!`; mutates `model` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/structure/mass_properties.jl`
- [[vehicle.mass_properties_get_inertia_tensor|get_inertia_tensor]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:86-86`

**Downstream**

- `callees` → [[vehicle.mass_properties_get_spacecraft_mass|get_spacecraft_mass]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:96-96`
- `callees` → [[vehx.structure_assembly_graph_traverse_bodies|traverse_bodies]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:91-91`
<!-- vulcan:connections:end -->

## Limitations

Unlike `update_inertia_tensor!`, this function does not grow `model.inertia_tensors` when `root_index` exceeds its current length, so setting the tensor for a freshly added root raises a `BoundsError`. No check enforces that the supplied matrix is symmetric or positive definite, so a physically impossible inertia can be installed and will only surface later as an integration failure. Any subsequent `update_inertia_tensor!` call on the same root silently overwrites the manual value.

## Provenance
Mapped from `src/vehicle/structure/mass_properties.jl` line 90.
