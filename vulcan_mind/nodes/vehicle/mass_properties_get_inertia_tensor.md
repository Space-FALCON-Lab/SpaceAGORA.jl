---
id: vehicle.mass_properties_get_inertia_tensor
label: get_inertia_tensor
kind: function
source:
  file: src/vehicle/structure/mass_properties.jl
  symbol: get_inertia_tensor
  lines:
  - 71
  - 71
inputs:
- id: model
  type: SpacecraftModel
  units: n/a
  required: true
  description: Positional argument `model`.
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
  description: Return value of `get_inertia_tensor`. Returns `model.inertia_tensor`.
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

# get_inertia_tensor

## Purpose

`get_inertia_tensor` reads back an already-computed inertia tensor from a `SpacecraftModel` rather than recomputing one. The single-argument method returns `model.inertia_tensor`; the two-argument method returns `model.inertia_tensors[root_index]` for a specific assembly root.

## Design & Implementation

This is a pure accessor pairing with `update_inertia_tensor!`, which is what actually populates `model.inertia_tensors`. The indexed method asserts `root_index <= length(model.inertia_tensors)` with the message "Root index out of bounds" before indexing. Neither method copies, so the caller receives the stored `SMatrix{3,3,Float64}` value directly.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SpacecraftModel | n/a | yes | Positional argument `model`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `get_inertia_tensor`. Returns `model.inertia_tensor`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/structure/mass_properties.jl`
- [[vehicle.mass_properties_update_inertia_tensor|update_inertia_tensor]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:67-67`

**Downstream**

- `callees` → [[vehicle.mass_properties_set_inertia_tensor_bang|set_inertia_tensor!]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:86-86`
<!-- vulcan:connections:end -->

## Limitations

The returned tensor is only as current as the last `update_inertia_tensor!` call; nothing invalidates it after a link is rotated or repositioned, so stale values are possible. The bounds `@assert` covers the upper end only and would not catch a non-positive `root_index`, and assertions can be disabled at optimisation level 3. The two accessors read different fields (`inertia_tensor` versus `inertia_tensors`), so they can disagree if only one is maintained.

## Provenance
Mapped from `src/vehicle/structure/mass_properties.jl` line 71.
