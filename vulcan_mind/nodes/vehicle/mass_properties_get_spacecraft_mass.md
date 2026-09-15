---
id: vehicle.mass_properties_get_spacecraft_mass
label: get_spacecraft_mass
kind: function
source:
  file: src/vehicle/structure/mass_properties.jl
  symbol: get_spacecraft_mass
  lines:
  - 101
  - 101
inputs:
- id: model
  type: SpacecraftModel
  units: n/a
  required: true
  description: Positional argument `model`.
- id: dry
  type: Any
  units: n/a
  required: false
  description: Keyword argument `dry` (default `false`).
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
  description: 'Return value of `get_spacecraft_mass`. Returns `length(model.roots)
    == 1 ? masses[1] : masses`.'
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

# get_spacecraft_mass

## Purpose

`get_spacecraft_mass` totals the mass of a spacecraft assembly in kilograms, with an optional `dry` keyword that excludes propellant. Three methods cover the whole model, a single assembly reached from one `Link`, and an explicit `bodies` vector paired with a `root_index`.

## Design & Implementation

The innermost method sums `b.m` over `bodies` and, unless `dry=true`, adds `model.prop_mass[root_index]`. The `Link` method first calls `traverse_bodies(model, body)` to collect the connected assembly and its root index, then forwards. The whole-model method loops over `model.roots`, traverses each, and pushes each total into a `Float64[]`; it then returns `masses[1]` when there is exactly one root and the full vector otherwise.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `model` | SpacecraftModel | n/a | yes | Positional argument `model`. |
| in | `dry` | Any | n/a | no | Keyword argument `dry` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `get_spacecraft_mass`. Returns `length(model.roots) == 1 ? masses[1] : masses`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/structure/mass_properties.jl`
- [[vehicle.mass_properties_set_inertia_tensor_bang|set_inertia_tensor!]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:96-96`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:105-105`
- `callees` → [[vehx.structure_assembly_graph_traverse_bodies|traverse_bodies]] · `callers` · call · `src/vehicle/structure/mass_properties.jl:104-104`
<!-- vulcan:connections:end -->

## Limitations

The whole-model method has a union return type, `Float64` or `Vector{Float64}` depending on the number of roots, which forces callers to branch and defeats type inference. It also indexes `masses[1]` without checking that `model.roots` is non-empty. The inner sum allocates an intermediate array from the comprehension `[b.m for b in bodies]` on every call, which matters when it runs inside a dynamics right-hand side.

## Provenance
Mapped from `src/vehicle/structure/mass_properties.jl` line 101.
