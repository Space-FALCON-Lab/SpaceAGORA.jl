---
id: vehx.structure_assembly_graph_traverse_bodies
label: traverse_bodies
kind: function
source:
  file: src/vehicle/structure/assembly_graph.jl
  symbol: traverse_bodies
  lines:
  - 8
  - 31
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Vehicle package namespace that loads this file and exports the symbol.
- id: model
  type: SpacecraftModel
  units: n/a
  required: true
  description: Model supplying the joint list that defines connectivity between links.
- id: body
  type: Link
  units: n/a
  required: true
  description: Seed link from which the connected component is explored.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: bodies
  type: Vector{Link}
  units: n/a
  description: All links reachable from the seed through joints, including the seed
    itself.
- id: root_index
  type: Int
  units: n/a
  description: Index of the root link found within the component, used to select per-assembly
    quantities.
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

# traverse_bodies

## Purpose
`traverse_bodies` answers the question every structural property helper must ask first: which links belong to the same physical assembly as this one, and which root anchors them. Centre of mass, composite inertia, reference area, projected areas and vehicle length are all sums over a connected component, so each of those helpers begins by calling this traversal and then reduces over the returned vector.

## Model & Assumptions
The assembly is modelled as an undirected graph whose vertices are links and whose edges are joints. A joint is symmetric, so reachability is explored in both directions regardless of which link was recorded as the first or second endpoint. Exactly one vertex in a component is expected to carry the root flag, and its position within the model root list identifies the assembly for indexed quantities such as propellant mass or a stored inertia tensor.

## Design & Implementation
The implementation is a breadth-first search using an explicit `Link[]` queue and a `Set{Link}` of visited vertices, seeded with the starting body and marked visited before the loop so a self-joint cannot cause a repeat. Each iteration pops the front vertex, records the root index when the vertex is flagged as root, and then scans the entire joint list for edges incident on it. Both orientations are tested with an explicit membership check before enqueueing. After the queue drains, an assertion demands that a root was found, and the visited set is materialised with `collect`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `model` | SpacecraftModel | n/a | yes | Model supplying the joint list that defines connectivity between links. |
| in | `body` | Link | n/a | yes | Seed link from which the connected component is explored. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `bodies` | Vector{Link} | n/a | — | All links reachable from the seed through joints, including the seed itself. |
| out | `root_index` | Int | n/a | — | Index of the root link found within the component, used to select per-assembly quantities. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[dynamics.aerodynamic_wrench_models_calcforcetorque|calcForceTorque]] · `callees` → `callers` · call · `src/dynamics/coupled/aerodynamic_wrench_models.jl:593-593`
- [[gnc.closed_form_solution_closed_form_calculation|closed_form_calculation]] · `callees` → `callers` · call · `src/gnc/guidance/aerobraking/common/closed_form_solution.jl:156-156`
- [[vehicle.geometry_properties_get_sa_area|get_SA_area]] · `callees` → `callers` · call · `src/vehicle/structure/geometry_properties.jl:63-63`
- [[vehicle.geometry_properties_get_sc_area|get_SC_area]] · `callees` → `callers` · call · `src/vehicle/structure/geometry_properties.jl:88-88`
- [[vehicle.geometry_properties_get_spacecraft_length|get_spacecraft_length]] · `callees` → `callers` · call · `src/vehicle/structure/geometry_properties.jl:48-48`
- [[vehicle.mass_properties_get_spacecraft_mass|get_spacecraft_mass]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:104-104`
- [[vehicle.mass_properties_set_inertia_tensor_bang|set_inertia_tensor!]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:91-91`
- [[vehx.structure_geometry_properties_get_spacecraft_reference_area|get_spacecraft_reference_area]] · `callees` → `callers` · call · `src/vehicle/structure/geometry_properties.jl:9-9`
- [[vehx.structure_mass_properties_update_inertia_tensor_bang|update_inertia_tensor!]] · `callees` → `callers` · call · `src/vehicle/structure/mass_properties.jl:33-33`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/vehicle/structure/assembly_graph.jl:12-12`
<!-- vulcan:connections:end -->

## Limitations
The inner scan walks every joint for every popped vertex, so the cost is quadratic in the size of the model rather than linear in the number of edges; for the small assemblies this simulator handles that is acceptable but it does not scale. Returning a collected `Set` means the output order is unspecified, so any consumer that assumes a stable link ordering across calls is fragile. The root lookup uses `model.roots`, a plural field that the current `SpacecraftModel` no longer defines, so this path only works for model variants that still carry a root vector. A component with no root trips an assertion instead of returning a diagnosable value.

## Provenance
Mapped from `src/vehicle/structure/assembly_graph.jl:8-31`, included by the `Structure` module and called from both `mass_properties.jl` and `geometry_properties.jl`.
