---
id: vehx.structure_geometry_properties_get_spacecraft_reference_area
label: get_spacecraft_reference_area
kind: function
source:
  file: src/vehicle/structure/geometry_properties.jl
  symbol: get_spacecraft_reference_area
  lines:
  - 6
  - 13
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
  description: Vehicle whose assemblies are traversed to accumulate facet areas.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: area
  type: Union{Float64,Vector{Float64}}
  units: m^2
  description: Total reference area, scalar for a single assembly or one entry per
    assembly.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- vehicle
- structure
- geometry
charts:
- vehx
origin: agent
---

# get_spacecraft_reference_area

## Purpose
`get_spacecraft_reference_area` reports the aerodynamic and solar reference area of the vehicle. Drag and solar radiation pressure effectors need an area to scale their coefficients against, and mission analysis reports quote area-to-mass ratio, so this is the single point at which geometry is condensed into one number per assembly.

## Theory & Math
The drag acceleration used downstream is $\mathbf{a}_d = -\tfrac{1}{2}\rho\,\frac{C_d A}{m}\,\|\mathbf{v}_{rel}\|\,\mathbf{v}_{rel}$, so the returned $A = \sum_k A_k$ enters linearly through the ballistic coefficient $B = m/(C_d A)$.

## Model & Assumptions
Area is accumulated from the facet list attached to each link in an assembly. The method shown here is the whole-model entry point: it iterates the model roots, expands each into its connected component with `traverse_bodies`, and delegates the actual summation to the method that takes a vector of links. The reference area is therefore a geometric total rather than a projected area along a specific relative velocity direction; projection against a flow direction is handled separately by the normal and tangent vector helpers in the same file.

## Design & Implementation
The body builds a `Float64[]` accumulator, pushes one entry per root, and then returns either the bare scalar or the vector depending on whether exactly one assembly was found. That conditional return keeps the common single-vehicle case ergonomic for callers that expect a number, while still supporting a multi-assembly model. Two further methods in the file cover the case where the caller already has a specific body or an explicit vector of links, and all three converge on the same summation, so the geometry convention is defined in exactly one place.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Vehicle package namespace that loads this file and exports the symbol. |
| in | `model` | SpacecraftModel | n/a | yes | Vehicle whose assemblies are traversed to accumulate facet areas. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `area` | Union{Float64,Vector{Float64}} | m^2 | — | Total reference area, scalar for a single assembly or one entry per assembly. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- *none*

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/vehicle/structure/geometry_properties.jl:10-10`
- `callees` → [[vehx.structure_assembly_graph_traverse_bodies|traverse_bodies]] · `callers` · call · `src/vehicle/structure/geometry_properties.jl:9-9`
<!-- vulcan:connections:end -->

## Limitations
The union return type is not type stable, which forces callers into a dynamic dispatch or a runtime branch and blocks inlining into a hot effector loop. The method depends on `model.roots`, which the current `SpacecraftModel` struct does not expose, so it applies to model variants that retain a root vector. Facets that are shadowed by other structure are double-counted because no occlusion test is performed, and articulated appendages contribute their full area regardless of hinge angle.

## Provenance
Mapped from `src/vehicle/structure/geometry_properties.jl:6-13`, with the per-body and per-vector methods immediately below it and the projection helpers `get_normal_vector` and `get_tangent_vector` at the end of the same file.
