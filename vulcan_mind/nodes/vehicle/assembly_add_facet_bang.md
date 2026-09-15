---
id: vehicle.assembly_add_facet_bang
label: add_facet!
kind: function
source:
  file: src/vehicle/spacecraft/assembly.jl
  symbol: add_facet!
  lines:
  - 58
  - 58
inputs:
- id: link
  type: Link
  units: n/a
  required: true
  description: Positional argument `link`.
- id: facet
  type: Facet
  units: n/a
  required: true
  description: Positional argument `facet`.
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
  description: Return value of `add_facet!`; mutates `link` in place.
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

# add_facet!

## Purpose

`add_facet!` attaches solar-radiation-pressure surfaces to a `Link`. Each `Facet` describes one flat area of the vehicle body used by the SRP force model, and this is the only supported way to populate a link's `SRP_facets` collection during assembly.

## Design & Implementation

Two methods share the name. `add_facet!(link::Link, facet::Facet)` appends a single facet with `push!(link.SRP_facets, facet)`. `add_facet!(link::Link, facets::Vector{Facet})` splats the whole vector into one `push!` call, `push!(link.SRP_facets, facets...)`, so a batch of facets is added in the order given. Both mutate `link` in place and return the underlying `SRP_facets` vector; neither touches the owning `SpacecraftModel`, so no facet counter is maintained.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `link` | Link | n/a | yes | Positional argument `link`. |
| in | `facet` | Facet | n/a | yes | Positional argument `facet`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `add_facet!`; mutates `link` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.vehicle|Robotics]] · `api` → `module_api` · call · `src/vehicle/spacecraft/assembly.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/vehicle/spacecraft/assembly.jl:59-59`
<!-- vulcan:connections:end -->

## Limitations

There is no de-duplication, so calling twice with the same facet double-counts its area in the SRP sum. Facet normals, areas and optical coefficients are taken as supplied without validation, and the vector method builds no intermediate copy, so the splat cost grows with batch size. Mutation is unsynchronised and unsafe from several tasks at once.

## Provenance
Mapped from `src/vehicle/spacecraft/assembly.jl` line 58.
