---
id: gncz.rpo_reference_geometry_rporeferencegeometry
label: RPOReferenceGeometry
kind: struct
source:
  file: src/gnc/navigation/rpo_nav/reference_geometry/rpo_reference_geometry.jl
  symbol: RPOReferenceGeometry
  lines:
  - 2
  - 5
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: NavigationHooks namespace where the station and chaser geometry types
    are already declared by the earlier includes.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: geometry
  type: RPOReferenceGeometry
  units: n/a
  description: Paired station and chaser geometry passed as one argument to every
    clearance, normal, and path statistics query.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncz
origin: agent
---

# RPOReferenceGeometry

## Purpose
`RPOReferenceGeometry` binds the target station geometry and the chaser geometry into the single value that the proximity-operations distance layer takes as its geometric context. Every clearance and surface query accepts this pair rather than two arguments, which guarantees that a station and the chaser inflating its keepout are never mismatched at a call site.

## Model & Assumptions
The struct is immutable and holds the station geometry, which carries the point cloud, its KD-tree, the keepout radius, and a name, together with the chaser box half extents. The convenience constructor takes a station positionally and accepts the chaser as a keyword defaulting to a freshly built default CubeSat, so a scenario that only cares about the target can omit the chaser entirely and still get a usable geometry.

## Design & Implementation
Because the type is declared after the station and chaser files are included, its fields are concretely typed rather than abstract, so field access inside the distance loops compiles to a direct load. Holding both bodies in one immutable value also makes the geometry cheap to capture inside a planner closure that scores many candidate paths against the same target.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | NavigationHooks namespace where the station and chaser geometry types are already declared by the earlier includes. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `geometry` | RPOReferenceGeometry | n/a | — | Paired station and chaser geometry passed as one argument to every clearance, normal, and path statistics query. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.replanning_rpo_geometry_with_replanning_spheres|rpo_geometry_with_replanning_spheres]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr_planning/replanning.jl:168-168`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/navigation/rpo_nav/reference_geometry/rpo_reference_geometry.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The pairing is static, so a scenario with several chasers or several targets needs one geometry value per pair and cannot express a shared station with per-vehicle inflation in a single object. No relative pose is stored, which means the whole distance layer implicitly assumes that query points have already been expressed in the station body frame.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/reference_geometry/rpo_reference_geometry.jl:1-11`.
