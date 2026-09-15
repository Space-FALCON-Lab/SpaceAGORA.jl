---
id: gncz.cubesat_geometry_rpocubesatgeometry
label: RPOCubeSatGeometry
kind: struct
source:
  file: src/gnc/navigation/rpo_nav/reference_geometry/cubesat_geometry.jl
  symbol: RPOCubeSatGeometry
  lines:
  - 2
  - 4
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: NavigationHooks namespace where the chaser geometry type is declared
    and exported to the clearance routines.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: chaser_geometry
  type: RPOCubeSatGeometry
  units: m
  description: Immutable chaser body half extents used to inflate the station keepout
    distance during clearance checks.
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

# RPOCubeSatGeometry

## Purpose
`RPOCubeSatGeometry` is the chaser side of the proximity-operations collision model. It reduces the approaching spacecraft to a box described by its half extents in the body frame, which is the only chaser property the clearance computation consumes.

## Model & Assumptions
The struct stores a single static three-vector of half extents rather than full dimensions, because clearance needs a radius rather than a span. The keyword constructor takes full dimensions in metres, defaulting to a three-unit CubeSat form factor of one tenth by one tenth by three tenths of a metre, validates that every dimension is finite and strictly positive, and halves them before storing. Rejecting a zero dimension at construction prevents a degenerate body that would report optimistic clearance everywhere.

## Design & Implementation
The type is immutable and holds a statically sized vector, so instances are stack allocated and can be embedded in the combined reference geometry without indirection. Clearance uses the largest of the three half extents, treating the chaser as a circumscribing sphere, which means the stored triple is more information than the current consumer needs and leaves room for an orientation-aware check later.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | NavigationHooks namespace where the chaser geometry type is declared and exported to the clearance routines. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `chaser_geometry` | RPOCubeSatGeometry | m | — | Immutable chaser body half extents used to inflate the station keepout distance during clearance checks. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/navigation/rpo_nav/reference_geometry/cubesat_geometry.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
A box abstraction cannot represent deployed appendages such as solar arrays or booms, and the sphere reduction used downstream makes the model conservative in proportion to the aspect ratio of the box. No attitude is stored, so the geometry cannot tighten clearance when the long axis is known to point away from the station.

## Provenance
Mapped from `src/gnc/navigation/rpo_nav/reference_geometry/cubesat_geometry.jl:1-13`.
