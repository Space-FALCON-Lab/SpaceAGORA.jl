---
id: gnc.pso_path_planning_rpo_pso_project_to_segment
label: rpo_pso_project_to_segment
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_project_to_segment
  lines:
  - 21
  - 21
inputs:
- id: q
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
- id: a
  type: Any
  units: n/a
  required: true
  description: Positional argument `a`.
- id: b
  type: Any
  units: n/a
  required: true
  description: Positional argument `b`.
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
  description: Return value of `rpo_pso_project_to_segment`. Returns `av + α * ab`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gnc
origin: agent
---

# rpo_pso_project_to_segment

## Purpose
Projects a point onto a segment, used by culling to pull a reseeded waypoint toward the local chord of the current best path.

## Design & Implementation
Converts all three inputs to `SVector`, returns `a` if the segment is degenerate, otherwise computes the clamped projection parameter from the dot product and returns the point on the segment. This duplicates `rpo_refinement_project_to_segment` line for line; the two exist so the swarm and the refinement passes can evolve independently.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `a` | Any | n/a | yes | Positional argument `a`. |
| in | `b` | Any | n/a | yes | Positional argument `b`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_pso_project_to_segment`. Returns `av + α * ab`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_rpo_pso_cull_swarm_bang|rpo_pso_cull_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:84-84`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Being a duplicate, a fix applied to one copy does not reach the other; the clamp keeps results on the segment but the caller still has to test the result against the keep-out region.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 21.
