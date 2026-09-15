---
id: gnc.pso_refinement_rpo_refinement_project_to_segment
label: rpo_refinement_project_to_segment
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refinement_project_to_segment
  lines:
  - 152
  - 152
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
  description: Return value of `rpo_refinement_project_to_segment`. Returns `av +
    α * ab`.
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

# rpo_refinement_project_to_segment

## Purpose
Finds the closest point on a segment to a given point, used to generate candidate positions when pulling Bezier handles toward local chords.

## Design & Implementation
Converts all three inputs to `SVector`, forms the segment direction `ab`, and returns `a` outright when the segment is degenerate. Otherwise the projection parameter is the dot product of the offset with `ab` over the squared length, clamped to the unit interval so the result stays on the segment rather than its infinite extension.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `a` | Any | n/a | yes | Positional argument `a`. |
| in | `b` | Any | n/a | yes | Positional argument `b`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_refinement_project_to_segment`. Returns `av + α * ab`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_refinement_rpo_refine_tighten_handles|rpo_refine_tighten_handles]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:202-202`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
None of the geometric significance is validated; if the chord itself passes through the keep-out region the projected target is inside the station and the subsequent cost test is what rejects it.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 152.
