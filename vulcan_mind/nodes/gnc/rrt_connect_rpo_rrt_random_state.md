---
id: gnc.rrt_connect_rpo_rrt_random_state
label: rpo_rrt_random_state
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_random_state
  lines:
  - 54
  - 54
inputs:
- id: rng
  type: Any
  units: n/a
  required: true
  description: Positional argument `rng`.
- id: lo
  type: Any
  units: n/a
  required: true
  description: Positional argument `lo`.
- id: hi
  type: Any
  units: n/a
  required: true
  description: Positional argument `hi`.
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
  type: SVector
  units: n/a
  description: Return value of `rpo_rrt_random_state`. Returns `SVector{3, Float64}(`.
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

# rpo_rrt_random_state

## Purpose
Draws a uniformly random RTN position inside an axis-aligned bounding box for tree exploration.

## Design & Implementation
Takes an RNG and index-able bounds `lo`, `hi` (three components each, metres). Returns `SVector{3,Float64}(lo[k] + rand(rng)*(hi[k]-lo[k]))` for k = 1..3 using three independent `rand(rng)` draws in [0,1). The bounds come from `rpo_pso_bounds(start, goal, cfg)` in the planner entry points.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rng` | Any | n/a | yes | Positional argument `rng`. |
| in | `lo` | Any | n/a | yes | Positional argument `lo`. |
| in | `hi` | Any | n/a | yes | Positional argument `hi`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SVector | n/a | — | Return value of `rpo_rrt_random_state`. Returns `SVector{3, Float64}(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:546-546`
- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:366-366`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
No check that `hi[k] >= lo[k]`; inverted bounds still produce values but outside the intended box. Sampling is uniform in the box with no obstacle-aware or informed biasing, so narrow passages are found slowly. Only the first three entries of `lo`/`hi` are used.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 54.
