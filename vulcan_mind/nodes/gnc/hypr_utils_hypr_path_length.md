---
id: gnc.hypr_utils_hypr_path_length
label: hypr_path_length
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_path_length
  lines:
  - 17
  - 17
inputs:
- id: points
  type: Any
  units: n/a
  required: true
  description: Positional argument `points`.
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
  description: Return value of `hypr_path_length`. Returns `total`.
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

# hypr_path_length

## Purpose
Computes the polyline arc length of a waypoint matrix, the base measure every planner's length cost and sample-count heuristic uses.

## Design & Implementation
Converts to `Matrix{Float64}`, returns zero for fewer than two columns, and accumulates `norm(pts[:, j+1] - pts[:, j])` over consecutive columns under `@inbounds`. Working on columns rather than rows matches the three-by-N convention used throughout the planners.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `points` | Any | n/a | yes | Positional argument `points`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `hypr_path_length`. Returns `total`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_path_length|rpo_path_length]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:3-3`
- [[gnc.swarm_and_retiming__robot_arm_path_length|_robot_arm_path_length]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:118-118`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Each column slice allocates a small vector, so the loop allocates twice per segment; for the short polygons planners use this is negligible, but it is not suited to dense sampled paths in a hot loop.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 17.
