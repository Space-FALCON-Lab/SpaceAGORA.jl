---
id: gnc.hypr_utils_hypr_rrt_nearest_index
label: hypr_rrt_nearest_index
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_rrt_nearest_index
  lines:
  - 142
  - 142
inputs:
- id: tree
  type: Any
  units: n/a
  required: true
  description: Positional argument `tree`.
- id: q
  type: Any
  units: n/a
  required: true
  description: Positional argument `q`.
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
  description: Return value of `hypr_rrt_nearest_index`. Returns `best_idx`.
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

# hypr_rrt_nearest_index

## Purpose
Finds the tree node closest to a query state, the core query of every RRT extension step.

## Design & Implementation
A linear scan over `tree.nodes` comparing squared Euclidean distances, tracking the best index and value, with `@inbounds`. Squared distance avoids a square root per node. Returns the index, defaulting to one for an empty comparison.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | Any | n/a | yes | Positional argument `tree`. |
| in | `q` | Any | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `hypr_rrt_nearest_index`. Returns `best_idx`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_nearest_index|rpo_rrt_nearest_index]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:45-45`
- [[gnc.rrt_warmstart__robot_arm_rrt_nearest_index|_robot_arm_rrt_nearest_index]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:13-13`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
O(n) per query with no k-d tree or grid, so the total cost of building a tree of `n` nodes is O(n²); for the few-thousand-node trees the planners use this is acceptable but it does not scale.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 142.
