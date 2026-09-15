---
id: gnc.rrt_warmstart__robot_arm_rrt_nearest_index
label: _robot_arm_rrt_nearest_index
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_nearest_index
  lines:
  - 12
  - 12
inputs:
- id: tree
  type: RobotArmRRTConnectTree
  units: n/a
  required: true
  description: Positional argument `tree`.
- id: q
  type: AbstractVector{<:Real}
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
  description: Return value of `_robot_arm_rrt_nearest_index`. Returns `hypr_rrt_nearest_index(tree,
    Float64.(collect(q)))`.
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

# _robot_arm_rrt_nearest_index

## Purpose
`_robot_arm_rrt_nearest_index` returns the index of the tree node closest in joint space to a query configuration. `_robot_arm_rrt_extend!` uses it to choose which node to grow from when extending toward a random or biased sample.

## Design & Implementation
Signature `(tree::RobotArmRRTConnectTree, q::AbstractVector{<:Real})`. It converts `q` with `Float64.(collect(q))` and forwards to the shared `hypr_rrt_nearest_index(tree, q)`, which performs a linear scan over `tree.nodes` computing `sum(abs2, node - q)` and returning the index with the smallest squared distance (initialised to index 1 with `Inf`). The result is an `Int` in `1:length(tree.nodes)`. The function allocates the converted query and, inside the scan, one temporary difference vector per node.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | RobotArmRRTConnectTree | n/a | yes | Positional argument `tree`. |
| in | `q` | AbstractVector{<:Real} | n/a | yes | Positional argument `q`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_rrt_nearest_index`. Returns `hypr_rrt_nearest_index(tree, Float64.(collect(q)))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_warmstart__robot_arm_rrt_extend_bang|_robot_arm_rrt_extend!]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:68-68`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_rrt_nearest_index|hypr_rrt_nearest_index]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:13-13`
<!-- vulcan:connections:end -->

## Limitations
O(n) per query with a temporary allocation per node, which dominates runtime for large trees. The metric is plain Euclidean over joint angles, ignoring joint limits, wrap-around for revolute joints and the kinematic significance of each joint. Ties resolve to the lowest index. A query whose length differs from the node dimension raises a `DimensionMismatch` inside the subtraction.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 12.
