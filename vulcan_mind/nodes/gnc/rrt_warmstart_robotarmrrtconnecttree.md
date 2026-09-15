---
id: gnc.rrt_warmstart_robotarmrrtconnecttree
label: RobotArmRRTConnectTree
kind: struct
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: RobotArmRRTConnectTree
  lines:
  - 2
  - 2
inputs:
- id: nodes
  type: Vector{Vector{Float64}}
  units: n/a
  required: true
  description: Field `nodes`.
- id: parents
  type: Vector{Int}
  units: n/a
  required: true
  description: Field `parents`.
- id: costs
  type: Vector{Float64}
  units: n/a
  required: true
  description: Field `costs`.
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
  type: RobotArmRRTConnectTree
  units: n/a
  description: Constructed `RobotArmRRTConnectTree`.
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

# RobotArmRRTConnectTree

## Purpose
`RobotArmRRTConnectTree` is the growable tree container used by the robot-arm RRT-Connect warm start: it stores joint-space configurations, the parent index of each node, and the accumulated path cost from the root. Two instances (rooted at the start and goal configurations) are grown alternately by `_robot_arm_rrt_connect_warmstart_path`.

## Design & Implementation
An immutable struct with three parallel vectors: `nodes::Vector{Vector{Float64}}` (joint angles in radians), `parents::Vector{Int}` (`0` for the root) and `costs::Vector{Float64}` (cumulative Euclidean joint-space length). The convenience constructor `RobotArmRRTConnectTree(root)` converts `root` with `Float64.(collect(root))` and initialises the vectors to `[root]`, `[0]`, `[0.0]`. The struct itself is immutable but its vectors are appended to in place by `_robot_arm_rrt_extend!`, so a tree is effectively mutable. Generic HYPR helpers `hypr_rrt_nearest_index`, `hypr_rrt_tree_path` and `hypr_rrt_join_paths` operate on any object exposing these three fields.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `nodes` | Vector{Vector{Float64}} | n/a | yes | Field `nodes`. |
| in | `parents` | Vector{Int} | n/a | yes | Field `parents`. |
| in | `costs` | Vector{Float64} | n/a | yes | Field `costs`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RobotArmRRTConnectTree | n/a | — | Constructed `RobotArmRRTConnectTree`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path|_robot_arm_rrt_connect_warmstart_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:244-244`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Nearest-neighbour lookup is a linear scan over `nodes`, so tree growth is O(n^2) overall; there is no spatial index. The three vectors must be kept the same length, but nothing enforces that invariant beyond the discipline of the extend function. Node vectors are stored by reference, so external mutation of a pushed configuration corrupts the tree. Joint-space distance is unweighted Euclidean, treating every joint equally regardless of link length.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 2.
