---
id: gnc.rrt_connect_rpo_rrt_extend_bang
label: rpo_rrt_extend!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_extend!
  lines:
  - 154
  - 154
inputs:
- id: tree
  type: RPORRTConnectTree
  units: n/a
  required: true
  description: Positional argument `tree`.
- id: q_target
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `q_target`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: settings
  type: RPORRTConnectSettings
  units: n/a
  required: true
  description: Positional argument `settings`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: true
  description: Keyword argument `safe_distance_m`.
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
  description: Return value of `rpo_rrt_extend!`; mutates `tree` in place. Returns
    `:trapped, nearest_idx` or `status, length(tree.nodes)`.
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

# rpo_rrt_extend!

## Purpose
Grows the tree by one collision-checked step toward `q_target`, appending a node with its parent link and accumulated cost.

## Design & Implementation
Finds `nearest_idx` with `rpo_rrt_nearest_index`, steers with `rpo_rrt_steer(tree.nodes[nearest_idx], q_target, settings.step_size_m)`, and returns `(:trapped, nearest_idx)` if steering is trapped or if `rpo_rrt_segment_is_safe` rejects the new edge. Otherwise it mutates `tree` via three `push!` calls: the new node, `nearest_idx` as parent, and `costs[nearest_idx] + norm(q_new - nodes[nearest_idx])`. Returns `(status, length(tree.nodes))` where status is `:advanced` or `:reached`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | RPORRTConnectTree | n/a | yes | Positional argument `tree`. |
| in | `q_target` | SVector{3, Float64} | n/a | yes | Positional argument `q_target`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `settings` | RPORRTConnectSettings | n/a | yes | Positional argument `settings`. |
| in | `safe_distance_m` | Real | n/a | yes | Keyword argument `safe_distance_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_extend!`; mutates `tree` in place. Returns `:trapped, nearest_idx` or `status, length(tree.nodes)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_connect_bang|rpo_rrt_connect!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:191-191`
- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:370-370`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:173-173`
- `callees` → [[gnc.rrt_connect_rpo_rrt_nearest_index|rpo_rrt_nearest_index]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:161-161`
- `callees` → [[gnc.rrt_connect_rpo_rrt_segment_is_safe|rpo_rrt_segment_is_safe]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:164-164`
- `callees` → [[gnc.rrt_connect_rpo_rrt_steer|rpo_rrt_steer]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:162-162`
<!-- vulcan:connections:end -->

## Limitations
A `:trapped` return still reports the nearest index, so callers must not treat that index as newly added. The three `push!` operations are not atomic, so an exception between them would corrupt the tree invariants. No duplicate-node suppression: repeated sampling near an existing node adds near-coincident nodes.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 154.
