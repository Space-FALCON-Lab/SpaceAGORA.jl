---
id: gnc.rrt_connect_rpo_rrt_connect_bang
label: rpo_rrt_connect!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_connect!
  lines:
  - 180
  - 180
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
  description: Return value of `rpo_rrt_connect!`; mutates `tree` in place. Returns
    `status, idx`.
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

# rpo_rrt_connect!

## Purpose
Greedily extends a tree toward a fixed target point step after step until it either reaches the target, becomes trapped, or exhausts `connect_max_steps`.

## Design & Implementation
Loops `while status == :advanced && steps < max(1, settings.connect_max_steps)`, calling `rpo_rrt_extend!(tree, q_target, geometry, settings; safe_distance_m)` each iteration and counting steps. Returns the final `(status, idx)` pair; `:reached` signals a successful bidirectional join, and the caller reads `idx` as the tree index that touched the target.

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
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_connect!`; mutates `tree` in place. Returns `status, idx`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:380-380`

**Downstream**

- `callees` → [[gnc.rrt_connect_rpo_rrt_extend_bang|rpo_rrt_extend!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:191-191`
<!-- vulcan:connections:end -->

## Limitations
If the step cap is hit while still `:advanced`, the function returns `:advanced` rather than a distinct exhaustion status, so callers cannot distinguish the two. Each step re-runs a full nearest-neighbour scan even though the nearest node is almost always the one just added. `connect_max_steps` default of 10,000 with a 0.75 m step allows a 7.5 km chain, far beyond typical RPO bounds.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 180.
