---
id: gnc.hypr_utils_hypr_rrt_refresh_subtree_costs_bang
label: hypr_rrt_refresh_subtree_costs!
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_rrt_refresh_subtree_costs!
  lines:
  - 196
  - 196
inputs:
- id: tree
  type: Any
  units: n/a
  required: true
  description: Positional argument `tree`.
- id: parent_idx
  type: Integer
  units: n/a
  required: true
  description: Positional argument `parent_idx`.
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
  description: Return value of `hypr_rrt_refresh_subtree_costs!`; mutates `tree` in
    place. Returns `tree`.
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

# hypr_rrt_refresh_subtree_costs!

## Purpose
After a rewiring changes a node's parent, recomputes the accumulated path cost of every descendant so the tree's cost vector stays consistent.

## Design & Implementation
Breadth-first from `parent_idx` using a queue: for each popped parent it scans the whole parent vector for children, sets each child's cost to the parent's cost plus the edge length, and enqueues the child. Reads and writes through the layout-tolerant accessors, mutating the tree's cost vector in place, and returns the tree.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | Any | n/a | yes | Positional argument `tree`. |
| in | `parent_idx` | Integer | n/a | yes | Positional argument `parent_idx`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `hypr_rrt_refresh_subtree_costs!`; mutates `tree` in place. Returns `tree`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_refresh_subtree_costs_bang|rpo_rrt_refresh_subtree_costs!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:210-210`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- `callees` → [[gnc.hypr_utils__hypr_rrt_costs|_hypr_rrt_costs]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:199-199`
- `callees` → [[gnc.hypr_utils__hypr_rrt_parents|_hypr_rrt_parents]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:198-198`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:205-205`
<!-- vulcan:connections:end -->

## Limitations
Each queue pop scans every node to find children, so a refresh costs O(n × subtree size); a child-list index would make it linear in the subtree. It assumes edge cost is Euclidean length, so a tree using a different edge metric would be recosted incorrectly.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 196.
