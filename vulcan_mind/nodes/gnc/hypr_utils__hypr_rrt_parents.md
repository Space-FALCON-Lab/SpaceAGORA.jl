---
id: gnc.hypr_utils__hypr_rrt_parents
label: _hypr_rrt_parents
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: _hypr_rrt_parents
  lines:
  - 12
  - 12
inputs:
- id: tree
  type: Any
  units: n/a
  required: true
  description: Positional argument `tree`.
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
  description: 'Return value of `_hypr_rrt_parents`. Returns `hasproperty(tree, :parents)
    ? getproperty(tree, :parents) : getproperty(tree, :p`.'
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

# _hypr_rrt_parents

## Purpose
Reads the parent-index vector from an RRT tree regardless of whether it was built by the RPO planner or the robot-arm planner.

## Design & Implementation
A one-line accessor that returns `tree.parents` when that property exists and `tree.parent` otherwise, via `hasproperty` and `getproperty`. This is what lets the path-reconstruction and cost-refresh routines accept both tree layouts.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | Any | n/a | yes | Positional argument `tree`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_hypr_rrt_parents`. Returns `hasproperty(tree, :parents) ? getproperty(tree, :parents) : getproperty(tree, :p`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.hypr_utils_hypr_rrt_refresh_subtree_costs_bang|hypr_rrt_refresh_subtree_costs!]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:198-198`
- [[gnc.hypr_utils_hypr_rrt_tree_path|hypr_rrt_tree_path]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:179-179`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`hasproperty` on a struct is resolved at runtime, and the two names hide a genuine type split; a third tree type with yet another field name would fail with a `FieldError` here rather than a clear message.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 12.
