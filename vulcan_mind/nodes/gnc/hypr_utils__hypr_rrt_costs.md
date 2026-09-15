---
id: gnc.hypr_utils__hypr_rrt_costs
label: _hypr_rrt_costs
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: _hypr_rrt_costs
  lines:
  - 14
  - 14
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
  description: 'Return value of `_hypr_rrt_costs`. Returns `hasproperty(tree, :costs)
    ? getproperty(tree, :costs) : getproperty(tree, :cost)`.'
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

# _hypr_rrt_costs

## Purpose
Reads the accumulated-cost vector from an RRT tree under either of the two field names in use.

## Design & Implementation
Mirrors the parent accessor: returns `tree.costs` if present, else `tree.cost`. Used only by `hypr_rrt_refresh_subtree_costs!`, which must write into whichever vector the tree actually owns.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tree` | Any | n/a | yes | Positional argument `tree`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_hypr_rrt_costs`. Returns `hasproperty(tree, :costs) ? getproperty(tree, :costs) : getproperty(tree, :cost)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.hypr_utils_hypr_rrt_refresh_subtree_costs_bang|hypr_rrt_refresh_subtree_costs!]] · `callees` → `callers` · call · `src/gnc/hypr/hypr_utils.jl:199-199`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Returns the live vector rather than a copy, which is required for the in-place refresh but means a careless caller can corrupt the tree's cost bookkeeping.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 14.
