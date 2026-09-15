---
id: gnc.planner_comparison_rpo_comparison_cost_iteration_xvalues
label: rpo_comparison_cost_iteration_xvalues
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_cost_iteration_xvalues
  lines:
  - 1060
  - 1060
inputs:
- id: plan
  type: Any
  units: n/a
  required: true
  description: Positional argument `plan`.
- id: costs
  type: Any
  units: n/a
  required: true
  description: Positional argument `costs`.
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
  type: AbstractArray
  units: n/a
  description: Return value of `rpo_comparison_cost_iteration_xvalues`. Returns `[max(iterations,
    0)]` or `collect(1:n)`.
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

# rpo_comparison_cost_iteration_xvalues

## Purpose
Chooses x-axis values for a planner's cost-history trace: sequential iteration numbers when a history exists, or the planner's reported iteration count when only a single cost sample is available.

## Design & Implementation
`rpo_comparison_cost_iteration_xvalues(plan, costs)` returns `Int[]` for empty `costs`. It determines `iterations` from `plan.planner_iteration_count`, falling back to `plan.iterations`, then to `length(costs)`. For a single cost it returns `[max(iterations, 0)]` so the marker sits at the final iteration; otherwise `collect(1:n)`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `plan` | Any | n/a | yes | Positional argument `plan`. |
| in | `costs` | Any | n/a | yes | Positional argument `costs`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractArray | n/a | — | Return value of `rpo_comparison_cost_iteration_xvalues`. Returns `[max(iterations, 0)]` or `collect(1:n)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_comparison_cost_iteration_plot|rpo_comparison_cost_iteration_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1083-1083`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
For multi-sample histories the x values are `1:n` even if the planner's reported iteration count differs from the history length (for example when `costs` was filtered for finite values by the caller), so the axis can be offset from true iteration numbers. `Int(plan.planner_iteration_count)` throws for non-integral values.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 1060.
