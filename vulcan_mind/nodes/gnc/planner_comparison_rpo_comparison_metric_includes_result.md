---
id: gnc.planner_comparison_rpo_comparison_metric_includes_result
label: rpo_comparison_metric_includes_result
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: rpo_comparison_metric_includes_result
  lines:
  - 713
  - 713
inputs:
- id: row
  type: Any
  units: n/a
  required: true
  description: Positional argument `row`.
- id: metric
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `metric`.
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
  description: Return value of `rpo_comparison_metric_includes_result`. Returns `metric
    in (:success_pct, :keepout_violations) || row.success`.
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

# rpo_comparison_metric_includes_result

## Purpose
Inline predicate deciding whether a result row participates in a metric's aggregate: success rate and keep-out violations count every run, while all other metrics are computed over successful runs only, so failed trajectories do not distort fuel or clearance statistics.

## Design & Implementation
Defined as `@inline rpo_comparison_metric_includes_result(row, metric::Symbol) = metric in (:success_pct, :keepout_violations) || row.success`. Evaluated per row inside the summary plot loop; cheap tuple membership plus a Boolean field read.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `row` | Any | n/a | yes | Positional argument `row`. |
| in | `metric` | Symbol | n/a | yes | Positional argument `metric`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_comparison_metric_includes_result`. Returns `metric in (:success_pct, :keepout_violations) \|\| row.success`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_comparison_metric_summary_plot|rpo_comparison_metric_summary_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:762-762`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The inclusion rule is fixed in code; there is no option to show planner runtime or goal error for failed cases in the plot even though those values exist. Planners with zero successes contribute no points to seven of the nine subplots, which can make an empty category look like missing data rather than total failure.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 713.
