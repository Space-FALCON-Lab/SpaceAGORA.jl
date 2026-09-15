---
id: gnc.trajectory_optimizers_rpo_comparison_planner_label
label: rpo_comparison_planner_label
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: rpo_comparison_planner_label
  lines:
  - 49
  - 49
inputs:
- id: planner_type
  type: Any
  units: n/a
  required: true
  description: Positional argument `planner_type`.
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
  description: Return value of `rpo_comparison_planner_label`. Returns `uppercase(string(planner))`.
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

# rpo_comparison_planner_label

## Purpose
Returns the human-readable display string for a comparison planner, used in plot legends, CSV headers, and progress lines of the RPO planner comparison outputs.

## Design & Implementation
Calls `normalize_rpo_comparison_planner_type(planner_type)` first, so any alias accepted there is labelled consistently. The mapping is `:hypr` -> "HYPR", `:pso_unrefined` -> "PSO (unrefined)", `:rrt_connect` -> "RRT-Connect", `:rrt_connect_bezier` -> "RRT-Connect + Bezier", `:rrt_star` -> "RRT*", `:chomp` -> "CHOMP", `:stomp` -> "STOMP". The final `uppercase(string(planner))` fallback is unreachable in practice because normalisation throws on unknown names.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planner_type` | Any | n/a | yes | Positional argument `planner_type`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_comparison_planner_label`. Returns `uppercase(string(planner))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison__rpo_comparison_progress_line_bang|_rpo_comparison_progress_line!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:83-83`
- [[gnc.planner_comparison_rpo_comparison_cost_iteration_plot|rpo_comparison_cost_iteration_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1077-1077`
- [[gnc.planner_comparison_rpo_comparison_failed_paths_plot|rpo_comparison_failed_paths_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:992-992`
- [[gnc.planner_comparison_rpo_comparison_metric_summary_plot|rpo_comparison_metric_summary_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:760-760`
- [[gnc.planner_comparison_rpo_comparison_path_family_plot|rpo_comparison_path_family_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:868-868`
- [[gnc.planner_comparison_rpo_comparison_single_path_plot|rpo_comparison_single_path_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:926-926`
- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:281-281`
- [[gncy.planner_comparison_rpo_run_planner_comparison_batch|rpo_run_planner_comparison_batch]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:594-594`

**Downstream**

- `callees` → [[gnc.trajectory_optimizers_normalize_rpo_comparison_planner_type|normalize_rpo_comparison_planner_type]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:50-50`
<!-- vulcan:connections:end -->

## Limitations
Because normalisation throws `ArgumentError` for unknown identifiers, this function is not safe to call on untrusted input without a try/catch. Labels are hard-coded English strings; there is no localisation or override hook.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 49.
