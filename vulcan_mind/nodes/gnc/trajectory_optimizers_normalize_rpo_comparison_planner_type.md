---
id: gnc.trajectory_optimizers_normalize_rpo_comparison_planner_type
label: normalize_rpo_comparison_planner_type
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl
  symbol: normalize_rpo_comparison_planner_type
  lines:
  - 35
  - 35
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
  description: Return value of `normalize_rpo_comparison_planner_type`.
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

# normalize_rpo_comparison_planner_type

## Purpose
Canonicalises a user-supplied planner identifier (Symbol or string, any case, hyphens or underscores) into one of seven canonical symbols so the comparison harness can dispatch on a fixed vocabulary.

## Design & Implementation
Takes `planner_type`, lowercases it, replaces `-` with `_`, and converts to a `Symbol`. Alias groups map to canonical names: `:hypr`/`:pso` -> `:hypr`; `:pso_unrefined`, `:unrefined_pso`, `:hypr_unrefined`, `:unrefined_hypr` -> `:pso_unrefined`; `:rrt_connect`/`:rrtconnect` -> `:rrt_connect`; five Bezier spellings -> `:rrt_connect_bezier`; `:rrt_star`/`:rrtstar` -> `:rrt_star`; `:chomp` and `:stomp` pass through. Anything else throws `ArgumentError` listing the accepted names.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planner_type` | Any | n/a | yes | Positional argument `planner_type`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `normalize_rpo_comparison_planner_type`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_comparison_cost_iteration_plot|rpo_comparison_cost_iteration_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:1073-1073`
- [[gnc.planner_comparison_rpo_comparison_failed_paths_plot|rpo_comparison_failed_paths_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:991-991`
- [[gnc.planner_comparison_rpo_comparison_path_family_plot|rpo_comparison_path_family_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:866-866`
- [[gnc.planner_comparison_rpo_comparison_single_path_plot|rpo_comparison_single_path_plot]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:925-925`
- [[gnc.planner_comparison_rpo_plan_comparison_path|rpo_plan_comparison_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:268-268`
- [[gnc.trajectory_optimizers_rpo_comparison_planner_label|rpo_comparison_planner_label]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl:50-50`
- [[gncy.planner_comparison_rpo_run_planner_comparison_batch|rpo_run_planner_comparison_batch]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:547-547`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only hyphens are normalised; spaces or dots in the input are not stripped and will fall through to the `ArgumentError`. The function is pure and allocates a new string and symbol on every call, which matters only if invoked in a hot loop. Adding a planner requires editing both this function and `rpo_comparison_planner_label`.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/trajectory_optimizers.jl` line 35.
