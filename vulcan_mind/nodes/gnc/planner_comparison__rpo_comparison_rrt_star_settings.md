---
id: gnc.planner_comparison__rpo_comparison_rrt_star_settings
label: _rpo_comparison_rrt_star_settings
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: _rpo_comparison_rrt_star_settings
  lines:
  - 214
  - 214
inputs:
- id: cfg
  type: RPOPlannerComparisonConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: n_iters
  type: Integer
  units: n/a
  required: true
  description: Positional argument `n_iters`.
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
  type: RPORRTStarSettings
  units: n/a
  description: Return value of `_rpo_comparison_rrt_star_settings`. Returns `RPORRTStarSettings(`.
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

# _rpo_comparison_rrt_star_settings

## Purpose
Copies the batch's `RPORRTStarSettings` with a replaced iteration count so RRT* runs under the same budget-matching rules as the other comparison planners.

## Design & Implementation
`_rpo_comparison_rrt_star_settings(cfg, n_iters::Integer)` mirrors the RRT-Connect variant: builds a new `RPORRTStarSettings` with `n_iters = Int(n_iters)`, copying `step_size_m`, `goal_sample_rate`, `collision_sample_ds_m`, the six adaptive collision-sampling fields, `neighbor_radius_m`, and `shortcut_iters` from `cfg.rrt_star`. Only the `:rrt_star` branch of `rpo_plan_comparison_path` calls it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RPOPlannerComparisonConfig | n/a | yes | Positional argument `cfg`. |
| in | `n_iters` | Integer | n/a | yes | Positional argument `n_iters`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPORRTStarSettings | n/a | — | Return value of `_rpo_comparison_rrt_star_settings`. Returns `RPORRTStarSettings(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:343-343`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[gnc.rrt_connect_rporrtstarsettings|RPORRTStarSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:216-216`
<!-- vulcan:connections:end -->

## Limitations
Explicit field copying means future additions to `RPORRTStarSettings` are dropped unless this function is updated. `neighbor_radius_m` is passed through unchanged even though a reduced iteration count changes the tree density that radius was tuned for.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 214.
