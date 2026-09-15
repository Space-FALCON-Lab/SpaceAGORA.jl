---
id: gnc.planner_comparison__rpo_comparison_rrt_connect_settings
label: _rpo_comparison_rrt_connect_settings
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: _rpo_comparison_rrt_connect_settings
  lines:
  - 195
  - 195
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
  type: RPORRTConnectSettings
  units: n/a
  description: Return value of `_rpo_comparison_rrt_connect_settings`. Returns `RPORRTConnectSettings(`.
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

# _rpo_comparison_rrt_connect_settings

## Purpose
Copies the batch's `RPORRTConnectSettings` while substituting a caller-selected iteration cap, so RRT-Connect can be run under runtime- or HYPR-matched budgets without mutating the shared config.

## Design & Implementation
`_rpo_comparison_rrt_connect_settings(cfg::RPOPlannerComparisonConfig, n_iters::Integer)` reads `src = cfg.rrt_connect` and constructs a new `RPORRTConnectSettings` with `n_iters = Int(n_iters)` and the eleven remaining fields (`step_size_m`, `goal_sample_rate`, `collision_sample_ds_m`, the six adaptive collision-sampling fields, `connect_max_steps`, `shortcut_iters`) copied verbatim. Used by the `:rrt_connect` and `:rrt_connect_bezier` branches and by `_rpo_comparison_rrt_connect_seed_path`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RPOPlannerComparisonConfig | n/a | yes | Positional argument `cfg`. |
| in | `n_iters` | Integer | n/a | yes | Positional argument `n_iters`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPORRTConnectSettings | n/a | — | Return value of `_rpo_comparison_rrt_connect_settings`. Returns `RPORRTConnectSettings(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison__rpo_comparison_rrt_connect_seed_path|_rpo_comparison_rrt_connect_seed_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:243-243`
- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:301-301`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[gnc.rrt_connect_rporrtconnectsettings|RPORRTConnectSettings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:197-197`
<!-- vulcan:connections:end -->

## Limitations
The field list is hand-maintained; a new field added to `RPORRTConnectSettings` would silently revert to its default here. `Int(n_iters)` throws `InexactError` for non-integral inputs. No check that `n_iters > 0`.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 195.
