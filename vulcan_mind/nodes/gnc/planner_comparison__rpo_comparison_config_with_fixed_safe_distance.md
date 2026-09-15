---
id: gnc.planner_comparison__rpo_comparison_config_with_fixed_safe_distance
label: _rpo_comparison_config_with_fixed_safe_distance
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: _rpo_comparison_config_with_fixed_safe_distance
  lines:
  - 48
  - 48
inputs:
- id: cfg
  type: RPOPlannerComparisonConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  type: RPOPlannerComparisonConfig
  units: n/a
  description: Return value of `_rpo_comparison_config_with_fixed_safe_distance`.
    Returns `RPOPlannerComparisonConfig(`.
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

# _rpo_comparison_config_with_fixed_safe_distance

## Purpose
Returns a copy of an `RPOPlannerComparisonConfig` with the keep-out distance pinned to the published 0.5 m value in both the top-level `safe_distance_m` field and the embedded `pso_config`, guaranteeing every planner is scored against the same collision margin.

## Design & Implementation
Constructs a new `RPOPlannerComparisonConfig` by copying every field of `cfg` except two: `pso_config = rpo_pso_config(cfg.pso_config; safe_distance_m = RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M)` (which also re-syncs `sample_ds_m` to 0.5 m and re-validates) and `safe_distance_m = RPO_PLANNER_COMPARISON_SAFE_DISTANCE_M`. Called at the top of `rpo_plan_comparison_path` and `rpo_run_planner_comparison_batch`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RPOPlannerComparisonConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | RPOPlannerComparisonConfig | n/a | — | Return value of `_rpo_comparison_config_with_fixed_safe_distance`. Returns `RPOPlannerComparisonConfig(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_rpo_plan_comparison_path|rpo_plan_comparison_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:267-267`
- [[gncy.planner_comparison_rpo_run_planner_comparison_batch|rpo_run_planner_comparison_batch]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:545-545`

**Downstream**

- `callees` → [[gnc.planner_comparison_rpoplannercomparisonconfig|RPOPlannerComparisonConfig]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:49-49`
- `callees` → [[gnc.pso_parameters_rpo_pso_config|rpo_pso_config]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:51-51`
<!-- vulcan:connections:end -->

## Limitations
Silently overrides any caller-chosen keep-out distance, so the config field is effectively read-only for these entry points. Rebuilding `pso_config` through `rpo_pso_config` can throw `ArgumentError` if the base config is otherwise invalid. Allocates a full config copy on every call, including inside the per-path runner.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 48.
