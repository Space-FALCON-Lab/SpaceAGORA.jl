---
id: gnc.planner_comparison__rpo_comparison_rrt_connect_seed_path
label: _rpo_comparison_rrt_connect_seed_path
kind: function
source:
  file: src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl
  symbol: _rpo_comparison_rrt_connect_seed_path
  lines:
  - 233
  - 233
inputs:
- id: start_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `start_rtn`.
- id: goal_rtn
  type: Any
  units: n/a
  required: true
  description: Positional argument `goal_rtn`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: base_cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `base_cfg`.
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
- id: runtime_limit_s
  type: Real
  units: n/a
  required: true
  description: Keyword argument `runtime_limit_s`.
- id: rng
  type: Any
  units: n/a
  required: true
  description: Keyword argument `rng`.
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
  description: Return value of `_rpo_comparison_rrt_connect_seed_path`. Returns `plan.path`.
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

# _rpo_comparison_rrt_connect_seed_path

## Purpose
Produces the collision-free RRT-Connect seed path that warm-starts CHOMP and STOMP in the comparison, so those local optimizers begin from a feasible route rather than a straight line through the station.

## Design & Implementation
Signature `_rpo_comparison_rrt_connect_seed_path(start_rtn, goal_rtn, geometry, base_cfg::RPOPSOConfig, cfg::RPOPlannerComparisonConfig, n_iters; runtime_limit_s, rng)`. Builds settings via `_rpo_comparison_rrt_connect_settings(cfg, n_iters)`, calls `rpo_rrt_connect_plan_path(start_rtn, goal_rtn, geometry, base_cfg; safe_distance_m = cfg.safe_distance_m, settings, max_runtime_s = runtime_limit_s, rng)`, and returns only `plan.path`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `start_rtn` | Any | n/a | yes | Positional argument `start_rtn`. |
| in | `goal_rtn` | Any | n/a | yes | Positional argument `goal_rtn`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `base_cfg` | RPOPSOConfig | n/a | yes | Positional argument `base_cfg`. |
| in | `cfg` | RPOPlannerComparisonConfig | n/a | yes | Positional argument `cfg`. |
| in | `n_iters` | Integer | n/a | yes | Positional argument `n_iters`. |
| in | `runtime_limit_s` | Real | n/a | yes | Keyword argument `runtime_limit_s`. |
| in | `rng` | Any | n/a | yes | Keyword argument `rng`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_rpo_comparison_rrt_connect_seed_path`. Returns `plan.path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.planner_comparison_runtime_limited_iters|runtime_limited_iters]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:370-370`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl`

**Downstream**

- `callees` → [[gnc.planner_comparison__rpo_comparison_rrt_connect_settings|_rpo_comparison_rrt_connect_settings]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:243-243`
- `callees` → [[gncy.rrt_connect_rpo_rrt_connect_plan_path|rpo_rrt_connect_plan_path]] · `callers` · call · `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl:244-244`
<!-- vulcan:connections:end -->

## Limitations
The seed planning time is included in the CHOMP/STOMP `planner_compute_time` because the caller's `t0` is taken before this call, inflating those planners' reported runtime relative to HYPR. It also consumes `rng` state, so the optimizer's subsequent random draws differ from a run without seeding. If RRT-Connect fails, the returned path may be partial or straight and is passed on without a feasibility check.

## Provenance
Mapped from `src/gnc/guidance/rpo/comparison_methods/planner_comparison.jl` line 233.
