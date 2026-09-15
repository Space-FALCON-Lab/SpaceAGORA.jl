---
id: gnc.pso_path_planning_rpo_pso_early_stopping_feasible
label: rpo_pso_early_stopping_feasible
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_early_stopping_feasible
  lines:
  - 172
  - 172
inputs:
- id: components
  type: Any
  units: n/a
  required: true
  description: Positional argument `components`.
- id: cfg
  type: RPOPSOConfig
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
  type: Bool
  units: n/a
  description: Return value of `rpo_pso_early_stopping_feasible`.
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

# rpo_pso_early_stopping_feasible

## Purpose
Gates early stopping on feasibility so the planner never quits early while its best path still violates clearance.

## Design & Implementation
Returns true immediately if `early_stopping_require_feasible` is off. Otherwise returns false for a `nothing` components record and true only when its `violation_count` is zero.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `components` | Any | n/a | yes | Positional argument `components`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `rpo_pso_early_stopping_feasible`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:575-575`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:575-575`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Feasibility here means zero sampled violations, which inherits the sampling resolution's blind spots between samples.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 172.
