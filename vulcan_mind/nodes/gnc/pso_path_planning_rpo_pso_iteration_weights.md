---
id: gnc.pso_path_planning_rpo_pso_iteration_weights
label: rpo_pso_iteration_weights
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_iteration_weights
  lines:
  - 2
  - 2
inputs:
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: iter
  type: Int
  units: n/a
  required: true
  description: Positional argument `iter`.
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
  description: Return value of `rpo_pso_iteration_weights`. Returns `hypr_iteration_weights(`.
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

# rpo_pso_iteration_weights

## Purpose
Supplies the inertia and the two acceleration coefficients for a given PSO iteration, scheduled to shift from exploration toward exploitation as the run progresses.

## Design & Implementation
Forwards thirteen fields of the `RPOPSOConfig` — the schedule enable flag, iteration count, the base `w_inertia`, `c1` and `c2`, and the eight schedule shape parameters covering transition fraction, minimum inertia, end fractions and coefficient bounds — to the shared `hypr_iteration_weights` together with `iter`. Keeping the schedule logic in the HYPR core means the robot-arm and RPO planners taper identically.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `iter` | Int | n/a | yes | Positional argument `iter`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_pso_iteration_weights`. Returns `hypr_iteration_weights(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:588-588`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:588-588`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_iteration_weights|hypr_iteration_weights]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:3-3`
<!-- vulcan:connections:end -->

## Limitations
The schedule is a function of iteration index only, not of observed convergence, so a run that stalls early still keeps exploring until the scheduled transition.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 2.
