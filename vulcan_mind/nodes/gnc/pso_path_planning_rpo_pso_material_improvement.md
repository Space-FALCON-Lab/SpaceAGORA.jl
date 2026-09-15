---
id: gnc.pso_path_planning_rpo_pso_material_improvement
label: rpo_pso_material_improvement
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_material_improvement
  lines:
  - 162
  - 162
inputs:
- id: new_cost
  type: Real
  units: n/a
  required: true
  description: Positional argument `new_cost`.
- id: reference_cost
  type: Real
  units: n/a
  required: true
  description: Positional argument `reference_cost`.
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
  description: Return value of `rpo_pso_material_improvement`.
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

# rpo_pso_material_improvement

## Purpose
Tests whether a new best cost improves on a reference by enough to reset the early-stopping patience counter.

## Design & Implementation
Forwards `new_cost`, `reference_cost` and the config's `early_stopping_min_abs_improvement` and `early_stopping_min_rel_improvement` to `hypr_material_improvement`, which requires both thresholds to be exceeded. Returns `Bool`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `new_cost` | Real | n/a | yes | Positional argument `new_cost`. |
| in | `reference_cost` | Real | n/a | yes | Positional argument `reference_cost`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `rpo_pso_material_improvement`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:569-569`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:569-569`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_material_improvement|hypr_material_improvement]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:163-163`
<!-- vulcan:connections:end -->

## Limitations
With an infinite reference, as at the start of a run, any finite cost counts as material, so the first evaluated iteration always resets patience.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 162.
