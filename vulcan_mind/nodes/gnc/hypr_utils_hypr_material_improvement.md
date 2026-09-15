---
id: gnc.hypr_utils_hypr_material_improvement
label: hypr_material_improvement
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_material_improvement
  lines:
  - 115
  - 115
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
- id: min_abs_improvement
  type: Real
  units: n/a
  required: true
  description: Positional argument `min_abs_improvement`.
- id: min_rel_improvement
  type: Real
  units: n/a
  required: true
  description: Positional argument `min_rel_improvement`.
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
  description: Return value of `hypr_material_improvement`.
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

# hypr_material_improvement

## Purpose
Decides whether a cost decrease is large enough to count as progress for early-stopping and stagnation logic.

## Design & Implementation
Returns false for a non-finite new cost and true for a non-finite reference, so the first finite result always registers. Otherwise the improvement `reference - new` must exceed the larger of the absolute threshold and the relative threshold times `max(|reference|, 1)`, the floor of one preventing a vanishing relative threshold near zero cost.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `new_cost` | Real | n/a | yes | Positional argument `new_cost`. |
| in | `reference_cost` | Real | n/a | yes | Positional argument `reference_cost`. |
| in | `min_abs_improvement` | Real | n/a | yes | Positional argument `min_abs_improvement`. |
| in | `min_rel_improvement` | Real | n/a | yes | Positional argument `min_rel_improvement`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `hypr_material_improvement`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_rpo_pso_material_improvement|rpo_pso_material_improvement]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:163-163`
- [[gnc.swarm_and_retiming__robot_arm_hypr_material_improvement|_robot_arm_hypr_material_improvement]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:55-55`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:116-116`
<!-- vulcan:connections:end -->

## Limitations
Using the larger of the two thresholds means both must be satisfied, which is stricter than the source's docstring wording of absolute or relative; a caller tuning only one threshold can be surprised by the other still gating.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 115.
