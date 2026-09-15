---
id: gnc.pso_path_planning_rpo_pso_protected_particle_mask
label: rpo_pso_protected_particle_mask
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: rpo_pso_protected_particle_mask
  lines:
  - 157
  - 157
inputs:
- id: costs
  type: Any
  units: n/a
  required: true
  description: Positional argument `costs`.
- id: elite_fraction
  type: Any
  units: n/a
  required: true
  description: Positional argument `elite_fraction`.
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
  description: Return value of `rpo_pso_protected_particle_mask`. Returns `hypr_protected_particle_mask(costs,
    elite_fraction)`.
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

# rpo_pso_protected_particle_mask

## Purpose
Marks the elite fraction of particles that stagnation learning must not disturb.

## Design & Implementation
A thin forward to `hypr_protected_particle_mask` with the current costs and `elite_fraction`, which returns a boolean vector selecting the lowest-cost particles. Wrapping it keeps the RPO planner's calls self-contained.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `costs` | Any | n/a | yes | Positional argument `costs`. |
| in | `elite_fraction` | Any | n/a | yes | Positional argument `elite_fraction`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_pso_protected_particle_mask`. Returns `hypr_protected_particle_mask(costs, elite_fraction)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:444-444`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:444-444`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_protected_particle_mask|hypr_protected_particle_mask]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:158-158`
<!-- vulcan:connections:end -->

## Limitations
Protection is decided from the current-iteration cost vector rather than personal bests, so a normally strong particle having one bad iteration can lose protection and be perturbed.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 157.
