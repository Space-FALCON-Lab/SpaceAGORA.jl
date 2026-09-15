---
id: gnc.hypr_utils_hypr_protected_particle_mask
label: hypr_protected_particle_mask
kind: function
source:
  file: src/gnc/hypr/hypr_utils.jl
  symbol: hypr_protected_particle_mask
  lines:
  - 129
  - 129
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
  description: Return value of `hypr_protected_particle_mask`. Returns `mask`.
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

# hypr_protected_particle_mask

## Purpose
Marks the best fraction of a swarm as protected so culling and stagnation learning leave elites untouched.

## Design & Implementation
Returns an all-false mask for an empty swarm. Otherwise it computes `elite_count` as the ceiling of the clamped fraction times `n`, clamped to between one and `n`, sorts particles by cost with `sortperm` and sets the mask for the first `elite_count`. At least one particle is always protected.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `costs` | Any | n/a | yes | Positional argument `costs`. |
| in | `elite_fraction` | Any | n/a | yes | Positional argument `elite_fraction`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `hypr_protected_particle_mask`. Returns `mask`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_rpo_pso_protected_particle_mask|rpo_pso_protected_particle_mask]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:158-158`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/hypr/hypr_utils.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/hypr/hypr_utils.jl:133-133`
<!-- vulcan:connections:end -->

## Limitations
The full `sortperm` is O(n log n) per call even though only the top few are needed; and infinite costs sort last, so with a swarm where every cost is `Inf` the protected set is arbitrary.

## Provenance
Mapped from `src/gnc/hypr/hypr_utils.jl` line 129.
