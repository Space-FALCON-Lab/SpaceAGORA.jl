---
id: gnc.pso_refinement_rpo_try_accept_refinement
label: rpo_try_accept_refinement
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_try_accept_refinement
  lines:
  - 164
  - 164
inputs:
- id: candidate
  type: Any
  units: n/a
  required: true
  description: Positional argument `candidate`.
- id: geometry
  type: Any
  units: n/a
  required: true
  description: Positional argument `geometry`.
- id: cfg
  type: RPOPSOConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: current_components
  type: Any
  units: n/a
  required: true
  description: Positional argument `current_components`.
- id: safe_distance_m
  type: Real
  units: n/a
  required: false
  description: Keyword argument `safe_distance_m` (default `0.0`).
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
  description: Return value of `rpo_try_accept_refinement`. Returns `clamped, comps,
    true` or `nothing, current_components, false`.
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

# rpo_try_accept_refinement

## Purpose
The single gate through which every refinement candidate passes: clamp it, cost it, and keep it only if it beats the current path.

## Design & Implementation
Clamps the candidate into the search box, evaluates `rpo_normalized_path_cost_components` with the geometry and safety distance, and calls `rpo_refinement_better`. On success it returns the clamped path, its components and `true`; otherwise `nothing`, the unchanged current components and `false`. Returning the components alongside the path means callers never re-cost an accepted candidate.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `candidate` | Any | n/a | yes | Positional argument `candidate`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `current_components` | Any | n/a | yes | Positional argument `current_components`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_try_accept_refinement`. Returns `clamped, comps, true` or `nothing, current_components, false`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_refinement_rpo_refine_lower_degree|rpo_refine_lower_degree]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:247-247`
- [[gnc.pso_refinement_rpo_refine_shortcut_refit|rpo_refine_shortcut_refit]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:186-186`
- [[gnc.pso_refinement_rpo_refine_tighten_handles|rpo_refine_tighten_handles]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:210-210`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`

**Downstream**

- `callees` → [[gnc.path_costs_rpo_normalized_path_cost_components|rpo_normalized_path_cost_components]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:166-166`
- `callees` → [[gnc.pso_refinement_rpo_refinement_better|rpo_refinement_better]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:167-167`
- `callees` → [[gnc.pso_refinement_rpo_refinement_clamp_path|rpo_refinement_clamp_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:165-165`
<!-- vulcan:connections:end -->

## Limitations
The cost evaluation resamples the whole candidate path, so a refinement pass that proposes many candidates — handle tightening tries up to nine per handle — spends most of its time here.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 164.
