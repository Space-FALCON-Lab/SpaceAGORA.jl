---
id: gnc.pso_refinement_rpo_refine_tighten_handles
label: rpo_refine_tighten_handles
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refine_tighten_handles
  lines:
  - 190
  - 190
inputs:
- id: path
  type: Any
  units: n/a
  required: true
  description: Positional argument `path`.
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
  description: Return value of `rpo_refine_tighten_handles`. Returns `current, current_components,
    improved`.
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

# rpo_refine_tighten_handles

## Purpose
One refinement strategy: move each interior Bezier handle part-way toward a straighter position and keep the move if the objective improves.

## Design & Implementation
For each interior column `j` it forms three targets — the projection of the handle onto the start-goal chord, its projection onto the segment between its neighbours, and the point at fraction `(j-1)/(n-1)` along the start-goal chord — and for each target tries blends at `λ` of 0.25, 0.5 and 0.75. The first accepted candidate updates `current` and its components and moves on to the next handle. The function returns the possibly updated path, components and whether anything was accepted.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `path` | Any | n/a | yes | Positional argument `path`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `current_components` | Any | n/a | yes | Positional argument `current_components`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_refine_tighten_handles`. Returns `current, current_components, improved`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:298-298`

**Downstream**

- `callees` → [[gnc.pso_refinement_rpo_refinement_project_to_segment|rpo_refinement_project_to_segment]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:202-202`
- `callees` → [[gnc.pso_refinement_rpo_try_accept_refinement|rpo_try_accept_refinement]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:210-210`
<!-- vulcan:connections:end -->

## Limitations
Handles are visited in order and each is settled before the next, so the result depends on column order and is not a joint optimisation; up to nine full cost evaluations per handle makes this expensive on long polygons.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 190.
