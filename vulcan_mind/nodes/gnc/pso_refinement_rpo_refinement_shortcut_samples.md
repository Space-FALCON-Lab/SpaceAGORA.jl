---
id: gnc.pso_refinement_rpo_refinement_shortcut_samples
label: rpo_refinement_shortcut_samples
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refinement_shortcut_samples
  lines:
  - 69
  - 69
inputs:
- id: samples
  type: Any
  units: n/a
  required: true
  description: Positional argument `samples`.
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
  description: Return value of `rpo_refinement_shortcut_samples`. Returns `current`.
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

# rpo_refinement_shortcut_samples

## Purpose
Removes unnecessary detours from a sampled path by replacing runs of points with a straight chord wherever the chord is safe.

## Design & Implementation
Returns immediately for two or fewer columns. For up to `refinement_waypoint_passes` passes it walks an index `i` and, for each, tries the farthest possible `j` first, working backwards to `i + 2`, accepting the first chord that is no longer than the path it replaces and passes `rpo_refinement_segment_is_safe`. Acceptance rebuilds the matrix keeping columns `1:i` and `j:end` and restarts from the same `i`; otherwise `i` advances. A pass with no change ends the loop early.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Any | n/a | yes | Positional argument `samples`. |
| in | `geometry` | Any | n/a | yes | Positional argument `geometry`. |
| in | `cfg` | RPOPSOConfig | n/a | yes | Positional argument `cfg`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_refinement_shortcut_samples`. Returns `current`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_refinement_rpo_refine_shortcut_refit|rpo_refine_shortcut_refit]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:183-183`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`

**Downstream**

- `callees` → [[gnc.path_sampling_rpo_path_length|rpo_path_length]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:79-79`
- `callees` → [[gnc.pso_refinement_rpo_refinement_segment_is_safe|rpo_refinement_segment_is_safe]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:82-82`
<!-- vulcan:connections:end -->

## Limitations
Each acceptance reallocates the whole matrix, and the farthest-first search is quadratic in the number of samples per pass, so dense sampling makes this the most expensive refinement step; the length test with a 1e-9 tolerance always passes for a straight chord, so it is effectively only the safety check that gates removal.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 69.
