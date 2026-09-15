---
id: gnc.pso_refinement_rpo_refine_shortcut_refit
label: rpo_refine_shortcut_refit
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refine_shortcut_refit
  lines:
  - 174
  - 174
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
  description: Return value of `rpo_refine_shortcut_refit`. Returns `rpo_try_accept_refinement(candidate,
    geometry, cfg, current_components; safe_dis`.
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

# rpo_refine_shortcut_refit

## Purpose
One refinement strategy: densely sample the current curve, cut out detours with straight shortcuts, and refit a Bezier of the same degree to what remains.

## Design & Implementation
Samples the path with `rpo_sample_path` at the refinement density, runs `rpo_refinement_shortcut_samples`, and returns early with no change if the shortcutting removed nothing. Otherwise it fits a control polygon with the same column count as the input through `rpo_fit_bezier_fixed_endpoints` and submits it to `rpo_try_accept_refinement`.

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
| out | `result` | Any | n/a | — | Return value of `rpo_refine_shortcut_refit`. Returns `rpo_try_accept_refinement(candidate, geometry, cfg, current_components; safe_dis`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncy.pso_refinement_rpo_post_refine_path|rpo_post_refine_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:284-284`

**Downstream**

- `callees` → [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:175-175`
- `callees` → [[gnc.pso_parameters_rpo_hypr_refinement_sampling_density_m|rpo_hypr_refinement_sampling_density_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:180-180`
- `callees` → [[gnc.pso_refinement_rpo_fit_bezier_fixed_endpoints|rpo_fit_bezier_fixed_endpoints]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:185-185`
- `callees` → [[gnc.pso_refinement_rpo_refinement_shortcut_samples|rpo_refinement_shortcut_samples]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:183-183`
- `callees` → [[gnc.pso_refinement_rpo_try_accept_refinement|rpo_try_accept_refinement]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:186-186`
<!-- vulcan:connections:end -->

## Limitations
Refitting to the original degree can reintroduce curvature the shortcut removed, so even a substantial shortcut may be rejected once refit; the strategy has no way to try an intermediate degree.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 174.
