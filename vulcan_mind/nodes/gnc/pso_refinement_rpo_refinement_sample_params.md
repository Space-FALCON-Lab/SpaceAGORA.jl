---
id: gnc.pso_refinement_rpo_refinement_sample_params
label: rpo_refinement_sample_params
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_refinement.jl
  symbol: rpo_refinement_sample_params
  lines:
  - 103
  - 103
inputs:
- id: samples
  type: Any
  units: n/a
  required: true
  description: Positional argument `samples`.
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
  description: Return value of `rpo_refinement_sample_params`. Returns `collect(range(0.0,
    1.0; length=n))` or `s ./ total`.
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

# rpo_refinement_sample_params

## Purpose
Assigns each sampled point a parameter in the unit interval proportional to its arc length along the path, so the Bezier fit places control points where the curve actually spends its length.

## Design & Implementation
Returns a zero vector for one or zero points. Otherwise it computes cumulative arc length with `rpo_arc_length_params` and divides by the total. If the total is at or below machine epsilon — every sample coincident — it falls back to `range(0, 1)` so the fit matrix is still well-posed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `samples` | Any | n/a | yes | Positional argument `samples`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_refinement_sample_params`. Returns `collect(range(0.0, 1.0; length=n))` or `s ./ total`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_refinement_rpo_fit_bezier_fixed_endpoints|rpo_fit_bezier_fixed_endpoints]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:128-128`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl`

**Downstream**

- `callees` → [[gnc.path_retiming_rpo_arc_length_params|rpo_arc_length_params]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_refinement.jl:107-107`
<!-- vulcan:connections:end -->

## Limitations
Chord-length parameterisation is an approximation to the true Bezier parameter, so the fit minimises error at slightly wrong parameters and can lag on strongly curved paths.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_refinement.jl` line 103.
