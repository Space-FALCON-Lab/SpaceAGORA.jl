---
id: gnc.path_sampling_rpo_adaptive_sampling_step_m
label: rpo_adaptive_sampling_step_m
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_adaptive_sampling_step_m
  lines:
  - 92
  - 92
inputs:
- id: clearance
  type: Real
  units: n/a
  required: true
  description: Positional argument `clearance`.
- id: min_ds
  type: Real
  units: n/a
  required: true
  description: Positional argument `min_ds`.
- id: max_ds
  type: Real
  units: n/a
  required: true
  description: Positional argument `max_ds`.
- id: far_clearance_m
  type: Real
  units: n/a
  required: true
  description: Positional argument `far_clearance_m`.
- id: power
  type: Real
  units: n/a
  required: true
  description: Positional argument `power`.
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
  description: Return value of `rpo_adaptive_sampling_step_m`. Returns `min(step,
    max_step, clearance_excess + min_step)`.
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

# rpo_adaptive_sampling_step_m

## Purpose
Chooses the local sample spacing at one point from how much clearance it has: dense near the station, sparse in free space.

## Theory & Math
$$
u = \operatorname{clamp}\left(\frac{\max(c - d_s, 0)}{d_{\text{far}}}, 0, 1\right),\qquad \Delta s = \min\left( \Delta s_{\min} + u^{p}\,(\Delta s_{\max} - \Delta s_{\min}),\; \Delta s_{\max},\; (c - d_s) + \Delta s_{\min} \right)
$$

with $c$ the clearance, $d_s$ the safety distance, $d_{\text{far}}$ the clearance at which spacing reaches its maximum and $p$ the blending power.

## Design & Implementation
Clamps `min_ds` at 1e-9, `max_ds` at no less than `min_ds`, and `far_clearance_m` at no less than `min_ds`. The excess clearance beyond `safe_distance_m` is normalised by the far distance into `u` in the unit interval, and the step blends from minimum to maximum as `u^power`. The result is then capped a second time at `clearance_excess + min_ds`, guaranteeing a sample can never step past the point where clearance would reach the safety distance.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `clearance` | Real | n/a | yes | Positional argument `clearance`. |
| in | `min_ds` | Real | n/a | yes | Positional argument `min_ds`. |
| in | `max_ds` | Real | n/a | yes | Positional argument `max_ds`. |
| in | `far_clearance_m` | Real | n/a | yes | Positional argument `far_clearance_m`. |
| in | `power` | Real | n/a | yes | Positional argument `power`. |
| in | `safe_distance_m` | Real | n/a | no | Keyword argument `safe_distance_m` (default `0.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_adaptive_sampling_step_m`. Returns `min(step, max_step, clearance_excess + min_step)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_adaptive_segment_samples|rpo_adaptive_segment_samples]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:135-135`
- [[gncy.path_sampling_rpo_sample_path_bezier_adaptive|rpo_sample_path_bezier_adaptive]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:219-219`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:100-100`
<!-- vulcan:connections:end -->

## Limitations
A `power` below one makes the step grow steeply as soon as the path leaves the safety band, which can undersample moderately close regions; nothing validates `power` is positive.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/path_sampling.jl` line 92.
