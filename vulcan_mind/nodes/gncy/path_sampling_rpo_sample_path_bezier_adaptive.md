---
id: gncy.path_sampling_rpo_sample_path_bezier_adaptive
label: rpo_sample_path_bezier_adaptive
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_sampling.jl
  symbol: rpo_sample_path_bezier_adaptive
  lines:
  - 197
  - 254
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: control_points
  type: Matrix
  units: m
  required: true
  description: Three-row matrix of Bezier control points defining the candidate RPO
    path in the RTN frame.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: samples
  type: Matrix
  units: m
  description: Three-row matrix of sampled path positions spaced finely near obstacles
    and coarsely in open space.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- gnc
charts:
- gncy
origin: agent
---

# rpo_sample_path_bezier_adaptive

## Purpose
`rpo_sample_path_bezier_adaptive` discretises a Bezier RPO path with spacing that tightens near obstacles and relaxes in open space. Adaptive sampling is what makes the cost evaluation affordable, because uniform sampling fine enough to resolve a close approach would be wastefully dense across the rest of the transfer.

## Theory & Math
For a Bezier curve $B(t)$ the parameter step for a desired spatial step $\Delta s$ is $\Delta t \approx \Delta s / \lVert B'(t) \rVert$, with $\lVert B'(t)\rVert$ approximated by $\lVert B(t+\varepsilon)-B(t)\rVert / \varepsilon$ at $\varepsilon = 10^{-4}$.

## Model & Assumptions
The step length at each point is a function of local clearance, interpolating between `min_ds` and `max_ds` according to `adaptive_sampling_far_clearance_m` and the shaping exponent `adaptive_sampling_power`. Because a Bezier curve is parameterised by a non-arc-length variable, the routine converts a desired spatial step into a parameter step by dividing by a finite-difference speed estimate. `rpo_bezier_speed_estimate` perturbs the parameter by 1e-4, falling back to a backward difference at the upper end, and returns zero when neither perturbation moves the parameter.

## Design & Implementation
The minimum step comes from `rpo_adaptive_sampling_min_ds_m` given the base step, geometry, and safe distance, and the maximum is the larger of `adaptive_sampling_max_ds_m` and that minimum. A length reference is the largest of the control-polygon length, the straight-line endpoint distance, and the minimum step. The main loop marches the parameter from zero toward one, halving the parameter step while the resulting chord exceeds 1.25 times the target spacing, with a floor of 1e-6 on the step. A step budget of `ceil(length_ref / min_ds) + 2` guards against non-termination. After the loop the endpoint is appended if the parameter fell short, and both the first and last columns are overwritten with the exact control-point endpoints so the sampled path starts and ends exactly where the planner intended.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `control_points` | Matrix | m | yes | Three-row matrix of Bezier control points defining the candidate RPO path in the RTN frame. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `samples` | Matrix | m | — | Three-row matrix of sampled path positions spaced finely near obstacles and coarsely in open space. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:275-275`

**Downstream**

- `callees` → [[gnc.clearance_rpo_clearance_distance_to_station|rpo_clearance_distance_to_station]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:218-218`
- `callees` → [[gnc.path_sampling_rpo_adaptive_sampling_min_ds_m|rpo_adaptive_sampling_min_ds_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:206-206`
- `callees` → [[gnc.path_sampling_rpo_adaptive_sampling_step_m|rpo_adaptive_sampling_step_m]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:219-219`
- `callees` → [[gnc.path_sampling_rpo_bezier_point_bang|rpo_bezier_point!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:213-213`
- `callees` → [[gnc.path_sampling_rpo_bezier_speed_estimate|rpo_bezier_speed_estimate]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:227-227`
- `callees` → [[gnc.path_sampling_rpo_path_length|rpo_path_length]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:208-208`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_sampling.jl:214-214`
<!-- vulcan:connections:end -->

## Limitations
The speed estimate is a two-point finite difference, so on curves with rapidly changing speed the initial parameter step can badly overshoot before the halving loop corrects it. Clamping the parameter step at 1e-6 means a curve with an extreme speed ratio can exhaust the step budget and terminate with a coarser tail than requested. Overwriting the endpoints hides any accumulated parameter drift rather than correcting the interior spacing that caused it.

## Provenance
Mapped from path_sampling.jl lines 187-254; include site observed at guidance_hooks.jl line 67.
