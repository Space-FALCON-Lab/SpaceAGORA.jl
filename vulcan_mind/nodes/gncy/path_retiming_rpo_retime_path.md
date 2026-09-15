---
id: gncy.path_retiming_rpo_retime_path
label: rpo_retime_path
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/path_retiming.jl
  symbol: rpo_retime_path
  lines:
  - 115
  - 297
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: geometric_path
  type: Matrix
  units: m
  required: true
  description: Three-row geometric RTN path to be converted into a time-parameterised
    reference trajectory.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: timed_reference
  type: Tuple
  units: n/a
  description: Resampled positions, arc-length parameters, and the time stamps assigned
    to each sample by the speed profile.
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

# rpo_retime_path

## Purpose
`rpo_retime_path` converts a purely geometric RPO path into a time-parameterised reference trajectory. It assigns a speed to every sample based on how much clearance and how much curvature that sample has, then integrates the speed profile to produce the time stamps the tracker follows.

## Theory & Math
The clearance-limited speed solves $d = v\tau + v^2/(2a)$ for $v$, giving $v_{clear} = -a\tau + \sqrt{(a\tau)^2 + 2ad}$. The curvature limit is $v_{curve} = \sqrt{a/\kappa}$, and the applied speed is $v = \max\left(v_{\min},\ \min(s\cdot\min(v_{clear}, v_{curve}),\ v_{\max})\right)$.

## Model & Assumptions
Two independent speed limits are combined. The clearance limit comes from a braking argument: with maximum deceleration `retime_a_max_mps2` and reaction time `retime_reaction_time_s`, the fastest safe speed at available clearance is the positive root of the stopping-distance relation. The curvature limit comes from lateral acceleration, capping speed at the square root of maximum acceleration over local curvature. The smaller of the two is scaled by `retime_speed_scale`, clipped at `retime_max_speed_mps`, and finally raised to `retime_min_speed_mps` if that is positive.

## Design & Implementation
Sampling uses `rpo_sample_path` at `cfg.sample_ds_m` with the configured curve type, after which `rpo_remove_near_duplicate_samples` strips points closer than `duplicate_tol_m` so arc-length parameterisation and curvature estimation do not divide by near-zero spacing. Arc lengths come from `rpo_arc_length_params` and curvature from `rpo_curvature_from_samples`. A degenerate zero-length path emits a warning and returns a single-point reference. Clearance at each sample is queried through `rpo_clearance_to_station`. Samples whose speed limit evaluates to zero are collected into `infeasible_idxs` and assigned `fallback_speed`, itself clipped by the maximum speed and floored at machine epsilon, so a batch run degrades rather than stalling.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `geometric_path` | Matrix | m | yes | Three-row geometric RTN path to be converted into a time-parameterised reference trajectory. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `timed_reference` | Tuple | n/a | — | Resampled positions, arc-length parameters, and the time stamps assigned to each sample by the speed profile. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.rpo_reference_trajectory_rpo_reference_from_path|rpo_reference_from_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/rpo_reference_trajectory.jl:3-3`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:147-147`
- `callees` → [[gnc.path_retiming_rpo_arc_length_params|rpo_arc_length_params]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:133-133`
- `callees` → [[gnc.path_retiming_rpo_curvature_from_samples|rpo_curvature_from_samples]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:141-141`
- `callees` → [[gnc.path_retiming_rpo_interpolate_along_path|rpo_interpolate_along_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:250-250`
- `callees` → [[gnc.path_retiming_rpo_remove_near_duplicate_samples|rpo_remove_near_duplicate_samples]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:131-131`
- `callees` → [[gnc.path_sampling_rpo_sample_path|rpo_sample_path]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:123-123`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:193-193`
- `callees` → [[gncz.clearance_rpo_clearance_to_station|rpo_clearance_to_station]] · `callers` · call · `src/gnc/guidance/rpo/hypr/path_retiming.jl:166-166`
<!-- vulcan:connections:end -->

## Limitations
The fallback speed keeps the retimer running through geometrically infeasible regions but does not make those regions collision-free, as the inline comment states explicitly. Curvature is estimated from discrete samples, so under-sampled tight turns understate curvature and overstate the allowed speed. The braking model assumes clearance can be traded directly for stopping distance along the path tangent and ignores the direction of the nearest obstacle. Zero or negative `retime_a_max_mps2` collapses the clearance limit to zero and forces the whole path onto the fallback speed.

## Provenance
Mapped from path_retiming.jl lines 115-297; include site observed at guidance_hooks.jl line 66.
