---
id: gnc.rrt_warmstart__robot_arm_rrt_segment_samples
label: _robot_arm_rrt_segment_samples
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_segment_samples
  lines:
  - 31
  - 31
inputs:
- id: q_from
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_from`.
- id: q_to
  type: Any
  units: n/a
  required: true
  description: Positional argument `q_to`.
- id: sample_ds_rad
  type: Real
  units: n/a
  required: true
  description: Positional argument `sample_ds_rad`.
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
  description: Return value of `_robot_arm_rrt_segment_samples`. Returns `samples`.
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

# _robot_arm_rrt_segment_samples

## Purpose
`_robot_arm_rrt_segment_samples` discretises a straight joint-space edge between two configurations into evenly spaced samples so that clearance and limit checks can be evaluated along the whole edge rather than only at its endpoints.

## Design & Implementation
Signature `(q_from, q_to, sample_ds_rad::Real)`. After converting both endpoints to `Vector{Float64}`, it computes `n_steps = max(1, ceil(Int, norm(q1 - q0) / max(sample_ds_rad, 1e-9)))` and allocates a `Matrix{Float64}` of size `n × (n_steps + 1)`. Column `i + 1` holds the linear interpolation `(1 - α) q0 + α q1` with `α = i / n_steps` for `i in 0:n_steps`, so the first and last columns are exactly the endpoints and consecutive samples are at most `sample_ds_rad` apart. The loop is `@inbounds` and writes with broadcast assignment.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_from` | Any | n/a | yes | Positional argument `q_from`. |
| in | `q_to` | Any | n/a | yes | Positional argument `q_to`. |
| in | `sample_ds_rad` | Real | n/a | yes | Positional argument `sample_ds_rad`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_rrt_segment_samples`. Returns `samples`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_warmstart__robot_arm_rrt_segment_is_safe|_robot_arm_rrt_segment_is_safe]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:54-54`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:35-35`
<!-- vulcan:connections:end -->

## Limitations
Sample spacing is measured in joint-space radians, so the corresponding end-effector displacement between samples depends on link lengths and can exceed obstacle radii for long arms, letting thin obstacles slip between samples. Very small `sample_ds_rad` values produce large matrices (edge length / ds columns) with no upper bound. Endpoints of differing length throw `DimensionMismatch` at the subtraction.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 31.
