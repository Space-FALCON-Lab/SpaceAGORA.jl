---
id: gnc.swarm_and_retiming__robot_arm_hypr_reference_times_from_scales
label: _robot_arm_hypr_reference_times_from_scales
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_reference_times_from_scales
  lines:
  - 142
  - 142
inputs:
- id: nominal_dt
  type: Float64
  units: n/a
  required: true
  description: Positional argument `nominal_dt`.
- id: scales
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `scales`.
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
  description: Return value of `_robot_arm_hypr_reference_times_from_scales`. Returns
    `t_ref`.
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

# _robot_arm_hypr_reference_times_from_scales

## Purpose
Converts a vector of per-segment time-stretch factors into a cumulative reference time grid for the retimed joint trajectory.

## Design & Implementation
Takes `nominal_dt::Float64` (seconds per unscaled segment) and `scales::AbstractVector{<:Real}`. Allocates `t_ref = zeros(length(scales) + 1)` and fills `t_ref[k + 1] = t_ref[k] + nominal_dt * scales[k]` under `@inbounds`, so `t_ref[1] = 0` and `t_ref[end]` is the total retimed duration. Returns the `Vector{Float64}`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `nominal_dt` | Float64 | n/a | yes | Positional argument `nominal_dt`. |
| in | `scales` | AbstractVector{<:Real} | n/a | yes | Positional argument `scales`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_hypr_reference_times_from_scales`. Returns `t_ref`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.swarm_and_retiming__robot_arm_hypr_retime_reference|_robot_arm_hypr_retime_reference]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:398-398`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:145-145`
<!-- vulcan:connections:end -->

## Limitations
Negative or zero scales are not rejected and would produce a non-monotonic grid, which later breaks `searchsortedfirst` in the cloth estimator. Floating-point accumulation error grows linearly with segment count; there is no compensated summation.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 142.
