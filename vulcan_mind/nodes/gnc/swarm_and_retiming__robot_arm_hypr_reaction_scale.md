---
id: gnc.swarm_and_retiming__robot_arm_hypr_reaction_scale
label: _robot_arm_hypr_reaction_scale
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_hypr_reaction_scale
  lines:
  - 132
  - 132
inputs:
- id: cfg
  type: RobotArmHYPRConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: demand_ratio
  type: Real
  units: n/a
  required: true
  description: Positional argument `demand_ratio`.
- id: k
  type: Int
  units: n/a
  required: true
  description: Positional argument `k`.
- id: n_seg
  type: Int
  units: n/a
  required: true
  description: Positional argument `n_seg`.
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
  description: Return value of `_robot_arm_hypr_reaction_scale`. Returns `clamp(scale,
    cfg.retime_min_scale, cfg.retime_max_scale)`.
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

# _robot_arm_hypr_reaction_scale

## Purpose
Computes a per-segment time-stretch factor that slows the arm in the middle of a motion in proportion to how much it is already being stretched, as a proxy for spacecraft reaction demand when no explicit wrench limits are configured.

## Theory & Math
$$s_k = \operatorname{clamp}\!\left(1 + g\,\tau\,\max(\rho_k, 0)\,4\phi_k(1-\phi_k),\; s_{\min},\; s_{\max}\right),\qquad \phi_k = \frac{k-1}{\max(n_{\mathrm{seg}}-1,\,1)}$$ where $g$ is `retime_reaction_gain`, $\tau$ is `retime_reaction_time_s` (s), $\rho_k$ is the caller-supplied demand ratio for segment $k$, and $s_{\min}, s_{\max}$ are the configured scale bounds.

## Design & Implementation
Returns 1.0 immediately when `cfg.retime_reaction_time_s <= 0`. Otherwise `phase = clamp((k - 1) / max(n_seg - 1, 1), 0, 1)` maps the segment index `k` to [0, 1], and `taper = 4 * phase * (1 - phase)` is a parabola equal to 0 at both ends and 1 at mid-motion. The scale is `1 + retime_reaction_gain * retime_reaction_time_s * max(demand_ratio, 0) * taper`, clamped to `[cfg.retime_min_scale, cfg.retime_max_scale]`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | RobotArmHYPRConfig | n/a | yes | Positional argument `cfg`. |
| in | `demand_ratio` | Real | n/a | yes | Positional argument `demand_ratio`. |
| in | `k` | Int | n/a | yes | Positional argument `k`. |
| in | `n_seg` | Int | n/a | yes | Positional argument `n_seg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_hypr_reaction_scale`. Returns `clamp(scale, cfg.retime_min_scale, cfg.retime_max_scale)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.swarm_and_retiming__robot_arm_hypr_retime_reference|_robot_arm_hypr_retime_reference]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:387-387`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl:137-137`
<!-- vulcan:connections:end -->

## Limitations
The parabolic taper is a heuristic with no dynamic basis; it is only used on the code path where both `retime_max_base_force_n` and `retime_max_base_torque_nm` are infinite. For `n_seg = 1` the phase is 0 and the scale is always 1. The gain and time constant multiply, so their individual units are not independently meaningful.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 132.
