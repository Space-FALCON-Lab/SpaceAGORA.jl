---
id: gnc.rrt_warmstart__robot_arm_rrt_random_state
label: _robot_arm_rrt_random_state
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl
  symbol: _robot_arm_rrt_random_state
  lines:
  - 17
  - 17
inputs:
- id: rng
  type: Any
  units: n/a
  required: true
  description: Positional argument `rng`.
- id: lo
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `lo`.
- id: hi
  type: AbstractVector{<:Real}
  units: n/a
  required: true
  description: Positional argument `hi`.
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
  description: Return value of `_robot_arm_rrt_random_state`. Returns `q`.
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

# _robot_arm_rrt_random_state

## Purpose
`_robot_arm_rrt_random_state` draws a uniformly random joint configuration inside the per-joint limits, providing the exploration samples for the RRT-Connect warm start when the goal-bias coin flip does not fire.

## Design & Implementation
Signature `(rng, lo::AbstractVector{<:Real}, hi::AbstractVector{<:Real})`. It allocates `q = zeros(length(lo))` and, in an `@inbounds` loop, sets `q[i] = lo[i] + rand(rng) * (hi[i] - lo[i])` so each joint is sampled independently and uniformly on `[lo[i], hi[i])` (radians). The supplied `rng` makes the sequence reproducible. `lo` and `hi` are built by the caller from `joint.lower_rad` and `joint.upper_rad` of `model.joints`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rng` | Any | n/a | yes | Positional argument `rng`. |
| in | `lo` | AbstractVector{<:Real} | n/a | yes | Positional argument `lo`. |
| in | `hi` | AbstractVector{<:Real} | n/a | yes | Positional argument `hi`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_rrt_random_state`. Returns `q`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.rrt_warmstart__robot_arm_rrt_connect_warmstart_path|_robot_arm_rrt_connect_warmstart_path]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl:257-257`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`hi` is assumed to be at least as long as `lo`; with `@inbounds` a shorter `hi` reads out of bounds rather than throwing. Limits with `hi < lo` silently sample outside the intended range. Sampling is uniform in joint space, which does not correspond to uniform coverage of the end-effector workspace. No rejection of self-colliding configurations is performed here; collision is checked only on edges.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/rrt_warmstart.jl` line 17.
