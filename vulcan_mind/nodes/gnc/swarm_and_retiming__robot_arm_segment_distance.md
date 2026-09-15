---
id: gnc.swarm_and_retiming__robot_arm_segment_distance
label: _robot_arm_segment_distance
kind: function
source:
  file: src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl
  symbol: _robot_arm_segment_distance
  lines:
  - 433
  - 433
inputs:
- id: p
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `p`.
- id: a
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `a`.
- id: b
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `b`.
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
  description: Return value of `_robot_arm_segment_distance`. Returns `norm(p - (a
    + t * ab))`.
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

# _robot_arm_segment_distance

## Purpose
Computes the shortest Euclidean distance from a point to a finite 3-D line segment, used for link-versus-obstacle clearance checks.

## Theory & Math
$$t = \operatorname{clamp}\!\left(\frac{(\mathbf{p}-\mathbf{a})\cdot(\mathbf{b}-\mathbf{a})}{\|\mathbf{b}-\mathbf{a}\|^2},\,0,\,1\right),\qquad d = \left\|\mathbf{p} - \big(\mathbf{a} + t(\mathbf{b}-\mathbf{a})\big)\right\|$$ where $\mathbf{p}$ is the query point and $\mathbf{a}$, $\mathbf{b}$ are the segment endpoints (m).

## Design & Implementation
Takes `p`, `a`, `b` as `SVector{3, Float64}`. With `ab = b - a` and `denom = dot(ab, ab)`, a degenerate segment (`denom <= eps(Float64)`) returns `norm(p - a)`. Otherwise the projection parameter `t = clamp(dot(p - a, ab) / denom, 0, 1)` locates the closest point `a + t * ab` and the function returns `norm(p - (a + t*ab))`. Fully allocation-free thanks to static vectors.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | SVector{3, Float64} | n/a | yes | Positional argument `p`. |
| in | `a` | SVector{3, Float64} | n/a | yes | Positional argument `a`. |
| in | `b` | SVector{3, Float64} | n/a | yes | Positional argument `b`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_robot_arm_segment_distance`. Returns `norm(p - (a + t * ab))`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gncz.clearance_robot_arm_clearance_stats_from_samples|robot_arm_clearance_stats_from_samples]] · `callees` → `callers` · call · `src/gnc/robotics/robot_arm_hypr/clearance.jl:27-27`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The degeneracy threshold `eps(Float64)` is on squared length, so segments shorter than about `1.5e-8` m are treated as points. Only `Float64` static vectors are accepted; other element types need conversion by the caller.

## Provenance
Mapped from `src/gnc/robotics/robot_arm_hypr/swarm_and_retiming.jl` line 433.
