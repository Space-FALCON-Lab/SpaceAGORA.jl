---
id: gnc.rrt_connect_rpo_rrt_steer
label: rpo_rrt_steer
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/rrt_connect.jl
  symbol: rpo_rrt_steer
  lines:
  - 63
  - 63
inputs:
- id: q_near
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `q_near`.
- id: q_target
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `q_target`.
- id: step_size_m
  type: Real
  units: n/a
  required: true
  description: Positional argument `step_size_m`.
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
  description: Return value of `rpo_rrt_steer`. Returns `hypr_rrt_steer(q_near, q_target,
    step_size_m)`.
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

# rpo_rrt_steer

## Purpose
Moves from `q_near` toward `q_target` by at most `step_size_m`, returning the new point and a status of `:reached`, `:advanced`, or `:trapped`.

## Design & Implementation
Delegates to `hypr_rrt_steer(q_near, q_target, step_size_m)`. The convention consumed by callers is: `:reached` when the target lies within one step, `:advanced` when a truncated step was taken, `:trapped` when no progress is possible (zero-length direction). `rpo_rrt_extend!` and `rpo_rrt_star_plan_path` skip the iteration on `:trapped`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `q_near` | SVector{3, Float64} | n/a | yes | Positional argument `q_near`. |
| in | `q_target` | SVector{3, Float64} | n/a | yes | Positional argument `q_target`. |
| in | `step_size_m` | Real | n/a | yes | Positional argument `step_size_m`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `rpo_rrt_steer`. Returns `hypr_rrt_steer(q_near, q_target, step_size_m)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.rrt_connect_rpo_rrt_extend_bang|rpo_rrt_extend!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:162-162`
- [[gnc.rrt_connect_rpo_rrt_star_plan_path|rpo_rrt_star_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:548-548`
- [[module.gnc|CommandTypes]] · `api` → `module_api` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl`

**Downstream**

- `callees` → [[gnc.hypr_utils_hypr_rrt_steer|hypr_rrt_steer]] · `callers` · call · `src/gnc/guidance/rpo/hypr/rrt_connect.jl:64-64`
<!-- vulcan:connections:end -->

## Limitations
The step size is accepted as any `Real` and is not validated for positivity here. Because the wrapper adds no behaviour, changes to the shared `hypr_rrt_steer` status vocabulary silently change RPO planner control flow. The returned point is not collision-checked; that is the caller's responsibility.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/rrt_connect.jl` line 63.
