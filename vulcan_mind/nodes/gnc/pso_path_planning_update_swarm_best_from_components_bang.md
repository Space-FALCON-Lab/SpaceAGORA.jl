---
id: gnc.pso_path_planning_update_swarm_best_from_components_bang
label: update_swarm_best_from_components!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: update_swarm_best_from_components!
  lines:
  - 380
  - 380
inputs:
- id: pidx
  type: Any
  units: n/a
  required: true
  description: Positional argument `pidx`.
- id: comps
  type: Any
  units: n/a
  required: true
  description: Positional argument `comps`.
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
  type: Nothing
  units: n/a
  description: Return value of `update_swarm_best_from_components!`; mutates `pidx`
    in place. Returns `nothing`.
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

# update_swarm_best_from_components!

## Purpose
Applies one particle's evaluation to its personal best and, if warranted, to the global best.

## Design & Implementation
A closure that, when `comps.total` beats `pbest_cost[pidx]`, copies the particle's position into `pbest`, records the total and obstacle cost, and resets its stagnation counter; otherwise it increments the counter. Independently, if the total beats `gbest_cost`, it copies the position into `gbest` and stores the components. Both updates are element-wise loops over `dim` to avoid allocating.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pidx` | Any | n/a | yes | Positional argument `pidx`. |
| in | `comps` | Any | n/a | yes | Positional argument `comps`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `update_swarm_best_from_components!`; mutates `pidx` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_evaluate_swarm_bang|evaluate_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:413-413`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:380-380`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Not thread-safe: the threaded evaluation path collects components first and applies this serially afterwards for that reason, and a future caller invoking it from inside `@threads` would race on `gbest`.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 380.
