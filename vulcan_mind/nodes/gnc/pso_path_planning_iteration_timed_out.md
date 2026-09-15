---
id: gnc.pso_path_planning_iteration_timed_out
label: iteration_timed_out
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: iteration_timed_out
  lines:
  - 259
  - 259
inputs:
- id: iter_start_ns
  type: UInt64
  units: n/a
  required: true
  description: Positional argument `iter_start_ns`.
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
  description: Return value of `iteration_timed_out`. Returns `isfinite(limit) &&
    limit > 0.0 && iteration_elapsed_s(iter_start_ns) >= limit`.
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

# iteration_timed_out

## Purpose
Checks whether the current iteration has exhausted its runtime allowance.

## Design & Implementation
A closure reading `cfg.iteration_runtime_limit_s`, returning true only when the limit is finite, positive and `iteration_elapsed_s` has reached it. It is polled between every phase of an iteration — after reexploration, stagnation learning, culling, early-stopping checks, the velocity update and evaluation — so a timeout is attributed to the phase that crossed it.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `iter_start_ns` | UInt64 | n/a | yes | Positional argument `iter_start_ns`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `iteration_timed_out`. Returns `isfinite(limit) && limit > 0.0 && iteration_elapsed_s(iter_start_ns) >= limit`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:533-533`
- [[gnc.pso_path_planning_evaluate_swarm_bang|evaluate_swarm!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:408-408`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:259-259`

**Downstream**

- `callees` → [[gnc.pso_path_planning_iteration_elapsed_s|iteration_elapsed_s]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:261-261`
<!-- vulcan:connections:end -->

## Limitations
It is a poll, not a preemption, so a single long phase such as a threaded swarm evaluation overruns the limit by up to one full phase.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 259.
