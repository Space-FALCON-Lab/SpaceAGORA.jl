---
id: gnc.pso_path_planning_record_iteration_timeout_bang
label: record_iteration_timeout!
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: record_iteration_timeout!
  lines:
  - 272
  - 272
inputs:
- id: iter
  type: Any
  units: n/a
  required: true
  description: Positional argument `iter`.
- id: phase
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `phase`.
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
  type: Nothing
  units: n/a
  description: Return value of `record_iteration_timeout!`; mutates `iter` in place.
    Returns `nothing`.
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

# record_iteration_timeout!

## Purpose
Records a runtime-budget breach with the phase that caused it, for diagnostics in the returned result.

## Design & Implementation
A closure that sets the captured `optimization_timed_out`, `iteration_timeout_iter` and `iteration_timeout_phase`, and pushes a named tuple of iteration, phase and elapsed seconds onto `iteration_timeout_events`. The phase symbol is one of `reexplore_evaluate`, `reexplore`, `iteration_start`, `stagnation_learning`, `cull`, `early_stopping`, `velocity_update` or `swarm_evaluate`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `iter` | Any | n/a | yes | Positional argument `iter`. |
| in | `phase` | Symbol | n/a | yes | Positional argument `phase`. |
| in | `iter_start_ns` | UInt64 | n/a | yes | Positional argument `iter_start_ns`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `record_iteration_timeout!`; mutates `iter` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_apply_stagnation_learning_bang|apply_stagnation_learning!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:528-528`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:272-272`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:276-276`
- `callees` → [[gnc.pso_path_planning_iteration_elapsed_s|iteration_elapsed_s]] · `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:279-279`
<!-- vulcan:connections:end -->

## Limitations
Because the loop breaks immediately after each call, the events vector holds at most one entry in practice; it is a vector for forward compatibility rather than current need.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 272.
