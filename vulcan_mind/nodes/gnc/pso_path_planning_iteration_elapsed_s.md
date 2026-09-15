---
id: gnc.pso_path_planning_iteration_elapsed_s
label: iteration_elapsed_s
kind: function
source:
  file: src/gnc/guidance/rpo/hypr/pso_path_planning.jl
  symbol: iteration_elapsed_s
  lines:
  - 254
  - 254
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
  description: Return value of `iteration_elapsed_s`. Returns `(time_ns() - iter_start_ns)
    / 1.0e9`.
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

# iteration_elapsed_s

## Purpose
Reports wall-clock seconds elapsed since an iteration began, for the per-iteration runtime budget.

## Design & Implementation
A closure inside `rpo_pso_plan_path` computing `(time_ns() - iter_start_ns) / 1e9` from a `UInt64` start stamp. Using `time_ns` rather than `time()` avoids a floating-point subtraction of two large epoch values.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `iter_start_ns` | UInt64 | n/a | yes | Positional argument `iter_start_ns`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `iteration_elapsed_s`. Returns `(time_ns() - iter_start_ns) / 1.0e9`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[gnc.pso_path_planning_iteration_timed_out|iteration_timed_out]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:261-261`
- [[gnc.pso_path_planning_record_iteration_timeout_bang|record_iteration_timeout!]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:279-279`
- [[gncy.pso_path_planning_rpo_pso_plan_path|rpo_pso_plan_path]] · `callees` → `callers` · call · `src/gnc/guidance/rpo/hypr/pso_path_planning.jl:254-254`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Wall clock includes time the thread spent descheduled, so under heavy contention an iteration can time out without having done its budgeted work.

## Provenance
Mapped from `src/gnc/guidance/rpo/hypr/pso_path_planning.jl` line 254.
