---
id: simulation.monte_carlo__run_monte_carlo_process
label: _run_monte_carlo_process
kind: function
source:
  file: src/simulation/campaigns/monte_carlo.jl
  symbol: _run_monte_carlo_process
  lines:
  - 165
  - 165
inputs:
- id: f
  type: Any
  units: n/a
  required: true
  description: Positional argument `f`.
- id: seeds
  type: Vector
  units: n/a
  required: true
  description: Positional argument `seeds`.
- id: spec
  type: MonteCarloSpec
  units: n/a
  required: true
  description: Positional argument `spec`.
- id: worker_ids
  type: Vector{Int}
  units: n/a
  required: true
  description: Positional argument `worker_ids`.
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
  description: Return value of `_run_monte_carlo_process`. Returns `completed`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _run_monte_carlo_process

## Purpose
The process-backend counterpart of the threaded runner, dispatching each sample to a Distributed worker while keeping the same queue, ordering and fail-fast semantics.

## Design & Implementation
Uses the identical channel-of-jobs structure but each dispatcher is an `@async` task, not `Threads.@spawn`, because it only blocks on IPC waiting for `remotecall_fetch` and should not occupy an OS thread. The user function is wrapped in a `run_sample` closure sent through a `CachingPool` over `worker_ids`, so the potentially large configuration-capturing closure is serialized to each worker once and referenced by identity thereafter. A `try`/`finally` calls `Distributed.clear!(pool)` afterwards so captured state does not outlive the campaign on the workers. Results land at their own index and are filtered and fail-fast-checked exactly as in the threaded path.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `f` | Any | n/a | yes | Positional argument `f`. |
| in | `seeds` | Vector | n/a | yes | Positional argument `seeds`. |
| in | `spec` | MonteCarloSpec | n/a | yes | Positional argument `spec`. |
| in | `worker_ids` | Vector{Int} | n/a | yes | Positional argument `worker_ids`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_run_monte_carlo_process`. Returns `completed`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/campaigns/monte_carlo.jl`
- [[simulation.adaptive_routing__run_campaign_with_route_env|_run_campaign_with_route_env]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:187-187`

**Downstream**

- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:214-214`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:214-214`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:214-214`
- `callees` → [[simulation.monte_carlo__run_monte_carlo_sample|_run_monte_carlo_sample]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:176-176`
- `callees` → [[simulation.monte_carlo__throw_first_monte_carlo_failure|_throw_first_monte_carlo_failure]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:203-203`
- `callees` → [[simulation.run_monte_carlo|run_monte_carlo]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:209-209`
- `callees` → [[simulation.run_simulation|run_simulation]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:223-223`
- `callees` → [[simx.engine_execution_run_simulation|run_simulation]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:223-223`
- `callees` → [[spaceagora.run_simulation|run_simulation]] · `callers` · call · `src/simulation/campaigns/monte_carlo.jl:223-223`
<!-- vulcan:connections:end -->

## Limitations
Every sample result, including a large `value`, crosses the process boundary by serialization, so returning full solutions from `f` is expensive here in a way it is not on threads. A worker that dies mid-sample surfaces as a thrown `ProcessExitedException` from `remotecall_fetch` rather than as a recorded failed sample.

## Provenance
Mapped from `src/simulation/campaigns/monte_carlo.jl` line 165.
