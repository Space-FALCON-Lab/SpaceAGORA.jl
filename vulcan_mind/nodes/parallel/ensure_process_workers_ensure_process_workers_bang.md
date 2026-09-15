---
id: parallel.ensure_process_workers_ensure_process_workers_bang
label: ensure_process_workers!
kind: function
source:
  file: src/parallel/process/worker_pool.jl
  symbol: ensure_process_workers!
  lines:
  - 193
  - 215
inputs:
- id: module_api
  type: Module
  units: n/a
  required: false
  description: Re-exported through the owning module's public surface.
- id: n
  type: Int
  units: workers
  required: true
  description: Desired number of active process workers selected by the route policy.
- id: pool
  type: ProcessPool
  units: n/a
  required: true
  description: Mutable worker-pool state to grow, reuse, or reconcile.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
tags:
- parallel
- processes
charts:
- parallel
origin: agent
outputs:
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
---

# ensure_process_workers!

## Purpose
`ensure_process_workers!` reconciles a process pool with the worker count requested by the selected route. It is the lifecycle bridge between route selection and campaign execution: the campaign asks for a count, and this function makes the pool’s active worker set meet that request before sampling starts.

## Theory & Math
The operation is an integer reconciliation problem. If `n_current` is the current worker count and `n_target` is the requested count, the function adds `max(n_target-n_current, 0)` workers and removes surplus workers according to the pool policy. The result is a pool whose active count is intended to equal `n_target`.

## Model & Assumptions
`n` must be nonnegative and meaningful for the host. The process pool must contain enough information to identify reusable workers and to update its lifecycle state. Worker initialization is assumed to load the project environment and any methods required by the campaign function.

## Design & Implementation
The implementation inspects the current pool, launches missing workers, updates the worker collection, and returns or mutates the pool used by the caller. Adaptive routing supplies `n` from the selected route, while `ProcessPool` supplies the mutable ownership record. Errors from worker launch or project setup propagate to the campaign boundary.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `n` | Int | workers | yes | Desired number of active process workers selected by the route policy. |
| in | `pool` | ProcessPool | n/a | yes | Mutable worker-pool state to grow, reuse, or reconcile. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.worker_pool__warm_gram_wrapper_bang|_warm_gram_wrapper!]] · `callees` → `callers` · feedback · `src/parallel/process/worker_pool.jl:166-166`
- [[simulation.adaptive_routing__run_campaign_with_route_env|_run_campaign_with_route_env]] · `callees` → `callers` · call · `src/simulation/campaigns/adaptive_routing.jl:184-184`

**Downstream**

- `callees` → [[parallel.worker_pool__bootstrap_process_worker_bang|_bootstrap_process_worker!]] · `callers` · call · `src/parallel/process/worker_pool.jl:200-200`
- `callees` → [[parallel.worker_pool__process_worker_exeflags|_process_worker_exeflags]] · `callers` · call · `src/parallel/process/worker_pool.jl:198-198`
<!-- vulcan:connections:end -->

## Limitations
The function cannot guarantee worker health after creation, and repeated reconciliation can be expensive if routes oscillate. Removing workers may discard warmed state. The process count is not a substitute for controlling native-library concurrency; calls into shared SPICE/GRAM state still require the runtime lock.

## Provenance
Mapped from `src/parallel/process/worker_pool.jl:165-215`.
