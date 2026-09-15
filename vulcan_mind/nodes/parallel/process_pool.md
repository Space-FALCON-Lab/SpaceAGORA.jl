---
id: parallel.process_pool
label: ProcessPool
kind: function
source:
  file: src/parallel/process/worker_pool.jl
  symbol: ProcessPool
  lines:
  - 19
  - 37
outputs:
- id: pool
  type: ProcessPool
  units: n/a
  description: Mutable process-pool state containing worker identities and lifecycle
    metadata.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
- processes
charts:
- parallel
origin: agent
inputs:
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
---

# ProcessPool

## Purpose
`ProcessPool` is the state container used by campaign code to track worker processes. It separates the desired pool configuration from the lifecycle operation that ensures workers exist, allowing repeated Monte Carlo calls to reuse a pool rather than creating processes for every sample.

## Theory & Math
Pool sizing is discrete: the target worker count is an integer `n`, and the pool state records the set of active worker identifiers. The performance model is dominated by startup and communication overhead, not by a physical equation. Reusing workers reduces fixed launch cost when the campaign performs many independent samples.

## Model & Assumptions
Workers are assumed to run compatible Julia code and to be able to deserialize the campaign closure and its configuration. The pool’s ownership and mutation must be coordinated when multiple campaign tasks share it. A process identifier in the record is not evidence that the worker remains healthy after creation.

## Design & Implementation
`worker_pool.jl` declares `ProcessPool` near the beginning of the process subsystem. `ensure_process_workers!` receives this record, compares the desired count with the current worker set, and adds or removes workers according to the route policy. Campaign code passes the pool through its process route so worker management remains outside sample generation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `pool` | ProcessPool | n/a | — | Mutable process-pool state containing worker identities and lifecycle metadata. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/process/worker_pool.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Worker startup can fail because of project loading, memory, scheduler, or serialization errors. The pool does not model thread safety inside native libraries; simulation runtime locks still apply. A long-lived pool can retain stale module state and memory, so callers need an explicit shutdown or refresh policy for large campaigns.

## Provenance
Mapped from `src/parallel/process/worker_pool.jl:1-37`.
