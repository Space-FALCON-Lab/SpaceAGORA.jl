---
id: parcore.parallel_process_parallelprocess
label: ParallelProcess
kind: struct
source:
  file: src/parallel/process/parallel_process.jl
  symbol: ParallelProcess
  lines:
  - 2
  - 10
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Worker pool implementation included by this module from the sibling
    process directory.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: process_api
  type: Module
  units: n/a
  description: 'Process-route surface: the ProcessPool record, the campaign pool accessor,
    worker provisioning and pool shutdown.'
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parcore
origin: agent
---

# ParallelProcess

## Purpose
`ParallelProcess` is the module wrapper for the process-based outer route. It includes the worker pool implementation and exports the four names the campaign layer needs to obtain, size and release a pool of worker processes.

## Model & Assumptions
The process route exists because worker processes do not share the coordinator's thread pool: each runs with a single thread and its own address space, so a process worker is bounded by physical parallelism rather than by the coordinator's `Threads.nthreads()`. That isolation is exactly what makes the route attractive for Monte Carlo campaigns, where samples are independent and shared mutable state would otherwise need locking.

## Design & Implementation
The module is a ten-line aggregator: it opens the module, includes the worker pool file and exports `ProcessPool`, `campaign_process_pool`, `ensure_process_workers!` and `shutdown_process_pool!`. The split between a mutable pool record and an `ensure` operation is deliberate — a campaign holds one pool across many sample batches and reconciles it toward a target worker count, rather than launching and tearing down processes per batch, because worker startup cost is dominated by loading the package in the new process.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Worker pool implementation included by this module from the sibling process directory. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `process_api` | Module | n/a | — | Process-route surface: the ProcessPool record, the campaign pool accessor, worker provisioning and pool shutdown. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/process/parallel_process.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Worker startup can fail on project resolution, memory pressure or serialisation of the campaign closure, and a process identifier recorded in the pool is not evidence that the worker is still healthy. Results must cross a process boundary, so anything the closure returns has to be serialisable and large returns pay a copy. A long-lived pool retains module state and memory in each worker, so campaigns that shift configuration substantially benefit from an explicit refresh rather than reuse.

## Provenance
Mapped from `src/parallel/process/parallel_process.jl:2-10`.
