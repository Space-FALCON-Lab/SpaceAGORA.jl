---
id: parallel.worker_pool_shutdown_process_pool_bang
label: shutdown_process_pool!
kind: function
source:
  file: src/parallel/process/worker_pool.jl
  symbol: shutdown_process_pool!
  lines:
  - 224
  - 224
inputs:
- id: pool
  type: ProcessPool
  units: n/a
  required: true
  description: Positional argument `pool`.
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
  description: Return value of `shutdown_process_pool!`; mutates `pool` in place.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- parallel
charts:
- parallel
origin: agent
---

# shutdown_process_pool!

## Purpose
Tears down every worker in a process pool and empties it, mainly so tests do not leak processes between cases.

## Design & Implementation
Takes `pool.lock` for the whole operation, returns early when `pool.workers` is already empty, then calls `rmprocs` on the whole worker vector at once and clears it with `empty!`. Holding the lock across both steps means a concurrent `ensure_process_workers!` cannot observe a pool whose worker ids have been removed from Distributed but are still listed. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pool` | ProcessPool | n/a | yes | Positional argument `pool`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `shutdown_process_pool!`; mutates `pool` in place. |
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
Campaign code deliberately leaves the process-global pool warm across calls, so calling this on that pool discards the one-time SpaceAGORA and GRAMSuite precompilation cost and the next campaign pays it again; `rmprocs` is called without a timeout, so a wedged worker blocks the shutdown.

## Provenance
Mapped from `src/parallel/process/worker_pool.jl` line 224.
