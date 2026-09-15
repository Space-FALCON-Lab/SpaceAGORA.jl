---
id: parallel.types__spinbarrierpool
label: _SpinBarrierPool
kind: struct
source:
  file: src/parallel/policy/types.jl
  symbol: _SpinBarrierPool
  lines:
  - 127
  - 127
inputs:
- id: workers
  type: Int
  units: n/a
  required: true
  description: Field `workers`.
- id: worker_gen
  type: Vector{Threads.Atomic{Int}}
  units: n/a
  required: true
  description: Field `worker_gen`.
- id: done_count
  type: Threads.Atomic{Int}
  units: n/a
  required: true
  description: Field `done_count`.
- id: stop
  type: Threads.Atomic{Bool}
  units: n/a
  required: true
  description: Field `stop`.
- id: request
  type: Base.RefValue{Any}
  units: n/a
  required: true
  description: Field `request`.
- id: errors
  type: Vector{Union{Nothing, Base.CapturedException}}
  units: n/a
  required: true
  description: Field `errors`.
- id: run_lock
  type: ReentrantLock
  units: n/a
  required: true
  description: Field `run_lock`.
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
  type: _SpinBarrierPool
  units: n/a
  description: Constructed `_SpinBarrierPool`.
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

# _SpinBarrierPool

## Purpose
Lowest-latency worker pool in the parallel policy layer. Worker tasks busy-poll a per-worker atomic generation counter instead of sleeping on a channel, cutting dispatch latency to tens of nanoseconds so sub-microsecond kernels can still be parallelised across 32 to 128+ threads at the cost of burning idle CPU.

## Design & Implementation
Plain `mutable struct` with `workers::Int`, `worker_gen::Vector{Threads.Atomic{Int}}`, `done_count::Threads.Atomic{Int}`, `stop::Threads.Atomic{Bool}`, a shared `request::Base.RefValue{Any}`, an `errors::Vector{Union{Nothing, Base.CapturedException}}` with one slot per worker, and `run_lock::ReentrantLock`. The outer constructor `_SpinBarrierPool(workers::Int)` zero-initialises every atomic and fills `errors` with `nothing`. `_create_spin_barrier_pool` clamps `workers` to `Threads.nthreads() - 1` and spawns `_spin_barrier_worker_loop_w` per worker. Dispatch writes `request[]`, bumps `worker_gen[w]` for `workers - 1` pool workers, runs the last worker slot on the coordinator thread itself, then spins on `done_count` before subtracting `pool_workers`. Each worker writes only its own `errors[worker_id]` before the `atomic_add!` on `done_count`, which release-publishes the slot to the coordinator. Instances are cached in `_spin_barrier_pools::Dict{Tuple{UInt, Symbol}, _SpinBarrierPool}` under `_spin_barrier_lock`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `workers` | Int | n/a | yes | Field `workers`. |
| in | `worker_gen` | Vector{Threads.Atomic{Int}} | n/a | yes | Field `worker_gen`. |
| in | `done_count` | Threads.Atomic{Int} | n/a | yes | Field `done_count`. |
| in | `stop` | Threads.Atomic{Bool} | n/a | yes | Field `stop`. |
| in | `request` | Base.RefValue{Any} | n/a | yes | Field `request`. |
| in | `errors` | Vector{Union{Nothing, Base.CapturedException}} | n/a | yes | Field `errors`. |
| in | `run_lock` | ReentrantLock | n/a | yes | Field `run_lock`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | _SpinBarrierPool | n/a | — | Constructed `_SpinBarrierPool`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/policy/types.jl`
- [[parallel.context__create_spin_barrier_pool|_create_spin_barrier_pool]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:201-201`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Idle workers spin at 100 percent CPU and only call `GC.safepoint()`, never `yield()`, so any other Julia task scheduled on those threads is starved; the pool therefore reserves at least one thread for the coordinator or it deadlocks. `request::RefValue{Any}` is untyped, so every field access in the worker loop is dynamically dispatched. `worker_gen` uses a monotonically increasing `Int` that would wrap after 2^63 dispatches. There is no timeout on the coordinator's `done_count` spin, so a worker that never returns hangs the dispatcher.

## Provenance
Mapped from `src/parallel/policy/types.jl` line 127.
