---
id: parallel.context__destroy_persistent_foreach_scope_bang
label: _destroy_persistent_foreach_scope!
kind: function
source:
  file: src/parallel/policy/context.jl
  symbol: _destroy_persistent_foreach_scope!
  lines:
  - 294
  - 294
inputs:
- id: scope_id
  type: UInt
  units: n/a
  required: true
  description: Positional argument `scope_id`.
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
  description: Return value of `_destroy_persistent_foreach_scope!`; mutates `scope_id`
    in place.
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

# _destroy_persistent_foreach_scope!

## Purpose
Tears down every persistent channel pool and spin-barrier pool that was created under a given policy scope, removing their dictionary entries and stopping their worker tasks; called from the `finally` of `with_policy_context`.

## Design & Implementation
Collects victims in two phases to keep lock hold times short. Under `_persistent_foreach_lock` it scans both `_persistent_foreach_pools` and `_persistent_foreach_worker_pools` for keys whose first element equals `scope_id`, pushes the pools into `channel_pools`, and `delete!`s the keys. Under `_spin_barrier_lock` it does the same for `_spin_barrier_pools` into `spin_pools`. After releasing both locks it calls `_shutdown_persistent_foreach_pool!` on each channel pool (which posts `:stop` to every request channel) and `_shutdown_spin_barrier_pool!` on each spin pool (which sets `stop` and bumps generations). Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `scope_id` | UInt | n/a | yes | Positional argument `scope_id`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_destroy_persistent_foreach_scope!`; mutates `scope_id` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parcore.context_with_policy_context|with_policy_context]] · `callees` → `callers` · call · `src/parallel/policy/context.jl:339-339`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/parallel/policy/context.jl:302-302`
- `callees` → [[parallel.context__shutdown_persistent_foreach_pool_bang|_shutdown_persistent_foreach_pool!]] · `callers` · call · `src/parallel/policy/context.jl:324-324`
- `callees` → [[parallel.context__shutdown_spin_barrier_pool_bang|_shutdown_spin_barrier_pool!]] · `callers` · call · `src/parallel/policy/context.jl:327-327`
<!-- vulcan:connections:end -->

## Limitations
Shutdown is asynchronous: worker tasks exit some time after this returns, so a pool's threads may still be busy briefly. If a pool is mid-batch when `:stop` is posted, `_shutdown_persistent_foreach_pool!` blocks on `run_lock` until that batch finishes, which is safe but can delay scope exit. A scope whose id collides with a reused `objectid` could remove another live scope's pools, though this requires the earlier context to have been collected without cleanup.

## Provenance
Mapped from `src/parallel/policy/context.jl` line 294.
