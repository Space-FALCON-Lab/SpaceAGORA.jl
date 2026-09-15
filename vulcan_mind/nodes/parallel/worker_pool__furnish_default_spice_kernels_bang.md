---
id: parallel.worker_pool__furnish_default_spice_kernels_bang
label: _furnish_default_spice_kernels!
kind: function
source:
  file: src/parallel/process/worker_pool.jl
  symbol: _furnish_default_spice_kernels!
  lines:
  - 102
  - 102
inputs:
- id: worker
  type: Int
  units: n/a
  required: true
  description: Positional argument `worker`.
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
  description: Return value of `_furnish_default_spice_kernels!`; mutates `worker`
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

# _furnish_default_spice_kernels!

## Purpose
Replays on a worker the SPICE kernel furnishing side effect that happened on the coordinator, which deserialization alone cannot reproduce.

## Design & Implementation
Evaluates `SpaceAGORA.SimulationModel.Earth("")` in the worker's `Main` through `remotecall_eval`. Constructing a default Earth furnishes the shared leapseconds and DE44x planetary ephemeris set as a side effect. `remotecall_eval` is used instead of a closure naming `SpaceAGORA` because this function lives in the `ParallelProcess` submodule, which does not see its parent's name — `ParallelProcess` is included before SpaceAGORA finishes defining its other submodules, so `..SpaceAGORA` is not resolvable at include time either. The call is wrapped in `try` and warns on failure, then calls `_warm_gram_wrapper!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `worker` | Int | n/a | yes | Positional argument `worker`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_furnish_default_spice_kernels!`; mutates `worker` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/process/worker_pool.jl`
- [[parallel.worker_pool__bootstrap_process_worker_bang|_bootstrap_process_worker!]] · `callees` → `callers` · call · `src/parallel/process/worker_pool.jl:83-83`

**Downstream**

- `callees` → [[environment.planets_earth|Earth]] · `callers` · call · `src/parallel/process/worker_pool.jl:114-114`
- `callees` → [[parallel.worker_pool__warm_gram_wrapper_bang|_warm_gram_wrapper!]] · `callers` · call · `src/parallel/process/worker_pool.jl:118-118`
<!-- vulcan:connections:end -->

## Limitations
Only the default kernel set is furnished. A campaign on a non-default kernel directory, or a non-Earth-primary mission whose kernels are outside that shared set, must furnish its own on the pool's workers before dispatching; this covers the common case and does not attempt to be general.

## Provenance
Mapped from `src/parallel/process/worker_pool.jl` line 102.
