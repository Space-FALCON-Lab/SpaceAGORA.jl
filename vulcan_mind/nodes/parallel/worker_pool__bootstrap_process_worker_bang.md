---
id: parallel.worker_pool__bootstrap_process_worker_bang
label: _bootstrap_process_worker!
kind: function
source:
  file: src/parallel/process/worker_pool.jl
  symbol: _bootstrap_process_worker!
  lines:
  - 57
  - 57
inputs:
- id: worker
  type: Int
  units: n/a
  required: true
  description: Positional argument `worker`.
- id: project_path
  type: String
  units: n/a
  required: true
  description: Positional argument `project_path`.
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
  description: Return value of `_bootstrap_process_worker!`; mutates `worker` in place.
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

# _bootstrap_process_worker!

## Purpose
Brings one freshly spawned worker process to the state a campaign dispatch assumes: SpaceAGORA loaded, GRAMSuite importable where available, and default SPICE kernels furnished.

## Design & Implementation
Loads SpaceAGORA through `Distributed.remotecall_eval` rather than a `remotecall_wait` closure, because any closure defined inside this package has `SpaceAGORA.ParallelProcess` as its home module and Distributed cannot deserialize it on a worker that has not loaded SpaceAGORA yet. It then pushes the vendored `data/GRAMSuite.jl` path onto the worker's `LOAD_PATH` when `find_package` finds nothing, and issues `import GRAMSuite` rather than `using`, because both packages export a differently typed `InitialTime` and `using` both would make the unqualified name ambiguous. The GRAMSuite step is wrapped in `try` and downgraded to a warning, since a campaign with no GRAM density model does not need it. Finally it calls `_furnish_default_spice_kernels!`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `worker` | Int | n/a | yes | Positional argument `worker`. |
| in | `project_path` | String | n/a | yes | Positional argument `project_path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_bootstrap_process_worker!`; mutates `worker` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[parallel.ensure_process_workers_ensure_process_workers_bang|ensure_process_workers!]] · `callees` → `callers` · call · `src/parallel/process/worker_pool.jl:200-200`

**Downstream**

- `callees` → [[parallel.worker_pool__furnish_default_spice_kernels_bang|_furnish_default_spice_kernels!]] · `callers` · call · `src/parallel/process/worker_pool.jl:83-83`
<!-- vulcan:connections:end -->

## Limitations
GRAMSuite loading is best-effort, so a worker that failed it is still added to the pool and only fails later when a GRAM density sample is dispatched to it; the warning names the worker but nothing records the degraded state on the pool itself.

## Provenance
Mapped from `src/parallel/process/worker_pool.jl` line 57.
