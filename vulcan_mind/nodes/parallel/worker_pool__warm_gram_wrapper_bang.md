---
id: parallel.worker_pool__warm_gram_wrapper_bang
label: _warm_gram_wrapper!
kind: function
source:
  file: src/parallel/process/worker_pool.jl
  symbol: _warm_gram_wrapper!
  lines:
  - 150
  - 150
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
  description: Return value of `_warm_gram_wrapper!`; mutates `worker` in place.
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

# _warm_gram_wrapper!

## Purpose
Forces GRAM's dynamically loaded native wrapper module to settle on a worker before any real campaign work reaches it, eliminating a world-age `MethodError` on the first GRAM density sample.

## Design & Implementation
Evaluates a quoted block in the worker that constructs a `GRAMAtmosphereModel` for Earth and then calls `GRAMSuite.point_density_state` on its `core` at 150 km altitude, zero latitude, zero longitude, zero elapsed time and wind disabled. Both halves matter. GRAMSuite `Base.include`s its wrapper on first model construction in a process, defining new types at a later world age; construction alone was measured to be insufficient, because the ephemeris-bypass hook that calls `GRAMSuite.GRAM.EphemerisStateC` is reached only once a density query actually happens. Doing this in its own remote call gives the definitions a fresh world age before any campaign closure runs. Failures warn and are swallowed.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `worker` | Int | n/a | yes | Positional argument `worker`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_warm_gram_wrapper!`; mutates `worker` in place. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.parallel|ParallelProfiles]] · `api` → `module_api` · call · `src/parallel/process/worker_pool.jl`
- [[parallel.worker_pool__furnish_default_spice_kernels_bang|_furnish_default_spice_kernels!]] · `callees` → `callers` · call · `src/parallel/process/worker_pool.jl:118-118`

**Downstream**

- `callees` → [[environment.density_models_gramatmospheremodel|GRAMAtmosphereModel]] · `callers` · call · `src/parallel/process/worker_pool.jl:153-153`
- `callees` → [[gnc.heat_rate_control_f|f]] · `callers` · call · `src/parallel/process/worker_pool.jl:186-186`
- `callees` → [[gnc.struct_load_control_f|f]] · `callers` · call · `src/parallel/process/worker_pool.jl:186-186`
- `callees` → [[gnc.tracking_executor_f|f]] · `callers` · call · `src/parallel/process/worker_pool.jl:186-186`
- `callees` → [[parallel.ensure_process_workers_ensure_process_workers_bang|ensure_process_workers!]] · `callers` · feedback · `src/parallel/process/worker_pool.jl:166-166`
<!-- vulcan:connections:end -->

## Limitations
The warmup is hard-coded to an Earth model at one sample point, so a code path specialised for another planet or a surrogate model can still hit its own first-call cost; being best-effort, a failed warmup leaves the worker in the pool and the original world-age error surfaces inside the caller's first real dispatch.

## Provenance
Mapped from `src/parallel/process/worker_pool.jl` line 150.
