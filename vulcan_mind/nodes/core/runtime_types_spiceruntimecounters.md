---
id: core.runtime_types_spiceruntimecounters
label: SpiceRuntimeCounters
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: SpiceRuntimeCounters
  lines:
  - 486
  - 486
inputs:
- id: nbody_spkpos_runtime_calls
  type: Base.Threads.Atomic{Int64}
  units: n/a
  required: false
  description: Field `nbody_spkpos_runtime_calls` (default `Base.Threads.Atomic{Int64}(0)`).
- id: nbody_spkpos_cache_build_calls
  type: Base.Threads.Atomic{Int64}
  units: n/a
  required: false
  description: Field `nbody_spkpos_cache_build_calls` (default `Base.Threads.Atomic{Int64}(0)`).
- id: srp_spkpos_runtime_calls
  type: Base.Threads.Atomic{Int64}
  units: n/a
  required: false
  description: Field `srp_spkpos_runtime_calls` (default `Base.Threads.Atomic{Int64}(0)`).
- id: srp_spkpos_cache_build_calls
  type: Base.Threads.Atomic{Int64}
  units: n/a
  required: false
  description: Field `srp_spkpos_cache_build_calls` (default `Base.Threads.Atomic{Int64}(0)`).
- id: planet_pxform_runtime_calls
  type: Base.Threads.Atomic{Int64}
  units: n/a
  required: false
  description: Field `planet_pxform_runtime_calls` (default `Base.Threads.Atomic{Int64}(0)`).
- id: planet_pxform_cache_build_calls
  type: Base.Threads.Atomic{Int64}
  units: n/a
  required: false
  description: Field `planet_pxform_cache_build_calls` (default `Base.Threads.Atomic{Int64}(0)`).
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
  type: SpiceRuntimeCounters
  units: n/a
  description: Constructed `SpiceRuntimeCounters` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# SpiceRuntimeCounters

## Purpose
Thread-safe counters recording how many SPICE position and orientation calls a run made at runtime versus during cache construction, for verifying that the caches are actually being hit.

## Design & Implementation
A `@kwdef struct` of six `Threads.Atomic{Int64}` fields — runtime and cache-build counts for N-body `spkpos`, SRP `spkpos` and planet `pxform` — each defaulting to zero. The atomics let threaded RHS evaluations increment without a lock.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `nbody_spkpos_runtime_calls` | Base.Threads.Atomic{Int64} | n/a | no | Field `nbody_spkpos_runtime_calls` (default `Base.Threads.Atomic{Int64}(0)`). |
| in | `nbody_spkpos_cache_build_calls` | Base.Threads.Atomic{Int64} | n/a | no | Field `nbody_spkpos_cache_build_calls` (default `Base.Threads.Atomic{Int64}(0)`). |
| in | `srp_spkpos_runtime_calls` | Base.Threads.Atomic{Int64} | n/a | no | Field `srp_spkpos_runtime_calls` (default `Base.Threads.Atomic{Int64}(0)`). |
| in | `srp_spkpos_cache_build_calls` | Base.Threads.Atomic{Int64} | n/a | no | Field `srp_spkpos_cache_build_calls` (default `Base.Threads.Atomic{Int64}(0)`). |
| in | `planet_pxform_runtime_calls` | Base.Threads.Atomic{Int64} | n/a | no | Field `planet_pxform_runtime_calls` (default `Base.Threads.Atomic{Int64}(0)`). |
| in | `planet_pxform_cache_build_calls` | Base.Threads.Atomic{Int64} | n/a | no | Field `planet_pxform_cache_build_calls` (default `Base.Threads.Atomic{Int64}(0)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SpiceRuntimeCounters | n/a | — | Constructed `SpiceRuntimeCounters` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[core.runtime_types_sharedbuffers|SharedBuffers]] · `callees` → `callers` · call · `src/core/types/runtime_types.jl:737-737`
- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Counters are per `SharedBuffers` instance, so a campaign of many runs must aggregate them itself; they count calls, not time, so a few expensive calls and many cheap ones look identical.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 486.
