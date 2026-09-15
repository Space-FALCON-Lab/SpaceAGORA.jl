---
id: simulation.setup__clear_ephemeris_reuse_cache_bang
label: _clear_ephemeris_reuse_cache!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _clear_ephemeris_reuse_cache!
  lines:
  - 337
  - 337
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
outputs:
- id: result
  type: Nothing
  units: n/a
  description: Return value of `_clear_ephemeris_reuse_cache!`. Returns `nothing`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- simulation
charts:
- simulation
origin: agent
---

# _clear_ephemeris_reuse_cache!

## Purpose
Empties all four process-global ephemeris reuse dictionaries, releasing memory between unrelated scenario batches or in tests that need a clean slate.

## Design & Implementation
Under `_EPHEMERIS_REUSE_LOCK` it calls `empty!` on `_SRP_EPHEMERIS_REUSE_CACHE`, `_NBODY_EPHEMERIS_REUSE_CACHE`, `_PLANET_FRAME_EPHEMERIS_REUSE_CACHE`, and `_NBODY_EPHEMERIS_PREWARMED_CACHE` in that order, then returns `nothing`. Takes no arguments. Because all four are cleared atomically under one lock, a concurrent lookup sees either the full old state or the empty state.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_clear_ephemeris_reuse_cache!`. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Clearing does not free caches still referenced by live `SimulationParams` objects, so memory drops only after those are garbage collected. It also discards prewarmed caches the operator registered explicitly, with no option to clear only the automatic ones.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 337.
