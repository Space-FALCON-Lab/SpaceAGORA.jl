---
id: simulation.setup__initialize_nbody_ephemeris_cache_buffer_bang
label: _initialize_nbody_ephemeris_cache_buffer!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_nbody_ephemeris_cache_buffer!
  lines:
  - 1407
  - 1407
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
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
  description: Return value of `_initialize_nbody_ephemeris_cache_buffer!`; mutates
    `p` in place. Returns `nothing`.
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

# _initialize_nbody_ephemeris_cache_buffer!

## Purpose
Clears the N-body ephemeris cache slot at run start so a previous run's table cannot be used with a different epoch.

## Design & Implementation
Sets `shared_buffers.nbody_ephemeris_cache[] = nothing`. Returns `nothing`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_nbody_ephemeris_cache_buffer!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:198-198`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Run before the builder; if the builder is skipped the slot stays `nothing` and effectors fall back to live queries.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1407.
