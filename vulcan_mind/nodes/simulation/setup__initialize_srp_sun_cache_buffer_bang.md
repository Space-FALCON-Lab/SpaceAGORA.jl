---
id: simulation.setup__initialize_srp_sun_cache_buffer_bang
label: _initialize_srp_sun_cache_buffer!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_srp_sun_cache_buffer!
  lines:
  - 1402
  - 1402
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
  description: Return value of `_initialize_srp_sun_cache_buffer!`; mutates `p` in
    place. Returns `nothing`.
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

# _initialize_srp_sun_cache_buffer!

## Purpose
Clears the SRP Sun-position cache slot at run start so a table built for a different epoch cannot feed the solar radiation pressure effector.

## Design & Implementation
Sets `shared_buffers.srp_sun_ephemeris_cache[] = nothing` and returns `nothing`. It precedes `_initialize_srp_sun_ephemeris_cache!`, which rebuilds or reuses a table when SRP caching is enabled and an SRP effector is active.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_srp_sun_cache_buffer!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:199-199`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
When the builder declines — no SRP effector, caching disabled, or too many samples — the slot remains `nothing` and every SRP evaluation performs a live `spkpos` under the SPICE lock.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1402.
