---
id: simulation.setup__srp_ephemeris_cache_enabled
label: _srp_ephemeris_cache_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _srp_ephemeris_cache_enabled
  lines:
  - 192
  - 192
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
  type: Bool
  units: n/a
  description: Return value of `_srp_ephemeris_cache_enabled`.
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

# _srp_ephemeris_cache_enabled

## Purpose
Master switch for pre-sampling the Sun position relative to the primary body into an `SRPSunEphemerisCache` so solar-radiation-pressure evaluations interpolate instead of calling SPICE on every RHS step.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_SRP_EPHEMERIS_CACHE", true)`, so the cache is on by default and accepts the `1/0`, `true/false`, `yes/no`, `on/off` spellings. Unrecognised text raises `ArgumentError` from the parser. Consulted by `_initialize_srp_sun_ephemeris_cache!` before any sampling work is done.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_srp_ephemeris_cache_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_srp_sun_ephemeris_cache_bang|_initialize_srp_sun_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1752-1752`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/setup.jl:193-193`
<!-- vulcan:connections:end -->

## Limitations
Reads the raw process `ENV` via `ParallelPolicy`, not the engine override layer used by `_engine_env_get`, so engine-level overrides of this variable are ignored. When disabled, SRP falls back to direct ephemeris queries with no warning about the performance cost.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 192.
