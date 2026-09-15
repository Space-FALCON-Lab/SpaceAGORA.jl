---
id: simulation.setup__planet_frame_cache_enabled
label: _planet_frame_cache_enabled
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _planet_frame_cache_enabled
  lines:
  - 216
  - 216
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
  description: Return value of `_planet_frame_cache_enabled`.
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

# _planet_frame_cache_enabled

## Purpose
Switch for the `PlanetFrameEphemerisCache`, which tabulates the primary body's inertial-to-body-fixed rotation so gravity-harmonics and atmosphere lookups avoid recomputing planet orientation each step.

## Design & Implementation
Returns `SimulationModel.ParallelPolicy.parse_bool_env("SPACEAGORA_PLANET_FRAME_CACHE", true)`. Enabled by default; unrecognised values throw `ArgumentError`. `_initialize_planet_frame_ephemeris_cache!` reads it once and skips all sampling when false.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_planet_frame_cache_enabled`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.simulation|RuntimeServices]] · `api` → `module_api` · call · `src/simulation/engine/setup.jl`
- [[simulation.setup__initialize_planet_frame_ephemeris_cache_bang|_initialize_planet_frame_ephemeris_cache!]] · `callees` → `callers` · call · `src/simulation/engine/setup.jl:1856-1856`

**Downstream**

- `callees` → [[parallel.env_config_parse_bool_env|parse_bool_env]] · `callers` · call · `src/simulation/engine/setup.jl:217-217`
<!-- vulcan:connections:end -->

## Limitations
Bypasses the engine override layer. Disabling the cache is transparent to correctness but can substantially slow high-degree harmonics runs; nothing warns about that trade.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 216.
