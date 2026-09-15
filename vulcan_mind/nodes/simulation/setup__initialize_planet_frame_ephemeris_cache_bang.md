---
id: simulation.setup__initialize_planet_frame_ephemeris_cache_bang
label: _initialize_planet_frame_ephemeris_cache!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_planet_frame_ephemeris_cache!
  lines:
  - 1855
  - 1855
inputs:
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: et_start
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et_start`.
- id: mission_end_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mission_end_s`.
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
  description: Return value of `_initialize_planet_frame_ephemeris_cache!`; mutates
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

# _initialize_planet_frame_ephemeris_cache!

## Purpose
Tabulates the planet's orientation quaternion over the mission so `l_pi` can be interpolated rather than recomputed from SPICE at every step.

## Design & Implementation
Returns early if disabled or the mission end is not positive; computes the sample count from `dt_s` and disables with a warning above the maximum. Checks the reuse cache keyed by planet, ephemerides model, epoch, end and step. Otherwise it samples `planet_frame_lpi` at each time, converts to a quaternion with `dcm_to_quaternion`, increments the `planet_pxform_cache_build_calls` counter for SPICE models, stores in the reuse cache if enabled, and installs the result. Prints a summary when verbose.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `et_start` | Float64 | n/a | yes | Positional argument `et_start`. |
| in | `mission_end_s` | Float64 | n/a | yes | Positional argument `mission_end_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_planet_frame_ephemeris_cache!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:228-228`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/engine/setup.jl:1899-1899`
- `callees` → [[core.quaternion_utils_dcm_to_quaternion|dcm_to_quaternion]] · `callers` · call · `src/simulation/engine/setup.jl:1886-1886`
- `callees` → [[core.runtime_types_planetframeephemeriscache|PlanetFrameEphemerisCache]] · `callers` · call · `src/simulation/engine/setup.jl:1892-1892`
- `callees` → [[environment.simple_ephemerides_planet_frame_lpi|planet_frame_lpi]] · `callers` · call · `src/simulation/engine/setup.jl:1885-1885`
- `callees` → [[simulation.setup__ephemeris_reuse_enabled|_ephemeris_reuse_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:1871-1871`
- `callees` → [[simulation.setup__ephemeris_reuse_lookup|_ephemeris_reuse_lookup]] · `callers` · call · `src/simulation/engine/setup.jl:1873-1873`
- `callees` → [[simulation.setup__ephemeris_reuse_max_entries|_ephemeris_reuse_max_entries]] · `callers` · call · `src/simulation/engine/setup.jl:1895-1895`
- `callees` → [[simulation.setup__ephemeris_reuse_store_bang|_ephemeris_reuse_store!]] · `callers` · call · `src/simulation/engine/setup.jl:1895-1895`
- `callees` → [[simulation.setup__planet_frame_cache_dt_s|_planet_frame_cache_dt_s]] · `callers` · call · `src/simulation/engine/setup.jl:1861-1861`
- `callees` → [[simulation.setup__planet_frame_cache_enabled|_planet_frame_cache_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:1856-1856`
- `callees` → [[simulation.setup__planet_frame_cache_max_samples|_planet_frame_cache_max_samples]] · `callers` · call · `src/simulation/engine/setup.jl:1863-1863`
- `callees` → [[simulation.setup__planet_frame_ephemeris_reuse_key|_planet_frame_ephemeris_reuse_key]] · `callers` · call · `src/simulation/engine/setup.jl:1872-1872`
<!-- vulcan:connections:end -->

## Limitations
Samples are taken serially under whatever lock `planet_frame_lpi` requires, so building a fine-grained table for a long mission can take seconds before the solve starts.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1855.
