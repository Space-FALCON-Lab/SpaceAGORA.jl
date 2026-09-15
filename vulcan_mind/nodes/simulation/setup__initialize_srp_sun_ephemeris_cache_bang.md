---
id: simulation.setup__initialize_srp_sun_ephemeris_cache_bang
label: _initialize_srp_sun_ephemeris_cache!
kind: function
source:
  file: src/simulation/engine/setup.jl
  symbol: _initialize_srp_sun_ephemeris_cache!
  lines:
  - 1751
  - 1751
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
  description: Return value of `_initialize_srp_sun_ephemeris_cache!`; mutates `p`
    in place. Returns `nothing`.
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

# _initialize_srp_sun_ephemeris_cache!

## Purpose
Tabulates the Sun's position relative to the primary over the mission for the solar radiation pressure effector.

## Design & Implementation
Returns early unless SRP caching is enabled, an SRP effector is active and the mission end is positive; computes the sample count and disables with a warning above the maximum. Checks the reuse cache keyed by primary, epoch, end and step. Otherwise it samples `_spice_position_j2000_m_unlocked` for the Sun at each time inside a single `SPICE_LOCK` acquisition, incrementing the build counter per sample, stores in the reuse cache if enabled, installs the result and prints a summary when verbose.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `et_start` | Float64 | n/a | yes | Positional argument `et_start`. |
| in | `mission_end_s` | Float64 | n/a | yes | Positional argument `mission_end_s`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Nothing | n/a | — | Return value of `_initialize_srp_sun_ephemeris_cache!`; mutates `p` in place. Returns `nothing`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simx.engine_execution_run_simulation|run_simulation]] · `callees` → `callers` · call · `src/simulation/engine/execution.jl:227-227`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/simulation/engine/setup.jl:1794-1794`
- `callees` → [[core.runtime_types_srpsunephemeriscache|SRPSunEphemerisCache]] · `callers` · call · `src/simulation/engine/setup.jl:1787-1787`
- `callees` → [[dynamics.perturbations__spice_query_name|_spice_query_name]] · `callers` · call · `src/simulation/engine/setup.jl:1766-1766`
- `callees` → [[environment.simple_ephemerides__spice_position_j2000_m_unlocked|_spice_position_j2000_m_unlocked]] · `callers` · call · `src/simulation/engine/setup.jl:1782-1782`
- `callees` → [[simulation.setup__ephemeris_reuse_enabled|_ephemeris_reuse_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:1767-1767`
- `callees` → [[simulation.setup__ephemeris_reuse_lookup|_ephemeris_reuse_lookup]] · `callers` · call · `src/simulation/engine/setup.jl:1769-1769`
- `callees` → [[simulation.setup__ephemeris_reuse_max_entries|_ephemeris_reuse_max_entries]] · `callers` · call · `src/simulation/engine/setup.jl:1790-1790`
- `callees` → [[simulation.setup__ephemeris_reuse_store_bang|_ephemeris_reuse_store!]] · `callers` · call · `src/simulation/engine/setup.jl:1790-1790`
- `callees` → [[simulation.setup__has_active_srp_effector|_has_active_srp_effector]] · `callers` · call · `src/simulation/engine/setup.jl:1753-1753`
- `callees` → [[simulation.setup__srp_ephemeris_cache_dt_s|_srp_ephemeris_cache_dt_s]] · `callers` · call · `src/simulation/engine/setup.jl:1758-1758`
- `callees` → [[simulation.setup__srp_ephemeris_cache_enabled|_srp_ephemeris_cache_enabled]] · `callers` · call · `src/simulation/engine/setup.jl:1752-1752`
- `callees` → [[simulation.setup__srp_ephemeris_cache_max_samples|_srp_ephemeris_cache_max_samples]] · `callers` · call · `src/simulation/engine/setup.jl:1760-1760`
- `callees` → [[simulation.setup__srp_ephemeris_reuse_key|_srp_ephemeris_reuse_key]] · `callers` · call · `src/simulation/engine/setup.jl:1768-1768`
<!-- vulcan:connections:end -->

## Limitations
Holding the SPICE lock for the whole build blocks every other SPICE user in the process for the duration.

## Provenance
Mapped from `src/simulation/engine/setup.jl` line 1751.
