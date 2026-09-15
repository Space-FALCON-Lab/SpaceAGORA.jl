---
id: simulation.effector_sampling_sample_solar_ephemeris
label: sample_solar_ephemeris
kind: function
source:
  file: src/simulation/engine/effector_sampling.jl
  symbol: sample_solar_ephemeris
  lines:
  - 162
  - 162
inputs:
- id: x
  type: Any
  units: n/a
  required: true
  description: Positional argument `x`.
- id: p
  type: Any
  units: n/a
  required: true
  description: Positional argument `p`.
- id: sat_idx
  type: Int
  units: n/a
  required: true
  description: Positional argument `sat_idx`.
- id: t
  type: Float64
  units: n/a
  required: true
  description: Positional argument `t`.
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
  type: SolarEphemerisSample
  units: n/a
  description: Return value of `sample_solar_ephemeris`.
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

# sample_solar_ephemeris

## Purpose
Obtains the Sun's position relative to the primary body in J2000 at time `t`, for solar radiation pressure and shadow effectors.

## Design & Implementation
Forms ephemeris time from the run's `et_start` plus `t`, resolves the primary's SPICE name, and consults the `srp_sun_ephemeris_cache` if one is installed; a cache miss or absent cache falls through to `_srp_sun_position_from_spice_j2000_m`, which honours the RHS memo flag and increments the `srp_spkpos_runtime_calls` counter. Returns a `SolarEphemerisSample`. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `x` | Any | n/a | yes | Positional argument `x`. |
| in | `p` | Any | n/a | yes | Positional argument `p`. |
| in | `sat_idx` | Int | n/a | yes | Positional argument `sat_idx`. |
| in | `t` | Float64 | n/a | yes | Positional argument `t`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SolarEphemerisSample | n/a | — | Return value of `sample_solar_ephemeris`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[simulation.dynamics_rhs__prefill_shared_body_samples_bang|_prefill_shared_body_samples!]] · `callees` → `callers` · call · `src/simulation/engine/dynamics_rhs.jl:1222-1222`
- [[simulation.effector_sampling__sample_reusable_solar|_sample_reusable_solar]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:279-279`
- [[simx.engine_effector_sampling_sample_environment|sample_environment]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:247-247`

**Downstream**

- `callees` → [[core.effector_sampling_solarephemerissample|SolarEphemerisSample]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:190-190`
- `callees` → [[dynamics.perturbations__spice_query_name|_spice_query_name]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:165-165`
- `callees` → [[dynamics.perturbations__srp_sun_position_from_cache_j2000_m|_srp_sun_position_from_cache_j2000_m]] · `callers` · call · `src/simulation/engine/effector_sampling.jl:171-171`
<!-- vulcan:connections:end -->

## Limitations
The cache-or-SPICE branch is duplicated inline rather than factored, so the two fallback calls must be kept identical by hand; and a cache that covers only part of the mission silently mixes cached and live positions with different interpolation error.

## Provenance
Mapped from `src/simulation/engine/effector_sampling.jl` line 162.
