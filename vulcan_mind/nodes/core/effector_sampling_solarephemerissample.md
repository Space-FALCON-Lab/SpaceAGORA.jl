---
id: core.effector_sampling_solarephemerissample
label: SolarEphemerisSample
kind: struct
source:
  file: src/core/types/effector_sampling.jl
  symbol: SolarEphemerisSample
  lines:
  - 87
  - 87
inputs:
- id: sun_pos_ii
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `sun_pos_ii`.
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
  description: Constructed `SolarEphemerisSample`.
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

# SolarEphemerisSample

## Purpose
Stage-consistent Sun position in the inertial frame for effectors that request `solar = true`, primarily solar radiation pressure and solar-panel thermal or power models.

## Design & Implementation
Single-field immutable struct: `sun_pos_ii::SVector{3,Float64}`, the Sun's position relative to the central body in metres, J2000 axes. The engine evaluates it once per stage (through SPICE or an analytic ephemeris) so multiple effectors share the lookup. Effectors compute the spacecraft-to-Sun vector as `sun_pos_ii - x.pos_ii`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `sun_pos_ii` | SVector{3, Float64} | n/a | yes | Field `sun_pos_ii`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SolarEphemerisSample | n/a | — | Constructed `SolarEphemerisSample`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/effector_sampling.jl`
- [[simulation.effector_sampling__sample_reusable_solar|_sample_reusable_solar]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:278-278`
- [[simulation.effector_sampling_sample_solar_ephemeris|sample_solar_ephemeris]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:190-190`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only position is provided, so effectors needing the Sun's velocity (aberration, Doppler) or the solar flux scaling must compute distance themselves and assume the 1 AU reference value. Eclipse or shadow state is not included and must be derived by the effector from the planet radius. Light-time correction is whatever the engine applied; the struct does not record it.

## Provenance
Mapped from `src/core/types/effector_sampling.jl` line 87.
