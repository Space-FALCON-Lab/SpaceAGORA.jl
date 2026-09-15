---
id: core.effector_sampling_atmospheresample
label: AtmosphereSample
kind: struct
source:
  file: src/core/types/effector_sampling.jl
  symbol: AtmosphereSample
  lines:
  - 76
  - 76
inputs:
- id: rho_kg_m3
  type: Float64
  units: n/a
  required: true
  description: Field `rho_kg_m3`.
- id: temperature_k
  type: Float64
  units: n/a
  required: true
  description: Field `temperature_k`.
- id: wind_pp
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Field `wind_pp`.
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
  type: AtmosphereSample
  units: n/a
  description: Constructed `AtmosphereSample`.
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

# AtmosphereSample

## Purpose
Stage-consistent atmospheric state at the spacecraft location, provided to effectors that request `atmosphere = true`, such as aerodynamic drag and lift, heating and structural-load models.

## Design & Implementation
Immutable struct with `rho_kg_m3::Float64` (density), `temperature_k::Float64` and `wind_pp::SVector{3,Float64}` (wind velocity in the planet-fixed frame, m/s). The engine fills it from whichever density model is active (constant, exponential, none or GRAM) so effectors are agnostic of the atmosphere source. Combining `wind_pp` with `PlanetFrameSample.vel_pp` yields the relative airspeed used for dynamic pressure.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `rho_kg_m3` | Float64 | n/a | yes | Field `rho_kg_m3`. |
| in | `temperature_k` | Float64 | n/a | yes | Field `temperature_k`. |
| in | `wind_pp` | SVector{3, Float64} | n/a | yes | Field `wind_pp`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AtmosphereSample | n/a | — | Constructed `AtmosphereSample`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/effector_sampling.jl`
- [[simulation.effector_sampling__sample_atmosphere_from_planet_frame|_sample_atmosphere_from_planet_frame]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:82-82`
- [[simulation.effector_sampling_sample_buffered_atmosphere|sample_buffered_atmosphere]] · `callees` → `callers` · call · `src/simulation/engine/effector_sampling.jl:148-148`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only bulk density, temperature and wind are carried; species composition, mean molecular mass and pressure are not included, so a free-molecular model needing gas constants must obtain them from the planet model. `rho_kg_m3 = 0` above the interface is not distinguished from a missing sample except by the enclosing `EnvironmentSample` field being `nothing`. No NaN or negativity checks are applied at construction.

## Provenance
Mapped from `src/core/types/effector_sampling.jl` line 76.
