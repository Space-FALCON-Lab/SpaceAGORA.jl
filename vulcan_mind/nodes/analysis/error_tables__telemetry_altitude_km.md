---
id: analysis.error_tables__telemetry_altitude_km
label: _telemetry_altitude_km
kind: function
source:
  file: src/analysis/verification/telemetry_verification/error_tables.jl
  symbol: _telemetry_altitude_km
  lines:
  - 61
  - 61
inputs:
- id: pos_m
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `pos_m`.
- id: vel_mps
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `vel_mps`.
- id: planet
  type: Any
  units: n/a
  required: true
  description: Positional argument `planet`.
- id: altitude_mode
  type: Symbol
  units: n/a
  required: true
  description: Positional argument `altitude_mode`.
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
  type: Float64
  units: n/a
  description: Return value of `_telemetry_altitude_km`.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- analysis
charts:
- analysis
origin: agent
---

# _telemetry_altitude_km

## Purpose
Converts an inertial position and velocity pair into the scalar altitude in kilometres that telemetry comparison scores against, under whichever altitude convention the scenario declared.

## Design & Implementation
Branches on `altitude_mode`. Under `:vacuum` it takes the norm of `pos_m` minus the planet's equatorial radius `Rp_e` and scales by 1e-3. Under `:oblate` it first rotates the state into the planet-fixed frame with `r_intor_p!`, then reads the geodetic altitude out of the first element returned by `rtolatlong`. Any other symbol raises `ArgumentError` naming the offending mode, so a typo in a scenario file fails loudly rather than silently scoring against the wrong convention. Declared `@inline` with a `::Float64` return because it runs once per telemetry sample.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `pos_m` | SVector{3, Float64} | n/a | yes | Positional argument `pos_m`. |
| in | `vel_mps` | SVector{3, Float64} | n/a | yes | Positional argument `vel_mps`. |
| in | `planet` | Any | n/a | yes | Positional argument `planet`. |
| in | `altitude_mode` | Symbol | n/a | yes | Positional argument `altitude_mode`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_telemetry_altitude_km`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:170-170`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/error_tables.jl`

**Downstream**

- `callees` → [[core.reference_system_r_intor_p_bang|r_intor_p!]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:70-70`
- `callees` → [[core.reference_system_rtolatlong|rtolatlong]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:71-71`
<!-- vulcan:connections:end -->

## Limitations
The vacuum branch treats the planet as a sphere of radius `Rp_e`, so at high latitude it disagrees with the oblate branch by the flattening term; mixing modes between a simulation and its reference telemetry produces a systematic bias this function cannot detect.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/error_tables.jl` line 61.
