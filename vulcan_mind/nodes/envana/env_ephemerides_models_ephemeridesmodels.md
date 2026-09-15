---
id: envana.env_ephemerides_models_ephemeridesmodels
label: EphemeridesModels
kind: struct
source:
  file: src/environment/ephemerides/ephemerides_models.jl
  symbol: EphemeridesModels
  lines:
  - 1
  - 16
inputs:
- id: module_api
  type: Module
  units: n/a
  required: true
  description: Parent environment namespace that includes the ephemerides module file.
- id: callers
  type: call
  units: n/a
  required: false
  description: Invocations of this symbol observed in mapped callers.
outputs:
- id: ephemerides_api
  type: Module
  units: n/a
  description: Exported ephemerides model types and frame/time query functions.
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- environment
charts:
- envana
origin: agent
---
# EphemeridesModels

## Purpose
`EphemeridesModels` is the namespace that owns planetary time conversion and body-fixed frame rotation for the simulator. It exports `SpiceEphemeridesModel` and `SimpleEphemeridesModel` together with the query functions `ephemerides_time_seconds`, `planet_frame_lpi`, `ephemerides_requires_spice`, and `ephemerides_cache_key`.

## Theory & Math
The two exported model types answer the same mathematical question by different means: given an epoch, produce the rotation matrix `R` from the local planet-inertial frame to the body-fixed frame. SPICE-backed evaluation obtains `R` from kernel-defined orientation; the simple model composes a rotation about the spin axis through the angle `theta = pm + omega_3 * (et - t_ref)`, with `pm` the prime-meridian angle at the reference epoch in radians, `omega_3` the spin rate in rad/s, and `et - t_ref` the elapsed ephemeris time in seconds. Ephemeris time itself is the TDB seconds past J2000, so `et = 0` corresponds to 2000-01-01T12:00:00 TDB.

## Model & Assumptions
The module assumes `SPICE_LOCK` from `RuntimeServices` serialises all native SPICE calls, because the CSPICE library is not thread-safe. `AbstractEphemeridesModel` from `AbstractTypes` is the dispatch root, so every concrete model must supply methods for the four exported query functions. Time is carried in Float64 seconds, which limits epoch resolution to roughly microseconds over a century-wide span.

## Design & Implementation
The file is a pure namespace: it pulls in `Dates`, `AstroTime`, `StaticArrays`, and `SPICE`, imports the lock and the abstract supertype, declares two `export` lines, and then `include`s `simple_ephemerides.jl`, which supplies the concrete struct definitions and their methods. Placing the include after the exports means the exported names are resolved once the included file defines them.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `module_api` | Module | n/a | yes | Parent environment namespace that includes the ephemerides module file. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `ephemerides_api` | Module | n/a | — | Exported ephemerides model types and frame/time query functions. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.environment|EnvironmentModels]] · `api` → `module_api` · call · `src/environment/ephemerides/ephemerides_models.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Only one implementation file is included, so any additional ephemerides backend must be added here explicitly. The module does not validate that SPICE kernels are furnished; `ephemerides_requires_spice` merely reports the requirement and leaves loading to the caller. Nothing in the module caches rotation matrices beyond the key produced by `ephemerides_cache_key`.

## Provenance
Read directly from `src/environment/ephemerides/ephemerides_models.jl:1-16`, including its export list and its single `include` of `simple_ephemerides.jl`.
