---
id: analysis.telemetry_loading__planet_fixed_to_j2000_state
label: _planet_fixed_to_j2000_state
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _planet_fixed_to_j2000_state
  lines:
  - 41
  - 41
inputs:
- id: planet_name
  type: String
  units: n/a
  required: true
  description: Positional argument `planet_name`.
- id: et
  type: Float64
  units: n/a
  required: true
  description: Positional argument `et`.
- id: r_m
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `r_m`.
- id: v_mps
  type: SVector{3, Float64}
  units: n/a
  required: true
  description: Positional argument `v_mps`.
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
  type: Any
  units: n/a
  description: Return value of `_planet_fixed_to_j2000_state`. Returns `_transform_state(from_frame,
    "J2000", et, r_m, v_mps)` or `_transform_state(fallback, "J2000", et, r_m, v_mps)`.
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

# _planet_fixed_to_j2000_state

## Purpose
Rotates a position and velocity given in the planet's body-fixed frame into the J2000 inertial frame at ephemeris time `et`, transparently retrying with the lower-precision Earth frame when the high-precision kernel is unavailable.

## Design & Implementation
Marked `@inline`. Obtains `from_frame` via `_planet_fixed_frame_name(planet_name)` and `fallback` via `_planet_fixed_frame_fallback_name`. Calls `_transform_state(from_frame, "J2000", et, r_m, v_mps)`, which builds a 6x6 `SMatrix` from `sxform` and multiplies the stacked `[r; v]` state so the velocity transformation includes the frame-rotation-rate term. On any exception it rethrows if `fallback === nothing`, otherwise repeats with the fallback frame. Inputs and outputs are `SVector{3,Float64}` in metres and metres per second.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `et` | Float64 | n/a | yes | Positional argument `et`. |
| in | `r_m` | SVector{3, Float64} | n/a | yes | Positional argument `r_m`. |
| in | `v_mps` | SVector{3, Float64} | n/a | yes | Positional argument `v_mps`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_planet_fixed_to_j2000_state`. Returns `_transform_state(from_frame, "J2000", et, r_m, v_mps)` or `_transform_state(fallback, "J2000", et, r_m, v_mps)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__initial_condition_from_time_aligned_telemetry|_initial_condition_from_time_aligned_telemetry]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:106-106`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- `callees` → [[analysis.telemetry_loading__planet_fixed_frame_fallback_name|_planet_fixed_frame_fallback_name]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:43-43`
- `callees` → [[analysis.telemetry_loading__planet_fixed_frame_name|_planet_fixed_frame_name]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:42-42`
- `callees` → [[envana.ana_telemetry_loading_transform_state|_transform_state]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:45-45`
<!-- vulcan:connections:end -->

## Limitations
The bare `catch` swallows every error type, including kernel-loading problems unrelated to the frame name and even `InterruptException`, before deciding whether to retry. Each call performs a fresh `sxform` lookup, so converting a long telemetry series costs one SPICE call per sample. The transformation is applied only to the state; accelerations or covariances are not rotated. Time `et` is assumed to be TDB seconds consistent with `_initial_time_et`.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 41.
