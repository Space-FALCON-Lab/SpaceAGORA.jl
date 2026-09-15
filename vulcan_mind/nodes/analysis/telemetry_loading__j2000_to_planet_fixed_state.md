---
id: analysis.telemetry_loading__j2000_to_planet_fixed_state
label: _j2000_to_planet_fixed_state
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _j2000_to_planet_fixed_state
  lines:
  - 52
  - 52
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
  description: Return value of `_j2000_to_planet_fixed_state`. Returns `_transform_state("J2000",
    to_frame, et, r_m, v_mps)` or `_transform_state("J2000", fallback, et, r_m, v_mps)`.
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

# _j2000_to_planet_fixed_state

## Purpose
Inverse of `_planet_fixed_to_j2000_state`: rotates an inertial J2000 position and velocity into the planet's body-fixed frame at ephemeris time `et`, using the same primary/fallback frame selection so Earth telemetry works with or without the `ITRF93` kernels.

## Design & Implementation
Marked `@inline`. Resolves `to_frame = _planet_fixed_frame_name(planet_name)` and `fallback = _planet_fixed_frame_fallback_name(planet_name)`, then calls `_transform_state("J2000", to_frame, et, r_m, v_mps)`. The 6x6 `sxform` matrix maps both position and velocity, adding the `-ω x r` contribution automatically. Failures rethrow when no fallback exists; otherwise the transform is repeated with the fallback frame. Returns a tuple `(r_pp, v_pp)` of `SVector{3,Float64}`.

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
| out | `result` | Any | n/a | — | Return value of `_j2000_to_planet_fixed_state`. Returns `_transform_state("J2000", to_frame, et, r_m, v_mps)` or `_transform_state("J2000", fallback, et, r_m, v_mps)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.error_tables__time_aligned_rows_errors|_time_aligned_rows_errors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:150-150`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- `callees` → [[analysis.telemetry_loading__planet_fixed_frame_fallback_name|_planet_fixed_frame_fallback_name]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:54-54`
- `callees` → [[analysis.telemetry_loading__planet_fixed_frame_name|_planet_fixed_frame_name]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:53-53`
- `callees` → [[envana.ana_telemetry_loading_transform_state|_transform_state]] · `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:56-56`
<!-- vulcan:connections:end -->

## Limitations
Shares the broad `catch` of its sibling, so a transient failure is masked by a silent precision downgrade for Earth. There is no caching of the rotation matrix across samples at the same `et`. The function does not check that `planet_name` matches the planet model used elsewhere in the verification pipeline, so a mismatch between configuration and telemetry provenance goes unnoticed.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 52.
