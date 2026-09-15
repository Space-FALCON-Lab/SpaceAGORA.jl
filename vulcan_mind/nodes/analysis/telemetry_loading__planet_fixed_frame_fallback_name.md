---
id: analysis.telemetry_loading__planet_fixed_frame_fallback_name
label: _planet_fixed_frame_fallback_name
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _planet_fixed_frame_fallback_name
  lines:
  - 18
  - 18
inputs:
- id: planet_name
  type: String
  units: n/a
  required: true
  description: Positional argument `planet_name`.
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
  type: Union{Nothing,
  units: n/a
  description: Return value of `_planet_fixed_frame_fallback_name`.
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

# _planet_fixed_frame_fallback_name

## Purpose
Supplies a secondary SPICE body-fixed frame to try when the primary frame from `_planet_fixed_frame_name` cannot be evaluated. Only Earth has a fallback (`IAU_EARTH`, which needs just the generic text PCK); all other bodies return `nothing`.

## Design & Implementation
Marked `@inline`; the same `lowercase(strip(planet_name)) == "earth"` test as the primary-frame function, returning `_EARTH_FALLBACK_BODY_FIXED_FRAME = "IAU_EARTH"` or `nothing`. The return type is `Union{Nothing, String}` so callers must test with `=== nothing` before use. It is consumed by `_planet_fixed_to_j2000_state` and `_j2000_to_planet_fixed_state` inside their `catch` blocks.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Union{Nothing, | n/a | — | Return value of `_planet_fixed_frame_fallback_name`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.telemetry_loading__j2000_to_planet_fixed_state|_j2000_to_planet_fixed_state]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:54-54`
- [[analysis.telemetry_loading__planet_fixed_to_j2000_state|_planet_fixed_to_j2000_state]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:43-43`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Falling back from `ITRF93` to `IAU_EARTH` silently degrades pointing accuracy by up to tens of metres at the surface (polar motion and precise UT1 are dropped), and nothing records that the fallback was used. Because only Earth has a fallback, any other missing kernel surfaces as an unhandled SPICE error. The comparison is string-based; a planet object with a canonical id would be more robust.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 18.
