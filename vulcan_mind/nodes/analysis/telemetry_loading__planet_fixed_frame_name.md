---
id: analysis.telemetry_loading__planet_fixed_frame_name
label: _planet_fixed_frame_name
kind: function
source:
  file: src/analysis/verification/telemetry_verification/telemetry_loading.jl
  symbol: _planet_fixed_frame_name
  lines:
  - 14
  - 14
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
  type: String
  units: n/a
  description: Return value of `_planet_fixed_frame_name`.
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

# _planet_fixed_frame_name

## Purpose
Maps a planet name to the SPICE body-fixed frame used for telemetry frame conversions, selecting the high-precision `ITRF93` frame for Earth and the generic `IAU_<PLANET>` frame for every other body.

## Design & Implementation
Marked `@inline`; normalises with `lowercase(strip(planet_name))`, compares to `"earth"`, and returns the constant `_EARTH_HIGH_PREC_BODY_FIXED_FRAME = "ITRF93"` on a match, else concatenates `"IAU_" * uppercase(strip(planet_name))`. The result feeds `sxform` in `_transform_state`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `planet_name` | String | n/a | yes | Positional argument `planet_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_planet_fixed_frame_name`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.telemetry_loading__j2000_to_planet_fixed_state|_j2000_to_planet_fixed_state]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:53-53`
- [[analysis.telemetry_loading__planet_fixed_to_j2000_state|_planet_fixed_to_j2000_state]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl:42-42`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/telemetry_loading.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
`ITRF93` requires the Earth high-precision binary PCK (`earth_*.bpc`) to be loaded; if it is not, `sxform` throws and the caller must fall back via `_planet_fixed_frame_fallback_name`. Non-Earth names are not validated, so `"Mars Express"` becomes `"IAU_MARS EXPRESS"`, a nonexistent frame that fails only at SPICE call time. Multi-word or aliased body names (for example `"Luna"`) are not mapped.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/telemetry_loading.jl` line 14.
