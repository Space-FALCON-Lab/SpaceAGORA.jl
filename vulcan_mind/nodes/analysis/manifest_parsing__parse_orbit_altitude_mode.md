---
id: analysis.manifest_parsing__parse_orbit_altitude_mode
label: _parse_orbit_altitude_mode
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_orbit_altitude_mode
  lines:
  - 142
  - 142
inputs:
- id: raw
  type: String
  units: n/a
  required: true
  description: Positional argument `raw`.
- id: context
  type: String
  units: n/a
  required: true
  description: Positional argument `context`.
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
  type: Symbol
  units: n/a
  description: Return value of `_parse_orbit_altitude_mode`.
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

# _parse_orbit_altitude_mode

## Purpose
Maps the manifest's altitude convention string onto the `:vacuum` or `:oblate` symbol the error tables dispatch on.

## Design & Implementation
Lowercases and strips, accepts `vacuum`, `spherical` or `rp_e` for `:vacuum` and `true`, `oblate` or `geodetic` for `:oblate`, and raises `ArgumentError` naming the context otherwise. `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | String | n/a | yes | Positional argument `raw`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_parse_orbit_altitude_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:574-574`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The canonical spelling in the error message is `vacuum|true`, which does not match the symbol names, a small source of confusion when reading logs.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 142.
