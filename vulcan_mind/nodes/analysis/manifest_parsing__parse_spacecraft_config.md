---
id: analysis.manifest_parsing__parse_spacecraft_config
label: _parse_spacecraft_config
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_spacecraft_config
  lines:
  - 471
  - 471
inputs:
- id: tbl
  type: Any
  units: n/a
  required: true
  description: Positional argument `tbl`.
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
  type: SpacecraftConfig
  units: n/a
  description: Return value of `_parse_spacecraft_config`.
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

# _parse_spacecraft_config

## Purpose
Parses the mandatory `spacecraft` table into a `SpacecraftConfig` for the three-body example vehicle.

## Design & Implementation
Requires bus and panel dimensions as three-vectors, bus mass, per-panel mass, panel y offset, propellant mass and integer id; reads `bus_ram_face` as a symbol defaulting to `legacy`; and parses the three optional attitude quaternions.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SpacecraftConfig | n/a | — | Return value of `_parse_spacecraft_config`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:542-542`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__optional_str|_optional_str]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:481-481`
- `callees` → [[analysis.manifest_parsing__parse_attitude_q|_parse_attitude_q]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:482-482`
- `callees` → [[analysis.manifest_parsing__parse_vec3|_parse_vec3]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:474-474`
- `callees` → [[analysis.manifest_parsing__require_float|_require_float]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:476-476`
- `callees` → [[analysis.manifest_parsing__require_int|_require_int]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:480-480`
- `callees` → [[analysis.manifest_parsing__require_table|_require_table]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:472-472`
- `callees` → [[analysis.types_spacecraftconfig|SpacecraftConfig]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:473-473`
<!-- vulcan:connections:end -->

## Limitations
`bus_ram_face` is converted to a symbol without validation here; the spacecraft builder rejects anything other than `legacy` or `frontal` later.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 471.
