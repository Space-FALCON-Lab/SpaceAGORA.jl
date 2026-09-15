---
id: analysis.manifest_parsing__parse_gravity_model
label: _parse_gravity_model
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_gravity_model
  lines:
  - 430
  - 430
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
  description: Return value of `_parse_gravity_model`.
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

# _parse_gravity_model

## Purpose
Maps the manifest's gravity model string onto the `:inverse_squared` or `:inverse_squared_j2` symbol that selects the base gravity effector.

## Design & Implementation
Lowercases and strips the input, accepts `inverse_squared`, `is2` or `point_mass` for the point-mass model and `inverse_squared_j2`, `is2_j2` or `j2` for the J2 model, and raises `ArgumentError` naming the value and context otherwise. `@inline` with a `::Symbol` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | String | n/a | yes | Positional argument `raw`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_parse_gravity_model`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:543-543`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Higher-order gravity is configured through the separate `gravity_harmonics_degree` and file keys, and when those are set this symbol is ignored entirely by the effector builder, which the parser does not warn about.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 430.
