---
id: analysis.reporting__orbit_altitude_mode
label: _orbit_altitude_mode
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _orbit_altitude_mode
  lines:
  - 120
  - 120
inputs:
- id: cfg
  type: OrbitEventsScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  description: Return value of `_orbit_altitude_mode`. Returns `String(cfg.orbit_altitude_mode)`.
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

# _orbit_altitude_mode

## Purpose
Produces the string recorded in reports describing how orbit altitude was defined for the scenario (for example radius-based versus surface-relative), or `"n/a"` when the comparison does not involve orbit events.

## Design & Implementation
`OrbitEventsScenarioConfig` always yields `String(cfg.orbit_altitude_mode)`. `TimeAlignedScenarioConfig` yields the same only when `cfg.comparison_mode == :orbit_events`; otherwise the literal `"n/a"`. The `String(...)` call converts a `Symbol` field to text for CSV serialisation. Both methods are `@inline`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_orbit_altitude_mode`. Returns `String(cfg.orbit_altitude_mode)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:391-391`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The `"n/a"` sentinel is a plain string, so consumers cannot distinguish an unset mode from a scenario genuinely named `n/a` without also checking `comparison_mode`. No validation that `orbit_altitude_mode` is one of the values the analysis code understands.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 120.
