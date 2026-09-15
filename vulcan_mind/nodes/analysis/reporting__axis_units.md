---
id: analysis.reporting__axis_units
label: _axis_units
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _axis_units
  lines:
  - 19
  - 19
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
  type: Any
  units: n/a
  description: Return value of `_axis_units`. Returns `cfg.units_x`.
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

# _axis_units

## Purpose
Returns the unit label for the independent (x) axis of a telemetry scenario so plots and CSV headers agree on whether time is in seconds, hours, days, or orbit number.

## Design & Implementation
Two `@inline` methods, one for `OrbitEventsScenarioConfig` and one for `TimeAlignedScenarioConfig`, both simply forwarding the `units_x` field of the config. The pair exists so callers can be written against the abstract config without a type test. No mutation, no exceptions, no defaulting: whatever string the scenario TOML supplied is returned verbatim.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_axis_units`. Returns `cfg.units_x`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:389-389`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The returned string is not validated against any known unit vocabulary, so an unusual `units_x` value propagates unchanged into column names and plot labels. No method exists for other `AbstractScenarioConfig` subtypes.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 19.
