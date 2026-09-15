---
id: analysis.reporting__value_units
label: _value_units
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _value_units
  lines:
  - 22
  - 22
inputs:
- id: cfg
  type: OrbitEventsScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: event
  type: String
  units: n/a
  required: true
  description: Positional argument `event`.
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
  description: Return value of `_value_units`. Returns `get(cfg.units_y, event, "km/s")`
    or `get(cfg.units_y, event, "km")`.
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

# _value_units

## Purpose
Determines the physical unit of the dependent value for a named comparison event, defaulting velocity-like events to kilometres per second and everything else to kilometres unless the scenario overrides them.

## Design & Implementation
Identical bodies for `OrbitEventsScenarioConfig` and `TimeAlignedScenarioConfig`. The `event::String` is scanned with `occursin` for the substrings `"speed"`, `"vx"`, `"vy"`, `"vz"`; a hit yields `get(cfg.units_y, event, "km/s")`, otherwise `get(cfg.units_y, event, "km")`. Explicit entries in `cfg.units_y` always win over the heuristic. Both methods are `@inline` and pure.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `event` | String | n/a | yes | Positional argument `event`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_value_units`. Returns `get(cfg.units_y, event, "km/s")` or `get(cfg.units_y, event, "km")`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:390-390`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The substring test is case-sensitive and can misfire: an event named `"convex"` contains `"vx"` and would default to km/s. Events measured in neither km nor km/s (angles, altitude in metres) must be listed in `units_y` explicitly or they are mislabelled. The two methods are copy-pasted rather than shared, so a fix to one must be repeated in the other.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 22.
