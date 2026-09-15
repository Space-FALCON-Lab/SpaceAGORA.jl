---
id: analysis.reporting__source_file
label: _source_file
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _source_file
  lines:
  - 117
  - 117
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
  description: 'Return value of `_source_file`. Returns `event == "peri" ? cfg.telemetry_peri_path
    : cfg.telemetry_apo_path`.'
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

# _source_file

## Purpose
Reports which telemetry input file an event's reference data came from, for provenance columns in the summary output.

## Design & Implementation
For `OrbitEventsScenarioConfig` the `event::String` selects between two paths: `"peri"` returns `cfg.telemetry_peri_path`, any other value returns `cfg.telemetry_apo_path`. For `TimeAlignedScenarioConfig` the event is ignored and `cfg.telemetry_path` is returned. Both are `@inline` and side-effect free.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `event` | String | n/a | yes | Positional argument `event`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_source_file`. Returns `event == "peri" ? cfg.telemetry_peri_path : cfg.telemetry_apo_path`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:388-388`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The orbit-events method treats every non-`"peri"` event (including `"peri_speed"`) as apoapsis, so a periapsis-speed event is attributed to the apoapsis file. Paths are returned as stored, with no existence check or normalisation.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 117.
