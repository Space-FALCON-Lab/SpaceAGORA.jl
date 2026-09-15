---
id: analysis.reporting__maneuver_count
label: _maneuver_count
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _maneuver_count
  lines:
  - 125
  - 125
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
  description: Return value of `_maneuver_count`. Returns `length(cfg.maneuver_orbit_numbers)`.
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

# _maneuver_count

## Purpose
Reports how many manoeuvres a scenario replays, used in status lines and summary metadata to indicate whether the comparison crossed burn events.

## Design & Implementation
For `OrbitEventsScenarioConfig` it returns `length(cfg.maneuver_orbit_numbers)`, the number of orbit indices flagged as containing a manoeuvre. The `TimeAlignedScenarioConfig` method discards its argument and returns the constant `0`, since time-aligned state comparison does not model manoeuvre replay. Both are `@inline` with no side effects.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_maneuver_count`. Returns `length(cfg.maneuver_orbit_numbers)`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:392-392`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Duplicate orbit numbers in `maneuver_orbit_numbers` are counted twice. A time-aligned scenario that actually includes burns still reports zero, which can be misleading when the underlying telemetry contains them.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 125.
