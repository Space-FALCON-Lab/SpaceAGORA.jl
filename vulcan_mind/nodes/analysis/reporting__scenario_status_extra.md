---
id: analysis.reporting__scenario_status_extra
label: _scenario_status_extra
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _scenario_status_extra
  lines:
  - 134
  - 134
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
  description: Return value of `_scenario_status_extra`.
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

# _scenario_status_extra

## Purpose
Builds the trailing fragment of a scenario status line that names the comparison mode, altitude mode and manoeuvre count, giving operators a one-glance description of how each scenario was evaluated.

## Design & Implementation
Returns a `String` beginning with a comma and space so it can be appended directly to a prefix. For `OrbitEventsScenarioConfig` the format is `", altitude_mode=<mode>, maneuvers=<n>"` using `String(cfg.orbit_altitude_mode)` and `length(cfg.maneuver_orbit_numbers)`. For `TimeAlignedScenarioConfig` it emits `", comparison_mode=orbit_events, altitude_mode=<mode>"` when `cfg.comparison_mode == :orbit_events`, otherwise `", comparison_mode=time_aligned_state"`. Both are `@inline` and allocate one string.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | String | n/a | — | Return value of `_scenario_status_extra`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:376-376`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The leading comma is baked in, so the caller must supply a non-empty prefix or the line starts with punctuation. The two methods report different key sets (only one includes `maneuvers`), making the status lines non-uniform across scenario types.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 134.
