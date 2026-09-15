---
id: analysis.scenario_builders__has_campaign_maneuvers
label: _has_campaign_maneuvers
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _has_campaign_maneuvers
  lines:
  - 394
  - 394
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
  type: Bool
  units: n/a
  description: Return value of `_has_campaign_maneuvers`.
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

# _has_campaign_maneuvers

## Purpose
Tests whether an orbit-events scenario replays propulsive manoeuvres, by checking for any configured manoeuvre orbit numbers.

## Design & Implementation
Returns `!isempty(cfg.maneuver_orbit_numbers)`. `@inline` with a `::Bool` return. It gates `_with_campaign_maneuvers` so scenarios without burns keep their empty guidance and control models.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_has_campaign_maneuvers`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__with_campaign_maneuvers|_with_campaign_maneuvers]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:399-399`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
It looks only at `maneuver_orbit_numbers`, not the campaign variant `maneuver_orbit_numbers_campaign` that the error tables consult, so the two can disagree about whether a scenario has manoeuvres.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 394.
