---
id: analysis.reporting__maneuver_replay_scale_mode
label: _maneuver_replay_scale_mode
kind: function
source:
  file: src/analysis/verification/telemetry_verification/reporting.jl
  symbol: _maneuver_replay_scale_mode
  lines:
  - 128
  - 128
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
  description: Return value of `_maneuver_replay_scale_mode`. Returns `cfg.maneuver_replay_scale_mode`.
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

# _maneuver_replay_scale_mode

## Purpose
Exposes the scenario's manoeuvre replay scaling strategy for inclusion in report metadata, so a reader can tell whether replayed burns were scaled by delta-v or by another rule.

## Design & Implementation
The `OrbitEventsScenarioConfig` method returns `cfg.maneuver_replay_scale_mode` unchanged (a `String` field). The `TimeAlignedScenarioConfig` method ignores its argument and returns the literal `"delta_v"` as the nominal default. Both are `@inline` accessors with no validation or mutation.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | OrbitEventsScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_maneuver_replay_scale_mode`. Returns `cfg.maneuver_replay_scale_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:393-393`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/reporting.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The hard-coded `"delta_v"` for time-aligned scenarios is a label only; nothing in that path actually performs delta-v scaling, so the metadata can overstate what was done. Valid mode strings are not enumerated here.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/reporting.jl` line 128.
