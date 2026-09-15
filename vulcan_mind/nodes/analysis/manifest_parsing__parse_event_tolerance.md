---
id: analysis.manifest_parsing__parse_event_tolerance
label: _parse_event_tolerance
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_event_tolerance
  lines:
  - 498
  - 498
inputs:
- id: ttbl
  type: Any
  units: n/a
  required: true
  description: Positional argument `ttbl`.
- id: event
  type: String
  units: n/a
  required: true
  description: Positional argument `event`.
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
  type: EventTolerance
  units: n/a
  description: Return value of `_parse_event_tolerance`.
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

# _parse_event_tolerance

## Purpose
Parses one event's tolerance sub-table into an `EventTolerance` named tuple.

## Design & Implementation
Requires the event's table, reads mandatory `max_abs_km` and `max_nmae`, and optional `max_rmse_km` defaulting to `Inf`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `ttbl` | Any | n/a | yes | Positional argument `ttbl`. |
| in | `event` | String | n/a | yes | Positional argument `event`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | EventTolerance | n/a | — | Return value of `_parse_event_tolerance`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__parse_tolerances|_parse_tolerances]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:511-511`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__optional_float|_optional_float]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:503-503`
- `callees` → [[analysis.manifest_parsing__require_float|_require_float]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:501-501`
- `callees` → [[analysis.manifest_parsing__require_table|_require_table]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:499-499`
<!-- vulcan:connections:end -->

## Limitations
Tolerances are in kilometres by convention even for the speed channels, whose docs say km/s; the naming is not adjusted per channel.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 498.
