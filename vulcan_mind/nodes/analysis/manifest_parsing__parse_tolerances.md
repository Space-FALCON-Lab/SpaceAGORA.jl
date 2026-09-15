---
id: analysis.manifest_parsing__parse_tolerances
label: _parse_tolerances
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_tolerances
  lines:
  - 507
  - 507
inputs:
- id: tbl
  type: Any
  units: n/a
  required: true
  description: Positional argument `tbl`.
- id: key
  type: String
  units: n/a
  required: true
  description: Positional argument `key`.
- id: events
  type: Vector{String}
  units: n/a
  required: true
  description: Positional argument `events`.
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
  type: Dict{String,
  units: n/a
  description: Return value of `_parse_tolerances`.
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

# _parse_tolerances

## Purpose
Parses a full tolerance table for every event, plus any explicitly provided derived speed-channel tolerances.

## Design & Implementation
Requires the table at `key`, parses each event's tolerance, and for each event also checks for a `<event>_speed` table, parsing it when present. Returns a `Dict` from event name to tolerance.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `events` | Vector{String} | n/a | yes | Positional argument `events`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Dict{String, | n/a | — | Return value of `_parse_tolerances`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:538-538`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_event_tolerance|_parse_event_tolerance]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:511-511`
- `callees` → [[analysis.manifest_parsing__require_table|_require_table]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:508-508`
<!-- vulcan:connections:end -->

## Limitations
The source comment says speed channels inherit the base event's tolerances when absent, but the code simply omits them, so a missing speed tolerance leads to a lookup failure downstream rather than inheritance.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 507.
