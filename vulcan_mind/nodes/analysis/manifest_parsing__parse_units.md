---
id: analysis.manifest_parsing__parse_units
label: _parse_units
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_units
  lines:
  - 488
  - 488
inputs:
- id: tbl
  type: Any
  units: n/a
  required: true
  description: Positional argument `tbl`.
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
  type: Tuple{String,
  units: n/a
  description: Return value of `_parse_units`.
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

# _parse_units

## Purpose
Parses the `units` table, requiring an x-axis unit and one y-axis unit per event.

## Design & Implementation
Requires the table, reads `x`, and builds a `Dict` from each event name to its required unit string. Returns the pair.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `events` | Vector{String} | n/a | yes | Positional argument `events`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Tuple{String, | n/a | — | Return value of `_parse_units`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:537-537`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__require_str|_require_str]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:490-490`
- `callees` → [[analysis.manifest_parsing__require_table|_require_table]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:489-489`
<!-- vulcan:connections:end -->

## Limitations
Units are carried as free strings for plot labels and are not used to convert data; a mismatch between the declared unit and the telemetry file's actual unit is invisible here.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 488.
