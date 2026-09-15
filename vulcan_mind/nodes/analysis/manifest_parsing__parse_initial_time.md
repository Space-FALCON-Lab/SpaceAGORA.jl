---
id: analysis.manifest_parsing__parse_initial_time
label: _parse_initial_time
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_initial_time
  lines:
  - 440
  - 440
inputs:
- id: tbl
  type: Any
  units: n/a
  required: true
  description: Positional argument `tbl`.
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
  type: InitialTime
  units: n/a
  description: Return value of `_parse_initial_time`.
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

# _parse_initial_time

## Purpose
Parses the six-field `initial_time` table into an `InitialTime`.

## Design & Implementation
Requires `year`, `month`, `day`, `hour` and `minute` as integers and `second` as a float, forwarding each to the keyword constructor.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | InitialTime | n/a | — | Return value of `_parse_initial_time`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:541-541`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__require_float|_require_float]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:447-447`
- `callees` → [[analysis.manifest_parsing__require_int|_require_int]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:442-442`
- `callees` → [[core.simulation_configuration_initialtime|InitialTime]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:441-441`
<!-- vulcan:connections:end -->

## Limitations
No calendar validation; a month of 13 is accepted here and fails only when converted to an epoch.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 440.
