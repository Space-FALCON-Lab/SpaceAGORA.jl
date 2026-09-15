---
id: analysis.manifest_parsing__require_float
label: _require_float
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _require_float
  lines:
  - 71
  - 71
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
  type: Float64
  units: n/a
  description: Return value of `_require_float`.
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

# _require_float

## Purpose
Fetches a mandatory numeric field as `Float64`, used for masses, radii, angles, tolerances and the entry-interface altitude.

## Design & Implementation
Wraps `_require_key` in `Float64(...)`, so a TOML integer literal such as `220` is accepted where a float is expected, which keeps hand-written manifests forgiving. `@inline` with a `::Float64` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Float64 | n/a | — | Return value of `_require_float`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:544-544`
- [[analysis.manifest_parsing__parse_event_tolerance|_parse_event_tolerance]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:501-501`
- [[analysis.manifest_parsing__parse_initial_time|_parse_initial_time]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:447-447`
- [[analysis.manifest_parsing__parse_spacecraft_config|_parse_spacecraft_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:476-476`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:72-72`
- `callees` → [[analysis.manifest_parsing__require_key|_require_key]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:72-72`
<!-- vulcan:connections:end -->

## Limitations
No range or finiteness validation, so negative or non-finite values pass through to whichever consumer reads them; each scenario builder that cares must check separately.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 71.
