---
id: analysis.manifest_parsing__optional_float64_vector
label: _optional_float64_vector
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _optional_float64_vector
  lines:
  - 133
  - 133
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
  type: Vector{Float64}
  units: n/a
  description: Return value of `_optional_float64_vector`.
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

# _optional_float64_vector

## Purpose
Fetches an optional array of floats as `Vector{Float64}`, used for the manoeuvre delta-v list and the flight apoapsis altitude list.

## Design & Implementation
Returns an empty `Float64` vector when the key is absent; otherwise requires the value to be an `AbstractVector` and converts each element through `Float64`, so integer literals in the manifest are accepted. Declared `@inline` with a `::Vector{Float64}` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `key` | String | n/a | yes | Positional argument `key`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{Float64} | n/a | — | Return value of `_optional_float64_vector`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__parse_maneuver_config|_parse_maneuver_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:221-221`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:139-139`
<!-- vulcan:connections:end -->

## Limitations
Its error message hard-codes `manifest scenario table` as the context instead of accepting a context argument, and no finiteness or sign check is applied — the manoeuvre parser validates entries after reading.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 133.
