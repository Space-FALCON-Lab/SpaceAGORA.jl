---
id: analysis.manifest_parsing__parse_time_aligned_comparison_mode
label: _parse_time_aligned_comparison_mode
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_time_aligned_comparison_mode
  lines:
  - 166
  - 166
inputs:
- id: raw
  type: String
  units: n/a
  required: true
  description: Positional argument `raw`.
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
  type: Symbol
  units: n/a
  description: Return value of `_parse_time_aligned_comparison_mode`.
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

# _parse_time_aligned_comparison_mode

## Purpose
Maps the comparison-mode string of a time-aligned scenario onto `:time_aligned_state` or `:orbit_events`, which selects how the error tables score the run.

## Design & Implementation
Lowercases and strips, accepts `time_aligned_state`, `time_aligned` or `state` for direct state comparison and `orbit_events`, `apo_peri` or `extrema` for apsis-based comparison, raising `ArgumentError` with the accepted canonical names otherwise. `@inline` with a `::Symbol` return.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `raw` | String | n/a | yes | Positional argument `raw`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_parse_time_aligned_comparison_mode`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:633-633`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
The default when the key is absent is supplied by the caller as `time_aligned_state`, so this function itself never sees a missing value and cannot report one.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 166.
