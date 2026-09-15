---
id: analysis.manifest_parsing__parse_atmosphere_truth_config
label: _parse_atmosphere_truth_config
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_atmosphere_truth_config
  lines:
  - 297
  - 297
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
  type: AtmosphereTruthConfig
  units: n/a
  description: Return value of `_parse_atmosphere_truth_config`.
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

# _parse_atmosphere_truth_config

## Purpose
Parses the optional `atmosphere_truth` table into the `AtmosphereTruthConfig` that selects and parameterises the scenario's density source.

## Design & Implementation
Returns defaults when absent. Requires `atmosphere_model` in `GRAM`, `tabulated_flight`, `nrlmsise00` or `tabulated_time`, and enforces that the tabulated file keys are present only with their matching model, with the flight sigma within ±3 and the time scale positive. It requires dataset, space-weather and solar-flux identifiers, reads the GRAM seed (default 1001), four perturbation scales, optional minimum step, the `gram_offline_surrogate` and `gram_global_lock` vocabularies, the static-grid and track-cache flags, and the Mars-specific tuples and scalars as optional values that stay `nothing` when absent. File paths are resolved through `_resolve_repo_path`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AtmosphereTruthConfig | n/a | — | Return value of `_parse_atmosphere_truth_config`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:576-576`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__optional_bool|_optional_bool]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:366-366`
- `callees` → [[analysis.manifest_parsing__optional_float|_optional_float]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:308-308`
- `callees` → [[analysis.manifest_parsing__optional_float_tuple|_optional_float_tuple]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:340-340`
- `callees` → [[analysis.manifest_parsing__optional_int|_optional_int]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:339-339`
- `callees` → [[analysis.manifest_parsing__optional_str|_optional_str]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:302-302`
- `callees` → [[analysis.manifest_parsing__require_str|_require_str]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:303-303`
- `callees` → [[analysis.manifest_parsing__require_table|_require_table]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:301-301`
- `callees` → [[analysis.manifest_parsing__resolve_repo_path|_resolve_repo_path]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:378-378`
- `callees` → [[analysis.types_atmospheretruthconfig|AtmosphereTruthConfig]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:299-299`
<!-- vulcan:connections:end -->

## Limitations
The `haskey ? _optional_x : nothing` pattern is repeated for six Mars fields where a single optional accessor returning `nothing` would be simpler; the three required identifier strings are not validated against any vocabulary.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 297.
