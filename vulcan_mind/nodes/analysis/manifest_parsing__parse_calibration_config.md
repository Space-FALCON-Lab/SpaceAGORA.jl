---
id: analysis.manifest_parsing__parse_calibration_config
label: _parse_calibration_config
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_calibration_config
  lines:
  - 388
  - 388
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
  type: CalibrationConfig
  units: n/a
  description: Return value of `_parse_calibration_config`.
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

# _parse_calibration_config

## Purpose
Parses the optional `calibration` table controlling the drag-scale, reflectivity and bias search a verification run may perform.

## Design & Implementation
Returns defaults when absent. Reads the profile list, validates `objective` as `mean_nmae`, `mean_rmse_km` or `max_nmae`, requires positive step counts, and requires the maximum of each search range to be at least its minimum. Defaults are three steps over drag scale 0.85 to 1.15 and reflectivity 1.15 to 1.45, bias fitting on with a 500 km cap, and calibration disabled.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | CalibrationConfig | n/a | — | Return value of `_parse_calibration_config`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:577-577`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__optional_bool|_optional_bool]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:409-409`
- `callees` → [[analysis.manifest_parsing__optional_float|_optional_float]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:402-402`
- `callees` → [[analysis.manifest_parsing__optional_int|_optional_int]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:398-398`
- `callees` → [[analysis.manifest_parsing__optional_str|_optional_str]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:394-394`
- `callees` → [[analysis.manifest_parsing__optional_symbol_vector|_optional_symbol_vector]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:393-393`
- `callees` → [[analysis.manifest_parsing__require_table|_require_table]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:392-392`
- `callees` → [[analysis.types_calibrationconfig|CalibrationConfig]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:390-390`
<!-- vulcan:connections:end -->

## Limitations
Ranges are validated for ordering but not for physical plausibility; a drag scale range including zero would be accepted.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 388.
