---
id: analysis.manifest_parsing__parse_maneuver_config
label: _parse_maneuver_config
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _parse_maneuver_config
  lines:
  - 205
  - 205
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
  type: Any
  units: n/a
  description: Return value of `_parse_maneuver_config`. Returns `(`.
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

# _parse_maneuver_config

## Purpose
Parses the optional `maneuvers` table of an orbit-events scenario into the orbit numbers, delta-v values, replay scaling and thruster parameters used to replay campaign burns.

## Design & Implementation
Returns a no-manoeuvre tuple with default rates of 30 s guidance and 10 s control when the table is absent. Otherwise it requires non-empty positive orbit numbers and a matching-length delta-v list, validates `replay_scale_mode` as `delta_v` or `flight_apoapsis_ratio` and, in the latter mode, requires a matching list of positive flight apoapsis altitudes. `orbit_number_offset` converts campaign-numbered orbits to epoch-relative ones; burns that shift below one are dropped with a printed count, while the unshifted list is kept as `orbit_numbers_campaign` for diagnostics. Altitudes are converted to metres; thrust defaults to 4 N and Isp to 220 s.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `tbl` | Any | n/a | yes | Positional argument `tbl`. |
| in | `context` | String | n/a | yes | Positional argument `context`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Any | n/a | — | Return value of `_parse_maneuver_config`. Returns `(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:575-575`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__optional_float64_vector|_optional_float64_vector]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:221-221`
- `callees` → [[analysis.manifest_parsing__optional_float|_optional_float]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:279-279`
- `callees` → [[analysis.manifest_parsing__optional_int64_vector|_optional_int64_vector]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:220-220`
- `callees` → [[analysis.manifest_parsing__optional_int|_optional_int]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:234-234`
- `callees` → [[analysis.manifest_parsing__optional_str|_optional_str]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:242-242`
- `callees` → [[analysis.manifest_parsing__require_table|_require_table]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:219-219`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:264-264`
<!-- vulcan:connections:end -->

## Limitations
The pre-epoch drop is reported by `println` rather than a structured warning; a dropped burn that was actually after the epoch because of a wrong offset is silently lost.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 205.
