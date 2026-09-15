---
id: analysis.manifest_parsing__load_scenarios_from_manifest
label: _load_scenarios_from_manifest
kind: function
source:
  file: src/analysis/verification/telemetry_verification/manifest_parsing.jl
  symbol: _load_scenarios_from_manifest
  lines:
  - 522
  - 522
inputs:
- id: manifest_path
  type: String
  units: n/a
  required: true
  description: Positional argument `manifest_path`.
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
  type: Vector{AbstractScenarioConfig}
  units: n/a
  description: Return value of `_load_scenarios_from_manifest`.
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

# _load_scenarios_from_manifest

## Purpose
Reads the TOML manifest and builds the typed scenario configuration list that drives the whole verification study.

## Design & Implementation
Parses the file, requires a `scenarios` array of tables, and for each entry reads the common fields — name, kind, planet, events, units, both tolerance profiles, initial time, spacecraft, gravity model, entry interface, optional harmonics, N-body, SRP, drag, incidence mode, wind, altitude mode, manoeuvres, atmosphere truth and calibration. It rejects attitude quaternions under any incidence mode other than `attitude`, since the historical `max_drag` path reads link quaternions and would silently change physics. It then constructs an `OrbitEventsScenarioConfig` or `TimeAlignedScenarioConfig` depending on `kind`, the latter additionally reading the telemetry column map, comparison mode, extrema separation, frames, offsets and truth mask. An empty result is an error.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `manifest_path` | String | n/a | yes | Positional argument `manifest_path`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Vector{AbstractScenarioConfig} | n/a | — | Return value of `_load_scenarios_from_manifest`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_verification|_run_verification]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:355-355`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__optional_bool|_optional_bool]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:550-550`
- `callees` → [[analysis.manifest_parsing__optional_float|_optional_float]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:551-551`
- `callees` → [[analysis.manifest_parsing__optional_float_tuple|_optional_float_tuple]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:602-602`
- `callees` → [[analysis.manifest_parsing__optional_int|_optional_int]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:545-545`
- `callees` → [[analysis.manifest_parsing__optional_str|_optional_str]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:547-547`
- `callees` → [[analysis.manifest_parsing__optional_str_vector|_optional_str_vector]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:549-549`
- `callees` → [[analysis.manifest_parsing__parse_atmosphere_truth_config|_parse_atmosphere_truth_config]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:576-576`
- `callees` → [[analysis.manifest_parsing__parse_calibration_config|_parse_calibration_config]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:577-577`
- `callees` → [[analysis.manifest_parsing__parse_element_frame|_parse_element_frame]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:601-601`
- `callees` → [[analysis.manifest_parsing__parse_gravity_model|_parse_gravity_model]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:543-543`
- `callees` → [[analysis.manifest_parsing__parse_ic_offset|_parse_ic_offset]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:699-699`
- `callees` → [[analysis.manifest_parsing__parse_initial_time|_parse_initial_time]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:541-541`
- `callees` → [[analysis.manifest_parsing__parse_maneuver_config|_parse_maneuver_config]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:575-575`
- `callees` → [[analysis.manifest_parsing__parse_orbit_altitude_mode|_parse_orbit_altitude_mode]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:574-574`
- `callees` → [[analysis.manifest_parsing__parse_reference_frame|_parse_reference_frame]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:686-686`
- `callees` → [[analysis.manifest_parsing__parse_spacecraft_config|_parse_spacecraft_config]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:542-542`
- `callees` → [[analysis.manifest_parsing__parse_time_aligned_comparison_mode|_parse_time_aligned_comparison_mode]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:633-633`
- `callees` → [[analysis.manifest_parsing__parse_tolerances|_parse_tolerances]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:538-538`
- `callees` → [[analysis.manifest_parsing__parse_truth_mask|_parse_truth_mask]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:701-701`
- `callees` → [[analysis.manifest_parsing__parse_units|_parse_units]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:537-537`
- `callees` → [[analysis.manifest_parsing__require_float|_require_float]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:544-544`
- `callees` → [[analysis.manifest_parsing__require_int|_require_int]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:540-540`
- `callees` → [[analysis.manifest_parsing__require_key|_require_key]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:524-524`
- `callees` → [[analysis.manifest_parsing__require_str|_require_str]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:533-533`
- `callees` → [[analysis.manifest_parsing__require_str_vector|_require_str_vector]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:536-536`
- `callees` → [[analysis.manifest_parsing__require_table|_require_table]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:541-541`
- `callees` → [[analysis.manifest_parsing__resolve_repo_path|_resolve_repo_path]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:548-548`
- `callees` → [[analysis.types_orbiteventsscenarioconfig|OrbitEventsScenarioConfig]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:580-580`
- `callees` → [[analysis.types_timealignedscenarioconfig|TimeAlignedScenarioConfig]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:641-641`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:580-580`
<!-- vulcan:connections:end -->

## Limitations
Two `error(...)` calls raise `ErrorException` while everything else raises `ArgumentError`, so a caller catching by type misses them; the function is roughly 170 lines because both scenario constructors receive twenty-plus shared keywords by name.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/manifest_parsing.jl` line 522.
