---
id: analysis.scenario_builders__telemetry_coefficients_normalized_for_scenario
label: _telemetry_coefficients_normalized_for_scenario
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _telemetry_coefficients_normalized_for_scenario
  lines:
  - 73
  - 73
inputs:
- id: scenario_name
  type: String
  units: n/a
  required: true
  description: Positional argument `scenario_name`.
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
  type: Bool
  units: n/a
  description: Return value of `_telemetry_coefficients_normalized_for_scenario`.
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

# _telemetry_coefficients_normalized_for_scenario

## Purpose
Decides whether a scenario's harmonics coefficient file is read as fully normalised, with a per-scenario environment override for unnormalised files.

## Design & Implementation
Reads `SPACEAGORA_TELEMETRY_HARMONICS_NORMALIZED_DEFAULT`, defaulting to true and treating `0`, `false`, `no` and `off` as false. If `SPACEAGORA_TELEMETRY_HARMONICS_UNNORMALIZED_SCENARIOS` names the scenario in its comma-separated list the result is false; otherwise the default applies.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `scenario_name` | String | n/a | yes | Positional argument `scenario_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Bool | n/a | — | Return value of `_telemetry_coefficients_normalized_for_scenario`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:138-138`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:83-83`
<!-- vulcan:connections:end -->

## Limitations
Normalisation is a property of the coefficient file, but it is configured per scenario name in the environment rather than declared in the manifest beside the file path, so the two can be mismatched without any check.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 73.
