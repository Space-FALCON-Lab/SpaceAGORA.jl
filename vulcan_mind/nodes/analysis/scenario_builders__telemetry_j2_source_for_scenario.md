---
id: analysis.scenario_builders__telemetry_j2_source_for_scenario
label: _telemetry_j2_source_for_scenario
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _telemetry_j2_source_for_scenario
  lines:
  - 58
  - 58
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
  type: Symbol
  units: n/a
  description: Return value of `_telemetry_j2_source_for_scenario`.
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

# _telemetry_j2_source_for_scenario

## Purpose
Decides whether a scenario's J2 term comes from the harmonics file's C20 coefficient or from the planet's own J2 constant, with an environment override for named scenarios.

## Design & Implementation
Reads `SPACEAGORA_TELEMETRY_J2_SOURCE_DEFAULT`, treating `planet` or `planet_j2` as `:planet_j2` and anything else as `:file_c20`. If `SPACEAGORA_TELEMETRY_J2_SOURCE_PLANET_SCENARIOS` is set it is split on commas into a lowercase set, and a scenario whose lowercased name is in the set returns `:planet_j2` regardless of the default.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `scenario_name` | String | n/a | yes | Positional argument `scenario_name`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | Symbol | n/a | — | Return value of `_telemetry_j2_source_for_scenario`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:139-139`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:68-68`
<!-- vulcan:connections:end -->

## Limitations
The set is rebuilt from the environment on every call; and an unrecognised default value silently maps to `:file_c20` rather than raising.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 58.
