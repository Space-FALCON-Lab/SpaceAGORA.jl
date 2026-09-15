---
id: analysis.scenario_builders__with_study_settings
label: _with_study_settings
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _with_study_settings
  lines:
  - 638
  - 638
inputs:
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: quick
  type: Bool
  units: n/a
  required: false
  description: Keyword argument `quick` (default `false`).
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
  type: SimulationConfiguration
  units: n/a
  description: Return value of `_with_study_settings`.
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

# _with_study_settings

## Purpose
Applies the verification study's integrator tolerances, step limits and output settings to a scenario configuration, in quick or full profile.

## Design & Implementation
Base tolerances are 1e-7 relative and 1e-9 absolute in full profile, five and fifty times looser in quick, each capped by the `STRICT_*` constants. Environment variables `SPACEAGORA_TELEMETRY_{RELTOL,ABSTOL}_{ORBIT,ATM}` may tighten but never loosen them. Maximum step sizes depend on the profile and on `_has_high_fidelity_effectors`, ranging from 60 to 240 seconds in orbit and 0.2 to 5 seconds in atmosphere, again cappable by environment. Output is forced to results and CSV on, plots off, and 2000 saved steps.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `quick` | Bool | n/a | no | Keyword argument `quick` (default `false`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationConfiguration | n/a | — | Return value of `_with_study_settings`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:162-162`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.manifest_parsing__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:647-647`
- `callees` → [[analysis.scenario_builders__has_high_fidelity_effectors|_has_high_fidelity_effectors]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:639-639`
- `callees` → [[core.simulation_configuration_integrationtolerances|IntegrationTolerances]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:687-687`
- `callees` → [[core.simulation_configuration_missionconfiguration|MissionConfiguration]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:672-672`
- `callees` → [[core.simulation_configuration_simulationsettings|SimulationSettings]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:663-663`
- `callees` → [[parcore.simulation_configuration_simulationconfiguration|SimulationConfiguration]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:661-661`
- `callees` → [[simulation.setup__parse_positive_float_env|_parse_positive_float_env]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:647-647`
<!-- vulcan:connections:end -->

## Limitations
Six environment variables silently interact with two profiles and a fidelity flag, so the effective tolerances of a given run are hard to know without reproducing this function's logic; they are not echoed anywhere.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 638.
