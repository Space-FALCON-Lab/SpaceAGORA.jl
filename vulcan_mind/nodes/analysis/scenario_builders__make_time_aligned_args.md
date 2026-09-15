---
id: analysis.scenario_builders__make_time_aligned_args
label: _make_time_aligned_args
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _make_time_aligned_args
  lines:
  - 589
  - 589
inputs:
- id: cfg
  type: TimeAlignedScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: mission_time_s
  type: Float64
  units: n/a
  required: true
  description: Positional argument `mission_time_s`.
- id: ic
  type: AbstractInitialCondition
  units: n/a
  required: true
  description: Positional argument `ic`.
- id: cd_scale
  type: Float64
  units: n/a
  required: false
  description: Keyword argument `cd_scale` (default `1.0`).
- id: cr_override
  type: Union{Nothing, Float64}
  units: n/a
  required: false
  description: Keyword argument `cr_override` (default `nothing`).
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
  description: Return value of `_make_time_aligned_args`.
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

# _make_time_aligned_args

## Purpose
The top-level builder for a time-aligned scenario, where the initial condition and mission duration come from the telemetry rather than from orbit geometry.

## Design & Implementation
Resolves the planet, builds the spacecraft from the supplied `ic`, the effectors with optional overrides and the density model, then calls `make_example_config` with the caller's `mission_time_s` and finally applies the wind flag. Unlike the orbit builder, no orbit-counted termination or manoeuvre layer is added.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | TimeAlignedScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `mission_time_s` | Float64 | n/a | yes | Positional argument `mission_time_s`. |
| in | `ic` | AbstractInitialCondition | n/a | yes | Positional argument `ic`. |
| in | `cd_scale` | Float64 | n/a | no | Keyword argument `cd_scale` (default `1.0`). |
| in | `cr_override` | Union{Nothing, Float64} | n/a | no | Keyword argument `cr_override` (default `nothing`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationConfiguration | n/a | — | Return value of `_make_time_aligned_args`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:246-246`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.scenario_builders__make_spacecraft|_make_spacecraft]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:602-602`
- `callees` → [[analysis.scenario_builders__planet_from_name|_planet_from_name]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:596-596`
- `callees` → [[analysis.scenario_builders__scenario_density_model|_scenario_density_model]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:610-610`
- `callees` → [[analysis.scenario_builders__scenario_dynamic_effectors|_scenario_dynamic_effectors]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:603-603`
- `callees` → [[analysis.scenario_builders__with_environment_wind|_with_environment_wind]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:624-624`
- `callees` → [[envana.ana_example_support_make_example_config|make_example_config]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:612-612`
<!-- vulcan:connections:end -->

## Limitations
The historical entry-interface workaround mentioned in the source comment is no longer needed, but the function still exposes `EI_km` from the manifest, so a stale manifest value continues to affect where the atmosphere partition begins.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 589.
