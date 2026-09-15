---
id: analysis.scenario_builders__scenario_density_model
label: _scenario_density_model
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _scenario_density_model
  lines:
  - 176
  - 176
inputs:
- id: cfg
  type: AbstractScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
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
  type: SimulationModel.NRLMSISE00AtmosphereModel
  units: n/a
  description: Return value of `_scenario_density_model`. Returns `_make_tabulated_flight_density_model(cfg.initial_time,
    cfg.atmosphere_truth)` or `_make_time_tabulated_density_model(cfg.atmosphere_truth)`
    or `SimulationModel.NRLMSISE00AtmosphereModel(use_space_indices=true)` or `_make_required_gram_density_model(cfg.planet_name,
    cfg.initial_time, cfg.atmosph`.
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

# _scenario_density_model

## Purpose
Chooses the atmosphere truth model a scenario is simulated against, from five named sources.

## Design & Implementation
Returns `NoAtmosphereModel` when drag is disabled. Otherwise it dispatches on `atmosphere_truth.atmosphere_model`: `tabulated_flight` and `tabulated_time` build table-backed models from files; `nrlmsise00` is accepted only for Earth and constructs `NRLMSISE00AtmosphereModel` with space indices; anything else falls through to the required GRAM model with its library-missing fallback.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | AbstractScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.NRLMSISE00AtmosphereModel | n/a | — | Return value of `_scenario_density_model`. Returns `_make_tabulated_flight_density_model(cfg.initial_time, cfg.atmosphere_truth)` or `_make_time_tabulated_density_model(cfg.atmosphere_truth)` or `SimulationModel.NRLMSISE00AtmosphereModel(use_space_indices=true)` or `_make_required_gram_density_model(cfg.planet_name, cfg.initial_time, cfg.atmosph`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__make_orbit_args|_make_orbit_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:567-567`
- [[analysis.scenario_builders__make_time_aligned_args|_make_time_aligned_args]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:610-610`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.scenario_builders__make_required_gram_density_model|_make_required_gram_density_model]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:191-191`
- `callees` → [[analysis.scenario_builders__make_tabulated_flight_density_model|_make_tabulated_flight_density_model]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:179-179`
- `callees` → [[analysis.scenario_builders__make_time_tabulated_density_model|_make_time_tabulated_density_model]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:182-182`
- `callees` → [[environment.density_models_noatmospheremodel|NoAtmosphereModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:177-177`
- `callees` → [[environment.density_models_nrlmsise00atmospheremodel|NRLMSISE00AtmosphereModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:189-189`
<!-- vulcan:connections:end -->

## Limitations
An unknown model name is not rejected but treated as GRAM, so a typo in the manifest produces a GRAM run rather than an error.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 176.
