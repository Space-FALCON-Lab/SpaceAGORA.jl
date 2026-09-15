---
id: analysis.scenario_builders__make_time_tabulated_density_model
label: _make_time_tabulated_density_model
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _make_time_tabulated_density_model
  lines:
  - 251
  - 251
inputs:
- id: truth
  type: AtmosphereTruthConfig
  units: n/a
  required: true
  description: Positional argument `truth`.
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
  type: SimulationModel.TimeTabulatedAtmosphereModel
  units: n/a
  description: Return value of `_make_time_tabulated_density_model`. Returns `SimulationModel.TimeTabulatedAtmosphereModel(`.
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

# _make_time_tabulated_density_model

## Purpose
Builds a density model from a plain time series of density, for scenarios whose atmosphere truth is a published rho-versus-time product.

## Design & Implementation
Reads CSV or Arrow depending on the file extension, requires `time_s` and `rho_kgm3` columns, sorts by time, prints the node count and span in hours, and constructs `TimeTabulatedAtmosphereModel` with the configured `scale` and `temperature_k`.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `truth` | AtmosphereTruthConfig | n/a | yes | Positional argument `truth`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.TimeTabulatedAtmosphereModel | n/a | — | Return value of `_make_time_tabulated_density_model`. Returns `SimulationModel.TimeTabulatedAtmosphereModel(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_density_model|_scenario_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:182-182`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:262-262`
- `callees` → [[environment.density_models_timetabulatedatmospheremodel|TimeTabulatedAtmosphereModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:265-265`
<!-- vulcan:connections:end -->

## Limitations
Times are scenario elapsed seconds, so the table's epoch must equal the scenario `initial_time` and nothing here can verify that; the temperature is a single constant for the whole table.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 251.
