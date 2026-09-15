---
id: analysis.scenario_builders__make_tabulated_flight_density_model
label: _make_tabulated_flight_density_model
kind: function
source:
  file: src/analysis/verification/telemetry_verification/scenario_builders.jl
  symbol: _make_tabulated_flight_density_model
  lines:
  - 199
  - 199
inputs:
- id: initial_time
  type: InitialTime
  units: n/a
  required: true
  description: Positional argument `initial_time`.
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
  type: SimulationModel.TabulatedFlightAtmosphereModel
  units: n/a
  description: Return value of `_make_tabulated_flight_density_model`. Returns `SimulationModel.TabulatedFlightAtmosphereModel(`.
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

# _make_tabulated_flight_density_model

## Purpose
Builds a density model from a flight-derived table of per-pass, per-leg, one-kilometre altitude bins, so a scenario can fly through the atmosphere the spacecraft actually measured.

## Design & Implementation
Loads the Arrow file into a `DataFrame` and requires the columns `P`, `leg`, `alt_km`, `rho_kgm3`, `sigma_kgm3` and `t_peri_utc`. It forms the epoch from the scenario `initial_time` to the millisecond, then for each pass id computes periapsis elapsed seconds from the first row's `t_peri_utc` and splits rows into inbound and outbound legs, storing altitude in metres, log-density, and sigma divided by density as the log-density sigma, each leg sorted by altitude. Passes are sorted by periapsis time, archive gaps are counted and printed, and a `TabulatedFlightAtmosphereModel` is constructed with the sigma scale and the literal constants 3.4 and 188.92.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `initial_time` | InitialTime | n/a | yes | Positional argument `initial_time`. |
| in | `truth` | AtmosphereTruthConfig | n/a | yes | Positional argument `truth`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SimulationModel.TabulatedFlightAtmosphereModel | n/a | — | Return value of `_make_tabulated_flight_density_model`. Returns `SimulationModel.TabulatedFlightAtmosphereModel(`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.scenario_builders__scenario_density_model|_scenario_density_model]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:179-179`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl`

**Downstream**

- `callees` → [[analysis.ic_fit_float64|Float64]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:212-212`
- `callees` → [[cli.spaceagora_cli_println|println]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:239-239`
- `callees` → [[environment.density_models_tabulatedflightatmospheremodel|TabulatedFlightAtmosphereModel]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:241-241`
- `callees` → [[gnc.pso_parameters_push_bang|push!]] · `callers` · call · `src/analysis/verification/telemetry_verification/scenario_builders.jl:221-221`
<!-- vulcan:connections:end -->

## Limitations
Rows with non-positive density are dropped silently; `t_peri_utc` is parsed from its first nineteen characters so any timezone suffix is ignored; and the two trailing constants are undocumented literals passed positionally.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/scenario_builders.jl` line 199.
