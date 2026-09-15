---
id: analysis.error_tables__time_aligned_rows_errors
label: _time_aligned_rows_errors
kind: function
source:
  file: src/analysis/verification/telemetry_verification/error_tables.jl
  symbol: _time_aligned_rows_errors
  lines:
  - 76
  - 76
inputs:
- id: cfg
  type: TimeAlignedScenarioConfig
  units: n/a
  required: true
  description: Positional argument `cfg`.
- id: args
  type: SimulationConfiguration
  units: n/a
  required: true
  description: Positional argument `args`.
- id: results_df
  type: DataFrame
  units: n/a
  required: true
  description: Positional argument `results_df`.
- id: telemetry
  type: Any
  units: n/a
  required: true
  description: Positional argument `telemetry`.
- id: bias_by_event
  type: Dict{String, Float64}
  units: n/a
  required: false
  description: Keyword argument `bias_by_event` (default `Dict{String, Float64}()`).
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
  type: AbstractArray
  units: n/a
  description: Return value of `_time_aligned_rows_errors`. Returns `[peri_summary,
    apo_summary, peri_speed_summary, apo_speed_summary], [peri_errors` or `summaries,
    errors`.
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

# _time_aligned_rows_errors

## Purpose
Produces the per-channel error summaries and residual tables for a time-aligned verification scenario, comparing a simulation results frame against loaded reference telemetry.

## Design & Implementation
Two distinct paths. When `cfg.comparison_mode` is `:orbit_events` it derives periapsis and apoapsis extrema from both telemetry and simulation, then scores altitude and speed at each apsis through `_compare_orbit_curve`, adding an apoapsis decay diagnostic once at least three apoapses exist on both sides. Otherwise it pulls position and velocity columns by alternative names through `_require_column`, optionally rotates the whole history into the planet-fixed frame sample by sample when `comparison_frame` is `:planet_fixed`, applies the per-event biases from `bias_by_event`, and scores altitude and the three position components through `_compare_time_series`. Velocity channels are appended only when all three telemetry velocity column names are configured.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `cfg` | TimeAlignedScenarioConfig | n/a | yes | Positional argument `cfg`. |
| in | `args` | SimulationConfiguration | n/a | yes | Positional argument `args`. |
| in | `results_df` | DataFrame | n/a | yes | Positional argument `results_df`. |
| in | `telemetry` | Any | n/a | yes | Positional argument `telemetry`. |
| in | `bias_by_event` | Dict{String, Float64} | n/a | no | Keyword argument `bias_by_event` (default `Dict{String, Float64}()`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AbstractArray | n/a | — | Return value of `_time_aligned_rows_errors`. Returns `[peri_summary, apo_summary, peri_speed_summary, apo_speed_summary], [peri_errors` or `summaries, errors`. |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.runner__run_single_scenario|_run_single_scenario]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/runner.jl:257-257`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/error_tables.jl`

**Downstream**

- `callees` → [[analysis.comparison_metrics__apo_decay_diagnostic|_apo_decay_diagnostic]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:112-112`
- `callees` → [[analysis.comparison_metrics__compare_time_series|_compare_time_series]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:179-179`
- `callees` → [[analysis.error_tables__telemetry_altitude_km|_telemetry_altitude_km]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:170-170`
- `callees` → [[analysis.telemetry_loading__extract_extrema_from_time_aligned_telemetry|_extract_extrema_from_time_aligned_telemetry]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:85-85`
- `callees` → [[analysis.telemetry_loading__extract_extrema_series|_extract_extrema_series]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:91-91`
- `callees` → [[analysis.telemetry_loading__initial_time_et|_initial_time_et]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:148-148`
- `callees` → [[analysis.telemetry_loading__j2000_to_planet_fixed_state|_j2000_to_planet_fixed_state]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:150-150`
- `callees` → [[analysis.telemetry_loading__require_column|_require_column]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:140-140`
- `callees` → [[analysis.telemetry_loading__to_float_vector|_to_float_vector]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:139-139`
- `callees` → [[envana.ana_comparison_metrics_compare_orbit_curve|_compare_orbit_curve]] · `callers` · call · `src/analysis/verification/telemetry_verification/error_tables.jl:97-97`
<!-- vulcan:connections:end -->

## Limitations
The planet-fixed rotation mutates `sim_x_m` through `sim_vz_mps` in place, so the extracted column vectors must not be aliases of the caller's DataFrame columns; biases are additive constants only, so a drifting offset is not representable.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/error_tables.jl` line 76.
