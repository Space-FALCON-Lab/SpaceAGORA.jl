---
id: analysis.types_atmospheretruthconfig
label: AtmosphereTruthConfig
kind: struct
source:
  file: src/analysis/verification/telemetry_verification/types.jl
  symbol: AtmosphereTruthConfig
  lines:
  - 6
  - 6
inputs:
- id: assumption_id
  type: String
  units: n/a
  required: false
  description: Field `assumption_id` (default `"gram_default"`).
- id: atmosphere_model
  type: String
  units: n/a
  required: false
  description: Field `atmosphere_model` (default `"GRAM"`).
- id: atmosphere_dataset
  type: String
  units: n/a
  required: false
  description: Field `atmosphere_dataset` (default `"default"`).
- id: space_weather_model
  type: String
  units: n/a
  required: false
  description: Field `space_weather_model` (default `"default"`).
- id: solar_flux_model
  type: String
  units: n/a
  required: false
  description: Field `solar_flux_model` (default `"default"`).
- id: gram_seed
  type: Int
  units: n/a
  required: false
  description: Field `gram_seed` (default `1001`).
- id: gram_perturbation_scales
  type: NTuple{4, Float64}
  units: n/a
  required: false
  description: Field `gram_perturbation_scales` (default `(0.0, 0.0, 0.0, 0.0)`).
- id: gram_min_relative_step_size
  type: Union{Nothing, Float64}
  units: n/a
  required: false
  description: Field `gram_min_relative_step_size` (default `nothing`).
- id: gram_offline_surrogate
  type: String
  units: n/a
  required: false
  description: Field `gram_offline_surrogate` (default `"off"`).
- id: gram_static_grid
  type: Bool
  units: n/a
  required: false
  description: Field `gram_static_grid` (default `false`).
- id: gram_track_cache
  type: Bool
  units: n/a
  required: false
  description: Field `gram_track_cache` (default `false`).
- id: gram_global_lock
  type: String
  units: n/a
  required: false
  description: Field `gram_global_lock` (default `"on"`).
- id: mars_map_year
  type: Union{Nothing, Int}
  units: n/a
  required: false
  description: Field `mars_map_year` (default `nothing`).
- id: mars_mgcm_dust_levels
  type: Union{Nothing, NTuple{3, Float64}}
  units: n/a
  required: false
  description: Field `mars_mgcm_dust_levels` (default `nothing`).
- id: mars_dust_storm
  type: Union{Nothing, NTuple{6, Float64}}
  units: n/a
  required: false
  description: Field `mars_dust_storm` (default `nothing`).
- id: mars_f107
  type: Union{Nothing, Float64}
  units: n/a
  required: false
  description: Field `mars_f107` (default `nothing`).
- id: mars_wind_scales
  type: Union{Nothing, NTuple{2, Float64}}
  units: n/a
  required: false
  description: Field `mars_wind_scales` (default `nothing`).
- id: mars_mola_heights
  type: Union{Nothing, Bool}
  units: n/a
  required: false
  description: Field `mars_mola_heights` (default `nothing`).
- id: mars_min_max
  type: Union{Nothing, Int}
  units: n/a
  required: false
  description: Field `mars_min_max` (default `nothing`).
- id: tabulated_flight_file
  type: String
  units: n/a
  required: false
  description: Field `tabulated_flight_file` (default `""`).
- id: tabulated_flight_sigma
  type: Float64
  units: n/a
  required: false
  description: Field `tabulated_flight_sigma` (default `0.0`).
- id: tabulated_time_file
  type: String
  units: n/a
  required: false
  description: Field `tabulated_time_file` (default `""`).
- id: tabulated_time_scale
  type: Float64
  units: n/a
  required: false
  description: Field `tabulated_time_scale` (default `1.0`).
- id: tabulated_time_temperature_k
  type: Float64
  units: n/a
  required: false
  description: Field `tabulated_time_temperature_k` (default `900.0`).
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
  type: AtmosphereTruthConfig
  units: n/a
  description: Constructed `AtmosphereTruthConfig` (keyword constructor via @kwdef).
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

# AtmosphereTruthConfig

## Purpose
Immutable keyword struct describing the atmospheric truth assumption used by a telemetry verification scenario: which model (`GRAM`, `tabulated_flight`, `tabulated_time`, or another family), GRAM seeding and perturbation settings, Mars-GRAM options, and file-based tabulated density inputs.

## Design & Implementation
Built with `Base.@kwdef` so every field has a default and the struct can be constructed from a manifest with partial overrides. Key fields: `assumption_id="gram_default"`, `atmosphere_model="GRAM"`, `gram_seed=1001`, `gram_perturbation_scales::NTuple{4,Float64}` (all zero), `gram_offline_surrogate="off"`, `gram_static_grid`, `gram_track_cache`, and `gram_global_lock="on"` as strings or Bools. Mars-specific options (`mars_map_year`, `mars_mgcm_dust_levels::NTuple{3}`, `mars_dust_storm::NTuple{6}`, `mars_f107`, `mars_wind_scales`, `mars_mola_heights`, `mars_min_max`) default to `nothing`. `tabulated_flight_file`/`tabulated_flight_sigma` and `tabulated_time_file`/`tabulated_time_scale`/`tabulated_time_temperature_k=900.0` support the flight-measured and time-tabulated density replay modes.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `assumption_id` | String | n/a | no | Field `assumption_id` (default `"gram_default"`). |
| in | `atmosphere_model` | String | n/a | no | Field `atmosphere_model` (default `"GRAM"`). |
| in | `atmosphere_dataset` | String | n/a | no | Field `atmosphere_dataset` (default `"default"`). |
| in | `space_weather_model` | String | n/a | no | Field `space_weather_model` (default `"default"`). |
| in | `solar_flux_model` | String | n/a | no | Field `solar_flux_model` (default `"default"`). |
| in | `gram_seed` | Int | n/a | no | Field `gram_seed` (default `1001`). |
| in | `gram_perturbation_scales` | NTuple{4, Float64} | n/a | no | Field `gram_perturbation_scales` (default `(0.0, 0.0, 0.0, 0.0)`). |
| in | `gram_min_relative_step_size` | Union{Nothing, Float64} | n/a | no | Field `gram_min_relative_step_size` (default `nothing`). |
| in | `gram_offline_surrogate` | String | n/a | no | Field `gram_offline_surrogate` (default `"off"`). |
| in | `gram_static_grid` | Bool | n/a | no | Field `gram_static_grid` (default `false`). |
| in | `gram_track_cache` | Bool | n/a | no | Field `gram_track_cache` (default `false`). |
| in | `gram_global_lock` | String | n/a | no | Field `gram_global_lock` (default `"on"`). |
| in | `mars_map_year` | Union{Nothing, Int} | n/a | no | Field `mars_map_year` (default `nothing`). |
| in | `mars_mgcm_dust_levels` | Union{Nothing, NTuple{3, Float64}} | n/a | no | Field `mars_mgcm_dust_levels` (default `nothing`). |
| in | `mars_dust_storm` | Union{Nothing, NTuple{6, Float64}} | n/a | no | Field `mars_dust_storm` (default `nothing`). |
| in | `mars_f107` | Union{Nothing, Float64} | n/a | no | Field `mars_f107` (default `nothing`). |
| in | `mars_wind_scales` | Union{Nothing, NTuple{2, Float64}} | n/a | no | Field `mars_wind_scales` (default `nothing`). |
| in | `mars_mola_heights` | Union{Nothing, Bool} | n/a | no | Field `mars_mola_heights` (default `nothing`). |
| in | `mars_min_max` | Union{Nothing, Int} | n/a | no | Field `mars_min_max` (default `nothing`). |
| in | `tabulated_flight_file` | String | n/a | no | Field `tabulated_flight_file` (default `""`). |
| in | `tabulated_flight_sigma` | Float64 | n/a | no | Field `tabulated_flight_sigma` (default `0.0`). |
| in | `tabulated_time_file` | String | n/a | no | Field `tabulated_time_file` (default `""`). |
| in | `tabulated_time_scale` | Float64 | n/a | no | Field `tabulated_time_scale` (default `1.0`). |
| in | `tabulated_time_temperature_k` | Float64 | n/a | no | Field `tabulated_time_temperature_k` (default `900.0`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | AtmosphereTruthConfig | n/a | — | Constructed `AtmosphereTruthConfig` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__parse_atmosphere_truth_config|_parse_atmosphere_truth_config]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:299-299`
- [[analysis.types_orbiteventsscenarioconfig|OrbitEventsScenarioConfig]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/types.jl:151-151`
- [[analysis.types_timealignedscenarioconfig|TimeAlignedScenarioConfig]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/types.jl:214-214`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`

**Downstream**

- *none*
<!-- vulcan:connections:end -->

## Limitations
Several tri-state options are encoded as strings (`gram_offline_surrogate`, `gram_global_lock`) with no enumeration check at construction, so typos are only caught by consumers. The `tabulated_*` file fields are empty strings rather than `nothing`, so presence must be tested with `isempty`. Units are implied by field names only (`tabulated_time_temperature_k` in kelvin, `gram_min_relative_step_size` dimensionless). Nothing enforces that `atmosphere_model` and the populated option group are mutually consistent.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/types.jl` line 6.
