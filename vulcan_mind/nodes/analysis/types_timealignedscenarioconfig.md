---
id: analysis.types_timealignedscenarioconfig
label: TimeAlignedScenarioConfig
kind: struct
source:
  file: src/analysis/verification/telemetry_verification/types.jl
  symbol: TimeAlignedScenarioConfig
  lines:
  - 156
  - 156
inputs:
- id: name
  type: String
  units: n/a
  required: true
  description: Field `name`.
- id: planet_name
  type: String
  units: n/a
  required: true
  description: Field `planet_name`.
- id: telemetry_path
  type: String
  units: n/a
  required: true
  description: Field `telemetry_path`.
- id: telemetry_time_col
  type: String
  units: n/a
  required: true
  description: Field `telemetry_time_col`.
- id: telemetry_altitude_col
  type: String
  units: n/a
  required: true
  description: Field `telemetry_altitude_col`.
- id: telemetry_x_col
  type: String
  units: n/a
  required: true
  description: Field `telemetry_x_col`.
- id: telemetry_y_col
  type: String
  units: n/a
  required: true
  description: Field `telemetry_y_col`.
- id: telemetry_z_col
  type: String
  units: n/a
  required: true
  description: Field `telemetry_z_col`.
- id: telemetry_sma_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_sma_col` (default `nothing`).
- id: telemetry_ecc_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_ecc_col` (default `nothing`).
- id: telemetry_inc_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_inc_col` (default `nothing`).
- id: telemetry_aop_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_aop_col` (default `nothing`).
- id: telemetry_raan_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_raan_col` (default `nothing`).
- id: telemetry_ta_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_ta_col` (default `nothing`).
- id: telemetry_x_ic_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_x_ic_col` (default `nothing`).
- id: telemetry_y_ic_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_y_ic_col` (default `nothing`).
- id: telemetry_z_ic_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_z_ic_col` (default `nothing`).
- id: telemetry_vx_ic_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_vx_ic_col` (default `nothing`).
- id: telemetry_vy_ic_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_vy_ic_col` (default `nothing`).
- id: telemetry_vz_ic_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_vz_ic_col` (default `nothing`).
- id: telemetry_vx_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_vx_col` (default `nothing`).
- id: telemetry_vy_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_vy_col` (default `nothing`).
- id: telemetry_vz_col
  type: Union{Nothing, String}
  units: n/a
  required: false
  description: Field `telemetry_vz_col` (default `nothing`).
- id: max_points_quick
  type: Int
  units: n/a
  required: true
  description: Field `max_points_quick`.
- id: max_points_full
  type: Int
  units: n/a
  required: true
  description: Field `max_points_full`.
- id: min_eval_points
  type: Int
  units: n/a
  required: true
  description: Field `min_eval_points`.
- id: units_x
  type: String
  units: n/a
  required: true
  description: Field `units_x`.
- id: units_y
  type: Dict{String, String}
  units: n/a
  required: true
  description: Field `units_y`.
- id: tolerances_quick
  type: Dict{String, EventTolerance}
  units: n/a
  required: true
  description: Field `tolerances_quick`.
- id: tolerances_full
  type: Dict{String, EventTolerance}
  units: n/a
  required: true
  description: Field `tolerances_full`.
- id: initial_time
  type: InitialTime
  units: n/a
  required: true
  description: Field `initial_time`.
- id: spacecraft
  type: SpacecraftConfig
  units: n/a
  required: true
  description: Field `spacecraft`.
- id: gravity_model
  type: Symbol
  units: n/a
  required: true
  description: Field `gravity_model`.
- id: gravity_harmonics_degree
  type: Int
  units: n/a
  required: false
  description: Field `gravity_harmonics_degree` (default `0`).
- id: gravity_harmonics_order
  type: Int
  units: n/a
  required: false
  description: Field `gravity_harmonics_order` (default `0`).
- id: gravity_harmonics_file
  type: String
  units: n/a
  required: false
  description: Field `gravity_harmonics_file` (default `""`).
- id: nbody_bodies
  type: Vector{String}
  units: n/a
  required: false
  description: Field `nbody_bodies` (default `String[]`).
- id: srp_enabled
  type: Bool
  units: n/a
  required: false
  description: Field `srp_enabled` (default `false`).
- id: srp_cr
  type: Float64
  units: n/a
  required: false
  description: Field `srp_cr` (default `1.3`).
- id: srp_area_m2
  type: Float64
  units: n/a
  required: false
  description: Field `srp_area_m2` (default `0.0`).
- id: drag_enabled
  type: Bool
  units: n/a
  required: false
  description: Field `drag_enabled` (default `true`).
- id: aero_fixed_attitude_incidence
  type: Symbol
  units: n/a
  required: false
  description: Field `aero_fixed_attitude_incidence` (default `:max_drag`).
- id: include_wind
  type: Bool
  units: n/a
  required: false
  description: Field `include_wind` (default `false`).
- id: orbit_altitude_mode
  type: Symbol
  units: n/a
  required: false
  description: Field `orbit_altitude_mode` (default `:vacuum`).
- id: cartesian_ic_frame
  type: Symbol
  units: n/a
  required: false
  description: Field `cartesian_ic_frame` (default `:inertial`).
- id: comparison_frame
  type: Symbol
  units: n/a
  required: false
  description: Field `comparison_frame` (default `:inertial`).
- id: comparison_mode
  type: Symbol
  units: n/a
  required: false
  description: Field `comparison_mode` (default `:time_aligned_state`).
- id: extrema_min_separation_s
  type: Float64
  units: n/a
  required: false
  description: Field `extrema_min_separation_s` (default `500.0`).
- id: atmosphere_truth
  type: AtmosphereTruthConfig
  units: n/a
  required: false
  description: Field `atmosphere_truth` (default `AtmosphereTruthConfig()`).
- id: calibration
  type: CalibrationConfig
  units: n/a
  required: false
  description: Field `calibration` (default `CalibrationConfig()`).
- id: EI_km
  type: Float64
  units: n/a
  required: true
  description: Field `EI_km`.
- id: ic_offset_m
  type: NTuple{3, Float64}
  units: n/a
  required: false
  description: Field `ic_offset_m` (default `(0.0, 0.0, 0.0)`).
- id: ic_offset_mps
  type: NTuple{3, Float64}
  units: n/a
  required: false
  description: Field `ic_offset_mps` (default `(0.0, 0.0, 0.0)`).
- id: truth_mask
  type: Symbol
  units: n/a
  required: false
  description: Field `truth_mask` (default `:none`).
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
  type: TimeAlignedScenarioConfig
  units: n/a
  description: Constructed `TimeAlignedScenarioConfig` (keyword constructor via @kwdef).
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

# TimeAlignedScenarioConfig

## Purpose
Description of a time-aligned state-comparison telemetry verification scenario, where the simulation is compared pointwise in time against a telemetry table of altitude, position, and optionally velocity and orbital elements.

## Design & Implementation
A `Base.@kwdef struct` subtype of `AbstractScenarioConfig`. Required fields name the telemetry CSV and its columns (`telemetry_path`, `telemetry_time_col`, `telemetry_altitude_col`, `telemetry_x/y/z_col`), the quick/full `max_points_*`, `min_eval_points`, units, per-profile `EventTolerance` dictionaries, `initial_time`, `spacecraft`, `gravity_model`, and `EI_km`. Optional element columns (`telemetry_sma_col` through `telemetry_ta_col`), six Cartesian IC columns (`telemetry_*_ic_col`) that switch the initial condition to `CartesianInitialCondition` when all present, and three true velocity columns (`telemetry_vx/vy/vz_col`) all default to `nothing`. Frame and mode switches: `cartesian_ic_frame`, `comparison_frame` (both `:inertial`), `comparison_mode=:time_aligned_state`, `extrema_min_separation_s=500.0`. `ic_offset_m` and `ic_offset_mps` add constant J2000 offsets used by IC fitting, and `truth_mask` (`:none`, `:nightside`, `:dayside`) screens GNSS samples by illumination.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Field `name`. |
| in | `planet_name` | String | n/a | yes | Field `planet_name`. |
| in | `telemetry_path` | String | n/a | yes | Field `telemetry_path`. |
| in | `telemetry_time_col` | String | n/a | yes | Field `telemetry_time_col`. |
| in | `telemetry_altitude_col` | String | n/a | yes | Field `telemetry_altitude_col`. |
| in | `telemetry_x_col` | String | n/a | yes | Field `telemetry_x_col`. |
| in | `telemetry_y_col` | String | n/a | yes | Field `telemetry_y_col`. |
| in | `telemetry_z_col` | String | n/a | yes | Field `telemetry_z_col`. |
| in | `telemetry_sma_col` | Union{Nothing, String} | n/a | no | Field `telemetry_sma_col` (default `nothing`). |
| in | `telemetry_ecc_col` | Union{Nothing, String} | n/a | no | Field `telemetry_ecc_col` (default `nothing`). |
| in | `telemetry_inc_col` | Union{Nothing, String} | n/a | no | Field `telemetry_inc_col` (default `nothing`). |
| in | `telemetry_aop_col` | Union{Nothing, String} | n/a | no | Field `telemetry_aop_col` (default `nothing`). |
| in | `telemetry_raan_col` | Union{Nothing, String} | n/a | no | Field `telemetry_raan_col` (default `nothing`). |
| in | `telemetry_ta_col` | Union{Nothing, String} | n/a | no | Field `telemetry_ta_col` (default `nothing`). |
| in | `telemetry_x_ic_col` | Union{Nothing, String} | n/a | no | Field `telemetry_x_ic_col` (default `nothing`). |
| in | `telemetry_y_ic_col` | Union{Nothing, String} | n/a | no | Field `telemetry_y_ic_col` (default `nothing`). |
| in | `telemetry_z_ic_col` | Union{Nothing, String} | n/a | no | Field `telemetry_z_ic_col` (default `nothing`). |
| in | `telemetry_vx_ic_col` | Union{Nothing, String} | n/a | no | Field `telemetry_vx_ic_col` (default `nothing`). |
| in | `telemetry_vy_ic_col` | Union{Nothing, String} | n/a | no | Field `telemetry_vy_ic_col` (default `nothing`). |
| in | `telemetry_vz_ic_col` | Union{Nothing, String} | n/a | no | Field `telemetry_vz_ic_col` (default `nothing`). |
| in | `telemetry_vx_col` | Union{Nothing, String} | n/a | no | Field `telemetry_vx_col` (default `nothing`). |
| in | `telemetry_vy_col` | Union{Nothing, String} | n/a | no | Field `telemetry_vy_col` (default `nothing`). |
| in | `telemetry_vz_col` | Union{Nothing, String} | n/a | no | Field `telemetry_vz_col` (default `nothing`). |
| in | `max_points_quick` | Int | n/a | yes | Field `max_points_quick`. |
| in | `max_points_full` | Int | n/a | yes | Field `max_points_full`. |
| in | `min_eval_points` | Int | n/a | yes | Field `min_eval_points`. |
| in | `units_x` | String | n/a | yes | Field `units_x`. |
| in | `units_y` | Dict{String, String} | n/a | yes | Field `units_y`. |
| in | `tolerances_quick` | Dict{String, EventTolerance} | n/a | yes | Field `tolerances_quick`. |
| in | `tolerances_full` | Dict{String, EventTolerance} | n/a | yes | Field `tolerances_full`. |
| in | `initial_time` | InitialTime | n/a | yes | Field `initial_time`. |
| in | `spacecraft` | SpacecraftConfig | n/a | yes | Field `spacecraft`. |
| in | `gravity_model` | Symbol | n/a | yes | Field `gravity_model`. |
| in | `gravity_harmonics_degree` | Int | n/a | no | Field `gravity_harmonics_degree` (default `0`). |
| in | `gravity_harmonics_order` | Int | n/a | no | Field `gravity_harmonics_order` (default `0`). |
| in | `gravity_harmonics_file` | String | n/a | no | Field `gravity_harmonics_file` (default `""`). |
| in | `nbody_bodies` | Vector{String} | n/a | no | Field `nbody_bodies` (default `String[]`). |
| in | `srp_enabled` | Bool | n/a | no | Field `srp_enabled` (default `false`). |
| in | `srp_cr` | Float64 | n/a | no | Field `srp_cr` (default `1.3`). |
| in | `srp_area_m2` | Float64 | n/a | no | Field `srp_area_m2` (default `0.0`). |
| in | `drag_enabled` | Bool | n/a | no | Field `drag_enabled` (default `true`). |
| in | `aero_fixed_attitude_incidence` | Symbol | n/a | no | Field `aero_fixed_attitude_incidence` (default `:max_drag`). |
| in | `include_wind` | Bool | n/a | no | Field `include_wind` (default `false`). |
| in | `orbit_altitude_mode` | Symbol | n/a | no | Field `orbit_altitude_mode` (default `:vacuum`). |
| in | `cartesian_ic_frame` | Symbol | n/a | no | Field `cartesian_ic_frame` (default `:inertial`). |
| in | `comparison_frame` | Symbol | n/a | no | Field `comparison_frame` (default `:inertial`). |
| in | `comparison_mode` | Symbol | n/a | no | Field `comparison_mode` (default `:time_aligned_state`). |
| in | `extrema_min_separation_s` | Float64 | n/a | no | Field `extrema_min_separation_s` (default `500.0`). |
| in | `atmosphere_truth` | AtmosphereTruthConfig | n/a | no | Field `atmosphere_truth` (default `AtmosphereTruthConfig()`). |
| in | `calibration` | CalibrationConfig | n/a | no | Field `calibration` (default `CalibrationConfig()`). |
| in | `EI_km` | Float64 | n/a | yes | Field `EI_km`. |
| in | `ic_offset_m` | NTuple{3, Float64} | n/a | no | Field `ic_offset_m` (default `(0.0, 0.0, 0.0)`). |
| in | `ic_offset_mps` | NTuple{3, Float64} | n/a | no | Field `ic_offset_mps` (default `(0.0, 0.0, 0.0)`). |
| in | `truth_mask` | Symbol | n/a | no | Field `truth_mask` (default `:none`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | TimeAlignedScenarioConfig | n/a | — | Constructed `TimeAlignedScenarioConfig` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:641-641`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`

**Downstream**

- `callees` → [[analysis.types_atmospheretruthconfig|AtmosphereTruthConfig]] · `callers` · call · `src/analysis/verification/telemetry_verification/types.jl:214-214`
- `callees` → [[analysis.types_calibrationconfig|CalibrationConfig]] · `callers` · call · `src/analysis/verification/telemetry_verification/types.jl:215-215`
<!-- vulcan:connections:end -->

## Limitations
The all-six-present rule for Cartesian IC columns is a convention enforced by the loader, not the struct; a partial set is silently ignored. Column names are strings with no existence check until the CSV is read. `truth_mask` and the frame symbols are unvalidated. The struct duplicates the force-model fields of `OrbitEventsScenarioConfig` rather than sharing a nested config, so the two must be kept in sync manually.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/types.jl` line 156.
