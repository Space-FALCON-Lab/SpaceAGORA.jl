---
id: analysis.types_orbiteventsscenarioconfig
label: OrbitEventsScenarioConfig
kind: struct
source:
  file: src/analysis/verification/telemetry_verification/types.jl
  symbol: OrbitEventsScenarioConfig
  lines:
  - 90
  - 90
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
- id: telemetry_peri_path
  type: String
  units: n/a
  required: true
  description: Field `telemetry_peri_path`.
- id: telemetry_apo_path
  type: String
  units: n/a
  required: true
  description: Field `telemetry_apo_path`.
- id: target_orbits_quick
  type: Int
  units: n/a
  required: true
  description: Field `target_orbits_quick`.
- id: target_orbits_full
  type: Int
  units: n/a
  required: true
  description: Field `target_orbits_full`.
- id: compare_points_quick
  type: Int
  units: n/a
  required: true
  description: Field `compare_points_quick`.
- id: compare_points_full
  type: Int
  units: n/a
  required: true
  description: Field `compare_points_full`.
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
- id: ra_m
  type: Float64
  units: n/a
  required: true
  description: Field `ra_m`.
- id: rp_altitude_m
  type: Float64
  units: n/a
  required: true
  description: Field `rp_altitude_m`.
- id: i_deg
  type: Float64
  units: n/a
  required: true
  description: Field `i_deg`.
- id: aop_deg
  type: Float64
  units: n/a
  required: true
  description: Field `aop_deg`.
- id: raan_deg
  type: Float64
  units: n/a
  required: true
  description: Field `raan_deg`.
- id: ta_deg
  type: Float64
  units: n/a
  required: true
  description: Field `ta_deg`.
- id: element_frame
  type: Symbol
  units: n/a
  required: false
  description: Field `element_frame` (default `:j2000`).
- id: initial_state_j2000_m
  type: Union{Nothing, NTuple{6, Float64}}
  units: n/a
  required: false
  description: Field `initial_state_j2000_m` (default `nothing`).
- id: epoch_orbit_offset
  type: Union{Nothing, Float64}
  units: n/a
  required: false
  description: Field `epoch_orbit_offset` (default `nothing`).
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
- id: maneuver_orbit_numbers
  type: Vector{Int64}
  units: n/a
  required: false
  description: Field `maneuver_orbit_numbers` (default `Int64[]`).
- id: maneuver_orbit_numbers_campaign
  type: Vector{Int64}
  units: n/a
  required: false
  description: Field `maneuver_orbit_numbers_campaign` (default `Int64[]`).
- id: maneuver_delta_v_mps
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `maneuver_delta_v_mps` (default `Float64[]`).
- id: maneuver_replay_scale_mode
  type: String
  units: n/a
  required: false
  description: Field `maneuver_replay_scale_mode` (default `"delta_v"`).
- id: maneuver_flight_apoapsis_alt_m
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `maneuver_flight_apoapsis_alt_m` (default `Float64[]`).
- id: maneuver_thrust_n
  type: Float64
  units: n/a
  required: false
  description: Field `maneuver_thrust_n` (default `0.0`).
- id: maneuver_isp_s
  type: Float64
  units: n/a
  required: false
  description: Field `maneuver_isp_s` (default `0.0`).
- id: maneuver_guidance_rate_s
  type: Float64
  units: n/a
  required: false
  description: Field `maneuver_guidance_rate_s` (default `30.0`).
- id: maneuver_control_rate_s
  type: Float64
  units: n/a
  required: false
  description: Field `maneuver_control_rate_s` (default `10.0`).
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
  type: OrbitEventsScenarioConfig
  units: n/a
  description: Constructed `OrbitEventsScenarioConfig` (keyword constructor via @kwdef).
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

# OrbitEventsScenarioConfig

## Purpose
Complete description of an apsis-event telemetry verification scenario: telemetry file paths for periapsis and apoapsis products, orbital-element or Cartesian initial condition, force-model switches, maneuver replay schedule, atmosphere truth, calibration settings, and per-profile tolerances.

## Design & Implementation
A `Base.@kwdef struct` subtype of `AbstractScenarioConfig`. Required fields (no defaults) include `name`, `planet_name`, `telemetry_peri_path`, `telemetry_apo_path`, quick/full `target_orbits_*` and `compare_points_*`, `min_eval_points`, `tolerances_quick`/`tolerances_full::Dict{String,EventTolerance}` (each an `EventTolerance` of `max_abs_km`, `max_nmae`, `max_rmse_km`), `initial_time::InitialTime`, the six elements `ra_m`, `rp_altitude_m`, `i_deg`, `aop_deg`, `raan_deg`, `ta_deg`, `spacecraft`, `gravity_model`, and `EI_km`. `initial_state_j2000_m::NTuple{6}` overrides the elements when present; `epoch_orbit_offset` aligns sim apsis numbering to the truth product. Maneuver replay is expressed through parallel vectors `maneuver_orbit_numbers`, `maneuver_delta_v_mps`, `maneuver_flight_apoapsis_alt_m`, plus `maneuver_replay_scale_mode` (`"delta_v"` or `"flight_apoapsis_ratio"`), thrust, Isp, and guidance/control rates.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `name` | String | n/a | yes | Field `name`. |
| in | `planet_name` | String | n/a | yes | Field `planet_name`. |
| in | `telemetry_peri_path` | String | n/a | yes | Field `telemetry_peri_path`. |
| in | `telemetry_apo_path` | String | n/a | yes | Field `telemetry_apo_path`. |
| in | `target_orbits_quick` | Int | n/a | yes | Field `target_orbits_quick`. |
| in | `target_orbits_full` | Int | n/a | yes | Field `target_orbits_full`. |
| in | `compare_points_quick` | Int | n/a | yes | Field `compare_points_quick`. |
| in | `compare_points_full` | Int | n/a | yes | Field `compare_points_full`. |
| in | `min_eval_points` | Int | n/a | yes | Field `min_eval_points`. |
| in | `units_x` | String | n/a | yes | Field `units_x`. |
| in | `units_y` | Dict{String, String} | n/a | yes | Field `units_y`. |
| in | `tolerances_quick` | Dict{String, EventTolerance} | n/a | yes | Field `tolerances_quick`. |
| in | `tolerances_full` | Dict{String, EventTolerance} | n/a | yes | Field `tolerances_full`. |
| in | `initial_time` | InitialTime | n/a | yes | Field `initial_time`. |
| in | `ra_m` | Float64 | n/a | yes | Field `ra_m`. |
| in | `rp_altitude_m` | Float64 | n/a | yes | Field `rp_altitude_m`. |
| in | `i_deg` | Float64 | n/a | yes | Field `i_deg`. |
| in | `aop_deg` | Float64 | n/a | yes | Field `aop_deg`. |
| in | `raan_deg` | Float64 | n/a | yes | Field `raan_deg`. |
| in | `ta_deg` | Float64 | n/a | yes | Field `ta_deg`. |
| in | `element_frame` | Symbol | n/a | no | Field `element_frame` (default `:j2000`). |
| in | `initial_state_j2000_m` | Union{Nothing, NTuple{6, Float64}} | n/a | no | Field `initial_state_j2000_m` (default `nothing`). |
| in | `epoch_orbit_offset` | Union{Nothing, Float64} | n/a | no | Field `epoch_orbit_offset` (default `nothing`). |
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
| in | `maneuver_orbit_numbers` | Vector{Int64} | n/a | no | Field `maneuver_orbit_numbers` (default `Int64[]`). |
| in | `maneuver_orbit_numbers_campaign` | Vector{Int64} | n/a | no | Field `maneuver_orbit_numbers_campaign` (default `Int64[]`). |
| in | `maneuver_delta_v_mps` | Vector{Float64} | n/a | no | Field `maneuver_delta_v_mps` (default `Float64[]`). |
| in | `maneuver_replay_scale_mode` | String | n/a | no | Field `maneuver_replay_scale_mode` (default `"delta_v"`). |
| in | `maneuver_flight_apoapsis_alt_m` | Vector{Float64} | n/a | no | Field `maneuver_flight_apoapsis_alt_m` (default `Float64[]`). |
| in | `maneuver_thrust_n` | Float64 | n/a | no | Field `maneuver_thrust_n` (default `0.0`). |
| in | `maneuver_isp_s` | Float64 | n/a | no | Field `maneuver_isp_s` (default `0.0`). |
| in | `maneuver_guidance_rate_s` | Float64 | n/a | no | Field `maneuver_guidance_rate_s` (default `30.0`). |
| in | `maneuver_control_rate_s` | Float64 | n/a | no | Field `maneuver_control_rate_s` (default `10.0`). |
| in | `atmosphere_truth` | AtmosphereTruthConfig | n/a | no | Field `atmosphere_truth` (default `AtmosphereTruthConfig()`). |
| in | `calibration` | CalibrationConfig | n/a | no | Field `calibration` (default `CalibrationConfig()`). |
| in | `EI_km` | Float64 | n/a | yes | Field `EI_km`. |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | OrbitEventsScenarioConfig | n/a | — | Constructed `OrbitEventsScenarioConfig` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[analysis.manifest_parsing__load_scenarios_from_manifest|_load_scenarios_from_manifest]] · `callees` → `callers` · call · `src/analysis/verification/telemetry_verification/manifest_parsing.jl:580-580`
- [[module.analysis|TelemetryVerification]] · `api` → `module_api` · call · `src/analysis/verification/telemetry_verification/types.jl`

**Downstream**

- `callees` → [[analysis.types_atmospheretruthconfig|AtmosphereTruthConfig]] · `callers` · call · `src/analysis/verification/telemetry_verification/types.jl:151-151`
- `callees` → [[analysis.types_calibrationconfig|CalibrationConfig]] · `callers` · call · `src/analysis/verification/telemetry_verification/types.jl:152-152`
<!-- vulcan:connections:end -->

## Limitations
Parallel maneuver vectors are not length-checked at construction; a mismatch surfaces only during replay. Mixed units coexist (`ra_m`, `rp_altitude_m` in metres; `EI_km` in kilometres; angles in degrees). `element_frame` and `orbit_altitude_mode` are unvalidated Symbols. `gravity_harmonics_degree/order` default to 0 while `gravity_model` is required, so an inconsistent combination is representable. The struct is large and immutable, so tweaking a single field requires rebuilding the whole object.

## Provenance
Mapped from `src/analysis/verification/telemetry_verification/types.jl` line 90.
