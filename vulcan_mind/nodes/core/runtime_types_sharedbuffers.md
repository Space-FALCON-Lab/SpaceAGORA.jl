---
id: core.runtime_types_sharedbuffers
label: SharedBuffers
kind: struct
source:
  file: src/core/types/runtime_types.jl
  symbol: SharedBuffers
  lines:
  - 711
  - 711
inputs:
- id: n_sats
  type: Int
  units: n/a
  required: true
  description: Field `n_sats`.
- id: densities
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `densities` (default `zeros(Float64, n_sats)`).
- id: temperatures
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `temperatures` (default `ones(Float64, n_sats)`).
- id: winds
  type: Vector{SVector{3,Float64}}
  units: n/a
  required: false
  description: Field `winds` (default `[SVector{3,Float64}(0.0, 0.0, 0.0) for _ in
    1:n_sats]`).
- id: density_sample_t
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `density_sample_t` (default `fill(NaN, n_sats)`).
- id: density_batch_altitudes
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `density_batch_altitudes` (default `zeros(Float64, n_sats)`).
- id: density_batch_latitudes
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `density_batch_latitudes` (default `zeros(Float64, n_sats)`).
- id: density_batch_longitudes
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `density_batch_longitudes` (default `zeros(Float64, n_sats)`).
- id: heat_rates
  type: Vector{Vector{Float64}}
  units: n/a
  required: false
  description: Field `heat_rates` (default `[Float64[] for _ in 1:n_sats]`).
- id: density_models
  type: Vector{_PerSatDensityModel}
  units: n/a
  required: false
  description: Field `density_models` (default `_PerSatDensityModel[]`).
- id: gram_density_cache
  type: Vector{Union{Nothing, GramTrackCache}}
  units: n/a
  required: false
  description: Field `gram_density_cache` (default `_typed_nothing_vector(GramTrackCache,
    n_sats)`).
- id: vacuum_gram_caches
  type: Vector{Union{Nothing, VacuumPredictedGRAMCache}}
  units: n/a
  required: false
  description: Field `vacuum_gram_caches` (default `_typed_nothing_vector(VacuumPredictedGRAMCache,
    n_sats)`).
- id: gram_isolated_pool_models
  type: Vector{GRAMAtmosphereModel}
  units: n/a
  required: false
  description: Field `gram_isolated_pool_models` (default `GRAMAtmosphereModel[]`).
- id: gram_isolated_pool_locks
  type: Vector{ReentrantLock}
  units: n/a
  required: false
  description: Field `gram_isolated_pool_locks` (default `ReentrantLock[]`).
- id: harmonics_workspaces
  type: Vector{Union{Nothing, _HarmonicsWorkspaceMap}}
  units: n/a
  required: false
  description: Field `harmonics_workspaces` (default `_typed_nothing_vector(_HarmonicsWorkspaceMap,
    n_sats)`).
- id: nbody_workspaces
  type: Vector{Union{Nothing, NBodyScratchWorkspace}}
  units: n/a
  required: false
  description: Field `nbody_workspaces` (default `_typed_nothing_vector(NBodyScratchWorkspace,
    n_sats)`).
- id: aero_workspaces
  type: Vector{Union{Nothing, AeroScratchWorkspace}}
  units: n/a
  required: false
  description: Field `aero_workspaces` (default `_typed_nothing_vector(AeroScratchWorkspace,
    n_sats)`).
- id: nbody_ephemeris_cache
  type: Base.RefValue{Union{Nothing, NBodyEphemerisCache}}
  units: n/a
  required: false
  description: Field `nbody_ephemeris_cache` (default `Ref{Union{Nothing, NBodyEphemerisCache}}(nothing)`).
- id: srp_sun_ephemeris_cache
  type: Base.RefValue{Union{Nothing, SRPSunEphemerisCache}}
  units: n/a
  required: false
  description: Field `srp_sun_ephemeris_cache` (default `Ref{Union{Nothing, SRPSunEphemerisCache}}(nothing)`).
- id: planet_frame_ephemeris_cache
  type: Base.RefValue{Union{Nothing, PlanetFrameEphemerisCache}}
  units: n/a
  required: false
  description: Field `planet_frame_ephemeris_cache` (default `Ref{Union{Nothing, PlanetFrameEphemerisCache}}(nothing)`).
- id: harmonics_lpi_lock
  type: ReentrantLock
  units: n/a
  required: false
  description: Field `harmonics_lpi_lock` (default `ReentrantLock()`).
- id: harmonics_lpi_key
  type: Base.RefValue{Any}
  units: n/a
  required: false
  description: Field `harmonics_lpi_key` (default `Ref{Any}(nothing)`).
- id: harmonics_lpi
  type: Base.RefValue{SMatrix{3,3,Float64,9}}
  units: n/a
  required: false
  description: Field `harmonics_lpi` (default `Ref(SMatrix{3,3,Float64,9}((1.0, 0.0,
    0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)))`).
- id: maneuver_commands
  type: Vector{PropulsiveManeuverCommand}
  units: n/a
  required: false
  description: Field `maneuver_commands` (default `[PropulsiveManeuverCommand() for
    _ in 1:n_sats]`).
- id: maneuver_burn_plans
  type: Vector{PropulsiveBurnPlan}
  units: n/a
  required: false
  description: Field `maneuver_burn_plans` (default `[PropulsiveBurnPlan() for _ in
    1:n_sats]`).
- id: spice_runtime_counters
  type: SpiceRuntimeCounters
  units: n/a
  required: false
  description: Field `spice_runtime_counters` (default `SpiceRuntimeCounters()`).
- id: spice_rhs_memo_enabled
  type: Base.RefValue{Bool}
  units: n/a
  required: false
  description: Field `spice_rhs_memo_enabled` (default `Ref(true)`).
- id: spice_rhs_memo
  type: SpiceRhsMemo
  units: n/a
  required: false
  description: Field `spice_rhs_memo` (default `SpiceRhsMemo()`).
- id: current_time
  type: Base.RefValue{Float64}
  units: n/a
  required: false
  description: Field `current_time` (default `Ref(0.0)`).
- id: et_start
  type: Base.RefValue{Float64}
  units: n/a
  required: false
  description: Field `et_start` (default `Ref(0.0)`).
- id: solve_segment_end_time
  type: Base.RefValue{Float64}
  units: n/a
  required: false
  description: Field `solve_segment_end_time` (default `Ref(NaN)`).
- id: debug_control
  type: Base.RefValue{Bool}
  units: n/a
  required: false
  description: Field `debug_control` (default `Ref(false)`).
- id: debug_initial_derivative
  type: Base.RefValue{Bool}
  units: n/a
  required: false
  description: Field `debug_initial_derivative` (default `Ref(false)`).
- id: effector_cost_ns_per_item
  type: Base.RefValue{Float64}
  units: n/a
  required: false
  description: Field `effector_cost_ns_per_item` (default `Ref(NaN)`).
- id: effector_cost_samples
  type: Base.RefValue{Int64}
  units: n/a
  required: false
  description: Field `effector_cost_samples` (default `Ref(Int64(0))`).
- id: rhs_effector_cost_ns
  type: Base.RefValue{Vector{Float64}}
  units: n/a
  required: false
  description: Field `rhs_effector_cost_ns` (default `Ref(Float64[])`).
- id: rhs_effector_cost_samples
  type: Base.RefValue{Vector{Int64}}
  units: n/a
  required: false
  description: Field `rhs_effector_cost_samples` (default `Ref(Int64[])`).
- id: rhs_flat_effector_partials
  type: Base.RefValue{Array{Float64, 3}}
  units: n/a
  required: false
  description: Field `rhs_flat_effector_partials` (default `Ref(Array{Float64, 3}(undef,
    0, 0, 0))`).
- id: rhs_flat_effector_totals
  type: Base.RefValue{Matrix{Float64}}
  units: n/a
  required: false
  description: Field `rhs_flat_effector_totals` (default `Ref(Matrix{Float64}(undef,
    0, 0))`).
- id: rhs_flat_state_samples
  type: Base.RefValue{Vector{Union{Nothing, StateSample}}}
  units: n/a
  required: false
  description: Field `rhs_flat_state_samples` (default `Ref(Vector{Union{Nothing,
    StateSample}}())`).
- id: rhs_flat_state_pos_ii
  type: Base.RefValue{Vector{SVector{3, Float64}}}
  units: n/a
  required: false
  description: Field `rhs_flat_state_pos_ii` (default `Ref(SVector{3, Float64}[])`).
- id: rhs_flat_state_vel_ii
  type: Base.RefValue{Vector{SVector{3, Float64}}}
  units: n/a
  required: false
  description: Field `rhs_flat_state_vel_ii` (default `Ref(SVector{3, Float64}[])`).
- id: rhs_flat_state_mass_kg
  type: Base.RefValue{Vector{Float64}}
  units: n/a
  required: false
  description: Field `rhs_flat_state_mass_kg` (default `Ref(Float64[])`).
- id: rhs_flat_state_q_ib
  type: Base.RefValue{Vector{SVector{4, Float64}}}
  units: n/a
  required: false
  description: Field `rhs_flat_state_q_ib` (default `Ref(SVector{4, Float64}[])`).
- id: rhs_flat_state_omega_body
  type: Base.RefValue{Vector{SVector{3, Float64}}}
  units: n/a
  required: false
  description: Field `rhs_flat_state_omega_body` (default `Ref(SVector{3, Float64}[])`).
- id: rhs_flat_planet_lpi
  type: Base.RefValue{SMatrix{3, 3, Float64, 9}}
  units: n/a
  required: false
  description: Field `rhs_flat_planet_lpi` (default `Ref(SMatrix{3,3,Float64,9}((1.0,
    0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)))`).
- id: rhs_flat_planet_pos_pp
  type: Base.RefValue{Vector{SVector{3, Float64}}}
  units: n/a
  required: false
  description: Field `rhs_flat_planet_pos_pp` (default `Ref(SVector{3, Float64}[])`).
- id: rhs_flat_planet_vel_pp
  type: Base.RefValue{Vector{SVector{3, Float64}}}
  units: n/a
  required: false
  description: Field `rhs_flat_planet_vel_pp` (default `Ref(SVector{3, Float64}[])`).
- id: rhs_flat_planet_alt_m
  type: Base.RefValue{Vector{Float64}}
  units: n/a
  required: false
  description: Field `rhs_flat_planet_alt_m` (default `Ref(Float64[])`).
- id: rhs_flat_planet_lat_rad
  type: Base.RefValue{Vector{Float64}}
  units: n/a
  required: false
  description: Field `rhs_flat_planet_lat_rad` (default `Ref(Float64[])`).
- id: rhs_flat_planet_lon_rad
  type: Base.RefValue{Vector{Float64}}
  units: n/a
  required: false
  description: Field `rhs_flat_planet_lon_rad` (default `Ref(Float64[])`).
- id: rhs_flat_solar_pos_ii
  type: Base.RefValue{SVector{3, Float64}}
  units: n/a
  required: false
  description: Field `rhs_flat_solar_pos_ii` (default `Ref(SVector{3, Float64}(0.0,
    0.0, 0.0))`).
- id: rhs_flat_solar_t
  type: Base.RefValue{Float64}
  units: n/a
  required: false
  description: Field `rhs_flat_solar_t` (default `Ref(NaN)`).
- id: rhs_flat_work_items
  type: Base.RefValue{Vector{Int}}
  units: n/a
  required: false
  description: Field `rhs_flat_work_items` (default `Ref(Int[])`).
- id: rhs_flat_packet_starts
  type: Base.RefValue{Vector{Int}}
  units: n/a
  required: false
  description: Field `rhs_flat_packet_starts` (default `Ref(Int[])`).
- id: rhs_flat_packet_ends
  type: Base.RefValue{Vector{Int}}
  units: n/a
  required: false
  description: Field `rhs_flat_packet_ends` (default `Ref(Int[])`).
- id: rhs_flat_packet_costs
  type: Base.RefValue{Vector{Float64}}
  units: n/a
  required: false
  description: Field `rhs_flat_packet_costs` (default `Ref(Float64[])`).
- id: rhs_flat_packet_elapsed_ns
  type: Base.RefValue{Vector{Int64}}
  units: n/a
  required: false
  description: Field `rhs_flat_packet_elapsed_ns` (default `Ref(Int64[])`).
- id: rhs_flat_packet_overhead_ema
  type: Base.RefValue{Float64}
  units: n/a
  required: false
  description: Field `rhs_flat_packet_overhead_ema` (default `Ref(NaN)`).
- id: rhs_flat_packet_overhead_samples
  type: Base.RefValue{Int64}
  units: n/a
  required: false
  description: Field `rhs_flat_packet_overhead_samples` (default `Ref(Int64(0))`).
- id: rhs_flat_packet_disabled
  type: Base.RefValue{Bool}
  units: n/a
  required: false
  description: Field `rhs_flat_packet_disabled` (default `Ref(false)`).
- id: rhs_planet_frame_prefilled
  type: Base.RefValue{Bool}
  units: n/a
  required: false
  description: Field `rhs_planet_frame_prefilled` (default `Ref(false)`).
- id: rhs_atmosphere_prefilled
  type: Base.RefValue{Bool}
  units: n/a
  required: false
  description: Field `rhs_atmosphere_prefilled` (default `Ref(false)`).
- id: rhs_solar_prefilled
  type: Base.RefValue{Bool}
  units: n/a
  required: false
  description: Field `rhs_solar_prefilled` (default `Ref(false)`).
- id: rhs_harmonics_batch_pool
  type: Base.RefValue{Any}
  units: n/a
  required: false
  description: Field `rhs_harmonics_batch_pool` (default `Ref{Any}(nothing)`).
- id: in_atmosphere
  type: Vector{Bool}
  units: n/a
  required: false
  description: Field `in_atmosphere` (default `fill(false, n_sats)`).
- id: in_atmosphere_sample_t
  type: Vector{Float64}
  units: n/a
  required: false
  description: Field `in_atmosphere_sample_t` (default `fill(NaN, n_sats)`).
- id: rhs_plan_override
  type: Base.RefValue{Union{Nothing, RhsExecutionPlan}}
  units: n/a
  required: false
  description: Field `rhs_plan_override` (default `Ref{Union{Nothing, RhsExecutionPlan}}(nothing)`).
- id: rhs_plan_step_cache
  type: Base.RefValue{Union{Nothing, RhsExecutionPlan}}
  units: n/a
  required: false
  description: Field `rhs_plan_step_cache` (default `Ref{Union{Nothing, RhsExecutionPlan}}(nothing)`).
- id: robot_arm_present
  type: Base.RefValue{Union{Nothing, Bool}}
  units: n/a
  required: false
  description: Field `robot_arm_present` (default `Ref{Union{Nothing, Bool}}(nothing)`).
- id: policy_env_config
  type: Base.RefValue{Union{Nothing, PolicyDecisionEnvConfig}}
  units: n/a
  required: false
  description: Field `policy_env_config` (default `Ref{Union{Nothing, PolicyDecisionEnvConfig}}(nothing)`).
- id: rhs_env_config
  type: Base.RefValue{Union{Nothing, RhsPlanEnvConfig}}
  units: n/a
  required: false
  description: Field `rhs_env_config` (default `Ref{Union{Nothing, RhsPlanEnvConfig}}(nothing)`).
- id: callback_env_config
  type: Base.RefValue{Union{Nothing, CallbackEnvConfig}}
  units: n/a
  required: false
  description: Field `callback_env_config` (default `Ref{Union{Nothing, CallbackEnvConfig}}(nothing)`).
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
  type: SharedBuffers
  units: n/a
  description: Constructed `SharedBuffers` (keyword constructor via @kwdef).
- id: callees
  type: call
  units: n/a
  description: Calls this symbol makes to other mapped symbols.
tags:
- core
charts:
- core
origin: agent
---

# SharedBuffers

## Purpose
The mutable state shared between callbacks and the integrator for one run: atmosphere samples, caches, scratch workspaces, ephemeris tables, SPICE counters and memo, timing, flat-RHS prefill buffers, and the run-scoped environment snapshots.

## Design & Implementation
A `@kwdef struct` keyed on a plain runtime `n_sats` field rather than a type parameter, because none of its buffers are statically sized and a type parameter forced a fresh specialisation of the entire RHS graph per satellite count. Per-satellite vectors are allocated from `n_sats`; scalar mutable state lives in `Ref` cells so the struct itself can stay immutable. Notable fields include the `rhs_plan_override` and `rhs_plan_step_cache` refs, typed concretely because a `Ref{Any}` there boxed every plan access, and the `in_atmosphere` flags with their `NaN`-until-staged timestamps.

## Interface (ICD)
<!-- vulcan:icd:begin -->
| Direction | Socket | Type | Units | Required | Description |
|---|---|---|---|---|---|
| in | `n_sats` | Int | n/a | yes | Field `n_sats`. |
| in | `densities` | Vector{Float64} | n/a | no | Field `densities` (default `zeros(Float64, n_sats)`). |
| in | `temperatures` | Vector{Float64} | n/a | no | Field `temperatures` (default `ones(Float64, n_sats)`). |
| in | `winds` | Vector{SVector{3,Float64}} | n/a | no | Field `winds` (default `[SVector{3,Float64}(0.0, 0.0, 0.0) for _ in 1:n_sats]`). |
| in | `density_sample_t` | Vector{Float64} | n/a | no | Field `density_sample_t` (default `fill(NaN, n_sats)`). |
| in | `density_batch_altitudes` | Vector{Float64} | n/a | no | Field `density_batch_altitudes` (default `zeros(Float64, n_sats)`). |
| in | `density_batch_latitudes` | Vector{Float64} | n/a | no | Field `density_batch_latitudes` (default `zeros(Float64, n_sats)`). |
| in | `density_batch_longitudes` | Vector{Float64} | n/a | no | Field `density_batch_longitudes` (default `zeros(Float64, n_sats)`). |
| in | `heat_rates` | Vector{Vector{Float64}} | n/a | no | Field `heat_rates` (default `[Float64[] for _ in 1:n_sats]`). |
| in | `density_models` | Vector{_PerSatDensityModel} | n/a | no | Field `density_models` (default `_PerSatDensityModel[]`). |
| in | `gram_density_cache` | Vector{Union{Nothing, GramTrackCache}} | n/a | no | Field `gram_density_cache` (default `_typed_nothing_vector(GramTrackCache, n_sats)`). |
| in | `vacuum_gram_caches` | Vector{Union{Nothing, VacuumPredictedGRAMCache}} | n/a | no | Field `vacuum_gram_caches` (default `_typed_nothing_vector(VacuumPredictedGRAMCache, n_sats)`). |
| in | `gram_isolated_pool_models` | Vector{GRAMAtmosphereModel} | n/a | no | Field `gram_isolated_pool_models` (default `GRAMAtmosphereModel[]`). |
| in | `gram_isolated_pool_locks` | Vector{ReentrantLock} | n/a | no | Field `gram_isolated_pool_locks` (default `ReentrantLock[]`). |
| in | `harmonics_workspaces` | Vector{Union{Nothing, _HarmonicsWorkspaceMap}} | n/a | no | Field `harmonics_workspaces` (default `_typed_nothing_vector(_HarmonicsWorkspaceMap, n_sats)`). |
| in | `nbody_workspaces` | Vector{Union{Nothing, NBodyScratchWorkspace}} | n/a | no | Field `nbody_workspaces` (default `_typed_nothing_vector(NBodyScratchWorkspace, n_sats)`). |
| in | `aero_workspaces` | Vector{Union{Nothing, AeroScratchWorkspace}} | n/a | no | Field `aero_workspaces` (default `_typed_nothing_vector(AeroScratchWorkspace, n_sats)`). |
| in | `nbody_ephemeris_cache` | Base.RefValue{Union{Nothing, NBodyEphemerisCache}} | n/a | no | Field `nbody_ephemeris_cache` (default `Ref{Union{Nothing, NBodyEphemerisCache}}(nothing)`). |
| in | `srp_sun_ephemeris_cache` | Base.RefValue{Union{Nothing, SRPSunEphemerisCache}} | n/a | no | Field `srp_sun_ephemeris_cache` (default `Ref{Union{Nothing, SRPSunEphemerisCache}}(nothing)`). |
| in | `planet_frame_ephemeris_cache` | Base.RefValue{Union{Nothing, PlanetFrameEphemerisCache}} | n/a | no | Field `planet_frame_ephemeris_cache` (default `Ref{Union{Nothing, PlanetFrameEphemerisCache}}(nothing)`). |
| in | `harmonics_lpi_lock` | ReentrantLock | n/a | no | Field `harmonics_lpi_lock` (default `ReentrantLock()`). |
| in | `harmonics_lpi_key` | Base.RefValue{Any} | n/a | no | Field `harmonics_lpi_key` (default `Ref{Any}(nothing)`). |
| in | `harmonics_lpi` | Base.RefValue{SMatrix{3,3,Float64,9}} | n/a | no | Field `harmonics_lpi` (default `Ref(SMatrix{3,3,Float64,9}((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)))`). |
| in | `maneuver_commands` | Vector{PropulsiveManeuverCommand} | n/a | no | Field `maneuver_commands` (default `[PropulsiveManeuverCommand() for _ in 1:n_sats]`). |
| in | `maneuver_burn_plans` | Vector{PropulsiveBurnPlan} | n/a | no | Field `maneuver_burn_plans` (default `[PropulsiveBurnPlan() for _ in 1:n_sats]`). |
| in | `spice_runtime_counters` | SpiceRuntimeCounters | n/a | no | Field `spice_runtime_counters` (default `SpiceRuntimeCounters()`). |
| in | `spice_rhs_memo_enabled` | Base.RefValue{Bool} | n/a | no | Field `spice_rhs_memo_enabled` (default `Ref(true)`). |
| in | `spice_rhs_memo` | SpiceRhsMemo | n/a | no | Field `spice_rhs_memo` (default `SpiceRhsMemo()`). |
| in | `current_time` | Base.RefValue{Float64} | n/a | no | Field `current_time` (default `Ref(0.0)`). |
| in | `et_start` | Base.RefValue{Float64} | n/a | no | Field `et_start` (default `Ref(0.0)`). |
| in | `solve_segment_end_time` | Base.RefValue{Float64} | n/a | no | Field `solve_segment_end_time` (default `Ref(NaN)`). |
| in | `debug_control` | Base.RefValue{Bool} | n/a | no | Field `debug_control` (default `Ref(false)`). |
| in | `debug_initial_derivative` | Base.RefValue{Bool} | n/a | no | Field `debug_initial_derivative` (default `Ref(false)`). |
| in | `effector_cost_ns_per_item` | Base.RefValue{Float64} | n/a | no | Field `effector_cost_ns_per_item` (default `Ref(NaN)`). |
| in | `effector_cost_samples` | Base.RefValue{Int64} | n/a | no | Field `effector_cost_samples` (default `Ref(Int64(0))`). |
| in | `rhs_effector_cost_ns` | Base.RefValue{Vector{Float64}} | n/a | no | Field `rhs_effector_cost_ns` (default `Ref(Float64[])`). |
| in | `rhs_effector_cost_samples` | Base.RefValue{Vector{Int64}} | n/a | no | Field `rhs_effector_cost_samples` (default `Ref(Int64[])`). |
| in | `rhs_flat_effector_partials` | Base.RefValue{Array{Float64, 3}} | n/a | no | Field `rhs_flat_effector_partials` (default `Ref(Array{Float64, 3}(undef, 0, 0, 0))`). |
| in | `rhs_flat_effector_totals` | Base.RefValue{Matrix{Float64}} | n/a | no | Field `rhs_flat_effector_totals` (default `Ref(Matrix{Float64}(undef, 0, 0))`). |
| in | `rhs_flat_state_samples` | Base.RefValue{Vector{Union{Nothing, StateSample}}} | n/a | no | Field `rhs_flat_state_samples` (default `Ref(Vector{Union{Nothing, StateSample}}())`). |
| in | `rhs_flat_state_pos_ii` | Base.RefValue{Vector{SVector{3, Float64}}} | n/a | no | Field `rhs_flat_state_pos_ii` (default `Ref(SVector{3, Float64}[])`). |
| in | `rhs_flat_state_vel_ii` | Base.RefValue{Vector{SVector{3, Float64}}} | n/a | no | Field `rhs_flat_state_vel_ii` (default `Ref(SVector{3, Float64}[])`). |
| in | `rhs_flat_state_mass_kg` | Base.RefValue{Vector{Float64}} | n/a | no | Field `rhs_flat_state_mass_kg` (default `Ref(Float64[])`). |
| in | `rhs_flat_state_q_ib` | Base.RefValue{Vector{SVector{4, Float64}}} | n/a | no | Field `rhs_flat_state_q_ib` (default `Ref(SVector{4, Float64}[])`). |
| in | `rhs_flat_state_omega_body` | Base.RefValue{Vector{SVector{3, Float64}}} | n/a | no | Field `rhs_flat_state_omega_body` (default `Ref(SVector{3, Float64}[])`). |
| in | `rhs_flat_planet_lpi` | Base.RefValue{SMatrix{3, 3, Float64, 9}} | n/a | no | Field `rhs_flat_planet_lpi` (default `Ref(SMatrix{3,3,Float64,9}((1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)))`). |
| in | `rhs_flat_planet_pos_pp` | Base.RefValue{Vector{SVector{3, Float64}}} | n/a | no | Field `rhs_flat_planet_pos_pp` (default `Ref(SVector{3, Float64}[])`). |
| in | `rhs_flat_planet_vel_pp` | Base.RefValue{Vector{SVector{3, Float64}}} | n/a | no | Field `rhs_flat_planet_vel_pp` (default `Ref(SVector{3, Float64}[])`). |
| in | `rhs_flat_planet_alt_m` | Base.RefValue{Vector{Float64}} | n/a | no | Field `rhs_flat_planet_alt_m` (default `Ref(Float64[])`). |
| in | `rhs_flat_planet_lat_rad` | Base.RefValue{Vector{Float64}} | n/a | no | Field `rhs_flat_planet_lat_rad` (default `Ref(Float64[])`). |
| in | `rhs_flat_planet_lon_rad` | Base.RefValue{Vector{Float64}} | n/a | no | Field `rhs_flat_planet_lon_rad` (default `Ref(Float64[])`). |
| in | `rhs_flat_solar_pos_ii` | Base.RefValue{SVector{3, Float64}} | n/a | no | Field `rhs_flat_solar_pos_ii` (default `Ref(SVector{3, Float64}(0.0, 0.0, 0.0))`). |
| in | `rhs_flat_solar_t` | Base.RefValue{Float64} | n/a | no | Field `rhs_flat_solar_t` (default `Ref(NaN)`). |
| in | `rhs_flat_work_items` | Base.RefValue{Vector{Int}} | n/a | no | Field `rhs_flat_work_items` (default `Ref(Int[])`). |
| in | `rhs_flat_packet_starts` | Base.RefValue{Vector{Int}} | n/a | no | Field `rhs_flat_packet_starts` (default `Ref(Int[])`). |
| in | `rhs_flat_packet_ends` | Base.RefValue{Vector{Int}} | n/a | no | Field `rhs_flat_packet_ends` (default `Ref(Int[])`). |
| in | `rhs_flat_packet_costs` | Base.RefValue{Vector{Float64}} | n/a | no | Field `rhs_flat_packet_costs` (default `Ref(Float64[])`). |
| in | `rhs_flat_packet_elapsed_ns` | Base.RefValue{Vector{Int64}} | n/a | no | Field `rhs_flat_packet_elapsed_ns` (default `Ref(Int64[])`). |
| in | `rhs_flat_packet_overhead_ema` | Base.RefValue{Float64} | n/a | no | Field `rhs_flat_packet_overhead_ema` (default `Ref(NaN)`). |
| in | `rhs_flat_packet_overhead_samples` | Base.RefValue{Int64} | n/a | no | Field `rhs_flat_packet_overhead_samples` (default `Ref(Int64(0))`). |
| in | `rhs_flat_packet_disabled` | Base.RefValue{Bool} | n/a | no | Field `rhs_flat_packet_disabled` (default `Ref(false)`). |
| in | `rhs_planet_frame_prefilled` | Base.RefValue{Bool} | n/a | no | Field `rhs_planet_frame_prefilled` (default `Ref(false)`). |
| in | `rhs_atmosphere_prefilled` | Base.RefValue{Bool} | n/a | no | Field `rhs_atmosphere_prefilled` (default `Ref(false)`). |
| in | `rhs_solar_prefilled` | Base.RefValue{Bool} | n/a | no | Field `rhs_solar_prefilled` (default `Ref(false)`). |
| in | `rhs_harmonics_batch_pool` | Base.RefValue{Any} | n/a | no | Field `rhs_harmonics_batch_pool` (default `Ref{Any}(nothing)`). |
| in | `in_atmosphere` | Vector{Bool} | n/a | no | Field `in_atmosphere` (default `fill(false, n_sats)`). |
| in | `in_atmosphere_sample_t` | Vector{Float64} | n/a | no | Field `in_atmosphere_sample_t` (default `fill(NaN, n_sats)`). |
| in | `rhs_plan_override` | Base.RefValue{Union{Nothing, RhsExecutionPlan}} | n/a | no | Field `rhs_plan_override` (default `Ref{Union{Nothing, RhsExecutionPlan}}(nothing)`). |
| in | `rhs_plan_step_cache` | Base.RefValue{Union{Nothing, RhsExecutionPlan}} | n/a | no | Field `rhs_plan_step_cache` (default `Ref{Union{Nothing, RhsExecutionPlan}}(nothing)`). |
| in | `robot_arm_present` | Base.RefValue{Union{Nothing, Bool}} | n/a | no | Field `robot_arm_present` (default `Ref{Union{Nothing, Bool}}(nothing)`). |
| in | `policy_env_config` | Base.RefValue{Union{Nothing, PolicyDecisionEnvConfig}} | n/a | no | Field `policy_env_config` (default `Ref{Union{Nothing, PolicyDecisionEnvConfig}}(nothing)`). |
| in | `rhs_env_config` | Base.RefValue{Union{Nothing, RhsPlanEnvConfig}} | n/a | no | Field `rhs_env_config` (default `Ref{Union{Nothing, RhsPlanEnvConfig}}(nothing)`). |
| in | `callback_env_config` | Base.RefValue{Union{Nothing, CallbackEnvConfig}} | n/a | no | Field `callback_env_config` (default `Ref{Union{Nothing, CallbackEnvConfig}}(nothing)`). |
| in | `module_api` | Module | n/a | no | Re-exported through the owning module's public surface. |
| in | `callers` | call | n/a | no | Invocations of this symbol observed in mapped callers. |
| out | `result` | SharedBuffers | n/a | — | Constructed `SharedBuffers` (keyword constructor via @kwdef). |
| out | `callees` | call | n/a | — | Calls this symbol makes to other mapped symbols. |
<!-- vulcan:icd:end -->

## Connections
<!-- vulcan:connections:begin -->
**Upstream**

- [[module.core|SimulationModel]] · `api` → `module_api` · call · `src/core/types/runtime_types.jl`

**Downstream**

- `callees` → [[core.runtime_types__typed_nothing_vector|_typed_nothing_vector]] · `callers` · call · `src/core/types/runtime_types.jl:722-722`
- `callees` → [[core.runtime_types_spicerhsmemo|SpiceRhsMemo]] · `callers` · call · `src/core/types/runtime_types.jl:739-739`
- `callees` → [[core.runtime_types_spiceruntimecounters|SpiceRuntimeCounters]] · `callers` · call · `src/core/types/runtime_types.jl:737-737`
- `callees` → [[gnc.command_types_propulsiveburnplan|PropulsiveBurnPlan]] · `callers` · call · `src/core/types/runtime_types.jl:736-736`
- `callees` → [[gnc.propulsive_maneuver_command|PropulsiveManeuverCommand]] · `callers` · call · `src/core/types/runtime_types.jl:735-735`
<!-- vulcan:connections:end -->

## Limitations
With over seventy fields it is the de facto global state of a run, and the ownership of each field — which callback writes it, which RHS path reads it — is documented only in comments; the `harmonics_lpi_key` and `rhs_harmonics_batch_pool` slots remain `Ref{Any}`.

## Provenance
Mapped from `src/core/types/runtime_types.jl` line 711.
