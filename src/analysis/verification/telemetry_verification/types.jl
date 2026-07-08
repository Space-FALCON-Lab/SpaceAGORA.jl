const EventTolerance = NamedTuple{
    (:max_abs_km, :max_nmae, :max_rmse_km),
    Tuple{Float64, Float64, Float64}
}

Base.@kwdef struct AtmosphereTruthConfig
    assumption_id::String = "gram_default"
    atmosphere_model::String = "GRAM"
    atmosphere_dataset::String = "default"
    space_weather_model::String = "default"
    solar_flux_model::String = "default"
    gram_seed::Int = 1001
    gram_perturbation_scales::NTuple{4, Float64} = (0.0, 0.0, 0.0, 0.0)
    gram_min_relative_step_size::Union{Nothing, Float64} = nothing
    gram_offline_surrogate::String = "off"
    gram_static_grid::Bool = false
    gram_track_cache::Bool = false
    gram_global_lock::String = "on"
    mars_map_year::Union{Nothing, Int} = nothing
    mars_mgcm_dust_levels::Union{Nothing, NTuple{3, Float64}} = nothing
    mars_dust_storm::Union{Nothing, NTuple{6, Float64}} = nothing
    mars_f107::Union{Nothing, Float64} = nothing
    mars_wind_scales::Union{Nothing, NTuple{2, Float64}} = nothing
    mars_mola_heights::Union{Nothing, Bool} = nothing
    mars_min_max::Union{Nothing, Int} = nothing
end

Base.@kwdef struct CalibrationConfig
    enabled::Bool = false
    profiles::Vector{Symbol} = Symbol[:full]
    search_on_quick_subset::Bool = true
    fit_cd_scale::Bool = true
    cd_scale_min::Float64 = 0.85
    cd_scale_max::Float64 = 1.15
    cd_scale_steps::Int = 3
    fit_cr::Bool = true
    cr_min::Float64 = 1.15
    cr_max::Float64 = 1.45
    cr_steps::Int = 3
    fit_bias::Bool = true
    bias_abs_max_km::Float64 = 500.0
    objective::String = "mean_nmae"
end

Base.@kwdef struct SpacecraftConfig
    bus_dims::NTuple{3, Float64}
    panel_dims::NTuple{3, Float64}
    bus_mass_kg::Float64
    panel_mass_each_kg::Float64
    panel_offset_y_m::Float64
    prop_mass_kg::Float64
    id::Int64
end

abstract type AbstractScenarioConfig end

Base.@kwdef struct OrbitEventsScenarioConfig <: AbstractScenarioConfig
    name::String
    planet_name::String
    telemetry_peri_path::String
    telemetry_apo_path::String
    target_orbits_quick::Int
    target_orbits_full::Int
    compare_points_quick::Int
    compare_points_full::Int
    min_eval_points::Int
    units_x::String
    units_y::Dict{String, String}
    tolerances_quick::Dict{String, EventTolerance}
    tolerances_full::Dict{String, EventTolerance}
    initial_time::InitialTime
    ra_m::Float64
    rp_altitude_m::Float64
    i_deg::Float64
    aop_deg::Float64
    raan_deg::Float64
    ta_deg::Float64
    element_frame::Symbol = :j2000
    spacecraft::SpacecraftConfig
    gravity_model::Symbol
    gravity_harmonics_degree::Int = 0
    gravity_harmonics_order::Int = 0
    gravity_harmonics_file::String = ""
    nbody_bodies::Vector{String} = String[]
    srp_enabled::Bool = false
    srp_cr::Float64 = 1.3
    srp_area_m2::Float64 = 0.0
    drag_enabled::Bool = true
    include_wind::Bool = false
    orbit_altitude_mode::Symbol = :vacuum
    maneuver_orbit_numbers::Vector{Int64} = Int64[]
    maneuver_delta_v_mps::Vector{Float64} = Float64[]
    maneuver_thrust_n::Float64 = 0.0
    maneuver_isp_s::Float64 = 0.0
    maneuver_guidance_rate_s::Float64 = 30.0
    maneuver_control_rate_s::Float64 = 10.0
    atmosphere_truth::AtmosphereTruthConfig = AtmosphereTruthConfig()
    calibration::CalibrationConfig = CalibrationConfig()
    EI_km::Float64
end

Base.@kwdef struct TimeAlignedScenarioConfig <: AbstractScenarioConfig
    name::String
    planet_name::String
    telemetry_path::String
    telemetry_time_col::String
    telemetry_altitude_col::String
    telemetry_x_col::String
    telemetry_y_col::String
    telemetry_z_col::String
    telemetry_sma_col::Union{Nothing, String} = nothing
    telemetry_ecc_col::Union{Nothing, String} = nothing
    telemetry_inc_col::Union{Nothing, String} = nothing
    telemetry_aop_col::Union{Nothing, String} = nothing
    telemetry_raan_col::Union{Nothing, String} = nothing
    telemetry_ta_col::Union{Nothing, String} = nothing
    # Optional: Cartesian IC columns. When all six are present, CartesianInitialCondition is used instead of Keplerian.
    telemetry_x_ic_col::Union{Nothing, String} = nothing
    telemetry_y_ic_col::Union{Nothing, String} = nothing
    telemetry_z_ic_col::Union{Nothing, String} = nothing
    telemetry_vx_ic_col::Union{Nothing, String} = nothing
    telemetry_vy_ic_col::Union{Nothing, String} = nothing
    telemetry_vz_ic_col::Union{Nothing, String} = nothing
    max_points_quick::Int
    max_points_full::Int
    min_eval_points::Int
    units_x::String
    units_y::Dict{String, String}
    tolerances_quick::Dict{String, EventTolerance}
    tolerances_full::Dict{String, EventTolerance}
    initial_time::InitialTime
    spacecraft::SpacecraftConfig
    gravity_model::Symbol
    gravity_harmonics_degree::Int = 0
    gravity_harmonics_order::Int = 0
    gravity_harmonics_file::String = ""
    nbody_bodies::Vector{String} = String[]
    srp_enabled::Bool = false
    srp_cr::Float64 = 1.3
    srp_area_m2::Float64 = 0.0
    drag_enabled::Bool = true
    include_wind::Bool = false
    orbit_altitude_mode::Symbol = :vacuum
    cartesian_ic_frame::Symbol = :inertial
    comparison_frame::Symbol = :inertial
    comparison_mode::Symbol = :time_aligned_state
    extrema_min_separation_s::Float64 = 500.0
    atmosphere_truth::AtmosphereTruthConfig = AtmosphereTruthConfig()
    calibration::CalibrationConfig = CalibrationConfig()
    EI_km::Float64
end

Base.@kwdef struct StudyConfig
    profile::Symbol
    out_summary::String
    out_errors::String
    manifest_path::String
    enforce::Bool
    generate_plots::Bool
end

"""
    VerificationRequest

Typed request for the telemetry verification study entrypoint. It controls the
profile, output paths, manifest path, enforcement mode, and optional plot
generation.
"""
Base.@kwdef struct VerificationRequest
    profile::Symbol = Symbol(lowercase(get(ENV, "SPACEAGORA_TELEMETRY_PROFILE", "quick")))
    out_summary::String = abspath(get(ENV, "SPACEAGORA_TELEMETRY_OUT_SUMMARY", joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_summary.csv")))
    out_errors::String = abspath(get(ENV, "SPACEAGORA_TELEMETRY_OUT_ERRORS", joinpath(DEFAULT_OUTPUT_DIR, "telemetry_orbit_accuracy_errors.csv")))
    manifest_path::String = abspath(get(ENV, "SPACEAGORA_TELEMETRY_MANIFEST", DEFAULT_MANIFEST_PATH))
    enforce::Bool = false
    generate_plots::Bool = _safe_parse_bool(get(ENV, "SPACEAGORA_TELEMETRY_PLOTS", "1"), true)
end

"""
    VerificationResult

Structured result returned by telemetry verification runs, including the summary
table, pointwise error table, emitted artifact paths, and total runtime.
"""
Base.@kwdef struct VerificationResult
    summary::DataFrame
    errors::DataFrame
    summary_path::String
    errors_path::String
    plots_dir::String
    profile::Symbol
    enforce::Bool
    total_runtime_s::Float64
end
