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
    # Flight-measured atmosphere mode: atmosphere_model = "tabulated_flight"
    # flies the sim on per-pass measured density profiles (see
    # TabulatedFlightAtmosphereModel). Used by the Odyssey nightly as a
    # digital-twin regression sentinel (PI decision, 2026-07) and by the
    # certification envelope; always visible in the summary via the
    # atmosphere_model column, so rows are never ambiguous about their mode.
    tabulated_flight_file::String = ""
    tabulated_flight_sigma::Float64 = 0.0
    # Time-tabulated atmosphere mode: atmosphere_model = "tabulated_time" flies
    # the sim on an externally supplied rho(t) table in scenario elapsed time
    # (see TimeTabulatedAtmosphereModel) — the orbital-arc counterpart of
    # tabulated_flight, used for flight-inferred along-track density replay
    # (digital-twin runs) and assimilated-product sampling.
    tabulated_time_file::String = ""
    tabulated_time_scale::Float64 = 1.0
    tabulated_time_temperature_k::Float64 = 900.0
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
    # :legacy = historical bus ref_area dims[1]*dims[3]; :frontal = the flow-normal
    # dims[2]*dims[3] face matching the Hart free-molecular coefficient
    # normalization (see make_three_body_spacecraft).
    bus_ram_face::Symbol = :legacy
    # Optional per-link attitude quaternions (x, y, z, w scalar-last),
    # normalized at parse time. The bus quaternion is expressed in the
    # flow-aligned reference frame (see the AerodynamicCoefficientfM
    # docstring); panel quaternions are relative to the BUS frame — the
    # kinematics convention, and the natural way to express a fixed panel
    # cant. `nothing` leaves the link at identity, bit-identical to the
    # pre-capability spacecraft. The manifest layer requires
    # aero_fixed_attitude_incidence = "attitude" whenever any of these is
    # set: the historical :max_drag path reads non-root quaternions, so the
    # combination would silently change default-mode physics.
    bus_attitude_q::Union{Nothing, NTuple{4, Float64}} = nothing
    panel_attitude_q_left::Union{Nothing, NTuple{4, Float64}} = nothing
    panel_attitude_q_right::Union{Nothing, NTuple{4, Float64}} = nothing
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
    # Optional exact J2000 Cartesian state (x, y, z in m; vx, vy, vz in m/s)
    # relative to the central body at initial_time, e.g. taken directly from the
    # mission NAV kernel. When present it overrides the element-based initial
    # condition above (the elements stay in the manifest as documentation).
    initial_state_j2000_m::Union{Nothing, NTuple{6, Float64}} = nothing
    # Campaign orbit number of the scenario epoch (the truth product's numbering
    # origin may predate the epoch, e.g. counting from orbit insertion). When
    # set, sim apsis events are placed at epoch_orbit_offset + k with unit step
    # (one apsis per orbit) and scoring is masked to the simulated span.
    epoch_orbit_offset::Union{Nothing, Float64} = nothing
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
    # Fixed-attitude incidence mode forwarded to AerodynamicCoefficientfM
    # (:max_drag | :attitude | :tumbling_average); :max_drag preserves the
    # historical accounting bit-identically.
    aero_fixed_attitude_incidence::Symbol = :max_drag
    include_wind::Bool = false
    orbit_altitude_mode::Symbol = :vacuum
    maneuver_orbit_numbers::Vector{Int64} = Int64[]
    maneuver_orbit_numbers_campaign::Vector{Int64} = Int64[]
    maneuver_delta_v_mps::Vector{Float64} = Float64[]
    # Diagnostic replay scaling: "delta_v" (benchmark default, replay flight
    # dv verbatim) or "flight_apoapsis_ratio" (scale each burn by flight/sim
    # apoapsis radius so it delivers the flight's periapsis change; injects
    # flight truth — diagnostics only, recorded in the summary).
    maneuver_replay_scale_mode::String = "delta_v"
    maneuver_flight_apoapsis_alt_m::Vector{Float64} = Float64[]
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
    # Optional: real per-sample velocity ground truth. When all three are present,
    # state_vx_time/state_vy_time/state_vz_time error rows are computed against
    # these telemetry columns instead of being omitted (the default vx_kmps/vy_kmps/
    # vz_kmps fallback in _load_time_aligned_telemetry is numerically differentiated
    # from position and is not compared against directly).
    telemetry_vx_col::Union{Nothing, String} = nothing
    telemetry_vy_col::Union{Nothing, String} = nothing
    telemetry_vz_col::Union{Nothing, String} = nothing
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
    # Fixed-attitude incidence mode forwarded to AerodynamicCoefficientfM
    # (:max_drag | :attitude | :tumbling_average); :max_drag preserves the
    # historical accounting bit-identically.
    aero_fixed_attitude_incidence::Symbol = :max_drag
    include_wind::Bool = false
    orbit_altitude_mode::Symbol = :vacuum
    cartesian_ic_frame::Symbol = :inertial
    comparison_frame::Symbol = :inertial
    comparison_mode::Symbol = :time_aligned_state
    extrema_min_separation_s::Float64 = 500.0
    atmosphere_truth::AtmosphereTruthConfig = AtmosphereTruthConfig()
    calibration::CalibrationConfig = CalibrationConfig()
    EI_km::Float64
    # Constant offsets added to the Cartesian telemetry IC (J2000, after any
    # frame conversion). The knobs behind differential-correction IC fitting
    # (fit_initial_state) and offset sweeps; zero by default.
    ic_offset_m::NTuple{3, Float64} = (0.0, 0.0, 0.0)
    ic_offset_mps::NTuple{3, Float64} = (0.0, 0.0, 0.0)
    # Illumination screen for single-frequency GNSS truth: ionospheric group
    # delay biases the dayside radial by O(100 m) (measured -162 m mean on
    # CYGNSS FM4), so :nightside keeps only samples with the position vector
    # anti-sunward (cos(sun angle) < 0). :dayside keeps the complement.
    truth_mask::Symbol = :none
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
