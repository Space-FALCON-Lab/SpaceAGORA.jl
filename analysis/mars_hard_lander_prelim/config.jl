const ANALYSIS_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

struct VehicleGeometryConfig
    body_label::String
    deployed_label::String
    nose_radius_m::Float64
    base_radius_m::Float64
    deployed_drag_surface_diameter_m::Float64
    drag_skirt_height_m::Float64
    cone_half_angle_deg::Float64
    hypersonic_pressure_model::Symbol
    body_alpha_deg::Float64
    panel_alpha_deg::Float64
    panel_count::Int
    panel_aspect_ratio::Float64
    panel_lateral_deflection_deg::Float64
    panel_skin_friction_coefficient::Float64
    panel_zero_lift_drag_coefficient::Float64
end

function shield_informed_vehicle_geometry(; body_alpha_deg::Real=0.0)
    return VehicleGeometryConfig(
        "SHIELD-like 70-deg sphere-cone body",
        "SHIELD-like 70-deg sphere-cone body + symmetric broadside flat plates",
        0.08,
        0.25,
        0.0,
        0.0,
        70.0,
        :modified_newtonian,
        Float64(body_alpha_deg),
        90.0,
        2,
        4.0,
        30.0,
        0.0,
        0.0,
    )
end

function shield_published_surrogate_geometry(; body_alpha_deg::Real=0.0)
    return VehicleGeometryConfig(
        "Published-size SHIELD stowed sphere-cone surrogate",
        "Published-size SHIELD body + drag-skirt-equivalent symmetric broadside area",
        0.288,
        0.9,
        2.2,
        0.75,
        50.0,
        :modified_newtonian,
        Float64(body_alpha_deg),
        90.0,
        2,
        4.0,
        30.0,
        0.0,
        0.0,
    )
end

function default_vehicle_geometry()
    return shield_informed_vehicle_geometry()
end

function with_body_alpha(geometry::VehicleGeometryConfig, body_alpha_deg::Real)
    return VehicleGeometryConfig(
        geometry.body_label,
        geometry.deployed_label,
        geometry.nose_radius_m,
        geometry.base_radius_m,
        geometry.deployed_drag_surface_diameter_m,
        geometry.drag_skirt_height_m,
        geometry.cone_half_angle_deg,
        geometry.hypersonic_pressure_model,
        Float64(body_alpha_deg),
        geometry.panel_alpha_deg,
        geometry.panel_count,
        geometry.panel_aspect_ratio,
        geometry.panel_lateral_deflection_deg,
        geometry.panel_skin_friction_coefficient,
        geometry.panel_zero_lift_drag_coefficient,
    )
end

struct PrelimConfig
    repo_root::String
    analysis_name::String
    output_root::String
    trajectory_dir::String
    plot_dir::String
    atmosphere_mode::Symbol
    use_gram_surrogate::Bool
    gram_cache_altitude_m::Float64
    gram_cache_longitude_deg::Float64
    gram_cache_time_s::Float64
    initial_time
    lat_deg::Float64
    initial_elapsed_time_s::Float64
    planet
    geometry::VehicleGeometryConfig
    h0_m::Float64
    V0_mps::Float64
    gamma0_deg::Float64
    theta0_rad::Float64
    beta_high_targets::Vector{Float64}
    beta_ratios::Vector{Float64}
    h_jettison_grid_m::Vector{Float64}
    h_switch_grid_m::Vector{Float64}
    secondary_V0_mps::Vector{Float64}
    secondary_gamma0_deg::Vector{Float64}
    run_secondary_sweeps::Bool
    include_reverse_bangbang_sweep::Bool
    enable_crossrange_proxy::Bool
    crossrange_sigmas::Vector{Float64}
    enable_target_range_guidance::Bool
    target_range_fraction::Float64
    target_range_tolerance_km::Float64
    target_range_max_iterations::Int
    enable_monte_carlo::Bool
    monte_carlo_samples::Int
    monte_carlo_seed::Int
    monte_carlo_density_scale_sigma::Float64
    monte_carlo_V0_sigma_mps::Float64
    monte_carlo_gamma0_sigma_deg::Float64
    monte_carlo_beta_sigma_fraction::Float64
    monte_carlo_h0_sigma_m::Float64
    monte_carlo_crossrange_sigma::Float64
    monte_carlo_use_winds::Bool
    lateral_guidance_crossrange_scale_km::Float64
    lateral_guidance_crossrange_rate_scale_mps::Float64
    enable_crossrange_targeted_monte_carlo::Bool
    crossrange_target_panel_cant_deg::Float64
    crossrange_target_differential_fraction::Float64
    trim_panel_area_fraction_of_deployed::Float64
    panel_system_areal_density_kg_m2::Float64
    saveat_s::Float64
    solver_reltol::Float64
    solver_abstol::Float64
    max_time_s::Float64
    impact_stop_distance_m::Float64
    representative_beta_high::Float64
    representative_beta_ratio::Float64
    representative_h_jettison_m::Float64
    representative_switches_m::Vector{Float64}
    generate_plots::Bool
end

function default_config(;
    output_root::Union{Nothing, String}=nothing,
    atmosphere_mode::Symbol=:gram,
    use_gram_surrogate::Bool=false,
    run_secondary_sweeps::Bool=false,
    include_reverse_bangbang_sweep::Bool=false,
    enable_crossrange_proxy::Bool=true,
    enable_target_range_guidance::Bool=true,
    enable_monte_carlo::Bool=true,
    generate_plots::Bool=true,
    h_switch_grid_m::Union{Nothing, AbstractVector{<:Real}}=nothing,
    beta_high_targets::Union{Nothing, AbstractVector{<:Real}}=nothing,
    beta_ratios::Union{Nothing, AbstractVector{<:Real}}=nothing,
    h_jettison_grid_m::Union{Nothing, AbstractVector{<:Real}}=nothing,
    representative_switches_m::Union{Nothing, AbstractVector{<:Real}}=nothing,
    crossrange_sigmas::Union{Nothing, AbstractVector{<:Real}}=nothing,
    representative_beta_high::Real=150.0,
    representative_beta_ratio::Real=4.0,
    representative_h_jettison_m::Real=20e3,
    target_range_fraction::Real=0.5,
    target_range_tolerance_km::Real=0.25,
    target_range_max_iterations::Integer=20,
    monte_carlo_samples::Integer=100,
    monte_carlo_seed::Integer=42,
    monte_carlo_density_scale_sigma::Real=0.20,
    monte_carlo_V0_sigma_mps::Real=50.0,
    monte_carlo_gamma0_sigma_deg::Real=0.25,
    monte_carlo_beta_sigma_fraction::Real=0.05,
    monte_carlo_h0_sigma_m::Real=500.0,
    monte_carlo_crossrange_sigma::Real=0.05,
    monte_carlo_use_winds::Bool=true,
    lateral_guidance_crossrange_scale_km::Real=20.0,
    lateral_guidance_crossrange_rate_scale_mps::Real=100.0,
    enable_crossrange_targeted_monte_carlo::Bool=true,
    crossrange_target_panel_cant_deg::Real=45.0,
    crossrange_target_differential_fraction::Real=0.40,
    trim_panel_area_fraction_of_deployed::Real=0.0,
    panel_system_areal_density_kg_m2::Real=5.0,
    geometry::Union{Nothing, VehicleGeometryConfig}=nothing,
)
    analysis_name = "mars_hard_lander_prelim"
    output_root_resolved = output_root === nothing ?
        joinpath(ANALYSIS_REPO_ROOT, "results", analysis_name) :
        String(output_root)
    trajectory_dir = joinpath(output_root_resolved, "trajectories")
    plot_dir = joinpath(output_root_resolved, "plots")
    return PrelimConfig(
        ANALYSIS_REPO_ROOT,
        analysis_name,
        output_root_resolved,
        trajectory_dir,
        plot_dir,
        atmosphere_mode,
        use_gram_surrogate,
        250.0,
        0.25,
        5.0,
        SM.InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
        0.0,
        0.0,
        SM.Mars(),
        geometry === nothing ? default_vehicle_geometry() : geometry,
        125e3,
        5.5e3,
        -12.0,
        0.0,
        beta_high_targets === nothing ? [75.0, 150.0, 300.0] : Float64.(collect(beta_high_targets)),
        beta_ratios === nothing ? [1.0, 2.0, 4.0, 8.0] : Float64.(collect(beta_ratios)),
        h_jettison_grid_m === nothing ? [0.0, 5e3, 10e3, 20e3, 30e3, 40e3] : Float64.(collect(h_jettison_grid_m)),
        h_switch_grid_m === nothing ? collect(120e3:-2e3:0.0) : Float64.(collect(h_switch_grid_m)),
        [4.5e3, 5.5e3, 6.5e3],
        [-10.0, -12.0, -15.0],
        run_secondary_sweeps,
        include_reverse_bangbang_sweep,
        enable_crossrange_proxy,
        crossrange_sigmas === nothing ? [0.02, 0.05, 0.10, 0.20] : Float64.(collect(crossrange_sigmas)),
        enable_target_range_guidance,
        Float64(target_range_fraction),
        Float64(target_range_tolerance_km),
        Int(target_range_max_iterations),
        enable_monte_carlo,
        Int(monte_carlo_samples),
        Int(monte_carlo_seed),
        Float64(monte_carlo_density_scale_sigma),
        Float64(monte_carlo_V0_sigma_mps),
        Float64(monte_carlo_gamma0_sigma_deg),
        Float64(monte_carlo_beta_sigma_fraction),
        Float64(monte_carlo_h0_sigma_m),
        Float64(monte_carlo_crossrange_sigma),
        monte_carlo_use_winds,
        Float64(lateral_guidance_crossrange_scale_km),
        Float64(lateral_guidance_crossrange_rate_scale_mps),
        enable_crossrange_targeted_monte_carlo,
        Float64(crossrange_target_panel_cant_deg),
        Float64(crossrange_target_differential_fraction),
        Float64(trim_panel_area_fraction_of_deployed),
        Float64(panel_system_areal_density_kg_m2),
        0.5,
        1e-8,
        1e-8,
        2_000.0,
        0.5,
        Float64(representative_beta_high),
        Float64(representative_beta_ratio),
        Float64(representative_h_jettison_m),
        representative_switches_m === nothing ? [100e3, 80e3, 60e3, 40e3, 20e3] : Float64.(collect(representative_switches_m)),
        generate_plots,
    )
end

function with_geometry(config::PrelimConfig, geometry::VehicleGeometryConfig; output_root::Union{Nothing, String}=nothing)
    return default_config(
        output_root=output_root === nothing ? config.output_root : output_root,
        atmosphere_mode=config.atmosphere_mode,
        use_gram_surrogate=config.use_gram_surrogate,
        run_secondary_sweeps=config.run_secondary_sweeps,
        include_reverse_bangbang_sweep=config.include_reverse_bangbang_sweep,
        enable_crossrange_proxy=config.enable_crossrange_proxy,
        enable_target_range_guidance=config.enable_target_range_guidance,
        enable_monte_carlo=config.enable_monte_carlo,
        generate_plots=config.generate_plots,
        h_switch_grid_m=config.h_switch_grid_m,
        beta_high_targets=config.beta_high_targets,
        beta_ratios=config.beta_ratios,
        h_jettison_grid_m=config.h_jettison_grid_m,
        representative_switches_m=config.representative_switches_m,
        crossrange_sigmas=config.crossrange_sigmas,
        representative_beta_high=config.representative_beta_high,
        representative_beta_ratio=config.representative_beta_ratio,
        representative_h_jettison_m=config.representative_h_jettison_m,
        target_range_fraction=config.target_range_fraction,
        target_range_tolerance_km=config.target_range_tolerance_km,
        target_range_max_iterations=config.target_range_max_iterations,
        monte_carlo_samples=config.monte_carlo_samples,
        monte_carlo_seed=config.monte_carlo_seed,
        monte_carlo_density_scale_sigma=config.monte_carlo_density_scale_sigma,
        monte_carlo_V0_sigma_mps=config.monte_carlo_V0_sigma_mps,
        monte_carlo_gamma0_sigma_deg=config.monte_carlo_gamma0_sigma_deg,
        monte_carlo_beta_sigma_fraction=config.monte_carlo_beta_sigma_fraction,
        monte_carlo_h0_sigma_m=config.monte_carlo_h0_sigma_m,
        monte_carlo_crossrange_sigma=config.monte_carlo_crossrange_sigma,
        monte_carlo_use_winds=config.monte_carlo_use_winds,
        lateral_guidance_crossrange_scale_km=config.lateral_guidance_crossrange_scale_km,
        lateral_guidance_crossrange_rate_scale_mps=config.lateral_guidance_crossrange_rate_scale_mps,
        enable_crossrange_targeted_monte_carlo=config.enable_crossrange_targeted_monte_carlo,
        crossrange_target_panel_cant_deg=config.crossrange_target_panel_cant_deg,
        crossrange_target_differential_fraction=config.crossrange_target_differential_fraction,
        trim_panel_area_fraction_of_deployed=config.trim_panel_area_fraction_of_deployed,
        panel_system_areal_density_kg_m2=config.panel_system_areal_density_kg_m2,
        geometry=geometry,
    )
end

function with_trim_fraction(config::PrelimConfig, trim_panel_area_fraction_of_deployed::Real; output_root::Union{Nothing, String}=nothing)
    return default_config(
        output_root=output_root === nothing ? config.output_root : output_root,
        atmosphere_mode=config.atmosphere_mode,
        use_gram_surrogate=config.use_gram_surrogate,
        run_secondary_sweeps=config.run_secondary_sweeps,
        include_reverse_bangbang_sweep=config.include_reverse_bangbang_sweep,
        enable_crossrange_proxy=config.enable_crossrange_proxy,
        enable_target_range_guidance=config.enable_target_range_guidance,
        enable_monte_carlo=config.enable_monte_carlo,
        generate_plots=config.generate_plots,
        h_switch_grid_m=config.h_switch_grid_m,
        beta_high_targets=config.beta_high_targets,
        beta_ratios=config.beta_ratios,
        h_jettison_grid_m=config.h_jettison_grid_m,
        representative_switches_m=config.representative_switches_m,
        crossrange_sigmas=config.crossrange_sigmas,
        representative_beta_high=config.representative_beta_high,
        representative_beta_ratio=config.representative_beta_ratio,
        representative_h_jettison_m=config.representative_h_jettison_m,
        target_range_fraction=config.target_range_fraction,
        target_range_tolerance_km=config.target_range_tolerance_km,
        target_range_max_iterations=config.target_range_max_iterations,
        monte_carlo_samples=config.monte_carlo_samples,
        monte_carlo_seed=config.monte_carlo_seed,
        monte_carlo_density_scale_sigma=config.monte_carlo_density_scale_sigma,
        monte_carlo_V0_sigma_mps=config.monte_carlo_V0_sigma_mps,
        monte_carlo_gamma0_sigma_deg=config.monte_carlo_gamma0_sigma_deg,
        monte_carlo_beta_sigma_fraction=config.monte_carlo_beta_sigma_fraction,
        monte_carlo_h0_sigma_m=config.monte_carlo_h0_sigma_m,
        monte_carlo_crossrange_sigma=config.monte_carlo_crossrange_sigma,
        monte_carlo_use_winds=config.monte_carlo_use_winds,
        lateral_guidance_crossrange_scale_km=config.lateral_guidance_crossrange_scale_km,
        lateral_guidance_crossrange_rate_scale_mps=config.lateral_guidance_crossrange_rate_scale_mps,
        enable_crossrange_targeted_monte_carlo=config.enable_crossrange_targeted_monte_carlo,
        crossrange_target_panel_cant_deg=config.crossrange_target_panel_cant_deg,
        crossrange_target_differential_fraction=config.crossrange_target_differential_fraction,
        trim_panel_area_fraction_of_deployed=trim_panel_area_fraction_of_deployed,
        panel_system_areal_density_kg_m2=config.panel_system_areal_density_kg_m2,
        geometry=config.geometry,
    )
end
