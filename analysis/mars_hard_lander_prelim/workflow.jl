struct PrelimResults
    config::PrelimConfig
    summary_df::DataFrame
    authority_df::DataFrame
    trade_metrics_df::DataFrame
    local_effectiveness_df::DataFrame
    crossrange_df::DataFrame
    crossrange_sensitivity_df::DataFrame
    alpha_sensitivity_df::DataFrame
    target_guidance_df::DataFrame
    monte_carlo_samples_df::DataFrame
    monte_carlo_summary_df::DataFrame
    output_root::String
    top_line::NamedTuple
end

function _prepare_output_dirs(config::PrelimConfig)
    rm(config.output_root; recursive=true, force=true)
    mkpath(config.output_root)
    mkpath(config.trajectory_dir)
    mkpath(config.plot_dir)
    return nothing
end

function _write_trajectory_csvs(config::PrelimConfig, representative::Dict{String, DataFrame})
    for (name, df) in representative
        CSV.write(joinpath(config.trajectory_dir, "$(name).csv"), df)
    end
    return nothing
end

function _representative_nominal_subset(summary_df::DataFrame, config::PrelimConfig)
    return filter(
        row -> row.policy == "bang_bang" &&
            _approx_eq(row.target_beta_high_kg_m2, config.representative_beta_high) &&
            _approx_eq(row.target_beta_ratio, config.representative_beta_ratio) &&
            _approx_eq(row.h_jettison_km, config.representative_h_jettison_m / 1e3),
        summary_df,
    )
end

function _select_top_line(summary_df::DataFrame, authority_df::DataFrame, config::PrelimConfig)
    nominal_authority = filter(
        row -> _approx_eq(row.target_beta_high_kg_m2, config.representative_beta_high) &&
            _approx_eq(row.target_beta_ratio, config.representative_beta_ratio) &&
            _approx_eq(row.h_jettison_km, config.representative_h_jettison_m / 1e3),
        authority_df,
    )
    nominal_row = nrow(nominal_authority) > 0 ? nominal_authority[argmax(nominal_authority.downrange_authority_km), :] : authority_df[argmax(authority_df.downrange_authority_km), :]
    nominal_summary = filter(
        row -> row.policy == "bang_bang" &&
            _approx_eq(row.target_beta_high_kg_m2, nominal_row.target_beta_high_kg_m2) &&
            _approx_eq(row.target_beta_ratio, nominal_row.target_beta_ratio) &&
            _approx_eq(row.h_jettison_km, nominal_row.h_jettison_km),
        summary_df,
    )
    min_idx = argmin(nominal_summary.impact_downrange_km)
    max_idx = argmax(nominal_summary.impact_downrange_km)
    return (
        beta_high_kg_m2=nominal_row.target_beta_high_kg_m2,
        beta_ratio=nominal_row.target_beta_ratio,
        h_jettison_km=nominal_row.h_jettison_km,
        downrange_authority_km=nominal_row.downrange_authority_km,
        min_impact_velocity_mps=nominal_summary.impact_velocity_mps[min_idx],
        max_impact_velocity_mps=nominal_summary.impact_velocity_mps[max_idx],
        min_peak_total_decel_earth_g=nominal_summary.peak_total_decel_earth_g[min_idx],
        max_peak_total_decel_earth_g=nominal_summary.peak_total_decel_earth_g[max_idx],
    )
end

function _write_readme(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    top_line::NamedTuple,
    authority_df::DataFrame,
    crossrange_df::DataFrame,
    crossrange_sensitivity_df::DataFrame,
    alpha_sensitivity_df::DataFrame,
    target_guidance_df::DataFrame,
    monte_carlo_summary_df::DataFrame,
)
    best_idx = argmax(authority_df.downrange_authority_km)
    best = authority_df[best_idx, :]
    shield_mode = occursin("SHIELD", config.geometry.body_label) || occursin("SHIELD", config.geometry.deployed_label)
    shield_deployed_area_m2 = deployed_drag_skirt_equivalent_area(config)
    switch_band_min_km = minimum(config.h_switch_grid_m) / 1e3
    switch_band_max_km = maximum(config.h_switch_grid_m) / 1e3
    open(joinpath(config.output_root, "README_prelim.md"), "w") do io
        println(io, shield_mode ? "# Published-Size SHIELD Reanalysis" : "# Mars Hard-Lander Preliminary Analysis")
        println(io)
        println(io, "## Included")
        println(io, "- Mars spherical gravity with a planar 3DOF point-mass entry model")
        println(io, "- MarsGRAM-backed density/temperature when available, with exponential fallback for no-GRAM runs")
        println(io, "- Geometry-derived continuum reference-shape aerodynamics for a sphere-cone body plus symmetric flat-plate deployed surfaces")
        println(io, "- Discrete aero-state switching: body-only, deployed-until-jettison, and bang-bang by switch altitude")
        if shield_mode
            println(io, "- Deterministic sweeps over subsonic deployment switch altitude, with stowed and deployed reference states derived from the published-size SHIELD surrogate geometry")
        else
            println(io, "- Deterministic sweeps over body-only ballistic coefficient, ballistic-coefficient ratio, jettison altitude, and switch altitude")
        end
        println(io)
        println(io, "## Not Included")
        println(io, "- Full 6DOF attitude dynamics")
        println(io, "- Transition-flow truth modeling or DSMC/CFD validation")
        println(io, "- Terrain-relative navigation, parachutes, retropropulsion, or powered descent")
        println(io, "- High-fidelity panel geometry, hinge loads, or structural dynamics")
        println(io)
        println(io, "## Top-Line Results")
        println(io, @sprintf(
            "- Representative case: β_high = %.0f kg/m², β_ratio = %.1f, h_jettison = %.0f km -> downrange authority %.2f km, impact velocity %.0f-%.0f m/s, peak aero load %.2f-%.2f g.",
            top_line.beta_high_kg_m2,
            top_line.beta_ratio,
            top_line.h_jettison_km,
            top_line.downrange_authority_km,
            top_line.min_impact_velocity_mps,
            top_line.max_impact_velocity_mps,
            top_line.min_peak_total_decel_earth_g,
            top_line.max_peak_total_decel_earth_g,
        ))
        println(io, @sprintf(
            "- Best authority found in this sweep: %.2f km at β_high = %.0f kg/m², β_ratio = %.1f, h_jettison = %.0f km.",
            best.downrange_authority_km,
            best.target_beta_high_kg_m2,
            best.target_beta_ratio,
            best.h_jettison_km,
        ))
        if !isempty(alpha_sensitivity_df)
            alpha_zero = filter(row -> _approx_eq(row.alpha_body_deg, 0.0), alpha_sensitivity_df)
            alpha_best = alpha_sensitivity_df[argmax(alpha_sensitivity_df.downrange_authority_km), :]
            if nrow(alpha_zero) > 0
                println(io, @sprintf(
                    "- Body AoA sensitivity (representative case): authority is %.2f km at α_body = 0 deg and peaks at %.2f km near α_body = %.0f deg.",
                    alpha_zero.downrange_authority_km[1],
                    alpha_best.downrange_authority_km,
                    alpha_best.alpha_body_deg,
                ))
            end
        end
        if !isempty(crossrange_df)
            best_cross_idx = argmax(crossrange_df.impact_crossrange_km)
            println(io, @sprintf(
                "- Crossrange proxy (representative bang-bang case): σ = %.2f gives %.2f km lateral displacement estimate.",
                crossrange_df.sigma[best_cross_idx],
                crossrange_df.impact_crossrange_km[best_cross_idx],
            ))
        end
        if !isempty(crossrange_sensitivity_df)
            best_cross_sens_idx = argmax(crossrange_sensitivity_df.abs_impact_crossrange_km)
            best_cross_sens = crossrange_sensitivity_df[best_cross_sens_idx, :]
            println(io, @sprintf(
                "- Differential panel crossrange sensitivity: best deterministic |crossrange| = %.2f km at α_body = %.0f deg, cant = %.0f deg, and differential fraction %.2f (%s).",
                best_cross_sens.abs_impact_crossrange_km,
                best_cross_sens.alpha_body_deg,
                best_cross_sens.panel_cant_deg,
                best_cross_sens.differential_deployment_fraction,
                best_cross_sens.differential_deployment_label,
            ))
        end
        if !isempty(target_guidance_df)
            println(io, @sprintf(
                "- Target-range switch solve: target %.2f km reached with h_switch = %.2f km and range error %.3f km (%s).",
                target_guidance_df.target_range_km[1],
                target_guidance_df.h_switch_km[1],
                target_guidance_df.range_error_km[1],
                target_guidance_df.status[1],
            ))
        end
        if !isempty(monte_carlo_summary_df)
            guided = filter(row -> row.policy == "guided_targeted_optimistic", monte_carlo_summary_df)
            body = filter(row -> row.policy == "body_only_monte_carlo", monte_carlo_summary_df)
            if nrow(guided) > 0 && nrow(body) > 0
                println(io, @sprintf(
                    "- Monte Carlo (N = %d, optimistic true-state switch solve): body-only σ_range = %.2f km vs guided σ_range = %.2f km, reduction factor %.2f.",
                    Int(guided.sample_count[1]),
                    body.std_impact_downrange_km[1],
                    guided.std_impact_downrange_km[1],
                    guided.std_reduction_factor_vs_body_only[1],
                ))
                println(io, @sprintf(
                    "- Guided 3D landing-footprint 95%% ellipse: major axis %.2f km, minor axis %.2f km.",
                    guided.ellipse_major_axis_km[1],
                    guided.ellipse_minor_axis_km[1],
                ))
            end
            crossrange_guided = filter(row -> row.policy == "guided_targeted_crossrange_bilateral", monte_carlo_summary_df)
            if nrow(crossrange_guided) > 0
                println(io, @sprintf(
                    "- Crossrange-aware Monte Carlo (cant %.0f deg, differential fraction %.2f): 95%% ellipse %.2f km x %.2f km, σ_range = %.2f km, σ_crossrange = %.2f km.",
                    config.crossrange_target_panel_cant_deg,
                    config.crossrange_target_differential_fraction,
                    crossrange_guided.ellipse_major_axis_km[1],
                    crossrange_guided.ellipse_minor_axis_km[1],
                    crossrange_guided.std_impact_downrange_km[1],
                    crossrange_guided.std_impact_crossrange_km[1],
                ))
            end
        end
        println(io)
        println(io, "## Initial Conditions")
        println(io, @sprintf(
            "- Nominal entry state: h0 = %.0f km, V0 = %.2f km/s, γ0 = %.1f deg, θ0 = %.1f deg, latitude = %.1f deg.",
            config.h0_m / 1e3,
            config.V0_mps / 1e3,
            config.gamma0_deg,
            rad2deg(config.theta0_rad),
            config.lat_deg,
        ))
        println(io)
        println(io, "## Parameter Interpretation")
        println(io, "- `β_high` in this workflow is the body-only reference ballistic coefficient at the calibration condition, so changing `β_high` mainly changes vehicle mass class for a fixed geometry family.")
        if shield_mode
            println(io, "- In this SHIELD rerun, `β_low` is derived from the deployed drag-skirt-equivalent geometry rather than prescribed directly. `β_ratio = β_high / β_low` therefore reflects the contrast between the stowed 1.8 m body and the deployed 2.2 m / 0.75 m skirt surrogate.")
        else
            println(io, "- `β_ratio = β_high / β_low` is the main control-authority lever. In this workflow it is achieved by adding deployed plate area, not by changing the body shape.")
        end
        println(io, "- `h_switch` is the control variable: in the SHIELD rerun it sets the passive skirt-deployment altitude, while any future rim-flap actuation would be a separate post-deployment guidance event.")
        println(io, "- `h_jettison` is an architecture / operations parameter: it sets how long the extra area is available before the surfaces are discarded.")
        println(io, "- `β` is not constant over the trajectory. The tables report target and achieved reference `β` values at the calibration condition because `C_D` changes with Mach and configuration during entry.")
        println(io)
        println(io, "## Aerodynamic Construction")
        println(io, "- Body: SHIELD-informed passive hard-lander surrogate using a `:sphere_cone` with nose radius $(config.geometry.nose_radius_m) m, base radius $(config.geometry.base_radius_m) m, and a $(config.geometry.cone_half_angle_deg)-deg nominal cone half-angle assumption.")
        println(io, "- Deployed state: the same body plus $(config.geometry.panel_count) symmetric flat plates represented with the continuum `:flat_plate` reference-shape model.")
        println(io, "- Coefficients are computed separately for the body and the deployed plate set, then combined as `(C_D A)_total = Σ C_D,i A_i` and `(C_L A)_total = Σ C_L,i A_i`.")
        println(io, "- Pressure model: `$(config.geometry.hypersonic_pressure_model)`. Nominal SHIELD-informed body angle of attack = $(config.geometry.body_alpha_deg) deg. Fixed plate effective angle = $(config.geometry.panel_alpha_deg) deg.")
        if shield_mode && shield_deployed_area_m2 > 0.0
            println(io, "- The deployed equivalent flat-plate area is derived from a frustum drag-skirt surrogate built from the published-size dimensions: stowed diameter $(2.0 * config.geometry.base_radius_m) m, deployed diameter $(config.geometry.deployed_drag_surface_diameter_m) m, and skirt height $(config.geometry.drag_skirt_height_m) m. That gives $(round(shield_deployed_area_m2; digits=3)) m² total equivalent deployed area before any differential crossrange commands.")
        else
            println(io, "- Current deployed plate areas per vehicle are set by the target `β_ratio`; per-plate area in each run is reported in `prelim_summary.csv`.")
        end
        println(io, "- Current implementation uses $(config.geometry.panel_count) explicit symmetric plates in the force build-up, each with equal area.")
        println(io, "- Differential panel deployment in the crossrange sensitivity study means the left and right panels no longer match: the favored side stays fully deployed and canted, while the opposite side is partially retracted and left closer to the neutral broadside state.")
        println(io, "- Access2Mars constrains the hard-lander mission class (passive, uncontrolled entry, ≥5 kg payload, large landing footprint) but does not provide a fully specified SHIELD aeroshell geometry. The $(config.geometry.cone_half_angle_deg)-deg sphere-cone remains a surrogate assumption, not a recovered SHIELD CAD truth model.")
        println(io)
        println(io, "## Assumptions")
        println(io, "- Atmosphere path used in this run: `$(atmosphere_label(adapter))`")
        println(io, "- Deterministic downrange sweeps use zero winds even when MarsGRAM is used for density and temperature, so the main authority plots isolate aerodynamic control rather than weather dispersions.")
        println(io, @sprintf("- SHIELD skirt deployment is only swept across the body-only subsonic band `h_switch = %.1f–%.1f km`; the current model does not deploy the skirt while the stowed vehicle is still supersonic.", switch_band_min_km, switch_band_max_km))
        println(io, "- The Monte Carlo landing ellipse uses a 3D point-mass propagation with MarsGRAM winds and differential left/right panel commands.")
        println(io, "- The crossrange-aware Monte Carlo uses a fixed asymmetric panel architecture with cant = $(config.crossrange_target_panel_cant_deg) deg and differential fraction = $(config.crossrange_target_differential_fraction), then chooses left or right sign per dispersed sample.")
        println(io, "- The differential panel crossrange sweep is deterministic and zero-wind; it is intended to quantify geometric lateral authority, not weather-driven footprint spread.")
        println(io, "- Nominal SHIELD-informed body angle of attack: $(config.geometry.body_alpha_deg) deg")
        println(io, "- Fixed deployed-panel effective angle: $(config.geometry.panel_alpha_deg) deg")
        println(io, "- Illustrative panel-system areal density used for mass-efficiency trades: $(config.panel_system_areal_density_kg_m2) kg/m²")
        println(io, "- Impact g-load proxy assumes a stopping distance of $(config.impact_stop_distance_m) m")
        println(io, "- Monte Carlo guidance uses an optimistic true-dispersed-state switch solve.")
        println(io, "- The 3D footprint model is still preliminary: it uses a point-mass vehicle with simple differential panel canting, not a bank-controlled or full 6DOF attitude simulation.")
        println(io)
        println(io, "## Recommended Slide Subset")
        if shield_mode
            println(io, @sprintf(
                "- Use the published-size SHIELD surrogate case directly: `β_high = %.1f kg/m²`, `β_ratio = %.2f`, fixed entry mass `120 kg`, and `h_jettison = %.0f km`.",
                config.representative_beta_high,
                config.representative_beta_ratio,
                config.representative_h_jettison_m / 1e3,
            ))
            println(io, @sprintf("- Show the subsonic deployment band `h_switch = %.1f–%.1f km` as the main authority sweep and keep the body-only / fully deployed endpoints visible for context.", switch_band_min_km, switch_band_max_km))
            println(io, "- Lead with downrange authority and local control-effectiveness. Present crossrange and the landing ellipse as preliminary differential-panel extensions, not as mature SHIELD guidance results.")
        else
            println(io, @sprintf(
                "- Use `β_high = %.0f kg/m²` as the nominal vehicle class and show the configured `β_ratio` sweep around the representative value `%.2f`.",
                config.representative_beta_high,
                config.representative_beta_ratio,
            ))
            println(io, @sprintf(
                "- Use `h_jettison = %.0f km` as the nominal case and show the configured switch band `%.0f–%.0f km` as the control-authority sweep.",
                config.representative_h_jettison_m / 1e3,
                minimum(config.h_switch_grid_m) / 1e3,
                maximum(config.h_switch_grid_m) / 1e3,
            ))
            println(io, "- Lead with downrange authority and local control-effectiveness. Use the landing ellipse as a preliminary 3D footprint result, but label it clearly as a point-mass / differential-panel first cut rather than a final guidance-validation figure.")
        end
        println(io)
        println(io, "## Recommended Next Step")
        println(io, "- Keep this 3DOF/3D MarsGRAM workflow as the main result engine for the proposal package, then refine the lateral panel model and wind-aware target selection before any framework-native 6DOF integration.")
        println(io)
        println(io, "## Output Files")
        println(io, "- `prelim_summary.csv`")
        println(io, "- `authority_summary.csv`")
        println(io, "- `authority_trade_metrics.csv`")
        println(io, "- `local_control_effectiveness.csv`")
        println(io, "- `alpha_body_sensitivity.csv`")
        println(io, "- `crossrange_sensitivity.csv`")
        println(io, "- `crossrange_sensitivity.md`")
        println(io, "- `target_guidance_summary.csv`")
        println(io, "- `monte_carlo_summary.csv`")
        println(io, "- `monte_carlo_samples.csv`")
        println(io, "- `proposal_trade_table.md`")
        println(io, "- `trajectories/`")
        println(io, "- `plots/`")
    end
    return nothing
end

function run_prelim(config::PrelimConfig=default_config())
    _prepare_output_dirs(config)
    adapter = build_atmosphere_adapter(config)
    @info "Running Mars hard-lander preliminary analysis" output_root=config.output_root atmosphere=atmosphere_label(adapter)

    summary_df, authority_df, representative_trajectories, aero_cases = run_deterministic_sweeps(config, adapter)
    trade_metrics_df = compute_authority_trade_metrics(config, summary_df, authority_df)

    CSV.write(joinpath(config.output_root, "prelim_summary.csv"), summary_df)
    CSV.write(joinpath(config.output_root, "authority_summary.csv"), authority_df)
    CSV.write(joinpath(config.output_root, "authority_trade_metrics.csv"), trade_metrics_df)
    _write_trajectory_csvs(config, representative_trajectories)

    nominal_bang_subset = _representative_nominal_subset(summary_df, config)
    body_traj_key = _named_trajectory_key(ControlPolicy("body_only", :body_only, NaN, NaN))
    body_traj = get(representative_trajectories, body_traj_key, DataFrame())
    local_effectiveness_df = isempty(nominal_bang_subset) || isempty(body_traj) ?
        DataFrame() :
        compute_local_effectiveness(nominal_bang_subset, body_traj)
    CSV.write(joinpath(config.output_root, "local_control_effectiveness.csv"), local_effectiveness_df)

    crossrange_df = DataFrame()
    crossrange_sensitivity_df = DataFrame()
    alpha_sensitivity_df = run_body_alpha_sensitivity(config, adapter)
    CSV.write(joinpath(config.output_root, "alpha_body_sensitivity.csv"), alpha_sensitivity_df)
    target_guidance_df = DataFrame()
    monte_carlo_samples_df = DataFrame()
    monte_carlo_summary_df = DataFrame()
    if config.enable_crossrange_proxy && nrow(nominal_bang_subset) > 0
        best_nominal_idx = argmax(nominal_bang_subset.impact_downrange_km)
        best_nominal = nominal_bang_subset[best_nominal_idx, :]
        nominal_case = aero_cases[(best_nominal.target_beta_high_kg_m2, best_nominal.target_beta_ratio)]
        bang_policy = ControlPolicy(
            "bang_bang",
            :bang_bang,
            best_nominal.h_switch_km * 1e3,
            best_nominal.h_jettison_km * 1e3,
        )
        bang_result = solve_entry_trajectory(config, adapter, nominal_case, bang_policy; save_trajectory=true)
        representative_trajectories["best_nominal_bang_bang"] = bang_result.trajectory
        CSV.write(joinpath(config.trajectory_dir, "best_nominal_bang_bang.csv"), bang_result.trajectory)
        crossrange_df = run_crossrange_proxy(bang_result.trajectory, config.crossrange_sigmas)
        CSV.write(joinpath(config.output_root, "crossrange_proxy_summary.csv"), crossrange_df)
    end

    nominal_target_case = nothing
    nominal_target_result = nothing
    if config.enable_target_range_guidance && !isempty(authority_df)
        target_case = select_target_guidance_case(authority_df, config)
        nominal_target_case = target_case
        nominal_aero_case = aero_cases[(target_case.beta_high_kg_m2, target_case.beta_ratio)]
        target_result = find_switch_altitude_for_target_range(
            config,
            adapter,
            nominal_aero_case,
            target_case.h_jettison_km * 1e3,
            target_case.target_range_km;
            save_trajectory=true,
        )
        nominal_target_result = target_result
        target_guidance_df = DataFrame([(
            beta_high_kg_m2=target_case.beta_high_kg_m2,
            beta_ratio=target_case.beta_ratio,
            h_jettison_km=target_case.h_jettison_km,
            min_reachable_range_km=target_case.min_range_km,
            max_reachable_range_km=target_case.max_range_km,
            target_range_km=target_result.target_range_km,
            achieved_range_km=target_result.achieved_range_km,
            range_error_km=target_result.range_error_km,
            h_switch_km=target_result.h_switch_m / 1e3,
            bracket_low_km=target_result.bracket_low_m / 1e3,
            bracket_high_km=target_result.bracket_high_m / 1e3,
            h_jettison_result_km=target_result.h_jettison_m / 1e3,
            success=target_result.success,
            status=String(target_result.status),
            iterations=target_result.iterations,
        )])
        CSV.write(joinpath(config.output_root, "target_guidance_summary.csv"), target_guidance_df)
        if !isempty(target_result.trajectory_result.trajectory)
            representative_trajectories["targeted_nominal"] = target_result.trajectory_result.trajectory
            CSV.write(joinpath(config.trajectory_dir, "targeted_nominal.csv"), target_result.trajectory_result.trajectory)
        end

        crossrange_sensitivity_df = run_crossrange_sensitivity(config, adapter)
        CSV.write(joinpath(config.output_root, "crossrange_sensitivity.csv"), crossrange_sensitivity_df)
        write_crossrange_sensitivity_note(config, crossrange_sensitivity_df)

        if config.enable_monte_carlo
            monte_carlo = run_small_monte_carlo(config, adapter, nominal_aero_case, target_result)
            monte_carlo_samples_df = monte_carlo.samples_df
            monte_carlo_summary_df = monte_carlo.summary_df
            CSV.write(joinpath(config.output_root, "monte_carlo_samples.csv"), monte_carlo_samples_df)
            CSV.write(joinpath(config.output_root, "monte_carlo_summary.csv"), monte_carlo_summary_df)
        end
    end

    if config.generate_plots
        save_representative_timeseries_plots(config, representative_trajectories)
        save_authority_heatmap(config, authority_df)
        save_terminal_metric_plot(config, summary_df, authority_df)
        if !isempty(local_effectiveness_df)
            save_local_effectiveness_plot(config, local_effectiveness_df)
        end
        if !isempty(crossrange_df)
            save_crossrange_proxy_plot(config, crossrange_df)
        end
        if !isempty(alpha_sensitivity_df)
            save_alpha_body_sensitivity_plot(config, alpha_sensitivity_df)
        end
        if !isempty(crossrange_sensitivity_df)
            save_crossrange_sensitivity_plot(config, crossrange_sensitivity_df)
        end
        if !isempty(monte_carlo_samples_df)
            save_landing_ellipse_plot(config, monte_carlo_samples_df, monte_carlo_summary_df)
        end
    end

    top_line = _select_top_line(summary_df, authority_df, config)
    _write_readme(
        config,
        adapter,
        top_line,
        authority_df,
        crossrange_df,
        crossrange_sensitivity_df,
        alpha_sensitivity_df,
        target_guidance_df,
        monte_carlo_summary_df,
    )
    write_proposal_trade_table(config, trade_metrics_df, target_guidance_df, monte_carlo_summary_df)

    return PrelimResults(
        config,
        summary_df,
        authority_df,
        trade_metrics_df,
        local_effectiveness_df,
        crossrange_df,
        crossrange_sensitivity_df,
        alpha_sensitivity_df,
        target_guidance_df,
        monte_carlo_samples_df,
        monte_carlo_summary_df,
        config.output_root,
        top_line,
    )
end

function main(; kwargs...)
    config = default_config(; kwargs...)
    return run_prelim(config)
end
