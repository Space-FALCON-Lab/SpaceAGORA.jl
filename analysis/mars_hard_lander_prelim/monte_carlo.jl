using Random

struct MonteCarloResults
    samples_df::DataFrame
    summary_df::DataFrame
end

function _clip_positive(x::Float64)
    return max(x, eps())
end

function _footprint_ellipse_metrics(x_km::Vector{Float64}, y_km::Vector{Float64})
    n = length(x_km)
    n == length(y_km) || throw(ArgumentError("Footprint vectors must have the same length."))
    if n <= 1
        return (major_axis_km=0.0, minor_axis_km=0.0, azimuth_deg=0.0)
    end
    x_centered = x_km .- mean(x_km)
    y_centered = y_km .- mean(y_km)
    cov_xx = sum(abs2, x_centered) / (n - 1)
    cov_yy = sum(abs2, y_centered) / (n - 1)
    cov_xy = sum(x_centered .* y_centered) / (n - 1)
    cov_mat = [cov_xx cov_xy; cov_xy cov_yy]
    eig = eigen(Symmetric(cov_mat))
    λ = max.(eig.values, 0.0)
    order = sortperm(λ; rev=true)
    λ_major = λ[order[1]]
    λ_minor = λ[order[2]]
    vec_major = eig.vectors[:, order[1]]
    scale = sqrt(5.991464547107979) # 95% confidence ellipse in 2D
    return (
        major_axis_km=2.0 * scale * sqrt(λ_major),
        minor_axis_km=2.0 * scale * sqrt(λ_minor),
        azimuth_deg=rad2deg(atan(vec_major[2], vec_major[1])),
    )
end

function _policy_summary_rows(samples_df::DataFrame)
    rows = NamedTuple[]
    body_std = NaN
    body_interval_width = NaN
    grouped = groupby(samples_df, :policy)
    for group in grouped
        x = Vector{Float64}(group.impact_downrange_km)
        y = Vector{Float64}(group.impact_crossrange_km)
        footprint = _footprint_ellipse_metrics(x, y)
        q_low = quantile(x, 0.025)
        q_high = quantile(x, 0.975)
        row = (
            policy=String(group.policy[1]),
            sample_count=nrow(group),
            mean_impact_downrange_km=mean(x),
            std_impact_downrange_km=std(x),
            p02_5_impact_downrange_km=q_low,
            p97_5_impact_downrange_km=q_high,
            interval_width_km=q_high - q_low,
            mean_impact_crossrange_km=mean(y),
            std_impact_crossrange_km=std(y),
            ellipse_major_axis_km=footprint.major_axis_km,
            ellipse_minor_axis_km=footprint.minor_axis_km,
            ellipse_azimuth_deg=footprint.azimuth_deg,
            successful_guidance_fraction=mean(group.guidance_success),
        )
        if row.policy == "body_only_monte_carlo"
            body_std = row.std_impact_downrange_km
            body_interval_width = row.interval_width_km
        end
        push!(rows, row)
    end
    summary_df = DataFrame(rows)
    if any(summary_df.policy .== "guided_targeted_optimistic")
        for i in 1:nrow(summary_df)
            if summary_df.policy[i] == "guided_targeted_optimistic" && isfinite(body_std) && body_std > 0.0
                summary_df[!, :std_reduction_factor_vs_body_only] = fill(NaN, nrow(summary_df))
                summary_df[!, :interval_reduction_factor_vs_body_only] = fill(NaN, nrow(summary_df))
                break
            end
        end
        if :std_reduction_factor_vs_body_only in propertynames(summary_df)
            for i in 1:nrow(summary_df)
                summary_df.std_reduction_factor_vs_body_only[i] = body_std / max(summary_df.std_impact_downrange_km[i], eps())
                summary_df.interval_reduction_factor_vs_body_only[i] = body_interval_width / max(summary_df.interval_width_km[i], eps())
            end
        end
    end
    return summary_df
end

function _crossrange_target_panel_command(config::PrelimConfig, favored_side_sign::Float64)
    return make_differential_panel_command(
        config;
        favored_side_sign=favored_side_sign,
        differential_fraction=config.crossrange_target_differential_fraction,
        cant_deg=config.crossrange_target_panel_cant_deg,
    )
end

function _trajectory_for_switch_3d_panel(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
    h_switch_m::Float64,
    h_jettison_m::Float64,
    panel_command::DifferentialPanelCommand;
    save_trajectory::Bool=false,
    h0_m::Real=config.h0_m,
    V0_mps::Real=config.V0_mps,
    gamma0_deg::Real=config.gamma0_deg,
    density_scale::Real=1.0,
    use_winds::Bool=config.monte_carlo_use_winds,
)
    return solve_entry_trajectory_3d(
        config,
        adapter,
        aero_case,
        ControlPolicy("crossrange_targeted", :bang_bang, h_switch_m, h_jettison_m);
        save_trajectory=save_trajectory,
        h0_m=h0_m,
        V0_mps=V0_mps,
        gamma0_deg=gamma0_deg,
        density_scale=density_scale,
        use_winds=use_winds,
        panel_command_override=panel_command,
    )
end

function _find_switch_altitude_for_target_range_3d_panel(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
    h_jettison_m::Float64,
    target_range_km::Float64,
    panel_command::DifferentialPanelCommand;
    h_switch_grid_m::AbstractVector{<:Real}=config.h_switch_grid_m,
    h0_m::Real=config.h0_m,
    V0_mps::Real=config.V0_mps,
    gamma0_deg::Real=config.gamma0_deg,
    density_scale::Real=1.0,
    tolerance_km::Real=config.target_range_tolerance_km,
    max_iterations::Integer=config.target_range_max_iterations,
    save_trajectory::Bool=false,
    use_winds::Bool=config.monte_carlo_use_winds,
)
    h_grid = sort(Float64.(collect(h_switch_grid_m)))
    coarse_ranges = Float64[]
    for h_switch_m in h_grid
        traj = _trajectory_for_switch_3d_panel(
            config,
            adapter,
            aero_case,
            h_switch_m,
            h_jettison_m,
            panel_command;
            save_trajectory=false,
            h0_m=h0_m,
            V0_mps=V0_mps,
            gamma0_deg=gamma0_deg,
            density_scale=density_scale,
            use_winds=use_winds,
        )
        push!(coarse_ranges, traj.summary.impact_downrange_km)
    end
    f_grid = coarse_ranges .- Float64(target_range_km)
    nearest_idx = argmin(abs.(f_grid))
    nearest_h_m = h_grid[nearest_idx]
    nearest_result = _trajectory_for_switch_3d_panel(
        config,
        adapter,
        aero_case,
        nearest_h_m,
        h_jettison_m,
        panel_command;
        save_trajectory=save_trajectory,
        h0_m=h0_m,
        V0_mps=V0_mps,
        gamma0_deg=gamma0_deg,
        density_scale=density_scale,
        use_winds=use_winds,
    )
    nearest_error_km = nearest_result.summary.impact_downrange_km - target_range_km
    abs(nearest_error_km) <= tolerance_km && return TargetGuidanceResult(
        true,
        :grid_hit,
        target_range_km,
        nearest_result.summary.impact_downrange_km,
        nearest_error_km,
        nearest_h_m,
        h_jettison_m,
        nearest_h_m,
        nearest_h_m,
        0,
        TrajectoryResult(nearest_result.summary, nearest_result.trajectory),
    )

    bracket = _find_root_bracket(h_grid, f_grid)
    if bracket === nothing
        status = target_range_km < minimum(coarse_ranges) ? :clamped_low : :clamped_high
        return TargetGuidanceResult(
            false,
            status,
            target_range_km,
            nearest_result.summary.impact_downrange_km,
            nearest_error_km,
            nearest_h_m,
            h_jettison_m,
            h_grid[1],
            h_grid[end],
            0,
            TrajectoryResult(nearest_result.summary, nearest_result.trajectory),
        )
    end

    i_lo, i_hi = bracket
    if i_lo == i_hi
        return TargetGuidanceResult(
            true,
            :grid_exact,
            target_range_km,
            nearest_result.summary.impact_downrange_km,
            nearest_error_km,
            h_grid[i_lo],
            h_jettison_m,
            h_grid[i_lo],
            h_grid[i_lo],
            0,
            TrajectoryResult(nearest_result.summary, nearest_result.trajectory),
        )
    end

    h_lo = h_grid[i_lo]
    h_hi = h_grid[i_hi]
    f_lo = f_grid[i_lo]
    f_hi = f_grid[i_hi]
    final_result = nearest_result
    final_h_m = nearest_h_m
    final_error_km = nearest_error_km
    for iter in 1:Int(max_iterations)
        h_mid = 0.5 * (h_lo + h_hi)
        traj_mid = _trajectory_for_switch_3d_panel(
            config,
            adapter,
            aero_case,
            h_mid,
            h_jettison_m,
            panel_command;
            save_trajectory=save_trajectory,
            h0_m=h0_m,
            V0_mps=V0_mps,
            gamma0_deg=gamma0_deg,
            density_scale=density_scale,
            use_winds=use_winds,
        )
        f_mid = traj_mid.summary.impact_downrange_km - target_range_km
        final_result = traj_mid
        final_h_m = h_mid
        final_error_km = f_mid
        abs(f_mid) <= tolerance_km && return TargetGuidanceResult(
            true,
            :converged,
            target_range_km,
            traj_mid.summary.impact_downrange_km,
            f_mid,
            h_mid,
            h_jettison_m,
            h_lo,
            h_hi,
            iter,
            TrajectoryResult(traj_mid.summary, traj_mid.trajectory),
        )
        if f_lo * f_mid <= 0.0
            h_hi = h_mid
            f_hi = f_mid
        else
            h_lo = h_mid
            f_lo = f_mid
        end
    end

    return TargetGuidanceResult(
        false,
        :max_iterations,
        target_range_km,
        final_result.summary.impact_downrange_km,
        final_error_km,
        final_h_m,
        h_jettison_m,
        h_lo,
        h_hi,
        Int(max_iterations),
        TrajectoryResult(final_result.summary, final_result.trajectory),
    )
end

function _crossrange_guidance_score(result_3d::Trajectory3DResult, target_range_km::Float64)
    range_error_km = result_3d.summary.impact_downrange_km - target_range_km
    return hypot(range_error_km, result_3d.summary.impact_crossrange_km)
end

function run_small_monte_carlo(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    nominal_aero_case::CalibratedAeroCase,
    target_result::TargetGuidanceResult,
)::MonteCarloResults
    rng = MersenneTwister(config.monte_carlo_seed)
    rows = NamedTuple[]
    for sample_id in 1:config.monte_carlo_samples
        density_scale = _clip_positive(1.0 + config.monte_carlo_density_scale_sigma * randn(rng))
        V0_mps = config.V0_mps + config.monte_carlo_V0_sigma_mps * randn(rng)
        gamma0_deg = config.gamma0_deg + config.monte_carlo_gamma0_sigma_deg * randn(rng)
        h0_m = config.h0_m + config.monte_carlo_h0_sigma_m * randn(rng)
        beta_scale = _clip_positive(1.0 + config.monte_carlo_beta_sigma_fraction * randn(rng))
        dispersed_case = scale_aero_case_mass(nominal_aero_case, beta_scale)
        body_result_3d = solve_entry_trajectory_3d(
            config,
            adapter,
            dispersed_case,
            ControlPolicy("body_only_monte_carlo", :body_only, NaN, NaN);
            save_trajectory=false,
            h0_m=h0_m,
            V0_mps=V0_mps,
            gamma0_deg=gamma0_deg,
            density_scale=density_scale,
            use_winds=config.monte_carlo_use_winds,
        )
        push!(rows, (
            sample_id=sample_id,
            policy="body_only_monte_carlo",
            guidance_success=true,
            guidance_status="body_only",
            density_scale=density_scale,
            V0_mps=V0_mps,
            gamma0_deg=gamma0_deg,
            h0_m=h0_m,
            beta_scale=beta_scale,
            h_switch_km=NaN,
            target_range_km=target_result.target_range_km,
            impact_downrange_km=body_result_3d.summary.impact_downrange_km,
            impact_crossrange_km=body_result_3d.summary.impact_crossrange_km,
            impact_velocity_mps=body_result_3d.summary.impact_velocity_mps,
            peak_total_decel_earth_g=body_result_3d.summary.peak_total_decel_earth_g,
            range_error_km=body_result_3d.summary.impact_downrange_km - target_result.target_range_km,
            chosen_side=missing,
            panel_cant_deg=NaN,
            differential_fraction=NaN,
        ))

        guided = find_switch_altitude_for_target_range(
            config,
            adapter,
            dispersed_case,
            target_result.h_jettison_m,
            target_result.target_range_km;
            h0_m=h0_m,
            V0_mps=V0_mps,
            gamma0_deg=gamma0_deg,
            density_scale=density_scale,
            save_trajectory=false,
        )
        guided_result_3d = solve_entry_trajectory_3d(
            config,
            adapter,
            dispersed_case,
            ControlPolicy("guided_targeted_optimistic", :bang_bang, guided.h_switch_m, target_result.h_jettison_m);
            save_trajectory=false,
            h0_m=h0_m,
            V0_mps=V0_mps,
            gamma0_deg=gamma0_deg,
            density_scale=density_scale,
            use_winds=config.monte_carlo_use_winds,
        )
        push!(rows, (
            sample_id=sample_id,
            policy="guided_targeted_optimistic",
            guidance_success=guided.success,
            guidance_status=String(guided.status),
            density_scale=density_scale,
            V0_mps=V0_mps,
            gamma0_deg=gamma0_deg,
            h0_m=h0_m,
            beta_scale=beta_scale,
            h_switch_km=guided.h_switch_m / 1e3,
            target_range_km=guided.target_range_km,
            impact_downrange_km=guided_result_3d.summary.impact_downrange_km,
            impact_crossrange_km=guided_result_3d.summary.impact_crossrange_km,
            impact_velocity_mps=guided_result_3d.summary.impact_velocity_mps,
            peak_total_decel_earth_g=guided_result_3d.summary.peak_total_decel_earth_g,
            range_error_km=guided_result_3d.summary.impact_downrange_km - guided.target_range_km,
            chosen_side=missing,
            panel_cant_deg=NaN,
            differential_fraction=NaN,
        ))

        if config.enable_crossrange_targeted_monte_carlo
            candidate_rows = NamedTuple[]
            for favored_side_sign in (-1.0, 1.0)
                panel_command = _crossrange_target_panel_command(config, favored_side_sign)
                crossrange_guided = _find_switch_altitude_for_target_range_3d_panel(
                    config,
                    adapter,
                    dispersed_case,
                    target_result.h_jettison_m,
                    target_result.target_range_km,
                    panel_command;
                    h0_m=h0_m,
                    V0_mps=V0_mps,
                    gamma0_deg=gamma0_deg,
                    density_scale=density_scale,
                    save_trajectory=false,
                    use_winds=config.monte_carlo_use_winds,
                )
                crossrange_result_3d = solve_entry_trajectory_3d(
                    config,
                    adapter,
                    dispersed_case,
                    ControlPolicy("guided_targeted_crossrange", :bang_bang, crossrange_guided.h_switch_m, target_result.h_jettison_m);
                    save_trajectory=false,
                    h0_m=h0_m,
                    V0_mps=V0_mps,
                    gamma0_deg=gamma0_deg,
                    density_scale=density_scale,
                    use_winds=config.monte_carlo_use_winds,
                    panel_command_override=panel_command,
                )
                push!(candidate_rows, (
                    favored_side_sign=favored_side_sign,
                    panel_command=panel_command,
                    guidance=crossrange_guided,
                    result_3d=crossrange_result_3d,
                    score=_crossrange_guidance_score(crossrange_result_3d, target_result.target_range_km),
                ))
            end
            chosen_idx = argmin(getproperty.(candidate_rows, :score))
            chosen = candidate_rows[chosen_idx]
            chosen_side_label = chosen.favored_side_sign < 0.0 ? "left" : "right"
            push!(rows, (
                sample_id=sample_id,
                policy="guided_targeted_crossrange_bilateral",
                guidance_success=chosen.guidance.success,
                guidance_status=String(chosen.guidance.status),
                density_scale=density_scale,
                V0_mps=V0_mps,
                gamma0_deg=gamma0_deg,
                h0_m=h0_m,
                beta_scale=beta_scale,
                h_switch_km=chosen.guidance.h_switch_m / 1e3,
                target_range_km=chosen.guidance.target_range_km,
                impact_downrange_km=chosen.result_3d.summary.impact_downrange_km,
                impact_crossrange_km=chosen.result_3d.summary.impact_crossrange_km,
                impact_velocity_mps=chosen.result_3d.summary.impact_velocity_mps,
                peak_total_decel_earth_g=chosen.result_3d.summary.peak_total_decel_earth_g,
                range_error_km=chosen.result_3d.summary.impact_downrange_km - chosen.guidance.target_range_km,
                chosen_side=chosen_side_label,
                panel_cant_deg=config.crossrange_target_panel_cant_deg,
                differential_fraction=config.crossrange_target_differential_fraction,
            ))
        end
    end
    samples_df = DataFrame(rows)
    summary_df = _policy_summary_rows(samples_df)
    return MonteCarloResults(samples_df, summary_df)
end
