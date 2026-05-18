function _prepare_plots_backend()
    if !haskey(Base.ENV, "GKSwstype")
        Base.ENV["GKSwstype"] = "100"
    end
    gr()
    default(
        size=(1100, 750),
        legend=:best,
        linewidth=2,
        left_margin=8Plots.mm,
        right_margin=8Plots.mm,
        top_margin=6Plots.mm,
        bottom_margin=8Plots.mm,
        guidefont=font(13),
        tickfont=font(11),
        titlefont=font(16),
        legendfont=font(10),
    )
    return nothing
end

function _trajectory_label(key::AbstractString)
    if startswith(key, "body_only")
        return "Body only"
    elseif startswith(key, "fixed_deployed")
        jet_match = Base.match(r"jet(\d+)", key)
        jet_km = jet_match === nothing ? "?" : string(parse(Int, jet_match.captures[1]))
        return "Always deployed until $(jet_km) km"
    elseif startswith(key, "bang_bang")
        sw_match = Base.match(r"sw(\d+)", key)
        jet_match = Base.match(r"jet(\d+)", key)
        sw_km = sw_match === nothing ? "?" : string(parse(Int, sw_match.captures[1]))
        jet_km = jet_match === nothing ? "?" : string(parse(Int, jet_match.captures[1]))
        return "Deploy below $(sw_km) km, jettison at $(jet_km) km"
    elseif startswith(key, "reverse_bang_bang")
        sw_match = Base.match(r"sw(\d+)", key)
        jet_match = Base.match(r"jet(\d+)", key)
        sw_km = sw_match === nothing ? "?" : string(parse(Int, sw_match.captures[1]))
        jet_km = jet_match === nothing ? "?" : string(parse(Int, jet_match.captures[1]))
        return "Deployed above $(sw_km) km, jettison at $(jet_km) km"
    end
    return replace(key, "_" => " ")
end

function save_representative_timeseries_plots(
    config::PrelimConfig,
    representative::Dict{String, DataFrame},
)
    _prepare_plots_backend()
    keys_sorted = sort(collect(keys(representative)))
    labels = _trajectory_label.(keys_sorted)

    p_alt = plot(xlabel="Downrange (km)", ylabel="Altitude (km)", title="Representative Trajectories: Altitude vs Downrange", legend=:outerright)
    p_vel = plot(xlabel="Velocity (m/s)", ylabel="Altitude (km)", title="Representative Trajectories: Velocity vs Altitude", legend=:outerright)
    p_q = plot(xlabel="Dynamic Pressure (Pa)", ylabel="Altitude (km)", title="Representative Trajectories: q vs Altitude", legend=:outerright)
    p_drag = plot(xlabel="Drag Acceleration (m/s²)", ylabel="Altitude (km)", title="Representative Trajectories: Drag Acceleration vs Altitude", legend=:outerright)

    for (plot_idx, key) in enumerate(keys_sorted)
        df = representative[key]
        label = labels[plot_idx]
        plot!(p_alt, df.downrange_km, df.altitude_km; label=label)
        plot!(p_vel, df.velocity_mps, df.altitude_km; label=label)
        plot!(p_q, df.q_pa, df.altitude_km; label=label)
        plot!(p_drag, df.drag_accel_mps2, df.altitude_km; label=label)
    end

    savefig(p_alt, joinpath(config.plot_dir, "altitude_vs_downrange.png"))
    savefig(p_vel, joinpath(config.plot_dir, "velocity_vs_altitude.png"))
    savefig(p_q, joinpath(config.plot_dir, "dynamic_pressure_vs_altitude.png"))
    savefig(p_drag, joinpath(config.plot_dir, "drag_accel_vs_altitude.png"))
end

function save_authority_heatmap(config::PrelimConfig, authority_df::DataFrame)
    _prepare_plots_backend()
    beta_values = sort(unique(authority_df.target_beta_high_kg_m2))
    layout_rows = 1
    layout_cols = length(beta_values)
    plot_list = Any[]
    for beta_high in beta_values
        subset = sort(filter(row -> row.target_beta_high_kg_m2 == beta_high, authority_df), [:h_jettison_km, :target_beta_ratio])
        h_jettison_values = sort(unique(subset.h_jettison_km))
        beta_ratio_values = sort(unique(subset.target_beta_ratio))
        z = fill(NaN, length(h_jettison_values), length(beta_ratio_values))
        for row in eachrow(subset)
            i = findfirst(==(row.h_jettison_km), h_jettison_values)
            j = findfirst(==(row.target_beta_ratio), beta_ratio_values)
            z[i, j] = row.downrange_authority_km
        end
        p = heatmap(
            beta_ratio_values,
            h_jettison_values,
            z;
            xlabel="β Ratio",
            ylabel="Jettison Altitude (km)",
            title=@sprintf("β High = %.0f kg/m²", beta_high),
            colorbar_title="km",
            xrotation=0,
            xticks=beta_ratio_values,
            yticks=h_jettison_values,
        )
        push!(plot_list, p)
    end
    combined = plot(plot_list...; layout=(layout_rows, layout_cols), size=(620 * layout_cols, 560))
    savefig(combined, joinpath(config.plot_dir, "downrange_authority_heatmap.png"))
end

function save_terminal_metric_plot(config::PrelimConfig, summary_df::DataFrame, authority_df::DataFrame)
    _prepare_plots_backend()
    rows = NamedTuple[]
    for beta_high in sort(unique(authority_df.target_beta_high_kg_m2))
        beta_subset = filter(row -> row.target_beta_high_kg_m2 == beta_high, authority_df)
        for beta_ratio in sort(unique(beta_subset.target_beta_ratio))
            ratio_subset = filter(row -> row.target_beta_ratio == beta_ratio, beta_subset)
            best_idx = argmax(ratio_subset.downrange_authority_km)
            best = ratio_subset[best_idx, :]
            summary_subset = filter(
                row -> row.policy == "bang_bang" &&
                    _approx_eq(row.target_beta_high_kg_m2, beta_high) &&
                    _approx_eq(row.target_beta_ratio, beta_ratio) &&
                    _approx_eq(row.h_jettison_km, best.h_jettison_km),
                summary_df,
            )
            min_idx = argmin(summary_subset.impact_downrange_km)
            max_idx = argmax(summary_subset.impact_downrange_km)
            for (endpoint, row_idx) in (("min_range", min_idx), ("max_range", max_idx))
                push!(rows, (
                    beta_high_kg_m2=beta_high,
                    beta_ratio=beta_ratio,
                    endpoint=endpoint,
                    impact_velocity_mps=summary_subset.impact_velocity_mps[row_idx],
                    peak_total_decel_earth_g=summary_subset.peak_total_decel_earth_g[row_idx],
                ))
            end
        end
    end
    metrics_df = DataFrame(rows)

    p_vel = plot(xlabel="β Ratio", ylabel="Impact Velocity (m/s)", title="Impact Velocity vs β Ratio", legend=:outerbottom)
    p_g = plot(xlabel="β Ratio", ylabel="Peak Aero Load (Earth-g)", title="Peak Aero Load vs β Ratio", legend=:outerbottom)
    for beta_high in sort(unique(metrics_df.beta_high_kg_m2))
        for endpoint in ("min_range", "max_range")
            subset = sort(
                filter(row -> row.beta_high_kg_m2 == beta_high && row.endpoint == endpoint, metrics_df),
                :beta_ratio,
            )
            label = @sprintf("β high %.0f, %s", beta_high, replace(endpoint, "_" => " "))
            plot!(p_vel, subset.beta_ratio, subset.impact_velocity_mps; marker=:circle, label=label)
            plot!(p_g, subset.beta_ratio, subset.peak_total_decel_earth_g; marker=:circle, label=label)
        end
    end
    combined = plot(p_vel, p_g; layout=(1, 2), size=(1450, 620))
    savefig(combined, joinpath(config.plot_dir, "terminal_velocity_and_peak_g_vs_beta_ratio.png"))
end

function save_local_effectiveness_plot(config::PrelimConfig, effectiveness_df::DataFrame)
    _prepare_plots_backend()
    p_line = plot(
        effectiveness_df.sensitivity_km_per_km,
        effectiveness_df.h_switch_km;
        xlabel="d(impact range) / d(h_switch) (km/km)",
        ylabel="Switch Altitude (km)",
        title="Local Control Effectiveness vs Switch Altitude",
        label="sensitivity",
    )
    p_scatter = scatter(
        effectiveness_df.velocity_mps ./ 1e3,
        effectiveness_df.h_switch_km;
        marker_z=effectiveness_df.sensitivity_abs_km_per_km,
        xlabel="Velocity (km/s)",
        ylabel="Altitude (km)",
        title="Where Control Authority Is Generated",
        colorbar_title="|sensitivity| (km/km)",
        label=false,
    )
    combined = plot(p_line, p_scatter; layout=(1, 2), size=(1450, 620))
    savefig(combined, joinpath(config.plot_dir, "local_control_effectiveness.png"))
end

function save_crossrange_proxy_plot(config::PrelimConfig, crossrange_df::DataFrame)
    isempty(crossrange_df) && return
    _prepare_plots_backend()
    p = bar(
        string.(crossrange_df.sigma),
        crossrange_df.impact_crossrange_km;
        xlabel="σ",
        ylabel="Impact Crossrange (km)",
        title="Crossrange Proxy vs Side-Force Fraction",
        legend=false,
    )
    savefig(p, joinpath(config.plot_dir, "crossrange_proxy_vs_sigma.png"))
end

function _ellipse_curve(center_x::Float64, center_y::Float64, major_km::Float64, minor_km::Float64, azimuth_deg::Float64; npts::Int=200)
    if !(major_km > 0.0) || !(minor_km >= 0.0)
        return Float64[], Float64[]
    end
    θ = range(0.0, 2π; length=npts)
    a = 0.5 * major_km
    b = 0.5 * minor_km
    ϕ = deg2rad(azimuth_deg)
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    x = similar(collect(θ))
    y = similar(collect(θ))
    for (i, t) in enumerate(θ)
        xr = a * cos(t)
        yr = b * sin(t)
        x[i] = center_x + xr * cosϕ - yr * sinϕ
        y[i] = center_y + xr * sinϕ + yr * cosϕ
    end
    return x, y
end

function save_landing_ellipse_plot(config::PrelimConfig, monte_carlo_samples_df::DataFrame, monte_carlo_summary_df::DataFrame)
    isempty(monte_carlo_samples_df) && return
    _prepare_plots_backend()
    p = plot(
        xlabel="Impact Downrange (km)",
        ylabel="Impact Crossrange (km)",
        title="Landing Footprint From 3D Monte Carlo",
        legend=:outerright,
    )
    for policy in ("body_only_monte_carlo", "guided_targeted_optimistic", "guided_targeted_crossrange_bilateral")
        samples = filter(row -> row.policy == policy, monte_carlo_samples_df)
        isempty(samples) && continue
        label = replace(policy, "_" => " ")
        scatter!(p, samples.impact_downrange_km, samples.impact_crossrange_km; label=label, markersize=3, alpha=0.7)
        summary = filter(row -> row.policy == policy, monte_carlo_summary_df)
        if nrow(summary) > 0
            ellipse = summary[1, :]
            if ellipse.ellipse_major_axis_km > 0.0
                x, y = _ellipse_curve(
                    mean(samples.impact_downrange_km),
                    mean(samples.impact_crossrange_km),
                    ellipse.ellipse_major_axis_km,
                    ellipse.ellipse_minor_axis_km,
                    ellipse.ellipse_azimuth_deg,
                )
                !isempty(x) && plot!(p, x, y; label="$(label) 95% ellipse", linestyle=:dash)
            end
        end
    end
    savefig(p, joinpath(config.plot_dir, "landing_ellipse_proxy.png"))
    savefig(p, joinpath(config.plot_dir, "landing_ellipse_3d.png"))

    p_centered = plot(
        xlabel="Relative Downrange (km)",
        ylabel="Crossrange (km)",
        title="Landing Footprint Centered On Mean Impact",
        legend=:outerright,
        aspect_ratio=:equal,
    )
    for policy in ("body_only_monte_carlo", "guided_targeted_optimistic", "guided_targeted_crossrange_bilateral")
        samples = filter(row -> row.policy == policy, monte_carlo_samples_df)
        isempty(samples) && continue
        label = replace(policy, "_" => " ")
        x_center = samples.impact_downrange_km .- mean(samples.impact_downrange_km)
        y_center = samples.impact_crossrange_km .- mean(samples.impact_crossrange_km)
        scatter!(p_centered, x_center, y_center; label=label, markersize=3, alpha=0.7)
        summary = filter(row -> row.policy == policy, monte_carlo_summary_df)
        if nrow(summary) > 0
            ellipse = summary[1, :]
            if ellipse.ellipse_major_axis_km > 0.0
                x, y = _ellipse_curve(
                    0.0,
                    0.0,
                    ellipse.ellipse_major_axis_km,
                    ellipse.ellipse_minor_axis_km,
                    ellipse.ellipse_azimuth_deg,
                )
                !isempty(x) && plot!(p_centered, x, y; label="$(label) 95% ellipse", linestyle=:dash)
            end
        end
    end
    savefig(p_centered, joinpath(config.plot_dir, "landing_ellipse_centered.png"))
end

function save_alpha_body_sensitivity_plot(config::PrelimConfig, alpha_df::DataFrame)
    isempty(alpha_df) && return
    _prepare_plots_backend()
    ordered = sort(alpha_df, :alpha_body_deg)
    p_authority = plot(
        ordered.alpha_body_deg,
        ordered.downrange_authority_km;
        xlabel="Body Angle of Attack (deg)",
        ylabel="Downrange Authority (km)",
        title="Representative Case: Authority vs Body AoA",
        marker=:circle,
        legend=false,
    )
    p_terminal = plot(
        ordered.alpha_body_deg,
        ordered.min_impact_velocity_mps;
        xlabel="Body Angle of Attack (deg)",
        ylabel="Impact Velocity (m/s)",
        title="Representative Case: Impact Velocity vs Body AoA",
        marker=:circle,
        label="min range",
        legend=:outerright,
    )
    plot!(p_terminal, ordered.alpha_body_deg, ordered.max_impact_velocity_mps; marker=:circle, label="max range")
    combined = plot(p_authority, p_terminal; layout=(1, 2), size=(1450, 620))
    savefig(combined, joinpath(config.plot_dir, "alpha_body_sensitivity.png"))
end

function save_crossrange_sensitivity_plot(config::PrelimConfig, crossrange_sensitivity_df::DataFrame)
    isempty(crossrange_sensitivity_df) && return
    _prepare_plots_backend()
    alpha_values = sort(unique(crossrange_sensitivity_df.alpha_body_deg))
    plots = Any[]
    for alpha_deg in alpha_values
        subset = filter(row -> _approx_eq(row.alpha_body_deg, alpha_deg), crossrange_sensitivity_df)
        p_cross = plot(
            xlabel="Panel Cant Angle (deg)",
            ylabel="|Impact Crossrange| (km)",
            title=@sprintf("α_body = %.0f deg: Crossrange", alpha_deg),
            legend=:outerright,
        )
        p_shift = plot(
            xlabel="Panel Cant Angle (deg)",
            ylabel="Downrange Shift vs Symmetric (km)",
            title=@sprintf("α_body = %.0f deg: Downrange Shift", alpha_deg),
            legend=:outerright,
        )
        for level in ["small", "medium", "large"]
            level_subset = sort(
                filter(row -> row.differential_deployment_label == level, subset),
                :panel_cant_deg,
            )
            isempty(level_subset) && continue
            plot!(
                p_cross,
                level_subset.panel_cant_deg,
                level_subset.abs_impact_crossrange_km;
                marker=:circle,
                label=level,
            )
            plot!(
                p_shift,
                level_subset.panel_cant_deg,
                level_subset.downrange_shift_vs_symmetric_km;
                marker=:circle,
                label=level,
            )
        end
        push!(plots, p_cross)
        push!(plots, p_shift)
    end
    combined = plot(plots...; layout=(length(alpha_values), 2), size=(1500, 520 * length(alpha_values)))
    savefig(combined, joinpath(config.plot_dir, "crossrange_sensitivity.png"))
end
