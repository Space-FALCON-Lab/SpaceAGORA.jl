@inline function _approx_eq(a::Real, b::Real; atol::Float64=1e-6)
    return isapprox(Float64(a), Float64(b); atol=atol, rtol=0.0)
end

function _summary_row(
    config::PrelimConfig,
    aero_case::CalibratedAeroCase,
    trajectory_result::TrajectoryResult,
    policy::ControlPolicy,
)
    s = trajectory_result.summary
    return (
        atmosphere_path=s.atmosphere_path,
        body_geometry=aero_case.body_label,
        deployed_geometry=aero_case.deployed_label,
        target_beta_high_kg_m2=aero_case.target_beta_high,
        target_beta_low_kg_m2=aero_case.target_beta_low,
        target_beta_ratio=aero_case.target_beta_ratio,
        achieved_beta_high_kg_m2=aero_case.achieved_beta_high,
        achieved_beta_low_kg_m2=aero_case.achieved_beta_low,
        achieved_beta_ratio=aero_case.achieved_beta_ratio,
        mass_kg=aero_case.mass_kg,
        panel_area_total_m2=aero_case.deployed_state.panel_area_total_m2,
        panel_area_each_m2=aero_case.deployed_state.panel_area_each_m2,
        reference_temperature_k=aero_case.reference_temperature_k,
        reference_mach=aero_case.reference_mach,
        entry_altitude_km=config.h0_m / 1e3,
        entry_velocity_km_s=config.V0_mps / 1e3,
        entry_flight_path_angle_deg=config.gamma0_deg,
        policy=policy.name,
        h_switch_km=isfinite(policy.h_switch_m) ? policy.h_switch_m / 1e3 : NaN,
        h_jettison_km=isfinite(policy.h_jettison_m) ? policy.h_jettison_m / 1e3 : NaN,
        impact_downrange_km=s.impact_downrange_km,
        flight_time_s=s.flight_time_s,
        impact_velocity_mps=s.impact_velocity_mps,
        impact_flight_path_angle_deg=s.impact_flight_path_angle_deg,
        peak_dynamic_pressure_pa=s.peak_dynamic_pressure_pa,
        peak_drag_accel_mps2=s.peak_drag_accel_mps2,
        peak_total_decel_earth_g=s.peak_total_decel_earth_g,
        max_heating_proxy=s.max_heating_proxy,
        impact_g_load_proxy_earth_g=s.impact_g_load_proxy_earth_g,
        integration_stopped_at_surface=s.integration_stopped_at_surface,
    )
end

function _named_trajectory_key(policy::ControlPolicy)
    switch_tag = isfinite(policy.h_switch_m) ? @sprintf("%03d", round(Int, policy.h_switch_m / 1e3)) : "nan"
    jettison_tag = isfinite(policy.h_jettison_m) ? @sprintf("%03d", round(Int, policy.h_jettison_m / 1e3)) : "nan"
    return "$(policy.kind)_sw$(switch_tag)_jet$(jettison_tag)"
end

function _should_save_representative(
    config::PrelimConfig,
    aero_case::CalibratedAeroCase,
    policy::ControlPolicy,
)
    if !_approx_eq(aero_case.target_beta_high, config.representative_beta_high) ||
       !_approx_eq(aero_case.target_beta_ratio, config.representative_beta_ratio)
        return false
    end
    if policy.kind === :body_only
        return true
    elseif policy.kind === :fixed_deployed
        return _approx_eq(policy.h_jettison_m, config.representative_h_jettison_m)
    elseif policy.kind === :bang_bang
        return _approx_eq(policy.h_jettison_m, config.representative_h_jettison_m) &&
            any(x -> _approx_eq(policy.h_switch_m, x), config.representative_switches_m)
    end
    return false
end

function compute_authority_table(summary_df::DataFrame)
    bang_df = filter(row -> row.policy == "bang_bang", summary_df)
    grouped = groupby(
        bang_df,
        [:target_beta_high_kg_m2, :target_beta_ratio, :target_beta_low_kg_m2, :achieved_beta_high_kg_m2, :achieved_beta_low_kg_m2, :achieved_beta_ratio, :h_jettison_km],
    )
    rows = NamedTuple[]
    for group in grouped
        ranges = group.impact_downrange_km
        max_idx = argmax(ranges)
        min_idx = argmin(ranges)
        push!(rows, (
            target_beta_high_kg_m2=group.target_beta_high_kg_m2[1],
            target_beta_low_kg_m2=group.target_beta_low_kg_m2[1],
            target_beta_ratio=group.target_beta_ratio[1],
            achieved_beta_high_kg_m2=group.achieved_beta_high_kg_m2[1],
            achieved_beta_low_kg_m2=group.achieved_beta_low_kg_m2[1],
            achieved_beta_ratio=group.achieved_beta_ratio[1],
            h_jettison_km=group.h_jettison_km[1],
            downrange_authority_km=maximum(ranges) - minimum(ranges),
            min_impact_downrange_km=minimum(ranges),
            max_impact_downrange_km=maximum(ranges),
            h_switch_min_km=group.h_switch_km[min_idx],
            h_switch_max_km=group.h_switch_km[max_idx],
        ))
    end
    return DataFrame(rows)
end

function compute_local_effectiveness(
    bang_summary_df::DataFrame,
    body_only_traj::DataFrame,
)
    sorted_df = sort(bang_summary_df, :h_switch_km)
    h_switch_km = Vector{Float64}(sorted_df.h_switch_km)
    range_km = Vector{Float64}(sorted_df.impact_downrange_km)
    n = length(h_switch_km)
    sensitivity = fill(0.0, n)
    for i in 1:n
        if i == 1
            sensitivity[i] = (range_km[i + 1] - range_km[i]) / (h_switch_km[i + 1] - h_switch_km[i])
        elseif i == n
            sensitivity[i] = (range_km[i] - range_km[i - 1]) / (h_switch_km[i] - h_switch_km[i - 1])
        else
            sensitivity[i] = (range_km[i + 1] - range_km[i - 1]) / (h_switch_km[i + 1] - h_switch_km[i - 1])
        end
    end

    body_sorted = sort(body_only_traj, :altitude_km, rev=true)
    velocity_at_switch = [
        _interp_descending(h, body_sorted.altitude_km, body_sorted.velocity_mps) for h in h_switch_km
    ]
    q_at_switch = [
        _interp_descending(h, body_sorted.altitude_km, body_sorted.q_pa) for h in h_switch_km
    ]

    return DataFrame(
        h_switch_km=h_switch_km,
        impact_downrange_km=range_km,
        sensitivity_km_per_km=sensitivity,
        sensitivity_abs_km_per_km=abs.(sensitivity),
        velocity_mps=velocity_at_switch,
        q_pa=q_at_switch,
    )
end

function run_deterministic_sweeps(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
)
    summary_rows = NamedTuple[]
    representative_trajectories = Dict{String, DataFrame}()
    aero_cases = Dict{Tuple{Float64, Float64}, CalibratedAeroCase}()

    total_cases = length(config.beta_high_targets) * length(config.beta_ratios)
    case_counter = 0
    for target_beta_high in config.beta_high_targets
        for target_beta_ratio in config.beta_ratios
            case_counter += 1
            @info "Calibrating aero case" case_counter total_cases target_beta_high target_beta_ratio
            aero_case = calibrate_aero_case(config, adapter, target_beta_high, target_beta_ratio)
            aero_cases[(target_beta_high, target_beta_ratio)] = aero_case

            body_policy = ControlPolicy("body_only", :body_only, NaN, NaN)
            body_result = solve_entry_trajectory(config, adapter, aero_case, body_policy; save_trajectory=_should_save_representative(config, aero_case, body_policy))
            push!(summary_rows, _summary_row(config, aero_case, body_result, body_policy))
            if !isempty(body_result.trajectory)
                representative_trajectories[_named_trajectory_key(body_policy)] = body_result.trajectory
            end

            for h_jettison_m in config.h_jettison_grid_m
                fixed_policy = ControlPolicy("fixed_deployed", :fixed_deployed, NaN, h_jettison_m)
                fixed_result = solve_entry_trajectory(config, adapter, aero_case, fixed_policy; save_trajectory=_should_save_representative(config, aero_case, fixed_policy))
                push!(summary_rows, _summary_row(config, aero_case, fixed_result, fixed_policy))
                if !isempty(fixed_result.trajectory)
                    representative_trajectories[_named_trajectory_key(fixed_policy)] = fixed_result.trajectory
                end

                for h_switch_m in config.h_switch_grid_m
                    bang_policy = ControlPolicy("bang_bang", :bang_bang, h_switch_m, h_jettison_m)
                    bang_result = solve_entry_trajectory(config, adapter, aero_case, bang_policy; save_trajectory=_should_save_representative(config, aero_case, bang_policy))
                    push!(summary_rows, _summary_row(config, aero_case, bang_result, bang_policy))
                    if !isempty(bang_result.trajectory)
                        representative_trajectories[_named_trajectory_key(bang_policy)] = bang_result.trajectory
                    end

                    if config.include_reverse_bangbang_sweep
                        reverse_policy = ControlPolicy("reverse_bang_bang", :reverse_bang_bang, h_switch_m, h_jettison_m)
                        reverse_result = solve_entry_trajectory(config, adapter, aero_case, reverse_policy; save_trajectory=false)
                        push!(summary_rows, _summary_row(config, aero_case, reverse_result, reverse_policy))
                    end
                end
            end
        end
    end

    summary_df = DataFrame(summary_rows)
    authority_df = compute_authority_table(summary_df)
    return summary_df, authority_df, representative_trajectories, aero_cases
end

