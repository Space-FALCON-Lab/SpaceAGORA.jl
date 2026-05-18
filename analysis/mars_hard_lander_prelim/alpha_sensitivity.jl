function run_body_alpha_sensitivity(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter;
    alpha_values_deg::AbstractVector{<:Real}=[0.0, 5.0, 10.0, 15.0, 20.0],
)
    rows = NamedTuple[]
    h_jettison_m = config.representative_h_jettison_m
    for alpha_deg in Float64.(collect(alpha_values_deg))
        cfg_alpha = with_geometry(config, with_body_alpha(config.geometry, alpha_deg))
        aero_case = calibrate_aero_case(
            cfg_alpha,
            adapter,
            cfg_alpha.representative_beta_high,
            cfg_alpha.representative_beta_ratio,
        )
        body_policy = ControlPolicy("alpha_body_sensitivity_body_only", :body_only, NaN, NaN)
        body_result = solve_entry_trajectory(cfg_alpha, adapter, aero_case, body_policy; save_trajectory=false)

        range_rows = NamedTuple[]
        for h_switch_m in cfg_alpha.h_switch_grid_m
            bang_policy = ControlPolicy("alpha_body_sensitivity", :bang_bang, h_switch_m, h_jettison_m)
            bang_result = solve_entry_trajectory(cfg_alpha, adapter, aero_case, bang_policy; save_trajectory=false)
            push!(range_rows, (
                h_switch_km=h_switch_m / 1e3,
                impact_downrange_km=bang_result.summary.impact_downrange_km,
                impact_velocity_mps=bang_result.summary.impact_velocity_mps,
                peak_total_decel_earth_g=bang_result.summary.peak_total_decel_earth_g,
            ))
        end
        range_df = DataFrame(range_rows)
        min_idx = argmin(range_df.impact_downrange_km)
        max_idx = argmax(range_df.impact_downrange_km)
        push!(rows, (
            alpha_body_deg=alpha_deg,
            target_beta_high_kg_m2=aero_case.target_beta_high,
            target_beta_ratio=aero_case.target_beta_ratio,
            h_jettison_km=h_jettison_m / 1e3,
            achieved_beta_high_kg_m2=aero_case.achieved_beta_high,
            achieved_beta_low_kg_m2=aero_case.achieved_beta_low,
            body_only_impact_downrange_km=body_result.summary.impact_downrange_km,
            downrange_authority_km=range_df.impact_downrange_km[max_idx] - range_df.impact_downrange_km[min_idx],
            min_range_h_switch_km=range_df.h_switch_km[min_idx],
            max_range_h_switch_km=range_df.h_switch_km[max_idx],
            min_impact_velocity_mps=range_df.impact_velocity_mps[min_idx],
            max_impact_velocity_mps=range_df.impact_velocity_mps[max_idx],
            min_peak_total_decel_earth_g=range_df.peak_total_decel_earth_g[min_idx],
            max_peak_total_decel_earth_g=range_df.peak_total_decel_earth_g[max_idx],
        ))
    end
    return DataFrame(rows)
end
