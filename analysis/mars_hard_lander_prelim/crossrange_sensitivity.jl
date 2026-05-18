function _crossrange_target_result(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
)
    h_jettison_m = config.representative_h_jettison_m
    impact_ranges_km = Float64[]
    for h_switch_m in config.h_switch_grid_m
        result = solve_entry_trajectory(
            config,
            adapter,
            aero_case,
            ControlPolicy("crossrange_target_bracket", :bang_bang, h_switch_m, h_jettison_m);
            save_trajectory=false,
        )
        push!(impact_ranges_km, result.summary.impact_downrange_km)
    end
    target_range_km = minimum(impact_ranges_km) +
        config.target_range_fraction * (maximum(impact_ranges_km) - minimum(impact_ranges_km))
    return find_switch_altitude_for_target_range(
        config,
        adapter,
        aero_case,
        h_jettison_m,
        target_range_km;
        save_trajectory=false,
    )
end

function run_crossrange_sensitivity(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter;
    alpha_body_values_deg::AbstractVector{<:Real}=[0.0, 5.0],
    panel_cant_angles_deg::AbstractVector{<:Real}=[0.0, 15.0, 30.0, 45.0],
    deployment_levels::AbstractVector{<:Tuple}=[
        ("small", 0.20),
        ("medium", 0.40),
        ("large", 0.60),
    ],
)
    rows = NamedTuple[]
    for alpha_body_deg in Float64.(collect(alpha_body_values_deg))
        cfg_alpha = with_geometry(config, with_body_alpha(config.geometry, alpha_body_deg))
        aero_case = calibrate_aero_case(
            cfg_alpha,
            adapter,
            cfg_alpha.representative_beta_high,
            cfg_alpha.representative_beta_ratio,
        )
        target_result = _crossrange_target_result(cfg_alpha, adapter, aero_case)
        policy = ControlPolicy(
            "crossrange_sensitivity",
            :bang_bang,
            target_result.h_switch_m,
            target_result.h_jettison_m,
        )
        symmetric_result = solve_entry_trajectory_3d(
            cfg_alpha,
            adapter,
            aero_case,
            policy;
            save_trajectory=false,
            use_winds=false,
            panel_command_override=neutral_panel_command(cfg_alpha),
        )
        for (deployment_label, deployment_fraction_raw) in deployment_levels
            deployment_fraction = Float64(deployment_fraction_raw)
            for cant_deg in Float64.(collect(panel_cant_angles_deg))
                command = make_differential_panel_command(
                    cfg_alpha;
                    favored_side_sign=1.0,
                    differential_fraction=deployment_fraction,
                    cant_deg=cant_deg,
                )
                result = solve_entry_trajectory_3d(
                    cfg_alpha,
                    adapter,
                    aero_case,
                    policy;
                    save_trajectory=false,
                    use_winds=false,
                    panel_command_override=command,
                )
                push!(rows, (
                    alpha_body_deg=alpha_body_deg,
                    target_beta_high_kg_m2=aero_case.target_beta_high,
                    target_beta_ratio=aero_case.target_beta_ratio,
                    h_jettison_km=target_result.h_jettison_m / 1e3,
                    h_switch_km=target_result.h_switch_m / 1e3,
                    target_range_km=target_result.target_range_km,
                    target_guidance_status=String(target_result.status),
                    panel_cant_deg=cant_deg,
                    differential_deployment_label=String(deployment_label),
                    differential_deployment_fraction=deployment_fraction,
                    favored_side="right",
                    left_area_scale=command.left_area_scale,
                    right_area_scale=command.right_area_scale,
                    left_alpha_deg=command.left_alpha_deg,
                    right_alpha_deg=command.right_alpha_deg,
                    left_cant_deg=command.left_cant_deg,
                    right_cant_deg=command.right_cant_deg,
                    symmetric_impact_downrange_km=symmetric_result.summary.impact_downrange_km,
                    symmetric_impact_crossrange_km=symmetric_result.summary.impact_crossrange_km,
                    impact_downrange_km=result.summary.impact_downrange_km,
                    impact_crossrange_km=result.summary.impact_crossrange_km,
                    abs_impact_crossrange_km=abs(result.summary.impact_crossrange_km),
                    downrange_shift_vs_symmetric_km=result.summary.impact_downrange_km - symmetric_result.summary.impact_downrange_km,
                    peak_side_accel_mps2=result.summary.peak_side_accel_mps2,
                    peak_side_accel_earth_g=result.summary.peak_side_accel_mps2 / _earth_g(),
                    impact_velocity_mps=result.summary.impact_velocity_mps,
                ))
            end
        end
    end
    return DataFrame(rows)
end

function write_crossrange_sensitivity_note(config::PrelimConfig, crossrange_sensitivity_df::DataFrame)
    isempty(crossrange_sensitivity_df) && return nothing
    best_idx = argmax(crossrange_sensitivity_df.abs_impact_crossrange_km)
    best = crossrange_sensitivity_df[best_idx, :]
    open(joinpath(config.output_root, "crossrange_sensitivity.md"), "w") do io
        println(io, "# Crossrange Sensitivity Note")
        println(io)
        println(io, "## What Differential Panel Deployment Means")
        println(io, "- Symmetric deployment means the left and right panels match, so their lateral force components cancel and the configuration behaves mostly like a drag-modulation system.")
        println(io, "- Differential deployment means the left and right panels intentionally do not match. One side can stay more deployed, more canted, or more deflected than the other.")
        println(io, "- In this study, the favored right panel stays fully deployed while the left panel is partially retracted by the differential-deployment fraction. The favored panel is also canted outward and mildly deflected away from the broadside drag state.")
        println(io, "- That asymmetry makes the left and right aerodynamic forces stop canceling, which produces a net side force and therefore crossrange motion.")
        println(io)
        println(io, "## Study Setup")
        println(io, "- Representative vehicle class only: β_high = $(config.representative_beta_high) kg/m², β_ratio = $(config.representative_beta_ratio), h_jettison = $(config.representative_h_jettison_m / 1e3) km.")
        println(io, "- Deterministic 3D point-mass propagation with zero winds to isolate geometric lateral authority.")
        println(io, "- Body angle-of-attack values swept: 0 deg and 5 deg.")
        println(io, "- Panel cant angles swept: 0, 15, 30, and 45 deg.")
        println(io, "- Differential-deployment fractions swept: small = 0.20, medium = 0.40, large = 0.60.")
        println(io)
        println(io, "## Best Deterministic Crossrange Case")
        println(io, @sprintf(
            "- Best case in this sweep: |crossrange| = %.2f km at α_body = %.0f deg, cant = %.0f deg, differential fraction = %.2f (%s), with downrange shift %.2f km and peak side acceleration %.3f m/s².",
            best.abs_impact_crossrange_km,
            best.alpha_body_deg,
            best.panel_cant_deg,
            best.differential_deployment_fraction,
            best.differential_deployment_label,
            best.downrange_shift_vs_symmetric_km,
            best.peak_side_accel_mps2,
        ))
        println(io)
        println(io, "## Interpretation")
        println(io, "- If cant = 0 deg, the configuration stays almost purely in-track because there is no meaningful lateral force component.")
        println(io, "- Larger cant and larger deployment mismatch both increase side force, but they also perturb the downrange solution.")
        println(io, "- This is still an open-loop, geometry-based sensitivity study. It is useful for showing when crossrange begins to matter, not for claiming a closed-loop 3D guidance capability yet.")
    end
    return nothing
end
