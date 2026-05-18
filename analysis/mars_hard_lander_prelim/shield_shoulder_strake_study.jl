struct ShoulderStrakeComponent
    name::String
    azimuth_deg::Float64
    area_m2::Float64
    position_m::SVector{3, Float64}
    radial_dir::SVector{3, Float64}
end

struct ShoulderStrakeAssumptions
    cg_axial_fraction_of_body_length::Float64
    static_margin_fraction_of_diameter::Float64
    control_deflection_deg::Float64
    trim_check_angle_deg::Float64
    cl_alpha_per_rad::Float64
    cd0::Float64
    induced_drag_factor::Float64
    cl_max::Float64
end

struct ShoulderStrakeTrajectoryCommand
    common_scale::Float64
    pitch_scale::Float64
    yaw_scale::Float64
end

function build_shield_shoulder_strakes(
    config::PrelimConfig;
    strake_area_fraction_of_stowed_ref::Real,
    strake_count::Integer=4,
)
    strake_count in (4, 8) || throw(ArgumentError("This first-cut SHIELD shoulder-strake study is defined for 4 or 8 surfaces."))
    total_area_m2 = Float64(strake_area_fraction_of_stowed_ref) * shield_stowed_reference_area_m2(config)
    each_area_m2 = total_area_m2 / strake_count
    x_arm_m = _stowed_body_flap_axial_arm_m(config)
    radius_m = config.geometry.base_radius_m
    layout = strake_count == 4 ? [
        ("top", 0.0),
        ("right", 90.0),
        ("bottom", 180.0),
        ("left", 270.0),
    ] : [
        ("top", 0.0),
        ("top_right", 45.0),
        ("right", 90.0),
        ("bottom_right", 135.0),
        ("bottom", 180.0),
        ("bottom_left", 225.0),
        ("left", 270.0),
        ("top_left", 315.0),
    ]
    strakes = ShoulderStrakeComponent[]
    for (name, azimuth_deg) in layout
        ϕ = deg2rad(azimuth_deg)
        radial_dir = SVector{3, Float64}(0.0, sin(ϕ), cos(ϕ))
        position_m = SVector{3, Float64}(x_arm_m, radius_m * sin(ϕ), radius_m * cos(ϕ))
        push!(strakes, ShoulderStrakeComponent(name, azimuth_deg, each_area_m2, position_m, radial_dir))
    end
    return (
        strake_count=strake_count,
        strake_area_total_m2=total_area_m2,
        strake_area_each_m2=each_area_m2,
        strakes=strakes,
        x_arm_m=x_arm_m,
    )
end

@inline function _shoulder_strake_cl_cd(
    deflection_deg::Real,
    assumptions::ShoulderStrakeAssumptions,
)
    α_rad = deg2rad(Float64(deflection_deg))
    CL = clamp(assumptions.cl_alpha_per_rad * α_rad, -assumptions.cl_max, assumptions.cl_max)
    CD = assumptions.cd0 + assumptions.induced_drag_factor * CL^2
    return CL, CD
end

function _shoulder_strake_force_moment_delta(
    strake::ShoulderStrakeComponent,
    deflection_deg::Real,
    assumptions::ShoulderStrakeAssumptions,
    cg_x_m::Real,
)
    CL, CD = _shoulder_strake_cl_cd(deflection_deg, assumptions)
    area = strake.area_m2
    Δforce_area = SVector{3, Float64}(CD * area, CL * area * strake.radial_dir[2], CL * area * strake.radial_dir[3])
    r_cg_m = strake.position_m - SVector{3, Float64}(Float64(cg_x_m), 0.0, 0.0)
    Δmoment_area = cross(r_cg_m, Δforce_area)
    return Δforce_area, Δmoment_area
end

function _shoulder_strake_command_deflections(layout, mode::String, deflection_deg::Real)
    δ = Float64(deflection_deg)
    if mode == "neutral"
        return Dict(strake.name => 0.0 for strake in layout.strakes)
    elseif mode == "common_drag"
        return Dict(strake.name => δ for strake in layout.strakes)
    elseif mode == "yaw_right"
        return Dict(strake.name => δ * sind(strake.azimuth_deg) for strake in layout.strakes)
    elseif mode == "pitch_up"
        return Dict(strake.name => δ * cosd(strake.azimuth_deg) for strake in layout.strakes)
    elseif startswith(mode, "single_")
        surf = replace(mode, "single_" => "")
        out = Dict(strake.name => 0.0 for strake in layout.strakes)
        haskey(out, surf) || throw(ArgumentError("Unsupported shoulder-strake name $(repr(surf)) for current layout."))
        out[surf] = δ
        return out
    end
    throw(ArgumentError("Unsupported shoulder-strake command mode $(repr(mode))."))
end

function _shoulder_strake_command_deflections(layout, command::ShoulderStrakeTrajectoryCommand, max_deflection_deg::Real)
    δ_common = clamp(command.common_scale, -1.0, 1.0) * Float64(max_deflection_deg)
    δ_pitch = clamp(command.pitch_scale, -1.0, 1.0) * Float64(max_deflection_deg)
    δ_yaw = clamp(command.yaw_scale, -1.0, 1.0) * Float64(max_deflection_deg)
    return Dict(
        strake.name => clamp(
            δ_common +
            δ_pitch * cosd(strake.azimuth_deg) +
            δ_yaw * sind(strake.azimuth_deg),
            -Float64(max_deflection_deg),
            Float64(max_deflection_deg),
        )
        for strake in layout.strakes
    )
end

function evaluate_shield_shoulder_strake_command(
    config::PrelimConfig,
    baseline::BodyFlapBaselineCondition,
    assumptions::ShoulderStrakeAssumptions;
    strake_area_fraction_of_stowed_ref::Real,
    mode::String,
    strake_count::Integer=4,
)
    layout = build_shield_shoulder_strakes(config; strake_area_fraction_of_stowed_ref=strake_area_fraction_of_stowed_ref, strake_count=strake_count)
    x_cg_m = assumptions.cg_axial_fraction_of_body_length * _stowed_body_flap_axial_arm_m(config)
    deflections = _shoulder_strake_command_deflections(layout, mode, assumptions.control_deflection_deg)

    Δforce_area = SVector{3, Float64}(0.0, 0.0, 0.0)
    Δmoment_area = SVector{3, Float64}(0.0, 0.0, 0.0)
    for strake in layout.strakes
        flap_Δforce_area, flap_Δmoment_area = _shoulder_strake_force_moment_delta(strake, deflections[strake.name], assumptions, x_cg_m)
        Δforce_area += flap_Δforce_area
        Δmoment_area += flap_Δmoment_area
    end
    outputs = _body_flap_control_outputs(Δforce_area, Δmoment_area)
    q_ref = baseline.representative_q_pa
    return (
        mode=mode,
        strake_count=layout.strake_count,
        strake_area_fraction_of_stowed_ref=Float64(strake_area_fraction_of_stowed_ref),
        strake_area_total_m2=layout.strake_area_total_m2,
        strake_area_each_m2=layout.strake_area_each_m2,
        cg_axial_fraction_of_body_length=assumptions.cg_axial_fraction_of_body_length,
        static_margin_fraction_of_diameter=assumptions.static_margin_fraction_of_diameter,
        control_deflection_deg=assumptions.control_deflection_deg,
        cl_alpha_per_rad=assumptions.cl_alpha_per_rad,
        cd0=assumptions.cd0,
        induced_drag_factor=assumptions.induced_drag_factor,
        cl_max=assumptions.cl_max,
        representative_altitude_km=baseline.representative_altitude_km,
        representative_mach=baseline.representative_mach,
        representative_q_pa=q_ref,
        representative_velocity_mps=baseline.representative_velocity_mps,
        outputs...,
        delta_drag_force_N=q_ref * outputs.delta_drag_area_m2,
        delta_side_force_N=q_ref * outputs.delta_side_area_m2,
        delta_vertical_force_N=q_ref * outputs.delta_vertical_area_m2,
        delta_roll_moment_Nm=q_ref * outputs.delta_roll_moment_area_m3,
        delta_pitch_up_moment_Nm=q_ref * outputs.delta_pitch_up_moment_area_m3,
        delta_yaw_right_moment_Nm=q_ref * outputs.delta_yaw_right_moment_area_m3,
    )
end

function _shoulder_strake_derivative_rows(
    config::PrelimConfig,
    baseline::BodyFlapBaselineCondition,
    assumptions::ShoulderStrakeAssumptions,
    strake_area_fraction::Float64,
    strake_count::Integer=4,
)
    rows = NamedTuple[]
    x_cg_m = assumptions.cg_axial_fraction_of_body_length * _stowed_body_flap_axial_arm_m(config)
    for mode in ("single_top", "single_right", "common_drag", "yaw_right", "pitch_up")
        δ = 1.0
        layout = build_shield_shoulder_strakes(config; strake_area_fraction_of_stowed_ref=strake_area_fraction, strake_count=strake_count)
        pos = _shoulder_strake_command_deflections(layout, mode, δ)
        neg = _shoulder_strake_command_deflections(layout, mode, -δ)
        Δforce_pos = SVector{3, Float64}(0.0, 0.0, 0.0)
        Δmoment_pos = SVector{3, Float64}(0.0, 0.0, 0.0)
        Δforce_neg = SVector{3, Float64}(0.0, 0.0, 0.0)
        Δmoment_neg = SVector{3, Float64}(0.0, 0.0, 0.0)
        for strake in layout.strakes
            fp, mp = _shoulder_strake_force_moment_delta(strake, pos[strake.name], assumptions, x_cg_m)
            fn, mn = _shoulder_strake_force_moment_delta(strake, neg[strake.name], assumptions, x_cg_m)
            Δforce_pos += fp
            Δmoment_pos += mp
            Δforce_neg += fn
            Δmoment_neg += mn
        end
        out_pos = _body_flap_control_outputs(Δforce_pos, Δmoment_pos)
        out_neg = _body_flap_control_outputs(Δforce_neg, Δmoment_neg)
        q_ref = baseline.representative_q_pa
        push!(rows, (
            mode=mode,
            strake_count=strake_count,
            strake_area_fraction_of_stowed_ref=strake_area_fraction,
            d_drag_area_m2_per_deg=0.5 * (out_pos.delta_drag_area_m2 - out_neg.delta_drag_area_m2),
            d_side_area_m2_per_deg=0.5 * (out_pos.delta_side_area_m2 - out_neg.delta_side_area_m2),
            d_vertical_area_m2_per_deg=0.5 * (out_pos.delta_vertical_area_m2 - out_neg.delta_vertical_area_m2),
            d_roll_moment_area_m3_per_deg=0.5 * (out_pos.delta_roll_moment_area_m3 - out_neg.delta_roll_moment_area_m3),
            d_pitch_up_moment_area_m3_per_deg=0.5 * (out_pos.delta_pitch_up_moment_area_m3 - out_neg.delta_pitch_up_moment_area_m3),
            d_yaw_right_moment_area_m3_per_deg=0.5 * (out_pos.delta_yaw_right_moment_area_m3 - out_neg.delta_yaw_right_moment_area_m3),
            d_side_force_N_per_deg=0.5 * q_ref * (out_pos.delta_side_area_m2 - out_neg.delta_side_area_m2),
            d_vertical_force_N_per_deg=0.5 * q_ref * (out_pos.delta_vertical_area_m2 - out_neg.delta_vertical_area_m2),
            d_pitch_up_moment_Nm_per_deg=0.5 * q_ref * (out_pos.delta_pitch_up_moment_area_m3 - out_neg.delta_pitch_up_moment_area_m3),
            d_yaw_right_moment_Nm_per_deg=0.5 * q_ref * (out_pos.delta_yaw_right_moment_area_m3 - out_neg.delta_yaw_right_moment_area_m3),
        ))
    end
    return DataFrame(rows)
end

function run_shield_shoulder_strake_stability_study(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter;
    entry_mass_kg::Real=120.0,
    strake_area_fractions_of_stowed_ref::AbstractVector{<:Real}=[0.01, 0.02, 0.05],
    strake_counts::AbstractVector{<:Integer}=[4],
    cg_axial_fractions_of_body_length::AbstractVector{<:Real}=[0.50, 0.60, 0.70],
    static_margin_fractions_of_diameter::AbstractVector{<:Real}=[0.05, 0.10, 0.15],
    trim_check_angle_deg::Real=5.0,
    control_deflection_deg::Real=15.0,
    cl_alpha_per_rad::Real=3.5,
    cd0::Real=0.05,
    induced_drag_factor::Real=0.15,
    cl_max::Real=1.2,
)
    mkpath(config.output_root)
    baseline_bundle = _representative_shield_body_flap_baseline_condition(config, adapter; entry_mass_kg=entry_mass_kg)
    baseline = baseline_bundle.condition
    targets = baseline_bundle.targets
    aero_case = baseline_bundle.aero_case
    surface_density_proxy_kg_m2 = shield_body_flap_areal_density_proxy_kg_m2(config, entry_mass_kg)
    rows = NamedTuple[]
    deriv_rows = NamedTuple[]
    for strake_count in Int.(collect(strake_counts))
        for cg_frac in Float64.(collect(cg_axial_fractions_of_body_length))
            for sm_frac in Float64.(collect(static_margin_fractions_of_diameter))
                assumptions = ShoulderStrakeAssumptions(
                    cg_frac,
                    sm_frac,
                    Float64(control_deflection_deg),
                    Float64(trim_check_angle_deg),
                    Float64(cl_alpha_per_rad),
                    Float64(cd0),
                    Float64(induced_drag_factor),
                    Float64(cl_max),
                )
                baseline_row = _baseline_stowed_body_stability_row(config, baseline, aero_case, targets, BodyFlapStabilityAssumptions(cg_frac, sm_frac, trim_check_angle_deg, control_deflection_deg))
                for strake_fraction in Float64.(collect(strake_area_fractions_of_stowed_ref))
                    surface_mass_proxy_kg = strake_fraction * shield_stowed_reference_area_m2(config) * surface_density_proxy_kg_m2
                    yaw_ctrl = evaluate_shield_shoulder_strake_command(config, baseline, assumptions; strake_area_fraction_of_stowed_ref=strake_fraction, mode="yaw_right", strake_count=strake_count)
                    pitch_ctrl = evaluate_shield_shoulder_strake_command(config, baseline, assumptions; strake_area_fraction_of_stowed_ref=strake_fraction, mode="pitch_up", strake_count=strake_count)
                    common_ctrl = evaluate_shield_shoulder_strake_command(config, baseline, assumptions; strake_area_fraction_of_stowed_ref=strake_fraction, mode="common_drag", strake_count=strake_count)
                    deriv_df = _shoulder_strake_derivative_rows(config, baseline, assumptions, strake_fraction, strake_count)
                    yaw_deriv = filter(row -> row.mode == "yaw_right", deriv_df)[1, :]
                    pitch_deriv = filter(row -> row.mode == "pitch_up", deriv_df)[1, :]
                    α_trim_cap_deg = abs(baseline_row.d_pitch_up_moment_Nm_per_deg) > 1e-12 ?
                        abs(pitch_ctrl.delta_pitch_up_moment_Nm / baseline_row.d_pitch_up_moment_Nm_per_deg) : Inf
                    β_trim_cap_deg = abs(baseline_row.d_yaw_right_moment_Nm_per_deg) > 1e-12 ?
                        abs(yaw_ctrl.delta_yaw_right_moment_Nm / baseline_row.d_yaw_right_moment_Nm_per_deg) : Inf
                    push!(rows, (
                        baseline_row...,
                        strake_count=strake_count,
                        strake_area_total_m2=pitch_ctrl.strake_area_total_m2,
                        strake_area_each_m2=pitch_ctrl.strake_area_each_m2,
                        strake_area_fraction_of_stowed_ref=strake_fraction,
                        strake_mass_proxy_kg=surface_mass_proxy_kg,
                        cl_alpha_per_rad=assumptions.cl_alpha_per_rad,
                        cd0=assumptions.cd0,
                        induced_drag_factor=assumptions.induced_drag_factor,
                        cl_max=assumptions.cl_max,
                        dCY_dδ_per_deg=yaw_deriv.d_side_area_m2_per_deg / baseline_row.stowed_reference_area_m2,
                        dCn_dδ_per_deg=yaw_deriv.d_yaw_right_moment_area_m3_per_deg / (baseline_row.stowed_reference_area_m2 * baseline_row.reference_length_m),
                        dCZ_dδ_per_deg=pitch_deriv.d_vertical_area_m2_per_deg / baseline_row.stowed_reference_area_m2,
                        dCm_dδ_per_deg=pitch_deriv.d_pitch_up_moment_area_m3_per_deg / (baseline_row.stowed_reference_area_m2 * baseline_row.reference_length_m),
                        d_side_force_N_per_deg_ctrl=yaw_deriv.d_side_force_N_per_deg,
                        d_yaw_right_moment_Nm_per_deg_ctrl=yaw_deriv.d_yaw_right_moment_Nm_per_deg,
                        d_vertical_force_N_per_deg_ctrl=pitch_deriv.d_vertical_force_N_per_deg,
                        d_pitch_up_moment_Nm_per_deg_ctrl=pitch_deriv.d_pitch_up_moment_Nm_per_deg,
                        common_drag_force_N=common_ctrl.delta_drag_force_N,
                        yaw_control_drag_coupling_N=yaw_ctrl.delta_drag_force_N,
                        pitch_control_drag_coupling_N=pitch_ctrl.delta_drag_force_N,
                        yaw_control_side_force_N=yaw_ctrl.delta_side_force_N,
                        yaw_control_moment_Nm=yaw_ctrl.delta_yaw_right_moment_Nm,
                        pitch_control_vertical_force_N=pitch_ctrl.delta_vertical_force_N,
                        pitch_control_moment_Nm=pitch_ctrl.delta_pitch_up_moment_Nm,
                        pitch_trim_cap_deg=α_trim_cap_deg,
                        yaw_trim_cap_deg=β_trim_cap_deg,
                        trim_plausible = baseline_row.stable_pitch && baseline_row.stable_yaw &&
                            α_trim_cap_deg >= assumptions.trim_check_angle_deg &&
                            β_trim_cap_deg >= assumptions.trim_check_angle_deg,
                    ))
                    for row in eachrow(deriv_df)
                        push!(deriv_rows, (
                            strake_count=strake_count,
                            cg_axial_fraction_of_body_length=cg_frac,
                            static_margin_fraction_of_diameter=sm_frac,
                            strake_area_fraction_of_stowed_ref=strake_fraction,
                            control_deflection_deg=assumptions.control_deflection_deg,
                            cl_alpha_per_rad=assumptions.cl_alpha_per_rad,
                            cd0=assumptions.cd0,
                            induced_drag_factor=assumptions.induced_drag_factor,
                            cl_max=assumptions.cl_max,
                            row...,
                        ))
                    end
                end
            end
        end
    end
    screen_df = DataFrame(rows)
    deriv_df = DataFrame(deriv_rows)
    csv_path = joinpath(config.output_root, "shield_shoulder_strake_stability_screen.csv")
    deriv_path = joinpath(config.output_root, "shield_shoulder_strake_derivatives.csv")
    note_path = joinpath(config.output_root, "shield_shoulder_strake_note.md")
    CSV.write(csv_path, screen_df)
    CSV.write(deriv_path, deriv_df)
    best = nrow(screen_df) > 0 ? screen_df[argmax(screen_df.pitch_trim_cap_deg .+ screen_df.yaw_trim_cap_deg), :] : nothing
    open(note_path, "w") do io
        println(io, "# SHIELD Shoulder-Petal / Strake Stability Screening")
        println(io)
        println(io, "This study tests a stronger pre-skirt guidance architecture than the earlier tiny body flaps: shoulder-mounted deployable petals / strakes with a higher-lift surrogate force model.")
        println(io)
        println(io, "## Frozen Baseline")
        println(io, @sprintf("- Same published-size SHIELD surrogate, entry mass %.1f kg", Float64(entry_mass_kg)))
        println(io, @sprintf("- Representative pre-deployment condition: altitude %.2f km, Mach %.2f, q = %.2f Pa", baseline.representative_altitude_km, baseline.representative_mach, baseline.representative_q_pa))
        println(io, "- Skirt deployment remains subsonic; this is strictly a pre-skirt control screen.")
        println(io)
        println(io, "## Shoulder-Strake Surrogate Assumptions")
        println(io, @sprintf("- Strake counts screened: %s", join(string.(sort(unique(screen_df.strake_count))), ", ")))
        println(io, @sprintf("- Lift slope: %.2f / rad", Float64(cl_alpha_per_rad)))
        println(io, @sprintf("- Zero-lift drag: %.3f", Float64(cd0)))
        println(io, @sprintf("- Induced-drag factor: %.3f", Float64(induced_drag_factor)))
        println(io, @sprintf("- CL cap: %.2f", Float64(cl_max)))
        println(io, "- Control surfaces are screened as stronger lifting strakes, not as the same flat-plate tabs used in the earlier body-flap study.")
        println(io)
        if best !== nothing
            println(io, "## Best Screened Case")
            println(io, @sprintf("- CG fraction: %.2f body lengths", best.cg_axial_fraction_of_body_length))
            println(io, @sprintf("- Strake count: %d", Int(best.strake_count)))
            println(io, @sprintf("- Static margin: %.0f%% of diameter", 100.0 * best.static_margin_fraction_of_diameter))
            println(io, @sprintf("- Total strake area: %.1f%% of stowed reference area", 100.0 * best.strake_area_fraction_of_stowed_ref))
            println(io, @sprintf("- Mass proxy: %.2f kg", best.strake_mass_proxy_kg))
            println(io, @sprintf("- Pitch trim cap: %.2f deg", best.pitch_trim_cap_deg))
            println(io, @sprintf("- Yaw trim cap: %.2f deg", best.yaw_trim_cap_deg))
            println(io, @sprintf("- Pitch control moment: %.2f N·m", best.pitch_control_moment_Nm))
            println(io, @sprintf("- Yaw control moment: %.2f N·m", best.yaw_control_moment_Nm))
            println(io, @sprintf("- Trim-plausible at %.1f deg requirement: %s", Float64(trim_check_angle_deg), best.trim_plausible ? "yes" : "no"))
        end
        println(io)
        println(io, "## Interpretation")
        if any(screen_df.trim_plausible)
            println(io, "- At least one screened shoulder-strake combination clears the trim gate. That makes it the first SHIELD-follow-on control architecture worth considering for a guided trajectory rerun.")
        else
            println(io, "- None of the screened shoulder-strake combinations clears the trim gate. That means even the stronger pre-skirt architecture does not yet justify a guided trajectory rerun.")
        end
        println(io, "- This is still a surrogate screen, not a thermal/structural closure result.")
    end
    return (
        baseline_condition=baseline,
        targets=targets,
        aero_case=aero_case,
        screen_df=screen_df,
        derivatives_df=deriv_df,
        csv_path=csv_path,
        derivatives_path=deriv_path,
        note_path=note_path,
    )
end

function _shoulder_strake_force_areas(
    config::PrelimConfig,
    assumptions::ShoulderStrakeAssumptions,
    strake_area_fraction_of_stowed_ref::Float64,
    command::ShoulderStrakeTrajectoryCommand,
    strake_count::Integer=4,
)
    layout = build_shield_shoulder_strakes(config; strake_area_fraction_of_stowed_ref=strake_area_fraction_of_stowed_ref, strake_count=strake_count)
    deflections = _shoulder_strake_command_deflections(layout, command, assumptions.control_deflection_deg)
    Δforce_area = SVector{3, Float64}(0.0, 0.0, 0.0)
    for strake in layout.strakes
        flap_Δforce_area, _ = _shoulder_strake_force_moment_delta(strake, deflections[strake.name], assumptions, 0.0)
        Δforce_area += flap_Δforce_area
    end
    return _body_flap_control_outputs(Δforce_area, SVector{3, Float64}(0.0, 0.0, 0.0))
end

function _sample_shoulder_strake_state(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
    policy::ControlPolicy,
    state::SVector{6, Float64},
    time_s::Float64,
    density_scale::Float64,
    use_winds::Bool,
    assumptions::ShoulderStrakeAssumptions,
    strake_area_fraction_of_stowed_ref::Float64,
    shoulder_command::ShoulderStrakeTrajectoryCommand,
    strake_count::Integer=4,
)
    r_vec = SVector{3, Float64}(state[1], state[2], state[3])
    v_vec = SVector{3, Float64}(state[4], state[5], state[6])
    r_norm = norm(r_vec)
    h_m = max(r_norm - config.planet.Rp_m, 0.0)
    east, north, up = _enu_basis(r_vec)
    lat, lon = _lat_lon(r_vec)

    v_e = dot(v_vec, east)
    v_n = dot(v_vec, north)
    v_u = dot(v_vec, up)
    rho_kg_m3, temperature_k, wind_enu = atmosphere_state(
        adapter,
        h_m,
        lat,
        lon,
        time_s;
        include_wind=use_winds,
    )
    rho_kg_m3 *= density_scale

    v_rel_enu = SVector{3, Float64}(v_e, v_n, v_u) - wind_enu
    v_rel_mag = max(norm(v_rel_enu), 1.0)
    mach = mach_number(adapter, v_rel_mag, temperature_k)
    deployed = deployed_active(policy, h_m)
    loads = aerodynamic_loads_3d(
        config,
        aero_case,
        deployed,
        mach,
        config.planet.γ,
        0.0,
        nothing,
        false,
    )

    if !deployed
        strake = _shoulder_strake_force_areas(config, assumptions, strake_area_fraction_of_stowed_ref, shoulder_command, strake_count)
        loads = (
            state_label = "stowed body + shoulder strakes",
            CDA_m2 = loads.CDA_m2 + strake.delta_drag_area_m2,
            CLA_m2 = loads.CLA_m2 + strake.delta_vertical_area_m2,
            CSA_m2 = loads.CSA_m2 + strake.delta_side_area_m2,
            beta_eff_kg_m2 = aero_case.mass_kg / max(loads.CDA_m2 + strake.delta_drag_area_m2, eps()),
        )
    end

    q_pa = 0.5 * rho_kg_m3 * v_rel_mag^2
    xw, yw, zw = _wind_frame(v_rel_enu)
    aero_enu = (q_pa / aero_case.mass_kg) * (
        loads.CDA_m2 * xw +
        loads.CLA_m2 * zw +
        loads.CSA_m2 * yw
    )
    aero_cart = aero_enu[1] * east + aero_enu[2] * north + aero_enu[3] * up
    gravity_cart = -(config.planet.μ / r_norm^3) * r_vec
    acc_cart = gravity_cart + aero_cart

    drag_accel_mps2 = q_pa * loads.CDA_m2 / aero_case.mass_kg
    lift_accel_mps2 = q_pa * abs(loads.CLA_m2) / aero_case.mass_kg
    side_accel_mps2 = q_pa * abs(loads.CSA_m2) / aero_case.mass_kg
    total_aero_g = sqrt(drag_accel_mps2^2 + lift_accel_mps2^2 + side_accel_mps2^2) / _earth_g()
    heating_proxy = sqrt(max(rho_kg_m3, 0.0)) * v_rel_mag^3
    horizontal_speed_mps = sqrt(v_e^2 + v_n^2)

    return (
        state=state,
        derivative=SVector{6, Float64}(v_vec[1], v_vec[2], v_vec[3], acc_cart[1], acc_cart[2], acc_cart[3]),
        altitude_km=h_m / 1e3,
        downrange_km=config.planet.Rp_m * lon / 1e3,
        crossrange_km=config.planet.Rp_m * lat / 1e3,
        heading_deg=rad2deg(atan(v_e, v_n)),
        velocity_mps=norm(v_vec),
        flight_path_angle_deg=rad2deg(atan(v_u, max(horizontal_speed_mps, 1.0))),
        rho_kg_m3=rho_kg_m3,
        temperature_k=temperature_k,
        q_pa=q_pa,
        drag_accel_mps2=drag_accel_mps2,
        lift_accel_mps2=lift_accel_mps2,
        side_accel_mps2=side_accel_mps2,
        total_aero_g=total_aero_g,
        heating_proxy=heating_proxy,
        mach=mach,
        beta_eff_kg_m2=loads.beta_eff_kg_m2,
        wind_east_mps=wind_enu[1],
        wind_north_mps=wind_enu[2],
        wind_up_mps=wind_enu[3],
        lateral_command=shoulder_command.yaw_scale,
        deployed=deployed,
        state_label=String(loads.state_label),
        kn=NaN,
    )
end

function _rk4_step_shoulder_strake(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
    policy::ControlPolicy,
    state::SVector{6, Float64},
    sample,
    time_s::Float64,
    dt_s::Float64,
    density_scale::Float64,
    use_winds::Bool,
    assumptions::ShoulderStrakeAssumptions,
    strake_area_fraction_of_stowed_ref::Float64,
    shoulder_command::ShoulderStrakeTrajectoryCommand,
    strake_count::Integer=4,
)
    k1 = sample.derivative
    k2 = _sample_shoulder_strake_state(config, adapter, aero_case, policy, state + 0.5 * dt_s * k1, time_s + 0.5 * dt_s, density_scale, use_winds, assumptions, strake_area_fraction_of_stowed_ref, shoulder_command, strake_count).derivative
    k3 = _sample_shoulder_strake_state(config, adapter, aero_case, policy, state + 0.5 * dt_s * k2, time_s + 0.5 * dt_s, density_scale, use_winds, assumptions, strake_area_fraction_of_stowed_ref, shoulder_command, strake_count).derivative
    k4 = _sample_shoulder_strake_state(config, adapter, aero_case, policy, state + dt_s * k3, time_s + dt_s, density_scale, use_winds, assumptions, strake_area_fraction_of_stowed_ref, shoulder_command, strake_count).derivative
    return state + (dt_s / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
end

function solve_shield_shoulder_strake_trajectory_3d(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
    policy::ControlPolicy,
    assumptions::ShoulderStrakeAssumptions;
    strake_area_fraction_of_stowed_ref::Real,
    shoulder_command::ShoulderStrakeTrajectoryCommand=ShoulderStrakeTrajectoryCommand(0.0, 0.0, 0.0),
    strake_count::Integer=4,
    save_trajectory::Bool=true,
    h0_m::Real=config.h0_m,
    V0_mps::Real=config.V0_mps,
    gamma0_deg::Real=config.gamma0_deg,
    heading_deg::Real=90.0,
    density_scale::Real=1.0,
    use_winds::Bool=false,
)
    r0 = SVector{3, Float64}(config.planet.Rp_m + Float64(h0_m), 0.0, 0.0)
    east0 = SVector{3, Float64}(0.0, 1.0, 0.0)
    north0 = SVector{3, Float64}(0.0, 0.0, 1.0)
    up0 = SVector{3, Float64}(1.0, 0.0, 0.0)
    γ0 = deg2rad(Float64(gamma0_deg))
    ψ0 = deg2rad(Float64(heading_deg))
    v_h = Float64(V0_mps) * cos(γ0)
    v_local0 = SVector{3, Float64}(v_h * sin(ψ0), v_h * cos(ψ0), Float64(V0_mps) * sin(γ0))
    v0 = v_local0[1] * east0 + v_local0[2] * north0 + v_local0[3] * up0
    current_state = SVector{6, Float64}(r0[1], r0[2], r0[3], v0[1], v0[2], v0[3])
    current_time_s = 0.0
    density_scale_float = Float64(density_scale)
    integration_dt_s = config.saveat_s

    time_s_vec = Float64[]; altitude_km = Float64[]; downrange_km = Float64[]; crossrange_km = Float64[]
    heading_deg_vec = Float64[]; velocity_mps = Float64[]; flight_path_angle_deg = Float64[]
    rho_kg_m3 = Float64[]; temperature_k = Float64[]; q_pa = Float64[]
    drag_accel_mps2 = Float64[]; lift_accel_mps2 = Float64[]; side_accel_mps2 = Float64[]
    total_aero_g = Float64[]; heating_proxy = Float64[]; mach_values = Float64[]
    beta_eff_kg_m2 = Float64[]; wind_east_mps = Float64[]; wind_north_mps = Float64[]; wind_up_mps = Float64[]
    lateral_command_vec = Float64[]; deployed_vec = Bool[]; state_label_vec = String[]; kn_values = Float64[]

    current_sample = _sample_shoulder_strake_state(config, adapter, aero_case, policy, current_state, current_time_s, density_scale_float, use_winds, assumptions, Float64(strake_area_fraction_of_stowed_ref), shoulder_command, strake_count)
    peak_dynamic_pressure_pa = current_sample.q_pa
    peak_drag_accel_mps2 = current_sample.drag_accel_mps2
    peak_side_accel_mps2 = current_sample.side_accel_mps2
    peak_total_decel_earth_g = current_sample.total_aero_g
    peak_heating_proxy = current_sample.heating_proxy
    final_sample = current_sample
    final_time_s = current_time_s

    if save_trajectory
        _push_trajectory_sample!(time_s_vec, altitude_km, downrange_km, crossrange_km, heading_deg_vec, velocity_mps, flight_path_angle_deg, rho_kg_m3, temperature_k, q_pa, drag_accel_mps2, lift_accel_mps2, side_accel_mps2, total_aero_g, heating_proxy, mach_values, beta_eff_kg_m2, wind_east_mps, wind_north_mps, wind_up_mps, lateral_command_vec, deployed_vec, state_label_vec, kn_values, current_time_s, current_sample)
    end

    integration_stopped_at_surface = false
    while current_time_s < config.max_time_s
        current_sample.altitude_km <= 0.0 && break
        dt_s = min(integration_dt_s, config.max_time_s - current_time_s)
        next_state = _rk4_step_shoulder_strake(config, adapter, aero_case, policy, current_state, current_sample, current_time_s, dt_s, density_scale_float, use_winds, assumptions, Float64(strake_area_fraction_of_stowed_ref), shoulder_command, strake_count)
        next_time_s = current_time_s + dt_s
        next_sample = _sample_shoulder_strake_state(config, adapter, aero_case, policy, next_state, next_time_s, density_scale_float, use_winds, assumptions, Float64(strake_area_fraction_of_stowed_ref), shoulder_command, strake_count)

        if next_sample.altitude_km <= 0.0
            h_prev_m = current_sample.altitude_km * 1e3
            h_next_m = next_sample.altitude_km * 1e3
            frac = clamp(h_prev_m / max(h_prev_m - h_next_m, eps()), 0.0, 1.0)
            impact_time_s = current_time_s + frac * dt_s
            impact_state = current_state + frac * (next_state - current_state)
            impact_sample = _sample_shoulder_strake_state(config, adapter, aero_case, policy, impact_state, impact_time_s, density_scale_float, use_winds, assumptions, Float64(strake_area_fraction_of_stowed_ref), shoulder_command, strake_count)
            peak_dynamic_pressure_pa = max(peak_dynamic_pressure_pa, impact_sample.q_pa)
            peak_drag_accel_mps2 = max(peak_drag_accel_mps2, impact_sample.drag_accel_mps2)
            peak_side_accel_mps2 = max(peak_side_accel_mps2, impact_sample.side_accel_mps2)
            peak_total_decel_earth_g = max(peak_total_decel_earth_g, impact_sample.total_aero_g)
            peak_heating_proxy = max(peak_heating_proxy, impact_sample.heating_proxy)
            final_sample = impact_sample
            final_time_s = impact_time_s
            if save_trajectory
                _push_trajectory_sample!(time_s_vec, altitude_km, downrange_km, crossrange_km, heading_deg_vec, velocity_mps, flight_path_angle_deg, rho_kg_m3, temperature_k, q_pa, drag_accel_mps2, lift_accel_mps2, side_accel_mps2, total_aero_g, heating_proxy, mach_values, beta_eff_kg_m2, wind_east_mps, wind_north_mps, wind_up_mps, lateral_command_vec, deployed_vec, state_label_vec, kn_values, impact_time_s, impact_sample)
            end
            integration_stopped_at_surface = true
            break
        end

        current_state = next_state
        current_time_s = next_time_s
        current_sample = next_sample
        peak_dynamic_pressure_pa = max(peak_dynamic_pressure_pa, current_sample.q_pa)
        peak_drag_accel_mps2 = max(peak_drag_accel_mps2, current_sample.drag_accel_mps2)
        peak_side_accel_mps2 = max(peak_side_accel_mps2, current_sample.side_accel_mps2)
        peak_total_decel_earth_g = max(peak_total_decel_earth_g, current_sample.total_aero_g)
        peak_heating_proxy = max(peak_heating_proxy, current_sample.heating_proxy)
        final_sample = current_sample
        final_time_s = current_time_s
        if save_trajectory
            _push_trajectory_sample!(time_s_vec, altitude_km, downrange_km, crossrange_km, heading_deg_vec, velocity_mps, flight_path_angle_deg, rho_kg_m3, temperature_k, q_pa, drag_accel_mps2, lift_accel_mps2, side_accel_mps2, total_aero_g, heating_proxy, mach_values, beta_eff_kg_m2, wind_east_mps, wind_north_mps, wind_up_mps, lateral_command_vec, deployed_vec, state_label_vec, kn_values, current_time_s, current_sample)
        end
    end

    trajectory_df = save_trajectory ? DataFrame(
        time_s=time_s_vec, altitude_km=altitude_km, downrange_km=downrange_km, crossrange_km=crossrange_km,
        heading_deg=heading_deg_vec, velocity_mps=velocity_mps, flight_path_angle_deg=flight_path_angle_deg,
        rho_kg_m3=rho_kg_m3, temperature_k=temperature_k, q_pa=q_pa, drag_accel_mps2=drag_accel_mps2,
        lift_accel_mps2=lift_accel_mps2, side_accel_mps2=side_accel_mps2, total_aero_g=total_aero_g,
        heating_proxy=heating_proxy, mach=mach_values, beta_eff_kg_m2=beta_eff_kg_m2,
        wind_east_mps=wind_east_mps, wind_north_mps=wind_north_mps, wind_up_mps=wind_up_mps,
        lateral_command=lateral_command_vec, kn=kn_values, deployed=deployed_vec, state_label=state_label_vec,
    ) : DataFrame()

    impact_velocity_mps = final_sample.velocity_mps
    impact_g_load_proxy = impact_velocity_mps^2 / (2.0 * config.impact_stop_distance_m * _earth_g())
    summary = (
        atmosphere_path=atmosphere_label(adapter),
        impact_downrange_km=final_sample.downrange_km,
        impact_crossrange_km=final_sample.crossrange_km,
        flight_time_s=final_time_s,
        impact_velocity_mps=impact_velocity_mps,
        impact_flight_path_angle_deg=final_sample.flight_path_angle_deg,
        impact_heading_deg=final_sample.heading_deg,
        peak_dynamic_pressure_pa=peak_dynamic_pressure_pa,
        peak_drag_accel_mps2=peak_drag_accel_mps2,
        peak_side_accel_mps2=peak_side_accel_mps2,
        peak_total_decel_earth_g=peak_total_decel_earth_g,
        max_heating_proxy=peak_heating_proxy,
        impact_g_load_proxy_earth_g=impact_g_load_proxy,
        integration_stopped_at_surface=integration_stopped_at_surface,
    )
    return Trajectory3DResult(summary, trajectory_df)
end

function run_shield_shoulder_strake_trajectory_study(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter;
    entry_mass_kg::Real=120.0,
    strake_count::Integer=4,
    cg_axial_fraction_of_body_length::Real=0.60,
    static_margin_fraction_of_diameter::Real=0.10,
    strake_area_fraction_of_stowed_ref::Real=0.05,
    control_deflection_deg::Real=15.0,
    trim_check_angle_deg::Real=5.0,
    cl_alpha_per_rad::Real=3.5,
    cd0::Real=0.05,
    induced_drag_factor::Real=0.15,
    cl_max::Real=1.2,
    pitch_command_grid::AbstractVector{<:Real}=collect(-1.0:0.1:1.0),
    yaw_command_grid::AbstractVector{<:Real}=[-1.0, 1.0],
)
    mkpath(config.output_root)
    baseline_bundle = _representative_shield_body_flap_baseline_condition(config, adapter; entry_mass_kg=entry_mass_kg)
    baseline = baseline_bundle.condition
    aero_case = baseline_bundle.aero_case
    targets = baseline_bundle.targets
    deploy_h_m = 0.5 * floor(baseline.subsonic_onset_km / 0.5) * 1e3
    assumptions = ShoulderStrakeAssumptions(
        Float64(cg_axial_fraction_of_body_length),
        Float64(static_margin_fraction_of_diameter),
        Float64(control_deflection_deg),
        Float64(trim_check_angle_deg),
        Float64(cl_alpha_per_rad),
        Float64(cd0),
        Float64(induced_drag_factor),
        Float64(cl_max),
    )
    policy = ControlPolicy("shoulder_strake_guided", :bang_bang, deploy_h_m, 0.0)

    body_only = solve_entry_trajectory_3d(config, adapter, aero_case, ControlPolicy("body_only", :body_only, NaN, NaN); save_trajectory=true, use_winds=false)
    passive_skirt = solve_entry_trajectory_3d(config, adapter, aero_case, policy; save_trajectory=true, use_winds=false)

    pitch_rows = NamedTuple[]
    pitch_results = Dict{Float64, Trajectory3DResult}()
    for cmd in Float64.(collect(pitch_command_grid))
        result = solve_shield_shoulder_strake_trajectory_3d(
            config, adapter, aero_case, policy, assumptions;
            strake_area_fraction_of_stowed_ref=strake_area_fraction_of_stowed_ref,
            shoulder_command=ShoulderStrakeTrajectoryCommand(0.0, cmd, 0.0),
            strake_count=strake_count,
            save_trajectory=true,
            use_winds=false,
        )
        pitch_results[cmd] = result
        push!(pitch_rows, (
            pitch_command_scale=cmd,
            pitch_deflection_deg=cmd * assumptions.control_deflection_deg,
            impact_downrange_km=result.summary.impact_downrange_km,
            impact_crossrange_km=result.summary.impact_crossrange_km,
            impact_velocity_mps=result.summary.impact_velocity_mps,
            peak_total_decel_earth_g=result.summary.peak_total_decel_earth_g,
        ))
    end
    pitch_df = sort(DataFrame(pitch_rows), :pitch_command_scale)
    min_row = pitch_df[argmin(pitch_df.impact_downrange_km), :]
    max_row = pitch_df[argmax(pitch_df.impact_downrange_km), :]
    target_range_km = min_row.impact_downrange_km + config.target_range_fraction * (max_row.impact_downrange_km - min_row.impact_downrange_km)
    target_idx = argmin(abs.(pitch_df.impact_downrange_km .- target_range_km))
    target_pitch_cmd = pitch_df.pitch_command_scale[target_idx]
    target_result = solve_shield_shoulder_strake_trajectory_3d(
        config, adapter, aero_case, policy, assumptions;
        strake_area_fraction_of_stowed_ref=strake_area_fraction_of_stowed_ref,
        shoulder_command=ShoulderStrakeTrajectoryCommand(0.0, target_pitch_cmd, 0.0),
        strake_count=strake_count,
        save_trajectory=true,
        use_winds=false,
    )

    yaw_rows = NamedTuple[]
    yaw_results = Dict{Float64, Trajectory3DResult}()
    for cmd in Float64.(collect(yaw_command_grid))
        result = solve_shield_shoulder_strake_trajectory_3d(
            config, adapter, aero_case, policy, assumptions;
            strake_area_fraction_of_stowed_ref=strake_area_fraction_of_stowed_ref,
            shoulder_command=ShoulderStrakeTrajectoryCommand(0.0, 0.0, cmd),
            strake_count=strake_count,
            save_trajectory=true,
            use_winds=false,
        )
        yaw_results[cmd] = result
        push!(yaw_rows, (
            yaw_command_scale=cmd,
            yaw_deflection_deg=cmd * assumptions.control_deflection_deg,
            impact_downrange_km=result.summary.impact_downrange_km,
            impact_crossrange_km=result.summary.impact_crossrange_km,
            impact_velocity_mps=result.summary.impact_velocity_mps,
            peak_total_decel_earth_g=result.summary.peak_total_decel_earth_g,
        ))
    end
    yaw_df = sort(DataFrame(yaw_rows), :yaw_command_scale)

    summary_rows = [
        (case="body_only", command_value=NaN, impact_downrange_km=body_only.summary.impact_downrange_km, impact_crossrange_km=body_only.summary.impact_crossrange_km, impact_velocity_mps=body_only.summary.impact_velocity_mps, peak_total_decel_earth_g=body_only.summary.peak_total_decel_earth_g),
        (case="passive_skirt_subsonic", command_value=NaN, impact_downrange_km=passive_skirt.summary.impact_downrange_km, impact_crossrange_km=passive_skirt.summary.impact_crossrange_km, impact_velocity_mps=passive_skirt.summary.impact_velocity_mps, peak_total_decel_earth_g=passive_skirt.summary.peak_total_decel_earth_g),
        (case="pitch_min_range", command_value=min_row.pitch_command_scale, impact_downrange_km=min_row.impact_downrange_km, impact_crossrange_km=min_row.impact_crossrange_km, impact_velocity_mps=min_row.impact_velocity_mps, peak_total_decel_earth_g=min_row.peak_total_decel_earth_g),
        (case="pitch_max_range", command_value=max_row.pitch_command_scale, impact_downrange_km=max_row.impact_downrange_km, impact_crossrange_km=max_row.impact_crossrange_km, impact_velocity_mps=max_row.impact_velocity_mps, peak_total_decel_earth_g=max_row.peak_total_decel_earth_g),
        (case="pitch_targeted", command_value=target_pitch_cmd, impact_downrange_km=target_result.summary.impact_downrange_km, impact_crossrange_km=target_result.summary.impact_crossrange_km, impact_velocity_mps=target_result.summary.impact_velocity_mps, peak_total_decel_earth_g=target_result.summary.peak_total_decel_earth_g),
    ]
    for row in eachrow(yaw_df)
        push!(summary_rows, (case=row.yaw_command_scale < 0 ? "yaw_left" : "yaw_right", command_value=row.yaw_command_scale, impact_downrange_km=row.impact_downrange_km, impact_crossrange_km=row.impact_crossrange_km, impact_velocity_mps=row.impact_velocity_mps, peak_total_decel_earth_g=row.peak_total_decel_earth_g))
    end
    summary_df = DataFrame(summary_rows)
    summary_path = joinpath(config.output_root, "shield_shoulder_strake_guided_summary.csv")
    pitch_path = joinpath(config.output_root, "shield_shoulder_strake_pitch_sweep.csv")
    yaw_path = joinpath(config.output_root, "shield_shoulder_strake_yaw_sweep.csv")
    note_path = joinpath(config.output_root, "shield_shoulder_strake_guided_note.md")
    CSV.write(summary_path, summary_df)
    CSV.write(pitch_path, pitch_df)
    CSV.write(yaw_path, yaw_df)
    mkpath(config.trajectory_dir)
    CSV.write(joinpath(config.trajectory_dir, "body_only.csv"), body_only.trajectory)
    CSV.write(joinpath(config.trajectory_dir, "passive_skirt_subsonic.csv"), passive_skirt.trajectory)
    for (cmd, result) in sort(collect(pitch_results); by=first)
        cmd_label = @sprintf("%+04.1f", cmd)
        cmd_label = replace(cmd_label, "." => "p")
        CSV.write(joinpath(config.trajectory_dir, "pitch_cmd_$(cmd_label).csv"), result.trajectory)
    end
    for (cmd, result) in sort(collect(yaw_results); by=first)
        cmd_label = @sprintf("%+04.1f", cmd)
        cmd_label = replace(cmd_label, "." => "p")
        CSV.write(joinpath(config.trajectory_dir, "yaw_cmd_$(cmd_label).csv"), result.trajectory)
    end
    CSV.write(joinpath(config.trajectory_dir, "pitch_min_range.csv"), pitch_results[min_row.pitch_command_scale].trajectory)
    CSV.write(joinpath(config.trajectory_dir, "pitch_max_range.csv"), pitch_results[max_row.pitch_command_scale].trajectory)
    CSV.write(joinpath(config.trajectory_dir, "pitch_targeted.csv"), target_result.trajectory)
    CSV.write(joinpath(config.trajectory_dir, "yaw_left.csv"), yaw_results[minimum(Float64.(collect(yaw_command_grid)))].trajectory)
    CSV.write(joinpath(config.trajectory_dir, "yaw_right.csv"), yaw_results[maximum(Float64.(collect(yaw_command_grid)))].trajectory)
    open(note_path, "w") do io
        println(io, "# SHIELD Shoulder-Strake Guided Trajectory Study")
        println(io)
        println(io, "This is the first deterministic trajectory rerun using a shoulder-petal / strake control surrogate that cleared the trim gate in the prior stability screen.")
        println(io)
        println(io, "## Nominal Guided Architecture")
        println(io, @sprintf("- CG fraction: %.2f body lengths", assumptions.cg_axial_fraction_of_body_length))
        println(io, @sprintf("- Static margin: %.0f%% of diameter", 100.0 * assumptions.static_margin_fraction_of_diameter))
        println(io, @sprintf("- Total shoulder-strake area: %.1f%% of stowed reference area", 100.0 * Float64(strake_area_fraction_of_stowed_ref)))
        println(io, @sprintf("- Max control deflection: %.0f deg", assumptions.control_deflection_deg))
        println(io, @sprintf("- Passive skirt deployment altitude: %.2f km", deploy_h_m / 1e3))
        println(io, "- Pitch control uses equal-and-opposite top/bottom strake deflections; yaw control uses equal-and-opposite left/right strake deflections.")
        println(io)
        println(io, "## Deterministic Results")
        println(io, @sprintf("- Body only impact range: %.2f km", body_only.summary.impact_downrange_km))
        println(io, @sprintf("- Passive subsonic skirt deployment impact range: %.2f km", passive_skirt.summary.impact_downrange_km))
        println(io, @sprintf("- Pitch-command range envelope: %.2f km to %.2f km", min_row.impact_downrange_km, max_row.impact_downrange_km))
        println(io, @sprintf("- Downrange authority from pitch command sweep: %.2f km", max_row.impact_downrange_km - min_row.impact_downrange_km))
        println(io, @sprintf("- Targeted pitch command: %.2f -> impact range %.2f km for nominal target %.2f km", target_pitch_cmd, target_result.summary.impact_downrange_km, target_range_km))
        println(io, @sprintf("- Yaw-left / yaw-right crossrange: %.2f km / %.2f km", yaw_results[minimum(Float64.(collect(yaw_command_grid)))].summary.impact_crossrange_km, yaw_results[maximum(Float64.(collect(yaw_command_grid)))].summary.impact_crossrange_km))
        println(io)
        println(io, "## Caveat")
        println(io, "This is still a point-mass trajectory rerun using trimmed pre-skirt force commands. It assumes the shoulder-strake architecture can hold the screened force state; it is not a 6DOF attitude simulation.")
    end
    return (
        baseline_condition=baseline,
        targets=targets,
        aero_case=aero_case,
        assumptions=assumptions,
        deploy_h_m=deploy_h_m,
        body_only=body_only,
        passive_skirt=passive_skirt,
        pitch_df=pitch_df,
        yaw_df=yaw_df,
        target_pitch_cmd=target_pitch_cmd,
        target_range_km=target_range_km,
        target_result=target_result,
        summary_df=summary_df,
        summary_path=summary_path,
        pitch_path=pitch_path,
        yaw_path=yaw_path,
        note_path=note_path,
    )
end
