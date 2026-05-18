struct BodyFlapComponent
    name::String
    azimuth_deg::Float64
    area_m2::Float64
    position_m::SVector{3, Float64}
    radial_dir::SVector{3, Float64}
end

struct BodyFlapBaselineCondition
    subsonic_onset_km::Float64
    representative_altitude_km::Float64
    representative_mach::Float64
    representative_q_pa::Float64
    representative_velocity_mps::Float64
end

struct BodyFlapStabilityAssumptions
    cg_axial_fraction_of_body_length::Float64
    static_margin_fraction_of_diameter::Float64
    trim_check_angle_deg::Float64
    control_deflection_deg::Float64
end

@inline function shield_stowed_reference_area_m2(config::PrelimConfig)
    return π * config.geometry.base_radius_m^2
end

@inline function shield_stowed_surface_areal_density_proxy_kg_m2(config::PrelimConfig, entry_mass_kg::Real)
    area_m2 = shield_stowed_reference_area_m2(config)
    area_m2 > 0.0 || throw(ArgumentError("Stowed reference area must be positive for SHIELD pre-separation control-surface mass proxy."))
    return Float64(entry_mass_kg) / area_m2
end

@inline function shield_body_flap_areal_density_proxy_kg_m2(config::PrelimConfig, entry_mass_kg::Real)
    return shield_stowed_surface_areal_density_proxy_kg_m2(config, entry_mass_kg)
end

@inline function _stowed_body_flap_axial_arm_m(config::PrelimConfig)
    return config.geometry.base_radius_m / tan(deg2rad(config.geometry.cone_half_angle_deg))
end

@inline function _stowed_body_reference_length_m(config::PrelimConfig)
    return 2.0 * config.geometry.base_radius_m
end

@inline function _body_flap_control_outputs(Δforce_area::SVector{3, Float64}, Δmoment_area::SVector{3, Float64})
    return (
        delta_drag_area_m2=Δforce_area[1],
        delta_side_area_m2=Δforce_area[2],
        delta_vertical_area_m2=Δforce_area[3],
        delta_roll_moment_area_m3=Δmoment_area[1],
        delta_pitch_up_moment_area_m3=-Δmoment_area[2],
        delta_yaw_right_moment_area_m3=Δmoment_area[3],
    )
end

function build_shield_body_flaps(
    config::PrelimConfig;
    flap_area_fraction_of_stowed_ref::Real,
    flap_count::Integer=4,
)
    flap_count == 4 || throw(ArgumentError("This first-cut SHIELD body-flap study is defined for 4 flaps."))
    total_flap_area_m2 = Float64(flap_area_fraction_of_stowed_ref) * shield_stowed_reference_area_m2(config)
    each_flap_area_m2 = total_flap_area_m2 / flap_count
    x_arm_m = _stowed_body_flap_axial_arm_m(config)
    radius_m = config.geometry.base_radius_m

    layout = [
        ("top", 0.0),
        ("right", 90.0),
        ("bottom", 180.0),
        ("left", 270.0),
    ]
    flaps = BodyFlapComponent[]
    for (name, azimuth_deg) in layout
        ϕ = deg2rad(azimuth_deg)
        radial_dir = SVector{3, Float64}(0.0, sin(ϕ), cos(ϕ))
        position_m = SVector{3, Float64}(x_arm_m, radius_m * sin(ϕ), radius_m * cos(ϕ))
        push!(flaps, BodyFlapComponent(name, azimuth_deg, each_flap_area_m2, position_m, radial_dir))
    end
    return (
        flap_count=flap_count,
        flap_area_total_m2=total_flap_area_m2,
        flap_area_each_m2=each_flap_area_m2,
        flaps=flaps,
        x_arm_m=x_arm_m,
    )
end

function _body_flap_force_moment_delta(
    config::PrelimConfig,
    flap::BodyFlapComponent,
    deflection_deg::Real,
    cg_x_m::Real=0.0,
)
    plate = FlatPlateReferenceBody(
        flap.area_m2,
        config.geometry.panel_skin_friction_coefficient,
        config.geometry.panel_zero_lift_drag_coefficient,
    )
    α_rad = deg2rad(Float64(deflection_deg))
    CL, CD = _flat_plate_cl_cd(α_rad, plate)
    area = flap.area_m2
    Δforce_area = SVector{3, Float64}(CD * area, CL * area * flap.radial_dir[2], CL * area * flap.radial_dir[3])
    r_cg_m = flap.position_m - SVector{3, Float64}(Float64(cg_x_m), 0.0, 0.0)
    Δmoment_area = cross(r_cg_m, Δforce_area)
    return Δforce_area, Δmoment_area
end

function _body_flap_command_deflections(mode::String, deflection_deg::Real)
    δ = Float64(deflection_deg)
    if mode == "neutral"
        return Dict("top" => 0.0, "right" => 0.0, "bottom" => 0.0, "left" => 0.0)
    elseif mode == "common_drag"
        return Dict("top" => δ, "right" => δ, "bottom" => δ, "left" => δ)
    elseif mode == "yaw_right"
        return Dict("top" => 0.0, "right" => δ, "bottom" => 0.0, "left" => -δ)
    elseif mode == "yaw_left"
        return Dict("top" => 0.0, "right" => -δ, "bottom" => 0.0, "left" => δ)
    elseif mode == "pitch_up"
        return Dict("top" => δ, "right" => 0.0, "bottom" => -δ, "left" => 0.0)
    elseif mode == "pitch_down"
        return Dict("top" => -δ, "right" => 0.0, "bottom" => δ, "left" => 0.0)
    elseif startswith(mode, "single_")
        flap_name = replace(mode, "single_" => "")
        return Dict("top" => 0.0, "right" => 0.0, "bottom" => 0.0, "left" => 0.0, flap_name => δ)
    end
    throw(ArgumentError("Unsupported body-flap command mode $(repr(mode))."))
end

function evaluate_shield_body_flap_command(
    config::PrelimConfig,
    baseline::BodyFlapBaselineCondition;
    flap_area_fraction_of_stowed_ref::Real,
    deflection_deg::Real,
    mode::String,
    cg_x_m::Real=0.0,
)
    layout = build_shield_body_flaps(config; flap_area_fraction_of_stowed_ref=flap_area_fraction_of_stowed_ref)
    deflections = _body_flap_command_deflections(mode, deflection_deg)

    Δforce_area = SVector{3, Float64}(0.0, 0.0, 0.0)
    Δmoment_area = SVector{3, Float64}(0.0, 0.0, 0.0)
    for flap in layout.flaps
        flap_Δforce_area, flap_Δmoment_area = _body_flap_force_moment_delta(config, flap, deflections[flap.name], cg_x_m)
        Δforce_area += flap_Δforce_area
        Δmoment_area += flap_Δmoment_area
    end

    outputs = _body_flap_control_outputs(Δforce_area, Δmoment_area)
    q_ref = baseline.representative_q_pa
    return (
        mode=mode,
        flap_area_fraction_of_stowed_ref=Float64(flap_area_fraction_of_stowed_ref),
        flap_area_total_m2=layout.flap_area_total_m2,
        flap_area_each_m2=layout.flap_area_each_m2,
        body_flap_axial_arm_m=layout.x_arm_m,
        deflection_deg=Float64(deflection_deg),
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

function _representative_shield_body_flap_baseline_condition(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter;
    entry_mass_kg::Real=120.0,
)
    targets = derive_fixed_mass_beta_targets_from_deployed_geometry(
        config,
        adapter;
        entry_mass_kg=entry_mass_kg,
    )
    aero_case = calibrate_aero_case(config, adapter, targets.beta_high_kg_m2, targets.beta_ratio)
    body_result = solve_entry_trajectory(
        config,
        adapter,
        aero_case,
        ControlPolicy("body_only", :body_only, NaN, NaN);
        save_trajectory=true,
    )
    subsonic_idx = findfirst(x -> x <= 1.0, body_result.trajectory.mach)
    subsonic_idx === nothing && throw(ArgumentError("Published-size SHIELD surrogate never becomes subsonic before impact."))
    supersonic_segment = filter(row -> row.mach > 1.0, body_result.trajectory)
    nrow(supersonic_segment) > 0 || throw(ArgumentError("Published-size SHIELD surrogate has no supersonic body-only segment for body-flap screening."))
    q_idx = argmax(supersonic_segment.q_pa)
    sample = supersonic_segment[q_idx, :]
    return (
        condition=BodyFlapBaselineCondition(
            body_result.trajectory.altitude_km[subsonic_idx],
            sample.altitude_km,
            sample.mach,
            sample.q_pa,
            sample.velocity_mps,
        ),
        targets=targets,
        aero_case=aero_case,
    )
end

function _body_flap_derivative_rows(config::PrelimConfig, baseline::BodyFlapBaselineCondition, flap_area_fraction::Float64; cg_x_m::Real=0.0)
    rows = NamedTuple[]
    for mode in ("single_top", "single_right", "single_bottom", "single_left", "common_drag", "yaw_right", "pitch_up")
        positive = evaluate_shield_body_flap_command(config, baseline; flap_area_fraction_of_stowed_ref=flap_area_fraction, deflection_deg=1.0, mode=mode, cg_x_m=cg_x_m)
        negative = evaluate_shield_body_flap_command(config, baseline; flap_area_fraction_of_stowed_ref=flap_area_fraction, deflection_deg=-1.0, mode=mode, cg_x_m=cg_x_m)
        push!(rows, (
            mode=mode,
            flap_area_fraction_of_stowed_ref=flap_area_fraction,
            d_drag_area_m2_per_deg=0.5 * (positive.delta_drag_area_m2 - negative.delta_drag_area_m2),
            d_side_area_m2_per_deg=0.5 * (positive.delta_side_area_m2 - negative.delta_side_area_m2),
            d_vertical_area_m2_per_deg=0.5 * (positive.delta_vertical_area_m2 - negative.delta_vertical_area_m2),
            d_roll_moment_area_m3_per_deg=0.5 * (positive.delta_roll_moment_area_m3 - negative.delta_roll_moment_area_m3),
            d_pitch_up_moment_area_m3_per_deg=0.5 * (positive.delta_pitch_up_moment_area_m3 - negative.delta_pitch_up_moment_area_m3),
            d_yaw_right_moment_area_m3_per_deg=0.5 * (positive.delta_yaw_right_moment_area_m3 - negative.delta_yaw_right_moment_area_m3),
            d_side_force_N_per_deg=0.5 * (positive.delta_side_force_N - negative.delta_side_force_N),
            d_vertical_force_N_per_deg=0.5 * (positive.delta_vertical_force_N - negative.delta_vertical_force_N),
            d_pitch_up_moment_Nm_per_deg=0.5 * (positive.delta_pitch_up_moment_Nm - negative.delta_pitch_up_moment_Nm),
            d_yaw_right_moment_Nm_per_deg=0.5 * (positive.delta_yaw_right_moment_Nm - negative.delta_yaw_right_moment_Nm),
        ))
    end
    return DataFrame(rows)
end

@inline function _body_flap_cg_location_m(config::PrelimConfig, assumptions::BodyFlapStabilityAssumptions)
    return assumptions.cg_axial_fraction_of_body_length * _stowed_body_flap_axial_arm_m(config)
end

@inline function _body_flap_cp_location_m(config::PrelimConfig, assumptions::BodyFlapStabilityAssumptions)
    return _body_flap_cg_location_m(config, assumptions) -
           assumptions.static_margin_fraction_of_diameter * _stowed_body_reference_length_m(config)
end

function _baseline_body_force_derivatives(config::PrelimConfig, aero_case::CalibratedAeroCase, representative_mach::Float64)
    body_link = aero_case.body_state.body_link
    γ = config.planet.γ
    δdeg = 1.0
    δrad = deg2rad(δdeg)
    cla_plus_m2, cda_plus_m2 = _component_cla_cda(body_link, δrad, representative_mach, γ)
    cla_minus_m2, cda_minus_m2 = _component_cla_cda(body_link, -δrad, representative_mach, γ)
    cla_zero_m2, cda_zero_m2 = _component_cla_cda(body_link, 0.0, representative_mach, γ)
    dcla_dalpha_m2_per_deg = (cla_plus_m2 - cla_minus_m2) / (2.0 * δdeg)
    dcda_dalpha_m2_per_deg = (cda_plus_m2 - cda_minus_m2) / (2.0 * δdeg)
    return (
        cla_zero_m2=cla_zero_m2,
        cda_zero_m2=cda_zero_m2,
        dcla_dalpha_m2_per_deg=dcla_dalpha_m2_per_deg,
        dcla_dalpha_m2_per_rad=dcla_dalpha_m2_per_deg * (180.0 / π),
        dcda_dalpha_m2_per_deg=dcda_dalpha_m2_per_deg,
        dcda_dalpha_m2_per_rad=dcda_dalpha_m2_per_deg * (180.0 / π),
    )
end

function _baseline_stowed_body_stability_row(
    config::PrelimConfig,
    baseline::BodyFlapBaselineCondition,
    aero_case::CalibratedAeroCase,
    targets,
    assumptions::BodyFlapStabilityAssumptions,
)
    body_force = _baseline_body_force_derivatives(config, aero_case, baseline.representative_mach)
    S_ref_m2 = shield_stowed_reference_area_m2(config)
    L_ref_m = _stowed_body_reference_length_m(config)
    x_cg_m = _body_flap_cg_location_m(config, assumptions)
    x_cp_m = _body_flap_cp_location_m(config, assumptions)
    restoring_arm_m = x_cp_m - x_cg_m

    d_pitch_up_moment_area_m3_per_deg = restoring_arm_m * body_force.dcla_dalpha_m2_per_deg
    d_pitch_up_moment_area_m3_per_rad = restoring_arm_m * body_force.dcla_dalpha_m2_per_rad
    d_yaw_right_moment_area_m3_per_deg = restoring_arm_m * body_force.dcla_dalpha_m2_per_deg
    d_yaw_right_moment_area_m3_per_rad = restoring_arm_m * body_force.dcla_dalpha_m2_per_rad
    pitch_restoring = d_pitch_up_moment_area_m3_per_deg * body_force.dcla_dalpha_m2_per_deg < 0.0
    yaw_restoring = d_yaw_right_moment_area_m3_per_deg * body_force.dcla_dalpha_m2_per_deg < 0.0

    q_ref = baseline.representative_q_pa
    return (
        stowed_beta_high_kg_m2=targets.beta_high_kg_m2,
        deployed_beta_low_kg_m2=targets.beta_low_kg_m2,
        representative_altitude_km=baseline.representative_altitude_km,
        representative_mach=baseline.representative_mach,
        representative_q_pa=q_ref,
        representative_velocity_mps=baseline.representative_velocity_mps,
        cg_axial_fraction_of_body_length=assumptions.cg_axial_fraction_of_body_length,
        static_margin_fraction_of_diameter=assumptions.static_margin_fraction_of_diameter,
        trim_check_angle_deg=assumptions.trim_check_angle_deg,
        control_deflection_deg=assumptions.control_deflection_deg,
        body_length_m=_stowed_body_flap_axial_arm_m(config),
        reference_length_m=L_ref_m,
        stowed_reference_area_m2=S_ref_m2,
        x_cg_m=x_cg_m,
        x_cp_m=x_cp_m,
        restoring_arm_m=restoring_arm_m,
        cla0_m2=body_force.cla_zero_m2,
        cda0_m2=body_force.cda_zero_m2,
        dcz_area_m2_per_deg=body_force.dcla_dalpha_m2_per_deg,
        dcz_area_m2_per_rad=body_force.dcla_dalpha_m2_per_rad,
        dcy_area_m2_per_deg=body_force.dcla_dalpha_m2_per_deg,
        dcy_area_m2_per_rad=body_force.dcla_dalpha_m2_per_rad,
        d_drag_area_m2_per_deg=body_force.dcda_dalpha_m2_per_deg,
        d_pitch_up_moment_area_m3_per_deg=d_pitch_up_moment_area_m3_per_deg,
        d_pitch_up_moment_area_m3_per_rad=d_pitch_up_moment_area_m3_per_rad,
        d_yaw_right_moment_area_m3_per_deg=d_yaw_right_moment_area_m3_per_deg,
        d_yaw_right_moment_area_m3_per_rad=d_yaw_right_moment_area_m3_per_rad,
        CZ_alpha_per_deg=body_force.dcla_dalpha_m2_per_deg / S_ref_m2,
        CZ_alpha_per_rad=body_force.dcla_dalpha_m2_per_rad / S_ref_m2,
        CY_beta_per_deg=body_force.dcla_dalpha_m2_per_deg / S_ref_m2,
        CY_beta_per_rad=body_force.dcla_dalpha_m2_per_rad / S_ref_m2,
        Cma_per_deg=d_pitch_up_moment_area_m3_per_deg / (S_ref_m2 * L_ref_m),
        Cma_per_rad=d_pitch_up_moment_area_m3_per_rad / (S_ref_m2 * L_ref_m),
        Cnb_per_deg=d_yaw_right_moment_area_m3_per_deg / (S_ref_m2 * L_ref_m),
        Cnb_per_rad=d_yaw_right_moment_area_m3_per_rad / (S_ref_m2 * L_ref_m),
        recovered_pitch_static_margin=-((d_pitch_up_moment_area_m3_per_rad / (S_ref_m2 * L_ref_m)) / (body_force.dcla_dalpha_m2_per_rad / S_ref_m2)),
        recovered_yaw_static_margin=-((d_yaw_right_moment_area_m3_per_rad / (S_ref_m2 * L_ref_m)) / (body_force.dcla_dalpha_m2_per_rad / S_ref_m2)),
        d_vertical_force_N_per_deg=q_ref * body_force.dcla_dalpha_m2_per_deg,
        d_side_force_N_per_deg=q_ref * body_force.dcla_dalpha_m2_per_deg,
        d_pitch_up_moment_Nm_per_deg=q_ref * d_pitch_up_moment_area_m3_per_deg,
        d_yaw_right_moment_Nm_per_deg=q_ref * d_yaw_right_moment_area_m3_per_deg,
        stable_pitch = pitch_restoring,
        stable_yaw = yaw_restoring,
    )
end

function _body_flap_stability_rows(
    config::PrelimConfig,
    baseline::BodyFlapBaselineCondition,
    aero_case::CalibratedAeroCase,
    targets,
    flap_density_proxy_kg_m2::Float64;
    flap_area_fractions_of_stowed_ref::AbstractVector{<:Real},
    assumptions_grid::AbstractVector{BodyFlapStabilityAssumptions},
)
    rows = NamedTuple[]
    for assumptions in assumptions_grid
        baseline_row = _baseline_stowed_body_stability_row(config, baseline, aero_case, targets, assumptions)
        x_cg_m = baseline_row.x_cg_m
        for flap_fraction in Float64.(collect(flap_area_fractions_of_stowed_ref))
            flap_mass_proxy_kg = flap_fraction * shield_stowed_reference_area_m2(config) * flap_density_proxy_kg_m2
            yaw_ctrl = evaluate_shield_body_flap_command(
                config,
                baseline;
                flap_area_fraction_of_stowed_ref=flap_fraction,
                deflection_deg=assumptions.control_deflection_deg,
                mode="yaw_right",
                cg_x_m=x_cg_m,
            )
            pitch_ctrl = evaluate_shield_body_flap_command(
                config,
                baseline;
                flap_area_fraction_of_stowed_ref=flap_fraction,
                deflection_deg=assumptions.control_deflection_deg,
                mode="pitch_up",
                cg_x_m=x_cg_m,
            )
            derivative_df = _body_flap_derivative_rows(config, baseline, flap_fraction; cg_x_m=x_cg_m)
            yaw_deriv = filter(row -> row.mode == "yaw_right", derivative_df)[1, :]
            pitch_deriv = filter(row -> row.mode == "pitch_up", derivative_df)[1, :]
            α_trim_cap_deg = abs(baseline_row.d_pitch_up_moment_Nm_per_deg) > 1e-12 ?
                abs(pitch_ctrl.delta_pitch_up_moment_Nm / baseline_row.d_pitch_up_moment_Nm_per_deg) : Inf
            β_trim_cap_deg = abs(baseline_row.d_yaw_right_moment_Nm_per_deg) > 1e-12 ?
                abs(yaw_ctrl.delta_yaw_right_moment_Nm / baseline_row.d_yaw_right_moment_Nm_per_deg) : Inf
            push!(rows, (
                baseline_row...,
                flap_area_fraction_of_stowed_ref=flap_fraction,
                flap_mass_proxy_kg=flap_mass_proxy_kg,
                dCY_dδ_per_deg=yaw_deriv.d_side_area_m2_per_deg / baseline_row.stowed_reference_area_m2,
                dCn_dδ_per_deg=yaw_deriv.d_yaw_right_moment_area_m3_per_deg / (baseline_row.stowed_reference_area_m2 * baseline_row.reference_length_m),
                dCZ_dδ_per_deg=pitch_deriv.d_vertical_area_m2_per_deg / baseline_row.stowed_reference_area_m2,
                dCm_dδ_per_deg=pitch_deriv.d_pitch_up_moment_area_m3_per_deg / (baseline_row.stowed_reference_area_m2 * baseline_row.reference_length_m),
                d_side_force_N_per_deg_ctrl=yaw_deriv.d_side_force_N_per_deg,
                d_yaw_right_moment_Nm_per_deg_ctrl=yaw_deriv.d_yaw_right_moment_Nm_per_deg,
                d_vertical_force_N_per_deg_ctrl=pitch_deriv.d_vertical_force_N_per_deg,
                d_pitch_up_moment_Nm_per_deg_ctrl=pitch_deriv.d_pitch_up_moment_Nm_per_deg,
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
        end
    end
    return DataFrame(rows)
end

function _write_shield_body_flap_stability_note(
    config::PrelimConfig,
    baseline::BodyFlapBaselineCondition,
    stability_df::DataFrame,
)
    outpath = joinpath(config.output_root, "shield_body_flap_stability_note.md")
    best_rows = filter(row -> row.trim_plausible, stability_df)
    nominal_rows = filter(
        row -> isapprox(row.static_margin_fraction_of_diameter, 0.10; atol=1e-9) &&
            isapprox(row.flap_area_fraction_of_stowed_ref, 0.05; atol=1e-9),
        stability_df,
    )
    nominal = nrow(nominal_rows) > 0 ? nominal_rows[1, :] : stability_df[1, :]
    open(outpath, "w") do io
        println(io, "# SHIELD Stowed-Body Flap Stability Screening")
        println(io)
        println(io, "This note adds a first-cut baseline moment and stability model on top of the stowed-body flap force/moment study. It is still a surrogate, but it is the first pass that asks whether the neutral body-flap architecture can plausibly sit on top of a restoring passive SHIELD body.")
        println(io)
        println(io, "## Representative Condition")
        println(io, @sprintf("- Altitude %.2f km, Mach %.2f, q = %.2f Pa, velocity %.2f m/s", baseline.representative_altitude_km, baseline.representative_mach, baseline.representative_q_pa, baseline.representative_velocity_mps))
        println(io)
        println(io, "## Stability Assumptions")
        println(io, "- Body flaps are assumed flush in the neutral state, so the passive stowed-body force model still defines the baseline.")
        println(io, "- CG location is assumed at a fixed fraction of the stowed cone body length.")
        println(io, "- Static stability is parameterized by an assumed axial static margin, expressed as a fraction of body diameter.")
        println(io, "- Because the surrogate body-force sign convention is inherited from the existing lift model, stability is judged here by **restoring opposition**: the incremental pitch/yaw moment must oppose the corresponding force perturbation.")
        println(io)
        println(io, "## Nominal Result")
        println(io, @sprintf("- Nominal assumption: CG at %.0f%% of body length, static margin %.0f%% of diameter", 100.0 * nominal.cg_axial_fraction_of_body_length, 100.0 * nominal.static_margin_fraction_of_diameter))
        println(io, @sprintf("- Baseline derivatives: `CZ_alpha = %.3f / rad`, `Cma = %.3f / rad`, `CY_beta = %.3f / rad`, `Cnb = %.3f / rad`", nominal.CZ_alpha_per_rad, nominal.Cma_per_rad, nominal.CY_beta_per_rad, nominal.Cnb_per_rad))
        println(io, @sprintf("- Recovered static margins from the derivative ratios: pitch %.3f diameters, yaw %.3f diameters", nominal.recovered_pitch_static_margin, nominal.recovered_yaw_static_margin))
        println(io, @sprintf("- For %.1f%% total flap area and %.0f deg control deflection: pitch trim cap ≈ %.2f deg, yaw trim cap ≈ %.2f deg", 100.0 * nominal.flap_area_fraction_of_stowed_ref, nominal.control_deflection_deg, nominal.pitch_trim_cap_deg, nominal.yaw_trim_cap_deg))
        println(io, @sprintf("- Corresponding control moments: pitch %.3f N·m, yaw %.3f N·m", nominal.pitch_control_moment_Nm, nominal.yaw_control_moment_Nm))
        println(io)
        println(io, "## Interpretation")
        if nrow(best_rows) > 0
            println(io, "- At least one screened combination is trim-plausible under the current surrogate assumptions.")
        else
            println(io, "- None of the screened combinations closes trim under the current surrogate assumptions.")
        end
        println(io, "- This is still not a 6DOF or aeroelastic proof. It only checks whether restoring baseline derivatives and control derivatives are of compatible sign and magnitude.")
        println(io, "- If the proposal continues, this is the right gate before new trajectory reruns: restoring baseline + enough flap moment to trim a few degrees of perturbation.")
        println(io)
        println(io, "## Output Files")
        println(io, "- `shield_body_flap_stability_screen.csv`")
        println(io, "- `shield_body_flap_stability_note.md`")
    end
    return outpath
end

function _write_shield_body_flap_note(
    config::PrelimConfig,
    baseline::BodyFlapBaselineCondition,
    targets,
    flap_rows::DataFrame,
    derivative_df::DataFrame,
    flap_density_proxy_kg_m2::Float64,
    csv_path::String,
    entry_mass_kg::Float64,
)
    outpath = joinpath(config.output_root, "shield_body_flap_note.md")
    neutral_rows = filter(row -> row.mode == "common_drag" && isapprox(row.deflection_deg, 10.0; atol=1e-9), flap_rows)
    yaw_rows = filter(row -> row.mode == "yaw_right" && isapprox(row.deflection_deg, 10.0; atol=1e-9), flap_rows)
    pitch_rows = filter(row -> row.mode == "pitch_up" && isapprox(row.deflection_deg, 10.0; atol=1e-9), flap_rows)
    open(outpath, "w") do io
        println(io, "# SHIELD Stowed-Body Flap Force And Moment Screening")
        println(io)
        println(io, "This study freezes the published-size SHIELD passive baseline, keeps skirt deployment subsonic, and asks whether stowed-body flaps are a more credible main guidance channel than skirt-deployment timing alone.")
        println(io)
        println(io, "## Frozen Passive Baseline")
        println(io, @sprintf("- Entry mass: %.1f kg", entry_mass_kg))
        println(io, @sprintf("- Stowed β_high = %.2f kg/m²", targets.beta_high_kg_m2))
        println(io, @sprintf("- Deployed β_low = %.2f kg/m²", targets.beta_low_kg_m2))
        println(io, @sprintf("- Body-only subsonic onset: %.2f km", baseline.subsonic_onset_km))
        println(io, @sprintf("- Representative pre-deployment body-flap sample: altitude %.2f km, Mach %.2f, q = %.2f Pa, velocity %.2f m/s", baseline.representative_altitude_km, baseline.representative_mach, baseline.representative_q_pa, baseline.representative_velocity_mps))
        println(io)
        println(io, "## Four-Flap Body Architecture Assumptions")
        println(io, "- Four equally spaced body flaps located around the stowed-body shoulder: top, right, bottom, and left.")
        println(io, "- Total controllable flap area is defined as a fraction of the stowed SHIELD frontal reference area, not the deployed skirt area.")
        println(io, @sprintf("- Mass proxy uses a SHIELD-derived stowed-area areal density: %.2f kg/m²", flap_density_proxy_kg_m2))
        println(io, "- Neutral deflection is zero additional control force or moment; this is a stowed-body guidance extension, not a modification to the passive deployed skirt.")
        println(io)
        println(io, "## What Was Evaluated")
        println(io, "- `common_drag`: all four body flaps deflect together")
        println(io, "- `yaw_right`: right flap out and left flap in, producing +y side force and yaw-right moment")
        println(io, "- `pitch_up`: top flap out and bottom flap in, producing +z force and pitch-up moment")
        println(io, "- Single-flap unit derivatives were computed to inspect sign and channel availability")
        println(io)
        println(io, "## Key Findings")
        if nrow(neutral_rows) > 0
            row = neutral_rows[1, :]
            println(io, @sprintf("- Common-mode 10 deg body-flap deflection at %.1f%% total flap area changes drag by %.3f N while keeping net side force %.3f N, pitch-up moment %.3f N·m, and yaw-right moment %.3f N·m balanced by symmetry.", 100.0 * row.flap_area_fraction_of_stowed_ref, row.delta_drag_force_N, row.delta_side_force_N, row.delta_pitch_up_moment_Nm, row.delta_yaw_right_moment_Nm))
        end
        if nrow(yaw_rows) > 0
            best = yaw_rows[argmax(abs.(yaw_rows.delta_yaw_right_moment_Nm)), :]
            println(io, @sprintf("- Yaw mode is available: at %.1f%% total flap area and 10 deg deflection, the model gives side force %.3f N and yaw-right moment %.3f N·m with drag coupling %.3f N.", 100.0 * best.flap_area_fraction_of_stowed_ref, best.delta_side_force_N, best.delta_yaw_right_moment_Nm, best.delta_drag_force_N))
        end
        if nrow(pitch_rows) > 0
            best = pitch_rows[argmax(abs.(pitch_rows.delta_pitch_up_moment_Nm)), :]
            println(io, @sprintf("- Pitch mode is available: at %.1f%% total flap area and 10 deg deflection, the model gives vertical force %.3f N and pitch-up moment %.3f N·m with drag coupling %.3f N.", 100.0 * best.flap_area_fraction_of_stowed_ref, best.delta_vertical_force_N, best.delta_pitch_up_moment_Nm, best.delta_drag_force_N))
        end
        no_roll = all(abs.(derivative_df.d_roll_moment_area_m3_per_deg) .< 1e-10)
        println(io, no_roll ?
            "- This simplified body-flap model does not create an independent roll channel, so the first control question is still pitch/yaw trim rather than full 3-axis control." :
            "- The derivative table shows some roll coupling, but it is not the dominant first-order channel in this simplified model.")
        println(io, "- The key architectural difference from the rim-flap study is timing: these body flaps are available before the skirt deploys, where trajectory sensitivity is expected to be much larger than the 6 km -> 0 km subsonic deployment band.")
        println(io, "- Static stability is still unresolved: this study only screens incremental control force and moment, not whether the stowed SHIELD body plus flaps remains trim-stable in flight.")
        println(io)
        println(io, "## Recommendation")
        println(io, "- If the proposal needs meaningful guidance authority, stowed-body flaps are a more plausible main control channel than skirt-deployment timing alone. The next gate is still force/moment/stability closure before any guided trajectory rerun.")
        println(io)
        println(io, "## Output Files")
        println(io, "- `$(basename(csv_path))`")
        println(io, "- `shield_body_flap_derivatives.csv`")
        println(io, "- `shield_body_flap_note.md`")
    end
    return outpath
end

function run_shield_body_flap_study(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter;
    entry_mass_kg::Real=120.0,
    flap_area_fractions_of_stowed_ref::AbstractVector{<:Real}=[0.01, 0.02, 0.05],
    deflection_grid_deg::AbstractVector{<:Real}=[5.0, 10.0, 15.0],
)
    mkpath(config.output_root)
    baseline_bundle = _representative_shield_body_flap_baseline_condition(config, adapter; entry_mass_kg=entry_mass_kg)
    baseline = baseline_bundle.condition
    targets = baseline_bundle.targets
    flap_density_proxy_kg_m2 = shield_body_flap_areal_density_proxy_kg_m2(config, entry_mass_kg)

    rows = NamedTuple[]
    for flap_fraction in Float64.(collect(flap_area_fractions_of_stowed_ref))
        flap_mass_proxy_kg = flap_fraction * shield_stowed_reference_area_m2(config) * flap_density_proxy_kg_m2
        for deflection_deg in Float64.(collect(deflection_grid_deg))
            for mode in ("common_drag", "yaw_right", "pitch_up")
                result = evaluate_shield_body_flap_command(
                    config,
                    baseline;
                    flap_area_fraction_of_stowed_ref=flap_fraction,
                    deflection_deg=deflection_deg,
                    mode=mode,
                )
                push!(rows, (
                    flap_mass_proxy_kg=flap_mass_proxy_kg,
                    result...,
                ))
            end
        end
    end
    flap_df = DataFrame(rows)
    derivatives_df = vcat([
        _body_flap_derivative_rows(config, baseline, flap_fraction)
        for flap_fraction in Float64.(collect(flap_area_fractions_of_stowed_ref))
    ]...)

    csv_path = joinpath(config.output_root, "shield_body_flap_force_moment_table.csv")
    deriv_path = joinpath(config.output_root, "shield_body_flap_derivatives.csv")
    CSV.write(csv_path, flap_df)
    CSV.write(deriv_path, derivatives_df)
    note_path = _write_shield_body_flap_note(config, baseline, targets, flap_df, derivatives_df, flap_density_proxy_kg_m2, csv_path, Float64(entry_mass_kg))

    return (
        baseline_condition=baseline,
        targets=targets,
        flap_density_proxy_kg_m2=flap_density_proxy_kg_m2,
        flap_df=flap_df,
        derivatives_df=derivatives_df,
        csv_path=csv_path,
        derivatives_path=deriv_path,
        note_path=note_path,
    )
end

function run_shield_body_flap_stability_study(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter;
    entry_mass_kg::Real=120.0,
    flap_area_fractions_of_stowed_ref::AbstractVector{<:Real}=[0.01, 0.02, 0.05],
    cg_axial_fraction_of_body_length::Real=0.60,
    static_margin_fractions_of_diameter::AbstractVector{<:Real}=[0.05, 0.10, 0.15],
    trim_check_angle_deg::Real=5.0,
    control_deflection_deg::Real=15.0,
)
    mkpath(config.output_root)
    baseline_bundle = _representative_shield_body_flap_baseline_condition(config, adapter; entry_mass_kg=entry_mass_kg)
    baseline = baseline_bundle.condition
    targets = baseline_bundle.targets
    aero_case = baseline_bundle.aero_case
    flap_density_proxy_kg_m2 = shield_body_flap_areal_density_proxy_kg_m2(config, entry_mass_kg)
    assumptions_grid = [
        BodyFlapStabilityAssumptions(
            Float64(cg_axial_fraction_of_body_length),
            Float64(sm),
            Float64(trim_check_angle_deg),
            Float64(control_deflection_deg),
        ) for sm in static_margin_fractions_of_diameter
    ]
    stability_df = _body_flap_stability_rows(
        config,
        baseline,
        aero_case,
        targets,
        flap_density_proxy_kg_m2;
        flap_area_fractions_of_stowed_ref=flap_area_fractions_of_stowed_ref,
        assumptions_grid=assumptions_grid,
    )
    csv_path = joinpath(config.output_root, "shield_body_flap_stability_screen.csv")
    CSV.write(csv_path, stability_df)
    note_path = _write_shield_body_flap_stability_note(config, baseline, stability_df)
    return (
        baseline_condition=baseline,
        targets=targets,
        aero_case=aero_case,
        flap_density_proxy_kg_m2=flap_density_proxy_kg_m2,
        stability_df=stability_df,
        csv_path=csv_path,
        note_path=note_path,
    )
end
