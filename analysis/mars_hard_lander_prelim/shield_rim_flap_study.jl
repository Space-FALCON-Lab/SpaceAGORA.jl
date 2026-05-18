struct RimFlapComponent
    name::String
    azimuth_deg::Float64
    area_m2::Float64
    position_m::SVector{3, Float64}
    radial_dir::SVector{3, Float64}
end

struct RimFlapBaselineCondition
    target_h_switch_km::Float64
    deployment_band_min_km::Float64
    deployment_band_max_km::Float64
    representative_altitude_km::Float64
    representative_mach::Float64
    representative_q_pa::Float64
    representative_velocity_mps::Float64
end

@inline function shield_surface_areal_density_proxy_kg_m2(config::PrelimConfig, entry_mass_kg::Real)
    area_m2 = deployed_drag_skirt_equivalent_area(config)
    area_m2 > 0.0 || throw(ArgumentError("Deployed drag-skirt equivalent area must be positive for SHIELD rim-flap study."))
    return Float64(entry_mass_kg) / area_m2
end

@inline function _rim_flap_control_outputs(Δforce_area::SVector{3, Float64}, Δmoment_area::SVector{3, Float64})
    return (
        delta_drag_area_m2=Δforce_area[1],
        delta_side_area_m2=Δforce_area[2],
        delta_vertical_area_m2=Δforce_area[3],
        delta_roll_moment_area_m3=Δmoment_area[1],
        delta_pitch_up_moment_area_m3=-Δmoment_area[2],
        delta_yaw_right_moment_area_m3=Δmoment_area[3],
    )
end

function build_shield_rim_flaps(
    config::PrelimConfig;
    flap_area_fraction_of_skirt::Real,
    flap_count::Integer=4,
)
    flap_count == 4 || throw(ArgumentError("This first-cut SHIELD rim-flap study is defined for 4 flaps."))
    total_flap_area_m2 = Float64(flap_area_fraction_of_skirt) * deployed_drag_skirt_equivalent_area(config)
    each_flap_area_m2 = total_flap_area_m2 / flap_count
    x_arm_m = config.geometry.drag_skirt_height_m
    deployed_radius_m = 0.5 * config.geometry.deployed_drag_surface_diameter_m

    layout = [
        ("top", 0.0),
        ("right", 90.0),
        ("bottom", 180.0),
        ("left", 270.0),
    ]
    flaps = RimFlapComponent[]
    for (name, azimuth_deg) in layout
        ϕ = deg2rad(azimuth_deg)
        radial_dir = SVector{3, Float64}(0.0, sin(ϕ), cos(ϕ))
        position_m = SVector{3, Float64}(x_arm_m, deployed_radius_m * sin(ϕ), deployed_radius_m * cos(ϕ))
        push!(flaps, RimFlapComponent(name, azimuth_deg, each_flap_area_m2, position_m, radial_dir))
    end
    return (
        flap_count=flap_count,
        flap_area_total_m2=total_flap_area_m2,
        flap_area_each_m2=each_flap_area_m2,
        flaps=flaps,
    )
end

function _rim_flap_force_moment_delta(
    config::PrelimConfig,
    flap::RimFlapComponent,
    deflection_deg::Real,
)
    plate = FlatPlateReferenceBody(
        flap.area_m2,
        config.geometry.panel_skin_friction_coefficient,
        config.geometry.panel_zero_lift_drag_coefficient,
    )
    α0_rad = deg2rad(config.geometry.panel_alpha_deg)
    α1_rad = deg2rad(config.geometry.panel_alpha_deg - Float64(deflection_deg))
    CL0, CD0 = _flat_plate_cl_cd(α0_rad, plate)
    CL1, CD1 = _flat_plate_cl_cd(α1_rad, plate)
    area = flap.area_m2

    force0_area = SVector{3, Float64}(CD0 * area, CL0 * area * flap.radial_dir[2], CL0 * area * flap.radial_dir[3])
    force1_area = SVector{3, Float64}(CD1 * area, CL1 * area * flap.radial_dir[2], CL1 * area * flap.radial_dir[3])
    Δforce_area = force1_area - force0_area
    Δmoment_area = cross(flap.position_m, Δforce_area)
    return Δforce_area, Δmoment_area
end

function _rim_flap_command_deflections(mode::String, deflection_deg::Real)
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
    throw(ArgumentError("Unsupported rim-flap command mode $(repr(mode))."))
end

function evaluate_shield_rim_flap_command(
    config::PrelimConfig,
    baseline::RimFlapBaselineCondition;
    flap_area_fraction_of_skirt::Real,
    deflection_deg::Real,
    mode::String,
)
    layout = build_shield_rim_flaps(config; flap_area_fraction_of_skirt=flap_area_fraction_of_skirt)
    deflections = _rim_flap_command_deflections(mode, deflection_deg)

    Δforce_area = SVector{3, Float64}(0.0, 0.0, 0.0)
    Δmoment_area = SVector{3, Float64}(0.0, 0.0, 0.0)
    for flap in layout.flaps
        flap_Δforce_area, flap_Δmoment_area = _rim_flap_force_moment_delta(config, flap, deflections[flap.name])
        Δforce_area += flap_Δforce_area
        Δmoment_area += flap_Δmoment_area
    end

    outputs = _rim_flap_control_outputs(Δforce_area, Δmoment_area)
    q_ref = baseline.representative_q_pa
    return (
        mode=mode,
        flap_area_fraction_of_skirt=Float64(flap_area_fraction_of_skirt),
        flap_area_total_m2=layout.flap_area_total_m2,
        flap_area_each_m2=layout.flap_area_each_m2,
        deflection_deg=Float64(deflection_deg),
        representative_altitude_km=baseline.representative_altitude_km,
        representative_mach=baseline.representative_mach,
        representative_q_pa=q_ref,
        outputs...,
        delta_drag_force_N=q_ref * outputs.delta_drag_area_m2,
        delta_side_force_N=q_ref * outputs.delta_side_area_m2,
        delta_vertical_force_N=q_ref * outputs.delta_vertical_area_m2,
        delta_roll_moment_Nm=q_ref * outputs.delta_roll_moment_area_m3,
        delta_pitch_up_moment_Nm=q_ref * outputs.delta_pitch_up_moment_area_m3,
        delta_yaw_right_moment_Nm=q_ref * outputs.delta_yaw_right_moment_area_m3,
    )
end

function _representative_shield_baseline_condition(
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
    switch_max_km = 0.5 * floor(body_result.trajectory.altitude_km[subsonic_idx] / 0.5)
    h_switch_grid_m = collect((switch_max_km * 1e3):-500.0:0.0)
    target_cfg = default_config(
        output_root=config.output_root,
        atmosphere_mode=config.atmosphere_mode,
        use_gram_surrogate=config.use_gram_surrogate,
        generate_plots=false,
        enable_crossrange_proxy=false,
        enable_target_range_guidance=true,
        enable_monte_carlo=false,
        beta_high_targets=[targets.beta_high_kg_m2],
        beta_ratios=[targets.beta_ratio],
        h_jettison_grid_m=[0.0],
        h_switch_grid_m=h_switch_grid_m,
        representative_beta_high=targets.beta_high_kg_m2,
        representative_beta_ratio=targets.beta_ratio,
        representative_h_jettison_m=0.0,
        representative_switches_m=collect((switch_max_km * 1e3):-1000.0:0.0),
        geometry=config.geometry,
    )
    target_guidance = find_switch_altitude_for_target_range(
        target_cfg,
        adapter,
        aero_case,
        0.0,
        begin
            always_deployed = solve_entry_trajectory(
                config,
                adapter,
                aero_case,
                ControlPolicy("fixed_deployed", :fixed_deployed, NaN, 0.0);
                save_trajectory=false,
            )
            body_only = solve_entry_trajectory(
                config,
                adapter,
                aero_case,
                ControlPolicy("body_only", :body_only, NaN, NaN);
                save_trajectory=false,
            )
            body_only.summary.impact_downrange_km + config.target_range_fraction * (always_deployed.summary.impact_downrange_km - body_only.summary.impact_downrange_km)
        end;
        h_switch_grid_m=h_switch_grid_m,
        save_trajectory=true,
    )
    deployed_segment = filter(row -> row.deployed, target_guidance.trajectory_result.trajectory)
    nrow(deployed_segment) > 0 || throw(ArgumentError("Targeted SHIELD trajectory never enters a deployed segment."))
    q_idx = argmax(deployed_segment.q_pa)
    sample = deployed_segment[q_idx, :]
    return (
        condition=RimFlapBaselineCondition(
            target_guidance.h_switch_m / 1e3,
            minimum(h_switch_grid_m) / 1e3,
            maximum(h_switch_grid_m) / 1e3,
            sample.altitude_km,
            sample.mach,
            sample.q_pa,
            sample.velocity_mps,
        ),
        target_guidance=target_guidance,
        aero_case=aero_case,
        targets=targets,
    )
end

function _rim_flap_derivative_rows(config::PrelimConfig, baseline::RimFlapBaselineCondition, flap_area_fraction::Float64)
    rows = NamedTuple[]
    for mode in ("single_top", "single_right", "single_bottom", "single_left", "common_drag", "yaw_right", "pitch_up")
        positive = evaluate_shield_rim_flap_command(config, baseline; flap_area_fraction_of_skirt=flap_area_fraction, deflection_deg=1.0, mode=mode)
        negative = evaluate_shield_rim_flap_command(config, baseline; flap_area_fraction_of_skirt=flap_area_fraction, deflection_deg=-1.0, mode=mode)
        push!(rows, (
            mode=mode,
            flap_area_fraction_of_skirt=flap_area_fraction,
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

function _write_shield_rim_flap_note(
    config::PrelimConfig,
    baseline::RimFlapBaselineCondition,
    targets,
    flap_rows::DataFrame,
    derivative_df::DataFrame,
    flap_density_proxy_kg_m2::Float64,
    csv_path::String,
    entry_mass_kg::Float64,
)
    outpath = joinpath(config.output_root, "shield_rim_flap_note.md")
    neutral_rows = filter(row -> row.mode == "common_drag" && isapprox(row.deflection_deg, 10.0; atol=1e-9), flap_rows)
    yaw_rows = filter(row -> row.mode == "yaw_right" && isapprox(row.deflection_deg, 10.0; atol=1e-9), flap_rows)
    pitch_rows = filter(row -> row.mode == "pitch_up" && isapprox(row.deflection_deg, 10.0; atol=1e-9), flap_rows)
    open(outpath, "w") do io
        println(io, "# SHIELD Four Rim-Flap Force And Moment Screening")
        println(io)
        println(io, "This study freezes the published-size SHIELD surrogate baseline and treats skirt deployment timing only as a bounded sensitivity. The main question here is whether four small rim flaps on the deployed skirt create a usable neutral state and useful incremental control channels before any new trajectory rerun.")
        println(io)
        println(io, "## Frozen Passive Baseline")
        println(io, @sprintf("- Entry mass: %.1f kg", entry_mass_kg))
        println(io, @sprintf("- Stowed β_high = %.2f kg/m²", targets.beta_high_kg_m2))
        println(io, @sprintf("- Deployed β_low = %.2f kg/m²", targets.beta_low_kg_m2))
        println(io, @sprintf("- Subsonic skirt-deployment band from body-only trajectory: %.1f to %.1f km", baseline.deployment_band_max_km, baseline.deployment_band_min_km))
        println(io, @sprintf("- Representative targeted skirt deployment switch: %.2f km", baseline.target_h_switch_km))
        println(io, @sprintf("- Representative deployed-condition sample used for force/moment scaling: altitude %.2f km, Mach %.2f, q = %.2f Pa, velocity %.2f m/s", baseline.representative_altitude_km, baseline.representative_mach, baseline.representative_q_pa, baseline.representative_velocity_mps))
        println(io)
        println(io, "## Four-Flap Architecture Assumptions")
        println(io, "- Four equally spaced rim flaps located at the deployed skirt tip: top, right, bottom, and left.")
        println(io, "- Flaps are treated as controllable subsets of the deployed skirt area, so neutral deflection is the deployed-skirt baseline rather than an added drag surface.")
        println(io, @sprintf("- Mass proxy uses the same SHIELD deployed-surface areal density: %.2f kg/m²", flap_density_proxy_kg_m2))
        println(io, "- Incremental flap control mass is therefore a conservative area-scaled proxy, not a hinge/actuator design result.")
        println(io)
        println(io, "## What Was Evaluated")
        println(io, "- `common_drag`: all four flaps deflect together")
        println(io, "- `yaw_right`: right flap deflects out and left flap in, producing +y side force and yaw-right moment")
        println(io, "- `pitch_up`: top flap deflects out and bottom flap in, producing +z force and pitch-up moment")
        println(io, "- Single-flap unit derivatives were also computed to inspect sign and channel availability")
        println(io)
        println(io, "## Key Findings")
        if nrow(neutral_rows) > 0
            row = neutral_rows[1, :]
            println(io, @sprintf("- Common-mode 10 deg deflection at %.1f%% flap area changes drag by %.3f N while keeping net side force %.3f N, pitch-up moment %.3f N·m, and yaw-right moment %.3f N·m essentially balanced by symmetry.", 100.0 * row.flap_area_fraction_of_skirt, row.delta_drag_force_N, row.delta_side_force_N, row.delta_pitch_up_moment_Nm, row.delta_yaw_right_moment_Nm))
        end
        if nrow(yaw_rows) > 0
            best = yaw_rows[argmax(abs.(yaw_rows.delta_yaw_right_moment_Nm)), :]
            println(io, @sprintf("- Yaw mode is available: at %.1f%% flap area and 10 deg deflection, the model gives side force %.3f N and yaw-right moment %.3f N·m with drag coupling %.3f N.", 100.0 * best.flap_area_fraction_of_skirt, best.delta_side_force_N, best.delta_yaw_right_moment_Nm, best.delta_drag_force_N))
        end
        if nrow(pitch_rows) > 0
            best = pitch_rows[argmax(abs.(pitch_rows.delta_pitch_up_moment_Nm)), :]
            println(io, @sprintf("- Pitch mode is available: at %.1f%% flap area and 10 deg deflection, the model gives vertical force %.3f N and pitch-up moment %.3f N·m with drag coupling %.3f N.", 100.0 * best.flap_area_fraction_of_skirt, best.delta_vertical_force_N, best.delta_pitch_up_moment_Nm, best.delta_drag_force_N))
        end
        no_roll = all(abs.(derivative_df.d_roll_moment_area_m3_per_deg) .< 1e-10)
        println(io, no_roll ?
            "- This simplified radial-force model does not produce an independent roll channel, which means the first useful control question is pitch/yaw trim rather than full 3-axis flap control." :
            "- The derivative table shows some roll coupling, but it is not the dominant first-order channel in this simplified model.")
        println(io, "- The neutral flap state is geometrically plausible in the surrogate because symmetric zero-deflection control adds no net side force or pitch/yaw moment by construction.")
        println(io, "- Actual static stability is still unresolved because the baseline deployed-skirt aerodynamic moments are not modeled here. This study only screens incremental flap control authority.")
        println(io)
        println(io, "## Recommendation")
        println(io, "- Do not rerun guided trajectories yet. The next gate is to add at least a first-cut deployed-skirt moment/stability model and then test whether these incremental flap channels sit on top of a plausible trimmed deployed baseline.")
        println(io)
        println(io, "## Output Files")
        println(io, "- `$(basename(csv_path))`")
        println(io, "- `shield_rim_flap_derivatives.csv`")
        println(io, "- `shield_rim_flap_note.md`")
    end
    return outpath
end

function run_shield_rim_flap_study(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter;
    entry_mass_kg::Real=120.0,
    flap_area_fractions_of_skirt::AbstractVector{<:Real}=[0.01, 0.02, 0.05],
    deflection_grid_deg::AbstractVector{<:Real}=[5.0, 10.0, 15.0],
)
    mkpath(config.output_root)
    baseline_bundle = _representative_shield_baseline_condition(config, adapter; entry_mass_kg=entry_mass_kg)
    baseline = baseline_bundle.condition
    targets = baseline_bundle.targets
    flap_density_proxy_kg_m2 = shield_surface_areal_density_proxy_kg_m2(config, entry_mass_kg)

    rows = NamedTuple[]
    for flap_fraction in Float64.(collect(flap_area_fractions_of_skirt))
        flap_mass_proxy_kg = flap_fraction * deployed_drag_skirt_equivalent_area(config) * flap_density_proxy_kg_m2
        for deflection_deg in Float64.(collect(deflection_grid_deg))
            for mode in ("common_drag", "yaw_right", "pitch_up")
                result = evaluate_shield_rim_flap_command(
                    config,
                    baseline;
                    flap_area_fraction_of_skirt=flap_fraction,
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
        _rim_flap_derivative_rows(config, baseline, flap_fraction)
        for flap_fraction in Float64.(collect(flap_area_fractions_of_skirt))
    ]...)

    csv_path = joinpath(config.output_root, "shield_rim_flap_force_moment_table.csv")
    deriv_path = joinpath(config.output_root, "shield_rim_flap_derivatives.csv")
    CSV.write(csv_path, flap_df)
    CSV.write(deriv_path, derivatives_df)
    note_path = _write_shield_rim_flap_note(config, baseline, targets, flap_df, derivatives_df, flap_density_proxy_kg_m2, csv_path, Float64(entry_mass_kg))

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
