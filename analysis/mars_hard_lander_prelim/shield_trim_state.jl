function shield_system_areal_density_proxy_kg_m2(config::PrelimConfig, entry_mass_kg::Real)
    area_m2 = deployed_drag_skirt_equivalent_area(config)
    area_m2 > 0.0 || throw(ArgumentError("Deployed drag-skirt equivalent area must be positive for SHIELD trim study."))
    return Float64(entry_mass_kg) / area_m2
end

function _mid_corridor_target_range(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter,
    aero_case::CalibratedAeroCase,
)
    ranges_km = Float64[]
    for h_switch_m in config.h_switch_grid_m
        result = solve_entry_trajectory(
            config,
            adapter,
            aero_case,
            ControlPolicy("trim_target_bracket", :bang_bang, h_switch_m, config.representative_h_jettison_m);
            save_trajectory=false,
        )
        push!(ranges_km, result.summary.impact_downrange_km)
    end
    return minimum(ranges_km) + config.target_range_fraction * (maximum(ranges_km) - minimum(ranges_km))
end

function _write_shield_state3_trim_note(
    config::PrelimConfig,
    rows::DataFrame,
    trim_area_fraction::Float64,
    trim_area_total_m2::Float64,
    trim_mass_kg::Float64,
    density_proxy_kg_m2::Float64,
)
    outpath = joinpath(config.output_root, "shield_state3_trim_note.md")
    open(outpath, "w") do io
        println(io, "# SHIELD State 3 Trim Study")
        println(io)
        println(io, "State 3 in this note means: deployed SHIELD drag-skirt surrogate plus added trim surfaces that activate only after skirt deployment.")
        println(io)
        println(io, "## Trim-Surface Assumptions")
        println(io, @sprintf("- Trim-area fraction of deployed skirt surrogate: %.1f%%", 100.0 * trim_area_fraction))
        println(io, @sprintf("- Total trim area: %.4f m² (%.4f m² per side)", trim_area_total_m2, trim_area_total_m2 / 2.0))
        println(io, @sprintf("- SHIELD areal-density proxy used for trim mass: %.2f kg/m²", density_proxy_kg_m2))
        println(io, @sprintf("- Added trim-system mass proxy: %.2f kg", trim_mass_kg))
        println(io)
        println(io, "## Interpretation")
        println(io, "- `baseline_shield_targeted` is the published-size SHIELD surrogate with no trim hardware added.")
        println(io, "- `state2_augmented_skirt_only` carries the added trim mass but keeps the trim surfaces inactive. That isolates the mass penalty of the added hardware.")
        println(io, "- `state3_symmetric_trim` activates the trim surfaces symmetrically after deployment. That isolates the extra drag effect without intentional crossrange steering.")
        println(io, "- `state3_asymmetric_trim` activates the same trim system asymmetrically in 3D, so any crossrange comes from the trim surfaces rather than from differentially steering the whole drag skirt.")
        println(io)
        println(io, "## Results")
        for row in eachrow(rows)
            if row.case == "state3_asymmetric_trim"
                println(io, @sprintf(
                    "- `%s`: impact downrange %.2f km, crossrange %.2f km, impact velocity %.2f m/s, peak side accel %.4f m/s².",
                    row.case,
                    row.impact_downrange_km,
                    row.impact_crossrange_km,
                    row.impact_velocity_mps,
                    row.peak_side_accel_mps2,
                ))
            else
                println(io, @sprintf(
                    "- `%s`: impact downrange %.2f km, impact velocity %.2f m/s, peak total decel %.2f g.",
                    row.case,
                    row.impact_downrange_km,
                    row.impact_velocity_mps,
                    row.peak_total_decel_earth_g,
                ))
            end
        end
    end
    return outpath
end

function run_shield_state3_trim_study(
    config::PrelimConfig,
    adapter::MarsAtmosphereAdapter;
    base_entry_mass_kg::Real=120.0,
    trim_panel_area_fraction_of_deployed::Real=0.02,
    favored_side_sign::Real=1.0,
    differential_fraction::Real=config.crossrange_target_differential_fraction,
    cant_deg::Real=config.crossrange_target_panel_cant_deg,
)
    mkpath(config.output_root)
    cfg_base = with_trim_fraction(config, 0.0)
    base_targets = derive_fixed_mass_beta_targets_from_deployed_geometry(
        cfg_base,
        adapter;
        entry_mass_kg=base_entry_mass_kg,
    )
    base_aero_case = calibrate_aero_case(cfg_base, adapter, base_targets.beta_high_kg_m2, base_targets.beta_ratio)

    density_proxy_kg_m2 = shield_system_areal_density_proxy_kg_m2(cfg_base, base_entry_mass_kg)
    trim_area_total_m2 = Float64(trim_panel_area_fraction_of_deployed) * base_aero_case.deployed_state.panel_area_total_m2
    trim_mass_kg = trim_area_total_m2 * density_proxy_kg_m2
    total_entry_mass_kg = Float64(base_entry_mass_kg) + trim_mass_kg

    cfg_trim = with_trim_fraction(config, trim_panel_area_fraction_of_deployed)
    trim_targets = derive_fixed_mass_beta_targets_from_deployed_geometry(
        cfg_trim,
        adapter;
        entry_mass_kg=total_entry_mass_kg,
    )
    trim_aero_case = calibrate_aero_case(cfg_trim, adapter, trim_targets.beta_high_kg_m2, trim_targets.beta_ratio)

    target_range_km = _mid_corridor_target_range(cfg_base, adapter, base_aero_case)
    target_result = find_switch_altitude_for_target_range(
        cfg_base,
        adapter,
        base_aero_case,
        cfg_base.representative_h_jettison_m,
        target_range_km;
        save_trajectory=true,
    )
    h_switch_m = target_result.h_switch_m

    baseline_policy = ControlPolicy("baseline_shield_targeted", :bang_bang, h_switch_m, cfg_base.representative_h_jettison_m)
    baseline_result = solve_entry_trajectory(cfg_base, adapter, base_aero_case, baseline_policy; save_trajectory=true)

    state2_policy = ControlPolicy("state2_augmented_skirt_only", :bang_bang, h_switch_m, cfg_trim.representative_h_jettison_m)
    state2_result = solve_entry_trajectory(cfg_trim, adapter, trim_aero_case, state2_policy; save_trajectory=true)

    state3_policy = ControlPolicy("state3_symmetric_trim", :bang_bang_trim, h_switch_m, cfg_trim.representative_h_jettison_m)
    state3_result = solve_entry_trajectory(cfg_trim, adapter, trim_aero_case, state3_policy; save_trajectory=true)

    symmetric_3d = solve_entry_trajectory_3d(
        cfg_trim,
        adapter,
        trim_aero_case,
        ControlPolicy("state3_trim_symmetric_3d", :bang_bang_trim, h_switch_m, cfg_trim.representative_h_jettison_m);
        save_trajectory=false,
        use_winds=false,
        panel_command_override=neutral_panel_command(cfg_trim),
    )
    asymmetric_cmd = make_differential_panel_command(
        cfg_trim;
        favored_side_sign=favored_side_sign,
        differential_fraction=differential_fraction,
        cant_deg=cant_deg,
    )
    asymmetric_3d = solve_entry_trajectory_3d(
        cfg_trim,
        adapter,
        trim_aero_case,
        ControlPolicy("state3_asymmetric_trim", :bang_bang_trim, h_switch_m, cfg_trim.representative_h_jettison_m);
        save_trajectory=false,
        use_winds=false,
        panel_command_override=asymmetric_cmd,
    )

    rows = DataFrame([
        (
            case="baseline_shield_targeted",
            entry_mass_kg=base_aero_case.mass_kg,
            trim_area_total_m2=0.0,
            trim_mass_proxy_kg=0.0,
            h_switch_km=h_switch_m / 1e3,
            impact_downrange_km=baseline_result.summary.impact_downrange_km,
            impact_crossrange_km=0.0,
            impact_velocity_mps=baseline_result.summary.impact_velocity_mps,
            peak_total_decel_earth_g=baseline_result.summary.peak_total_decel_earth_g,
            peak_side_accel_mps2=0.0,
        ),
        (
            case="state2_augmented_skirt_only",
            entry_mass_kg=trim_aero_case.mass_kg,
            trim_area_total_m2=trim_aero_case.deployed_state.trim_panel_area_total_m2,
            trim_mass_proxy_kg=trim_mass_kg,
            h_switch_km=h_switch_m / 1e3,
            impact_downrange_km=state2_result.summary.impact_downrange_km,
            impact_crossrange_km=0.0,
            impact_velocity_mps=state2_result.summary.impact_velocity_mps,
            peak_total_decel_earth_g=state2_result.summary.peak_total_decel_earth_g,
            peak_side_accel_mps2=0.0,
        ),
        (
            case="state3_symmetric_trim",
            entry_mass_kg=trim_aero_case.mass_kg,
            trim_area_total_m2=trim_aero_case.deployed_state.trim_panel_area_total_m2,
            trim_mass_proxy_kg=trim_mass_kg,
            h_switch_km=h_switch_m / 1e3,
            impact_downrange_km=state3_result.summary.impact_downrange_km,
            impact_crossrange_km=0.0,
            impact_velocity_mps=state3_result.summary.impact_velocity_mps,
            peak_total_decel_earth_g=state3_result.summary.peak_total_decel_earth_g,
            peak_side_accel_mps2=0.0,
        ),
        (
            case="state3_asymmetric_trim",
            entry_mass_kg=trim_aero_case.mass_kg,
            trim_area_total_m2=trim_aero_case.deployed_state.trim_panel_area_total_m2,
            trim_mass_proxy_kg=trim_mass_kg,
            h_switch_km=h_switch_m / 1e3,
            impact_downrange_km=asymmetric_3d.summary.impact_downrange_km,
            impact_crossrange_km=asymmetric_3d.summary.impact_crossrange_km,
            impact_velocity_mps=asymmetric_3d.summary.impact_velocity_mps,
            peak_total_decel_earth_g=asymmetric_3d.summary.peak_total_decel_earth_g,
            peak_side_accel_mps2=asymmetric_3d.summary.peak_side_accel_mps2,
        ),
    ])

    csv_path = joinpath(config.output_root, "shield_state3_trim_summary.csv")
    CSV.write(csv_path, rows)
    note_path = _write_shield_state3_trim_note(
        config,
        rows,
        Float64(trim_panel_area_fraction_of_deployed),
        trim_area_total_m2,
        trim_mass_kg,
        density_proxy_kg_m2,
    )

    return (
        rows=rows,
        csv_path=csv_path,
        note_path=note_path,
        base_entry_mass_kg=Float64(base_entry_mass_kg),
        total_entry_mass_kg=total_entry_mass_kg,
        trim_panel_area_fraction=Float64(trim_panel_area_fraction_of_deployed),
        trim_panel_area_total_m2=trim_area_total_m2,
        trim_panel_mass_proxy_kg=trim_mass_kg,
        shield_areal_density_proxy_kg_m2=density_proxy_kg_m2,
        h_switch_m=h_switch_m,
        baseline_target_range_km=target_range_km,
        symmetric_3d=symmetric_3d,
        asymmetric_3d=asymmetric_3d,
    )
end
