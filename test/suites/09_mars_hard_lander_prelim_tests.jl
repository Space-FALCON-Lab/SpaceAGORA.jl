if !isdefined(@__MODULE__, :MarsHardLanderPrelim)
    include(joinpath(REPO_ROOT, "analysis", "mars_hard_lander_prelim", "MarsHardLanderPrelim.jl"))
end

using .MarsHardLanderPrelim

@testset "Mars Hard-Lander Preliminary Analysis" begin
    @testset "Grant-Braun Sphere-Cone CN/CA" begin
        area = π * 1.25^2
        regular_body = MarsHardLanderPrelim.SphereConeReferenceBody(
            area,
            0.30,
            1.25,
            70.0,
            :regular_newtonian,
        )
        regular_alias_body = MarsHardLanderPrelim.SphereConeReferenceBody(
            area,
            0.30,
            1.25,
            70.0,
            :newtonian,
        )
        modified_body = MarsHardLanderPrelim.SphereConeReferenceBody(
            area,
            0.30,
            1.25,
            70.0,
            :modified_newtonian,
        )
        alpha_rad = deg2rad(10.0)
        gamma = 1.29

        CN_regular, CA_regular = MarsHardLanderPrelim._sphere_cone_cn_ca(
            alpha_rad,
            regular_body,
            5.0,
            gamma,
        )
        CN_alias, CA_alias = MarsHardLanderPrelim._sphere_cone_cn_ca(
            alpha_rad,
            regular_alias_body,
            20.0,
            gamma,
        )
        @test CN_regular ≈ 0.03987396864691749 rtol=1e-12
        @test CA_regular ≈ 1.7170715355728225 rtol=1e-12
        @test CN_alias ≈ CN_regular rtol=1e-14
        @test CA_alias ≈ CA_regular rtol=1e-14

        CN_regular_zero, CA_regular_zero = MarsHardLanderPrelim._sphere_cone_cn_ca(
            0.0,
            regular_body,
            5.0,
            gamma,
        )
        delta_rad = deg2rad(regular_body.cone_half_angle_deg)
        bluntness = regular_body.nose_radius / regular_body.base_radius
        @test CN_regular_zero ≈ 0.0 atol=1e-13
        @test CA_regular_zero ≈ 2.0 * sin(delta_rad)^2 + bluntness^2 * cos(delta_rad)^4 rtol=1e-14

        cp_max_m5 = MarsHardLanderPrelim._modified_newtonian_cp_max(5.0, gamma)
        @test cp_max_m5 ≈ 1.8443468538239867 rtol=1e-13
        CN_modified, CA_modified = MarsHardLanderPrelim._sphere_cone_cn_ca(
            alpha_rad,
            modified_body,
            5.0,
            gamma,
        )
        @test CN_modified / CN_regular ≈ cp_max_m5 / 2.0 rtol=1e-13
        @test CA_modified / CA_regular ≈ cp_max_m5 / 2.0 rtol=1e-13

        CL, CD = MarsHardLanderPrelim._sphere_cone_cl_cd(
            alpha_rad,
            modified_body,
            5.0,
            gamma,
        )
        @test CL ≈ CN_modified * cos(alpha_rad) - CA_modified * sin(alpha_rad) rtol=1e-14
        @test CD ≈ CA_modified * cos(alpha_rad) + CN_modified * sin(alpha_rad) rtol=1e-14

        cp_at_zero = MarsHardLanderPrelim._modified_newtonian_cp_max(0.0, gamma)
        cp_below_one = MarsHardLanderPrelim._modified_newtonian_cp_max(1.0 - 1e-8, gamma)
        cp_above_one = MarsHardLanderPrelim._modified_newtonian_cp_max(1.0 + 1e-8, gamma)
        @test cp_at_zero == 1.0
        @test cp_below_one ≈ cp_above_one rtol=1e-7

        sharp_cone = MarsHardLanderPrelim.SphereConeReferenceBody(
            π,
            0.0,
            1.0,
            45.0,
            :regular_newtonian,
        )
        CN_zero, CA_zero = MarsHardLanderPrelim._sphere_cone_cn_ca(0.0, sharp_cone, 10.0, 1.4)
        @test CN_zero ≈ 0.0 atol=1e-13
        @test CA_zero ≈ 1.0 atol=1e-14

        unsupported_body = MarsHardLanderPrelim.SphereConeReferenceBody(
            area,
            0.30,
            1.25,
            70.0,
            :unsupported,
        )
        @test_throws ArgumentError MarsHardLanderPrelim._sphere_cone_cn_ca(
            alpha_rad,
            unsupported_body,
            5.0,
            gamma,
        )
        @test_throws DomainError MarsHardLanderPrelim._modified_newtonian_cp_max(-1.0, gamma)
        @test_throws DomainError MarsHardLanderPrelim._modified_newtonian_cp_max(5.0, 1.0)
        CN_shadowed, CA_shadowed = MarsHardLanderPrelim._sphere_cone_cn_ca(
            deg2rad(71.0),
            regular_body,
            5.0,
            gamma,
        )
        @test isfinite(CN_shadowed)
        @test isfinite(CA_shadowed)
    end

    @testset "Atmosphere And Aero Calibration" begin
        cfg = default_config(
            atmosphere_mode=:exponential,
            generate_plots=false,
            enable_crossrange_proxy=false,
            enable_target_range_guidance=false,
            enable_monte_carlo=false,
            beta_high_targets=[150.0],
            beta_ratios=[2.0],
            h_jettison_grid_m=[20e3],
            h_switch_grid_m=[120e3, 60e3, 0.0],
            representative_beta_high=150.0,
            representative_beta_ratio=2.0,
            representative_h_jettison_m=20e3,
            representative_switches_m=[120e3, 60e3, 0.0],
        )
        adapter = build_atmosphere_adapter(cfg)
        rho_low, T_low = MarsHardLanderPrelim.density_temperature(adapter, 50e3, 0.0, 0.0)
        rho_high, T_high = MarsHardLanderPrelim.density_temperature(adapter, 100e3, 0.0, 0.0)
        @test rho_low > rho_high >= 0.0
        @test isfinite(T_low)
        @test isfinite(T_high)

        gram_cfg = default_config(
            atmosphere_mode=:gram,
            generate_plots=false,
            enable_crossrange_proxy=false,
            beta_high_targets=[150.0],
            beta_ratios=[2.0],
            h_jettison_grid_m=[20e3],
            h_switch_grid_m=[120e3, 60e3, 0.0],
            representative_beta_high=150.0,
            representative_beta_ratio=2.0,
            representative_h_jettison_m=20e3,
            representative_switches_m=[120e3, 60e3, 0.0],
        )
        gram_adapter = build_atmosphere_adapter(gram_cfg)
        rho_gram, T_gram = MarsHardLanderPrelim.density_temperature(gram_adapter, gram_cfg.h0_m, 0.0, 0.0)
        @test rho_gram >= 0.0
        @test isfinite(T_gram)

        aero_case = calibrate_aero_case(cfg, adapter, 150.0, 2.0)
        body_loads = MarsHardLanderPrelim.aerodynamic_loads(aero_case, false, aero_case.reference_mach, cfg.planet.γ)
        deployed_loads = MarsHardLanderPrelim.aerodynamic_loads(aero_case, true, aero_case.reference_mach, cfg.planet.γ)
        @test isfinite(body_loads.CDA_m2) && body_loads.CDA_m2 > 0.0
        @test isfinite(body_loads.CLA_m2)
        @test isfinite(deployed_loads.CDA_m2) && deployed_loads.CDA_m2 >= body_loads.CDA_m2
        @test isfinite(deployed_loads.CLA_m2)
        @test abs(aero_case.achieved_beta_high - 150.0) / 150.0 < 1e-10
        @test abs(aero_case.achieved_beta_ratio - 2.0) / 2.0 < 1e-10
    end

    @testset "Solver And Output Smoke" begin
        output_root = mktempdir()
        cfg = default_config(
            output_root=output_root,
            atmosphere_mode=:exponential,
            generate_plots=true,
            enable_crossrange_proxy=false,
            enable_target_range_guidance=true,
            enable_monte_carlo=true,
            monte_carlo_samples=12,
            beta_high_targets=[150.0],
            beta_ratios=[1.0, 2.0],
            h_jettison_grid_m=[0.0, 20e3],
            h_switch_grid_m=[120e3, 80e3, 40e3, 0.0],
            representative_beta_high=150.0,
            representative_beta_ratio=2.0,
            representative_h_jettison_m=20e3,
            representative_switches_m=[120e3, 80e3, 40e3, 0.0],
        )
        res = run_prelim(cfg)

        @test isfile(joinpath(output_root, "prelim_summary.csv"))
        @test isfile(joinpath(output_root, "authority_summary.csv"))
        @test isfile(joinpath(output_root, "authority_trade_metrics.csv"))
        @test isfile(joinpath(output_root, "local_control_effectiveness.csv"))
        @test isfile(joinpath(output_root, "alpha_body_sensitivity.csv"))
        @test isfile(joinpath(output_root, "crossrange_sensitivity.csv"))
        @test isfile(joinpath(output_root, "crossrange_sensitivity.md"))
        @test isfile(joinpath(output_root, "README_prelim.md"))
        @test isfile(joinpath(output_root, "proposal_trade_table.md"))
        @test isfile(joinpath(output_root, "target_guidance_summary.csv"))
        @test isfile(joinpath(output_root, "monte_carlo_summary.csv"))
        @test isfile(joinpath(output_root, "monte_carlo_samples.csv"))
        @test isfile(joinpath(output_root, "plots", "altitude_vs_downrange.png"))
        @test isfile(joinpath(output_root, "plots", "velocity_vs_altitude.png"))
        @test isfile(joinpath(output_root, "plots", "dynamic_pressure_vs_altitude.png"))
        @test isfile(joinpath(output_root, "plots", "drag_accel_vs_altitude.png"))
        @test isfile(joinpath(output_root, "plots", "downrange_authority_heatmap.png"))
        @test isfile(joinpath(output_root, "plots", "terminal_velocity_and_peak_g_vs_beta_ratio.png"))
        @test isfile(joinpath(output_root, "plots", "local_control_effectiveness.png"))
        @test isfile(joinpath(output_root, "plots", "alpha_body_sensitivity.png"))
        @test isfile(joinpath(output_root, "plots", "crossrange_sensitivity.png"))
        @test isfile(joinpath(output_root, "plots", "landing_ellipse_proxy.png"))
        @test isfile(joinpath(output_root, "plots", "landing_ellipse_centered.png"))

        @test size(res.summary_df, 1) > 0
        @test size(res.authority_df, 1) > 0
        @test size(res.trade_metrics_df, 1) > 0
        @test size(res.local_effectiveness_df, 1) > 0
        @test size(res.crossrange_sensitivity_df, 1) == 24
        @test size(res.alpha_sensitivity_df, 1) == 5
        @test size(res.target_guidance_df, 1) == 1
        @test size(res.monte_carlo_samples_df, 1) == 36
        @test size(res.monte_carlo_summary_df, 1) == 3
        @test all(res.summary_df.peak_dynamic_pressure_pa .>= 0.0)
        @test all(res.summary_df.peak_drag_accel_mps2 .>= 0.0)
        @test all(res.summary_df.integration_stopped_at_surface)

        body_only_rows = filter(row -> row.policy == "body_only" && MarsHardLanderPrelim._approx_eq(row.target_beta_high_kg_m2, 150.0), res.summary_df)
        fixed_rows = filter(
            row -> row.policy == "fixed_deployed" &&
                MarsHardLanderPrelim._approx_eq(row.target_beta_high_kg_m2, 150.0) &&
                MarsHardLanderPrelim._approx_eq(row.target_beta_ratio, 2.0) &&
                MarsHardLanderPrelim._approx_eq(row.h_jettison_km, 20.0),
            res.summary_df,
        )
        @test size(body_only_rows, 1) == 2
        @test size(fixed_rows, 1) == 1
        @test body_only_rows.impact_downrange_km[1] > fixed_rows.impact_downrange_km[1]

        ratio_one_authority = filter(row -> MarsHardLanderPrelim._approx_eq(row.target_beta_ratio, 1.0), res.authority_df)
        @test all(abs.(ratio_one_authority.downrange_authority_km) .< 1e-6)

        h = res.local_effectiveness_df.h_switch_km
        s = abs.(res.local_effectiveness_df.sensitivity_km_per_km)
        integral_km = sum(0.5 * (s[i] + s[i + 1]) * abs(h[i + 1] - h[i]) for i in 1:(length(h) - 1))
        span_km = maximum(res.local_effectiveness_df.impact_downrange_km) - minimum(res.local_effectiveness_df.impact_downrange_km)
        @test span_km == 0.0 || 0.5 <= integral_km / span_km <= 5.0
        @test abs(res.target_guidance_df.range_error_km[1]) <= cfg.target_range_tolerance_km
        @test "guided_targeted_optimistic" in res.monte_carlo_summary_df.policy
        @test "body_only_monte_carlo" in res.monte_carlo_summary_df.policy
        @test "guided_targeted_crossrange_bilateral" in res.monte_carlo_summary_df.policy
        @test all(res.monte_carlo_summary_df.ellipse_major_axis_km .>= res.monte_carlo_summary_df.ellipse_minor_axis_km)
        crossrange_guided_rows = filter(row -> row.policy == "guided_targeted_crossrange_bilateral", res.monte_carlo_samples_df)
        @test size(crossrange_guided_rows, 1) == 12
        @test all(side -> side in ("left", "right"), collect(skipmissing(crossrange_guided_rows.chosen_side)))
        @test MarsHardLanderPrelim._approx_eq(res.alpha_sensitivity_df.alpha_body_deg[1], 0.0)
        @test any(.!isnan.(res.trade_metrics_df.authority_per_added_area_km_per_m2))
        nominal_trade = filter(row -> MarsHardLanderPrelim._approx_eq(row.target_beta_high_kg_m2, 150.0) && MarsHardLanderPrelim._approx_eq(row.target_beta_ratio, 2.0), res.trade_metrics_df)
        @test size(nominal_trade, 1) > 0
        @test all(nominal_trade.authority_loss_vs_latest_jettison_km .>= -1e-9)
        large_cant0 = filter(
            row -> row.differential_deployment_label == "large" &&
                MarsHardLanderPrelim._approx_eq(row.panel_cant_deg, 0.0) &&
                MarsHardLanderPrelim._approx_eq(row.alpha_body_deg, 0.0),
            res.crossrange_sensitivity_df,
        )
        large_cant45 = filter(
            row -> row.differential_deployment_label == "large" &&
                MarsHardLanderPrelim._approx_eq(row.panel_cant_deg, 45.0) &&
                MarsHardLanderPrelim._approx_eq(row.alpha_body_deg, 0.0),
            res.crossrange_sensitivity_df,
        )
        @test size(large_cant0, 1) == 1
        @test size(large_cant45, 1) == 1
        @test abs(large_cant0.impact_crossrange_km[1]) <= abs(large_cant45.impact_crossrange_km[1]) + 1e-9
        @test all(res.crossrange_sensitivity_df.peak_side_accel_mps2 .>= -1e-12)
    end

    @testset "SHIELD Rim-Flap Static Screening" begin
        output_root = mktempdir()
        cfg = default_config(
            output_root=output_root,
            atmosphere_mode=:exponential,
            generate_plots=false,
            enable_crossrange_proxy=false,
            enable_target_range_guidance=false,
            enable_monte_carlo=false,
            geometry=shield_published_surrogate_geometry(),
        )
        adapter = build_atmosphere_adapter(cfg)
        res = run_shield_rim_flap_study(
            cfg,
            adapter;
            entry_mass_kg=120.0,
            flap_area_fractions_of_skirt=[0.02],
            deflection_grid_deg=[10.0],
        )

        @test isfile(res.csv_path)
        @test isfile(res.derivatives_path)
        @test isfile(res.note_path)
        @test res.baseline_condition.deployment_band_max_km > res.baseline_condition.deployment_band_min_km >= 0.0
        @test res.baseline_condition.target_h_switch_km >= res.baseline_condition.deployment_band_min_km
        @test res.baseline_condition.target_h_switch_km <= res.baseline_condition.deployment_band_max_km

        common = MarsHardLanderPrelim.evaluate_shield_rim_flap_command(
            cfg,
            res.baseline_condition;
            flap_area_fraction_of_skirt=0.02,
            deflection_deg=10.0,
            mode="common_drag",
        )
        @test abs(common.delta_side_area_m2) < 1e-10
        @test abs(common.delta_vertical_area_m2) < 1e-10
        @test abs(common.delta_pitch_up_moment_area_m3) < 1e-10
        @test abs(common.delta_yaw_right_moment_area_m3) < 1e-10

        yaw = MarsHardLanderPrelim.evaluate_shield_rim_flap_command(
            cfg,
            res.baseline_condition;
            flap_area_fraction_of_skirt=0.02,
            deflection_deg=10.0,
            mode="yaw_right",
        )
        pitch = MarsHardLanderPrelim.evaluate_shield_rim_flap_command(
            cfg,
            res.baseline_condition;
            flap_area_fraction_of_skirt=0.02,
            deflection_deg=10.0,
            mode="pitch_up",
        )
        @test yaw.delta_side_area_m2 > 0.0
        @test yaw.delta_yaw_right_moment_area_m3 > 0.0
        @test abs(yaw.delta_vertical_area_m2) < 1e-10
        @test pitch.delta_vertical_area_m2 > 0.0
        @test pitch.delta_pitch_up_moment_area_m3 > 0.0
        @test abs(pitch.delta_side_area_m2) < 1e-10
    end

    @testset "SHIELD Stowed-Body Flap Static Screening" begin
        output_root = mktempdir()
        cfg = default_config(
            output_root=output_root,
            atmosphere_mode=:exponential,
            generate_plots=false,
            enable_crossrange_proxy=false,
            enable_target_range_guidance=false,
            enable_monte_carlo=false,
            geometry=shield_published_surrogate_geometry(),
        )
        adapter = build_atmosphere_adapter(cfg)
        res = run_shield_body_flap_study(
            cfg,
            adapter;
            entry_mass_kg=120.0,
            flap_area_fractions_of_stowed_ref=[0.02],
            deflection_grid_deg=[10.0],
        )

        @test isfile(res.csv_path)
        @test isfile(res.derivatives_path)
        @test isfile(res.note_path)
        @test res.baseline_condition.subsonic_onset_km > 0.0
        @test res.baseline_condition.representative_altitude_km >= res.baseline_condition.subsonic_onset_km
        @test res.baseline_condition.representative_mach > 1.0

        common = MarsHardLanderPrelim.evaluate_shield_body_flap_command(
            cfg,
            res.baseline_condition;
            flap_area_fraction_of_stowed_ref=0.02,
            deflection_deg=10.0,
            mode="common_drag",
        )
        @test abs(common.delta_side_area_m2) < 1e-10
        @test abs(common.delta_vertical_area_m2) < 1e-10
        @test abs(common.delta_pitch_up_moment_area_m3) < 1e-10
        @test abs(common.delta_yaw_right_moment_area_m3) < 1e-10

        yaw = MarsHardLanderPrelim.evaluate_shield_body_flap_command(
            cfg,
            res.baseline_condition;
            flap_area_fraction_of_stowed_ref=0.02,
            deflection_deg=10.0,
            mode="yaw_right",
        )
        pitch = MarsHardLanderPrelim.evaluate_shield_body_flap_command(
            cfg,
            res.baseline_condition;
            flap_area_fraction_of_stowed_ref=0.02,
            deflection_deg=10.0,
            mode="pitch_up",
        )
        @test yaw.delta_side_area_m2 > 0.0
        @test yaw.delta_yaw_right_moment_area_m3 > 0.0
        @test abs(yaw.delta_vertical_area_m2) < 1e-10
        @test pitch.delta_vertical_area_m2 > 0.0
        @test pitch.delta_pitch_up_moment_area_m3 > 0.0
        @test abs(pitch.delta_side_area_m2) < 1e-10
    end

    @testset "SHIELD Stowed-Body Flap Stability Screening" begin
        output_root = mktempdir()
        cfg = default_config(
            output_root=output_root,
            atmosphere_mode=:exponential,
            generate_plots=false,
            enable_crossrange_proxy=false,
            enable_target_range_guidance=false,
            enable_monte_carlo=false,
            geometry=shield_published_surrogate_geometry(),
        )
        adapter = build_atmosphere_adapter(cfg)
        res = run_shield_body_flap_stability_study(
            cfg,
            adapter;
            entry_mass_kg=120.0,
            flap_area_fractions_of_stowed_ref=[0.05],
            cg_axial_fraction_of_body_length=0.60,
            static_margin_fractions_of_diameter=[0.10],
            trim_check_angle_deg=5.0,
            control_deflection_deg=15.0,
        )

        @test isfile(res.csv_path)
        @test isfile(res.note_path)
        @test size(res.stability_df, 1) == 1
        row = res.stability_df[1, :]
        @test row.stable_pitch
        @test row.stable_yaw
        @test row.recovered_pitch_static_margin > 0.0
        @test row.recovered_yaw_static_margin > 0.0
        @test row.dCm_dδ_per_deg > 0.0
        @test row.dCn_dδ_per_deg > 0.0
        @test row.pitch_trim_cap_deg > 0.0
        @test row.yaw_trim_cap_deg > 0.0
    end

    @testset "SHIELD Shoulder-Strake Stability Screening" begin
        output_root = mktempdir()
        cfg = default_config(
            output_root=output_root,
            atmosphere_mode=:exponential,
            generate_plots=false,
            enable_crossrange_proxy=false,
            enable_target_range_guidance=false,
            enable_monte_carlo=false,
            geometry=shield_published_surrogate_geometry(),
        )
        adapter = build_atmosphere_adapter(cfg)
        res = run_shield_shoulder_strake_stability_study(
            cfg,
            adapter;
            entry_mass_kg=120.0,
            strake_area_fractions_of_stowed_ref=[0.05],
            cg_axial_fractions_of_body_length=[0.50],
            static_margin_fractions_of_diameter=[0.10],
            trim_check_angle_deg=5.0,
            control_deflection_deg=15.0,
            cl_alpha_per_rad=3.5,
            cd0=0.05,
            induced_drag_factor=0.15,
            cl_max=1.2,
        )

        @test isfile(res.csv_path)
        @test isfile(res.derivatives_path)
        @test isfile(res.note_path)
        @test size(res.screen_df, 1) == 1
        row = res.screen_df[1, :]
        @test row.stable_pitch
        @test row.stable_yaw
        @test row.recovered_pitch_static_margin > 0.0
        @test row.yaw_control_moment_Nm > 0.0
        @test row.pitch_control_moment_Nm > 0.0
        @test row.pitch_trim_cap_deg > 0.0
        @test row.yaw_trim_cap_deg > 0.0
    end

    @testset "SHIELD Shoulder-Strake Guided Trajectory Screening" begin
        output_root = mktempdir()
        cfg = default_config(
            output_root=output_root,
            atmosphere_mode=:exponential,
            generate_plots=false,
            enable_crossrange_proxy=false,
            enable_target_range_guidance=false,
            enable_monte_carlo=false,
            geometry=shield_published_surrogate_geometry(),
        )
        adapter = build_atmosphere_adapter(cfg)
        res = run_shield_shoulder_strake_trajectory_study(
            cfg,
            adapter;
            entry_mass_kg=120.0,
            cg_axial_fraction_of_body_length=0.60,
            static_margin_fraction_of_diameter=0.10,
            strake_area_fraction_of_stowed_ref=0.05,
            control_deflection_deg=15.0,
            trim_check_angle_deg=5.0,
            cl_alpha_per_rad=3.5,
            cd0=0.05,
            induced_drag_factor=0.15,
            cl_max=1.2,
            pitch_command_grid=[-1.0, 0.0, 1.0],
            yaw_command_grid=[-1.0, 1.0],
        )

        @test isfile(res.summary_path)
        @test isfile(res.pitch_path)
        @test isfile(res.yaw_path)
        @test isfile(res.note_path)
        @test size(res.pitch_df, 1) == 3
        @test size(res.yaw_df, 1) == 2
        @test maximum(res.pitch_df.impact_downrange_km) > minimum(res.pitch_df.impact_downrange_km)
        @test all(abs.(res.yaw_df.impact_crossrange_km) .> 0.0)
    end
end
