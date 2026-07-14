@testset "Drag Dissipates Specific Orbital Energy" begin
    sc = make_spacecraft(
        ra_alt_m=220e3,
        rp_alt_m=100e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=180.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )

    df = run_case(args)
    eps = specific_energy(df, EARTH.μ)
    @test last(eps) < first(eps) - 1e5
end

@testset "Two-Body Kepler Invariants Stay Constant" begin
    ra_alt_m = 1_800e3
    rp_alt_m = 400e3
    sc = make_spacecraft(
        ra_alt_m=ra_alt_m,
        rp_alt_m=rp_alt_m,
        i_deg=28.0,
        ω_deg=15.0,
        Ω_deg=25.0,
        ν_deg=40.0
    )
    ra_m = EARTH.Rp_e + ra_alt_m
    rp_m = EARTH.Rp_e + rp_alt_m
    a0 = 0.5 * (ra_m + rp_m)
    period_s = 2pi * sqrt(a0^3 / EARTH.μ)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=2.5 * period_s,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-9,
            abstol_orbit=1e-9,
            dt_max_orbit=10.0
        )
    )

    df = run_case_silent(args)
    eps = specific_energy(df, EARTH.μ)
    hmag = specific_angular_momentum_magnitude(df)
    a_series = -EARTH.μ ./ (2.0 .* eps)

    eps0 = first(eps)
    h0 = first(hmag)
    a_rel_drift = maximum(abs.((a_series .- a_series[1]) ./ a0))
    eps_rel_drift = maximum(abs.((eps .- eps0) ./ eps0))
    h_rel_drift = maximum(abs.((hmag .- h0) ./ h0))

    @test eps_rel_drift < 1e-5
    @test h_rel_drift < 1e-5
    @test a_rel_drift < 1e-5
end

@testset "Gravity Backbone Improves Long-Horizon Two-Body Energy Drift" begin
    ra_alt_m = 1_400e3
    rp_alt_m = 400e3
    sc = make_spacecraft(
        ra_alt_m=ra_alt_m,
        rp_alt_m=rp_alt_m,
        i_deg=28.0,
        ω_deg=15.0,
        Ω_deg=25.0,
        ν_deg=40.0
    )
    a0 = 0.5 * ((EARTH.Rp_e + ra_alt_m) + (EARTH.Rp_e + rp_alt_m))
    period_s = 2pi * sqrt(a0^3 / EARTH.μ)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=15.0 * period_s,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-5,
            abstol_orbit=1e-5,
            dt_max_orbit=120.0
        )
    )

    df_tsit = DataFrame()
    df_backbone = DataFrame()
    tsit_time = @elapsed begin
        df_tsit = withenv("SPACEAGORA_SOLVER_MODE" => "tsit5") do
            run_case_silent(args)
        end
    end
    backbone_time = @elapsed begin
        df_backbone = withenv(
            "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
            "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "20.0"
        ) do
            run_case_silent(args)
        end
    end

    eps_tsit = specific_energy(df_tsit, EARTH.μ)
    eps_backbone = specific_energy(df_backbone, EARTH.μ)
    drift_tsit = maximum(abs.((eps_tsit .- first(eps_tsit)) ./ first(eps_tsit)))
    drift_backbone = maximum(abs.((eps_backbone .- first(eps_backbone)) ./ first(eps_backbone)))

    @info "gravity_backbone_benchmark" case="two_body_long_horizon" tsit5_seconds=tsit_time gravity_backbone_seconds=backbone_time tsit5_energy_drift=drift_tsit gravity_backbone_energy_drift=drift_backbone
    @test drift_backbone < drift_tsit
end

@testset "J2 Secular Rates Match Standard First-Order Drift" begin
    sc = make_spacecraft(
        ra_alt_m=2_000e3,
        rp_alt_m=400e3,
        i_deg=45.0,
        ω_deg=25.0,
        Ω_deg=40.0,
        ν_deg=20.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=86_400.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-9,
            abstol_orbit=1e-9,
            dt_max_orbit=10.0
        )
    )

    df = run_case_silent(args)
    times = Float64.(df.time)
    Ω_series = Vector{Float64}(undef, nrow(df))
    ω_series = Vector{Float64}(undef, nrow(df))

    @inbounds for idx in 1:nrow(df)
        pos = SVector{3, Float64}(
            Float64(df.sc1_pos_1[idx]),
            Float64(df.sc1_pos_2[idx]),
            Float64(df.sc1_pos_3[idx])
        )
        vel = SVector{3, Float64}(
            Float64(df.sc1_vel_1[idx]),
            Float64(df.sc1_vel_2[idx]),
            Float64(df.sc1_vel_3[idx])
        )
        oe = rvtoorbitalelement(pos, vel, EARTH)
        Ω_series[idx] = Float64(oe[4])
        ω_series[idx] = Float64(oe[5])
    end

    oe0 = rvtoorbitalelement(
        SVector{3, Float64}(Float64(df.sc1_pos_1[1]), Float64(df.sc1_pos_2[1]), Float64(df.sc1_pos_3[1])),
        SVector{3, Float64}(Float64(df.sc1_vel_1[1]), Float64(df.sc1_vel_2[1]), Float64(df.sc1_vel_3[1])),
        EARTH
    )
    Ωdot_expected, ωdot_expected = SimulationModel.GravityEffectors.j2_secular_rates(
        Float64(oe0[1]),
        Float64(oe0[2]),
        Float64(oe0[3]),
        EARTH
    )

    Ωdot_measured = linear_regression_slope(times, unwrap_angle_series(Ω_series))
    ωdot_measured = linear_regression_slope(times, unwrap_angle_series(ω_series))

    # Vallado / Montenbruck-Gill first-order secular rates describe the drift trend.
    # Compare slopes with moderate tolerance because the integrated elements are osculating, not mean.
    @test signbit(Ωdot_measured) == signbit(Ωdot_expected)
    @test signbit(ωdot_measured) == signbit(ωdot_expected)
    @test isapprox(Ωdot_measured, Ωdot_expected; rtol=0.10, atol=0.0)
    @test isapprox(ωdot_measured, ωdot_expected; rtol=0.15, atol=0.0)
end

@testset "Gravity Backbone J2 Preserves Secular Drift Trend" begin
    sc = make_spacecraft(
        ra_alt_m=2_000e3,
        rp_alt_m=400e3,
        i_deg=45.0,
        ω_deg=20.0,
        Ω_deg=15.0,
        ν_deg=40.0
    )
    a0 = 0.5 * ((EARTH.Rp_e + 2_000e3) + (EARTH.Rp_e + 400e3))
    period_s = 2pi * sqrt(a0^3 / EARTH.μ)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=6.0 * period_s,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredJ2GravityModel(),),
        keplerian=true,
        ephemerides_model=SimpleEphemeridesModel(),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            dt_max_orbit=30.0
        )
    )

    df = withenv(
        "SPACEAGORA_SOLVER_MODE" => "gravity_backbone_split",
        "SPACEAGORA_GRAVITY_BACKBONE_DT_S" => "5.0"
    ) do
        run_case_silent(args)
    end

    times = Float64.(df.time)
    Ω_series = Vector{Float64}(undef, nrow(df))
    ω_series = Vector{Float64}(undef, nrow(df))
    @inbounds for idx in 1:nrow(df)
        pos = SVector{3, Float64}(
            Float64(df.sc1_pos_1[idx]),
            Float64(df.sc1_pos_2[idx]),
            Float64(df.sc1_pos_3[idx])
        )
        vel = SVector{3, Float64}(
            Float64(df.sc1_vel_1[idx]),
            Float64(df.sc1_vel_2[idx]),
            Float64(df.sc1_vel_3[idx])
        )
        oe = rvtoorbitalelement(pos, vel, EARTH)
        Ω_series[idx] = Float64(oe[4])
        ω_series[idx] = Float64(oe[5])
    end

    oe0 = rvtoorbitalelement(
        SVector{3, Float64}(Float64(df.sc1_pos_1[1]), Float64(df.sc1_pos_2[1]), Float64(df.sc1_pos_3[1])),
        SVector{3, Float64}(Float64(df.sc1_vel_1[1]), Float64(df.sc1_vel_2[1]), Float64(df.sc1_vel_3[1])),
        EARTH
    )
    Ωdot_expected, ωdot_expected = SimulationModel.GravityEffectors.j2_secular_rates(
        Float64(oe0[1]),
        Float64(oe0[2]),
        Float64(oe0[3]),
        EARTH
    )

    Ωdot_measured = linear_regression_slope(times, unwrap_angle_series(Ω_series))
    ωdot_measured = linear_regression_slope(times, unwrap_angle_series(ω_series))
    @test signbit(Ωdot_measured) == signbit(Ωdot_expected)
    @test signbit(ωdot_measured) == signbit(ωdot_expected)
    @test isapprox(Ωdot_measured, Ωdot_expected; rtol=0.20, atol=0.0)
    @test isapprox(ωdot_measured, ωdot_expected; rtol=0.25, atol=0.0)
end

# The AGORA Earth golden regression test now lives at test/golden/agora_earth/
# and runs via test/golden/runtests.jl (see test/TEST_RESTRUCTURE_PLAN.md).

@testset "N-Body Gravity Typed API Smoke" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(), NBodyGravityModel(["Sun"], "Earth", SPICE_PATH)),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df = run_case(args)
    @test nrow(df) > 10
    @test all(isfinite, df.sc1_pos_1)
    @test all(isfinite, df.sc1_vel_1)
end

@testset "Two-Spacecraft Isolation vs Single-Craft Baselines" begin
    sc_a = make_spacecraft(ra_alt_m=520e3, rp_alt_m=430e3, ν_deg=150.0)
    sc_b = make_spacecraft(ra_alt_m=700e3, rp_alt_m=650e3, ν_deg=110.0)

    args_multi = build_config_multi(
        spacecraft=[sc_a, sc_b],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    args_a = build_config(
        spacecraft=make_spacecraft(ra_alt_m=520e3, rp_alt_m=430e3, ν_deg=150.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )
    args_b = build_config(
        spacecraft=make_spacecraft(ra_alt_m=700e3, rp_alt_m=650e3, ν_deg=110.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    df_multi = run_case_silent(args_multi)
    df_a = run_case_silent(args_a)
    df_b = run_case_silent(args_b)

    @test nrow(df_multi) > 10
    sample_idxs = round.(Int, range(1, nrow(df_multi), length=8))
    for idx in sample_idxs
        t = Float64(df_multi.time[idx])

        pa_m = SVector{3, Float64}(Float64(df_multi.sc1_pos_1[idx]), Float64(df_multi.sc1_pos_2[idx]), Float64(df_multi.sc1_pos_3[idx]))
        va_m = SVector{3, Float64}(Float64(df_multi.sc1_vel_1[idx]), Float64(df_multi.sc1_vel_2[idx]), Float64(df_multi.sc1_vel_3[idx]))
        pa_s = SVector{3, Float64}(
            interp_linear(df_a.time, df_a.sc1_pos_1, t),
            interp_linear(df_a.time, df_a.sc1_pos_2, t),
            interp_linear(df_a.time, df_a.sc1_pos_3, t)
        )
        va_s = SVector{3, Float64}(
            interp_linear(df_a.time, df_a.sc1_vel_1, t),
            interp_linear(df_a.time, df_a.sc1_vel_2, t),
            interp_linear(df_a.time, df_a.sc1_vel_3, t)
        )
        # Multi-satellite adaptive stepping can introduce modest trajectory differences vs single-body runs.
        @test norm(pa_m - pa_s) < 200.0
        @test norm(va_m - va_s) < 0.2

        pb_m = SVector{3, Float64}(Float64(df_multi.sc2_pos_1[idx]), Float64(df_multi.sc2_pos_2[idx]), Float64(df_multi.sc2_pos_3[idx]))
        vb_m = SVector{3, Float64}(Float64(df_multi.sc2_vel_1[idx]), Float64(df_multi.sc2_vel_2[idx]), Float64(df_multi.sc2_vel_3[idx]))
        pb_s = SVector{3, Float64}(
            interp_linear(df_b.time, df_b.sc1_pos_1, t),
            interp_linear(df_b.time, df_b.sc1_pos_2, t),
            interp_linear(df_b.time, df_b.sc1_pos_3, t)
        )
        vb_s = SVector{3, Float64}(
            interp_linear(df_b.time, df_b.sc1_vel_1, t),
            interp_linear(df_b.time, df_b.sc1_vel_2, t),
            interp_linear(df_b.time, df_b.sc1_vel_3, t)
        )
        @test norm(pb_m - pb_s) < 200.0
        @test norm(vb_m - vb_s) < 0.2
    end
end

@testset "Single-Link Drag Dissipates Specific Orbital Energy" begin
    sc = make_single_link_spacecraft(
        ra_alt_m=220e3,
        rp_alt_m=100e3,
        i_deg=35.0,
        ω_deg=40.0,
        Ω_deg=10.0,
        ν_deg=180.0
    )
    args = build_config(
        spacecraft=sc,
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=900.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )

    df = run_case(args)
    eps = specific_energy(df, EARTH.μ)
    @test last(eps) < first(eps) - 1e5
end
