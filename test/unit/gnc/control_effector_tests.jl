@testset "Control Burn Energy Sign (End-to-End)" begin
    sc0 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    args0 = build_config(
        spacecraft=sc0,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(),
        control_rates=Float64[],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df0 = run_case_silent(args0)
    eps0 = specific_energy(df0, EARTH.μ)
    Δeps0 = last(eps0) - first(eps0)
    @test abs(Δeps0) < 2e3

    scp = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    thr_pro = TimedTangentialThrusterModel(800.0, +1.0, 100.0, 101.0)
    argsp = build_config(
        spacecraft=scp,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thr_pro,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    dfp = run_case_silent(argsp)
    epsp = specific_energy(dfp, EARTH.μ)
    Δepsp = last(epsp) - first(epsp)
    @test Δepsp > 5e3

    scr = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    thr_retro = TimedTangentialThrusterModel(800.0, -1.0, 100.0, 101.0)
    argsr = build_config(
        spacecraft=scr,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=250.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thr_retro,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    dfr = run_case_silent(argsr)
    epsr = specific_energy(dfr, EARTH.μ)
    Δepsr = last(epsr) - first(epsr)
    @test Δepsr < -5e3
end

@testset "Control Mass Flow Coupling (End-to-End)" begin
    sc_no_control = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
    args_no_control = build_config(
        spacecraft=sc_no_control,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=180.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(),
        control_rates=Float64[],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df_no_control = run_case_silent(args_no_control)
    @test "sc1_mass" in names(df_no_control)
    mass_no_control = Vector{Float64}(df_no_control.sc1_mass)
    @test maximum(abs.(mass_no_control .- first(mass_no_control))) < 1e-9

    function run_burn_case(isp::Float64)
        sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=60.0)
        thruster = make_base_thruster_model(
            thrust=600.0,
            direction=0.0,
            Δv=0.0,
            start_burn_time=0.0,
            stop_burn_time=80.0,
            Isp=isp
        )
        args = build_config(
            spacecraft=sc,
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=180.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            control_effectors=(thruster,),
            control_rates=[1.0],
            keplerian=true,
            tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
        )
        return run_case_silent(args)
    end

    df_low_isp = run_burn_case(200.0)
    df_high_isp = run_burn_case(400.0)

    mass_low = Vector{Float64}(df_low_isp.sc1_mass)
    mass_high = Vector{Float64}(df_high_isp.sc1_mass)
    Δm_low = first(mass_low) - last(mass_low)
    Δm_high = first(mass_high) - last(mass_high)

    @test Δm_low > 5.0
    @test Δm_high > 2.0
    @test all(diff(mass_low) .<= 1e-7)
    @test all(diff(mass_high) .<= 1e-7)
    @test isapprox(Δm_low / Δm_high, 2.0; atol=0.08, rtol=0.0)
end

@testset "Control Callback Multi-Spacecraft Mapping" begin
    sc1 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
    sc2 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=500e3, ν_deg=120.0)
    shared_thruster = BaseThrusterModel(
        thrust=[800.0, 800.0],
        direction=[0.0, π],
        Δv=[20.0, -20.0],
        start_burn_time=[-1.0, -1.0],
        stop_burn_time=[-1.0, -1.0],
        Isp=[300.0, 300.0]
    )
    args = build_config_multi(
        spacecraft=[sc1, sc2],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=1000.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(shared_thruster,),
        control_rates=[1.0],
        keplerian=true,
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=5.0)
    )
    df = run_case_silent(args; isolate_state=false)

    for sat_idx in 1:2
        @test isfinite(shared_thruster.start_burn_time[sat_idx])
        @test isfinite(shared_thruster.stop_burn_time[sat_idx])
        @test shared_thruster.stop_burn_time[sat_idx] > shared_thruster.start_burn_time[sat_idx]
    end

    mass1 = Vector{Float64}(df.sc1_mass)
    mass2 = Vector{Float64}(df.sc2_mass)
    @test first(mass1) - last(mass1) > 0.1
    @test first(mass2) - last(mass2) > 0.1

    eps1 = 0.5 .* (df.sc1_vel_1.^2 .+ df.sc1_vel_2.^2 .+ df.sc1_vel_3.^2) .-
           EARTH.μ ./ sqrt.(df.sc1_pos_1.^2 .+ df.sc1_pos_2.^2 .+ df.sc1_pos_3.^2)
    eps2 = 0.5 .* (df.sc2_vel_1.^2 .+ df.sc2_vel_2.^2 .+ df.sc2_vel_3.^2) .-
           EARTH.μ ./ sqrt.(df.sc2_pos_1.^2 .+ df.sc2_pos_2.^2 .+ df.sc2_pos_3.^2)
    Δeps1 = last(eps1) - first(eps1)
    Δeps2 = last(eps2) - first(eps2)
    @test Δeps1 > 2e3
    @test Δeps2 < -2e3
end
