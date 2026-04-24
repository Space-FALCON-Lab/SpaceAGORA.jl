@testset "No-GRAM Onboarding Mode" begin
    @testset "Preset builders choose documented fallback models" begin
        env_default = make_no_gram_environment()
        @test env_default.planet isa Earth
        @test env_default.density_model isa NoAtmosphereModel
        @test env_default.ephemerides_model isa SimpleEphemeridesModel
        @test env_default.wind == false

        env_none = make_no_gram_environment(planet=:earth, atmosphere=:none, EI_km=120.0)
        @test env_none.planet isa Earth
        @test env_none.density_model isa NoAtmosphereModel
        @test env_none.ephemerides_model isa SimpleEphemeridesModel
        @test env_none.wind == false

        env_exp = make_no_gram_environment(planet=:mars, atmosphere=:exponential, EI_km=140.0)
        @test env_exp.planet isa Mars
        @test env_exp.density_model isa ExponentialAtmosphereModel
        @test env_exp.ephemerides_model isa SimpleEphemeridesModel
        @test env_exp.density_model.ρ_ref == env_exp.planet.ρ_ref
        @test env_exp.density_model.h_ref == env_exp.planet.h_ref
        @test env_exp.density_model.H == env_exp.planet.H
    end

    @testset "Simple ephemerides avoid kernel-dependent runtime" begin
        lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
            kclear()
        end

        earth_sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=400e3, ν_deg=175.0)
        earth_args = make_example_config(
            planet=Earth(),
            spacecraft=earth_sc,
            mission_time=120.0,
            initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(InverseSquaredJ2GravityModel(),),
            density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            orientation_sim=false,
            keplerian=true,
            EI_km=120.0,
            verbose=false,
            results=false,
            results_directory=joinpath(REPO_ROOT, "output", "test_no_gram_earth")
        )
        @test_nowarn run_simulation(earth_args)
        earth_et0 = ephemerides_time_seconds(earth_args.initial_time, earth_args.environment_model.ephemerides_model)
        earth_lpi0 = planet_frame_lpi(earth_args.environment_model.planet, earth_et0, earth_args.environment_model.ephemerides_model)
        @test !all(iszero, earth_lpi0)

        mars_sc = make_spacecraft(ra_alt_m=250e3, rp_alt_m=150e3, i_deg=30.0, ω_deg=10.0, Ω_deg=20.0, ν_deg=160.0)
        mars_args = make_example_config(
            planet=Mars(),
            spacecraft=mars_sc,
            mission_time=90.0,
            initial_time=InitialTime(year=2024, month=2, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(InverseSquaredGravityModel(),),
            density_model=ExponentialAtmosphereModel(Mars()),
            ephemerides_model=SimpleEphemeridesModel(),
            orientation_sim=false,
            keplerian=true,
            EI_km=140.0,
            verbose=false,
            results=false,
            results_directory=joinpath(REPO_ROOT, "output", "test_no_gram_mars")
        )
        @test_nowarn run_simulation(mars_args)
        mars_et0 = ephemerides_time_seconds(mars_args.initial_time, mars_args.environment_model.ephemerides_model)
        mars_lpi0 = planet_frame_lpi(mars_args.environment_model.planet, mars_et0, mars_args.environment_model.ephemerides_model)
        @test !all(iszero, mars_lpi0)
    end

    @testset "Earth starter-pack SPICE kernels fall back cleanly" begin
        lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
            kclear()
        end

        spice_path = joinpath(REPO_ROOT, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")
        earth = @test_nowarn Earth("", spice_path)
        et = ephemerides_time_seconds(
            InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
            SpiceEphemeridesModel()
        )
        l_pi = @test_nowarn planet_frame_lpi(earth, et, SpiceEphemeridesModel())
        @test !all(iszero, l_pi)
    end

    @testset "Simple ephemerides reject unsupported high-fidelity effectors" begin
        planet = Earth()
        spacecraft = make_spacecraft(ra_alt_m=500e3, rp_alt_m=400e3, ν_deg=175.0)
        args_srp = make_example_config(
            planet=planet,
            spacecraft=spacecraft,
            mission_time=60.0,
            initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
            dynamic_effectors=(SolarRadiationPressureModel(1.2, 10.0),),
            density_model=NoAtmosphereModel(),
            ephemerides_model=SimpleEphemeridesModel(),
            verbose=false,
            results=false
        )
        @test_throws ArgumentError run_simulation(args_srp)
    end
end
