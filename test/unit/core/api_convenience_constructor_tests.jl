@testset "API Convenience Constructors" begin
    @testset "Mission/Environment Config Validation" begin
        mc_default = MissionConfiguration()
        @test mc_default.mission_type == MissionTime

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_str = MissionConfiguration(mission_type="Time", mission_time=600.0, number_of_orbits=1, num_steps_to_save=10)
            @test mc_str.mission_type == MissionTime
            @test mc_str.mission_type == "Time"
            @test "Time" == mc_str.mission_type
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_sym = MissionConfiguration(mission_type=:orbits, mission_time=600.0, number_of_orbits=2, num_steps_to_save=10)
            @test mc_sym.mission_type == MissionOrbits
            @test mc_sym.mission_type == :orbits
            @test :orbits == mc_sym.mission_type
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            mc_enum = MissionConfiguration(mission_type=MissionOrbits, mission_time=600.0, number_of_orbits=3, num_steps_to_save=10)
            @test mc_enum.mission_type == MissionOrbits
            @test mc_enum.mission_type == "Orbits"
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "1") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
            @test_logs (:warn, r"mission_type=.*deprecated") MissionConfiguration(mission_type="Time")
            @test SimulationModel.SimConfig._deprecated_mission_type_input_warned[] == true
        end

        withenv("SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0") do
            SimulationModel.SimConfig._deprecated_mission_type_input_warned[] = false
            @test_logs MissionConfiguration(mission_type="Time")
            @test SimulationModel.SimConfig._deprecated_mission_type_input_warned[] == false
        end

        @test_throws ArgumentError MissionConfiguration(mission_type="invalid")
        @test_throws ArgumentError MissionConfiguration(mission_time=0.0)
        @test_throws ArgumentError MissionConfiguration(number_of_orbits=0)
        @test_throws ArgumentError MissionConfiguration(num_steps_to_save=0)

        env_ok = EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SpiceEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=8,
            topo_order=8,
            wind=false
        )
        @test env_ok.EI == 120.0
        @test env_ok.ephemerides_model isa SpiceEphemeridesModel

        env_no_gram_default = make_no_gram_environment()
        @test env_no_gram_default.ephemerides_model isa SimpleEphemeridesModel
        @test env_no_gram_default.density_model isa NoAtmosphereModel
        @test env_no_gram_default.planet isa Earth

        env_no_gram = make_no_gram_environment(planet=:mars, atmosphere=:exponential, EI_km=140.0)
        @test env_no_gram.ephemerides_model isa SimpleEphemeridesModel
        @test env_no_gram.density_model isa ExponentialAtmosphereModel

        exp_model = ExponentialAtmosphereModel(1.0e-4, 120.0e3, 12.0e3)
        @test exp_model.temperature_k == 200.0
        @test exp_model.valid_min_altitude_m == 120.0e3
        @test exp_model.valid_max_altitude_m == 180.0e3

        nrl_model = NRLMSISE00AtmosphereModel()
        @test nrl_model.f107a == 150.0
        @test nrl_model.f107 == 150.0
        @test nrl_model.ap == 4.0
        @test nrl_model.valid_min_altitude_m == 0.0
        @test nrl_model.valid_max_altitude_m == 1_000.0e3

        dt_nrl = DateTime(2024, 1, 1, 0, 0, 0)
        j2000_dt = DateTime(2000, 1, 1, 12, 0, 0)
        el_time_nrl = Dates.value(dt_nrl - j2000_dt) / 1000.0
        p_nrl = (
            args=(
                initial_time=InitialTime(year=2024, month=1, day=1, hour=0, minute=0, second=0.0),
                environment_model=(planet=EARTH,),
            ),
        )

        fixed_nrl = NRLMSISE00AtmosphereModel(f107a=120.0, f107=130.0, ap=6.0)
        expected_nrl_fixed = SatelliteToolbox.AtmosphericModels.nrlmsise00(
            dt_nrl, 400.0e3, 0.1, 0.2, 120.0, 130.0, 6.0
        )
        rho_nrl_fixed, T_nrl_fixed, wind_nrl_fixed = getDensity(fixed_nrl, 400.0e3, 0.1, 0.2, el_time_nrl, false)
        @test isapprox(rho_nrl_fixed, expected_nrl_fixed.total_density; atol=0.0, rtol=1e-12)
        @test isapprox(T_nrl_fixed, expected_nrl_fixed.temperature; atol=0.0, rtol=1e-12)
        @test wind_nrl_fixed == SVector{3, Float64}(0.0, 0.0, 0.0)

        rho_nrl_rel, T_nrl_rel, wind_nrl_rel = getDensity(fixed_nrl, 400.0e3, 0.1, 0.2, 0.0, false, p_nrl)
        @test isapprox(rho_nrl_rel, expected_nrl_fixed.total_density; atol=0.0, rtol=1e-12)
        @test isapprox(T_nrl_rel, expected_nrl_fixed.temperature; atol=0.0, rtol=1e-12)
        @test wind_nrl_rel == SVector{3, Float64}(0.0, 0.0, 0.0)

        provider_hits = Ref(0)
        provider_nrl = NRLMSISE00AtmosphereModel(
            index_provider=(instant, h, lat, lon) -> begin
                provider_hits[] += 1
                return (f107a=95.0, f107=105.0, ap=[8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0])
            end
        )
        expected_nrl_provider = SatelliteToolbox.AtmosphericModels.nrlmsise00(
            dt_nrl, 400.0e3, 0.1, 0.2, 95.0, 105.0, [8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0]
        )
        rho_nrl_provider, T_nrl_provider, wind_nrl_provider = getDensity(provider_nrl, 400.0e3, 0.1, 0.2, el_time_nrl, false)
        @test provider_hits[] == 1
        @test isapprox(rho_nrl_provider, expected_nrl_provider.total_density; atol=0.0, rtol=1e-12)
        @test isapprox(T_nrl_provider, expected_nrl_provider.temperature; atol=0.0, rtol=1e-12)
        @test wind_nrl_provider == SVector{3, Float64}(0.0, 0.0, 0.0)

        space_indices_nrl = NRLMSISE00AtmosphereModel(use_space_indices=true)
        @test space_indices_nrl.index_provider isa SimulationModel.EnvironmentModels.NRLMSISE00SpaceIndicesProvider
        low_altitude_indices = space_indices_nrl.index_provider(dt_nrl, 70.0e3, 0.1, 0.2)
        @test low_altitude_indices == (f107a=150.0, f107=150.0, ap=4.0)

        helper_calls = Tuple{Symbol, DateTime}[]
        helper_instant = DateTime(2024, 1, 2, 10, 30, 0)
        function fake_space_indices_lookup(index, instant::DateTime)
            index_sym = if index isa Val{:F10adj_avg_center81}
                :F10adj_avg_center81
            elseif index isa Val{:F10adj}
                :F10adj
            elseif index isa Val{:Ap_daily}
                :Ap_daily
            elseif index isa Val{:Ap}
                :Ap
            else
                error("Unexpected index lookup in test: $index")
            end
            push!(helper_calls, (index_sym, instant))
            if index_sym === :F10adj_avg_center81
                return 177.0
            elseif index_sym === :F10adj
                return 166.0
            elseif index_sym === :Ap_daily
                return 99
            end
            day = Date(instant)
            if day == Date(2024, 1, 2)
                return (21, 22, 23, 24, 25, 26, 27, 28)
            elseif day == Date(2024, 1, 1)
                return (11, 12, 13, 14, 15, 16, 17, 18)
            elseif day == Date(2023, 12, 31)
                return (1, 2, 3, 4, 5, 6, 7, 8)
            end
            error("Unexpected Ap lookup date in test: $day")
        end

        helper_indices = SimulationModel.EnvironmentModels._nrlmsise_space_indices_indices(
            fake_space_indices_lookup,
            helper_instant
        )
        @test helper_indices.f107a == 177.0
        @test helper_indices.f107 == 166.0
        @test helper_indices.ap == SVector{7, Float64}(99.0, 24.0, 23.0, 22.0, 21.0, 14.5, 4.5)
        @test (Symbol(:F10adj_avg_center81), helper_instant) in helper_calls
        @test (Symbol(:F10adj), helper_instant - Day(1)) in helper_calls
        @test SimulationModel.EnvironmentModels._nrlmsise_ap_slot_index(helper_instant) == 4

        @test_throws ArgumentError NRLMSISE00AtmosphereModel(ap=[1.0, 2.0])
        @test_throws ArgumentError NRLMSISE00AtmosphereModel(
            use_space_indices=true,
            index_provider=(instant) -> (f107a=1.0, f107=1.0, ap=1.0)
        )
        @test_throws ArgumentError NRLMSISE00AtmosphereModel(space_indices_force_download=true)
        @test_throws ArgumentError NRLMSISE00AtmosphereModel(valid_min_altitude_m=1.0, valid_max_altitude_m=0.0)

        ρ_mid = 1.225 * exp(-20.0e3 / 8.5e3)
        piecewise = PiecewiseExponentialAtmosphereModel(
            [0.0, 20.0e3, 60.0e3],
            [1.225, ρ_mid],
            [8.5e3, 12.0e3];
            temperature_k=210.0
        )
        @test piecewise.valid_min_altitude_m == 0.0
        @test piecewise.valid_max_altitude_m == 60.0e3

        ρ_piece_low, T_piece_low, wind_piece_low = getDensity(piecewise, 10.0e3, 0.0, 0.0, 0.0, false)
        @test isapprox(ρ_piece_low, 1.225 * exp(-10.0e3 / 8.5e3); atol=0.0, rtol=1e-12)
        @test T_piece_low == 210.0
        @test wind_piece_low == SVector{3, Float64}(0.0, 0.0, 0.0)

        ρ_piece_high, T_piece_high, wind_piece_high = getDensity(piecewise, 30.0e3, 0.0, 0.0, 0.0, false)
        @test isapprox(ρ_piece_high, ρ_mid * exp(-(30.0e3 - 20.0e3) / 12.0e3); atol=0.0, rtol=1e-12)
        @test T_piece_high == 210.0
        @test wind_piece_high == SVector{3, Float64}(0.0, 0.0, 0.0)

        rhos_piece = zeros(3)
        Ts_piece = zeros(3)
        winds_piece = [SVector{3, Float64}(0.0, 0.0, 0.0) for _ in 1:3]
        getDensityBatch!(
            rhos_piece,
            Ts_piece,
            winds_piece,
            piecewise,
            [10.0e3, 30.0e3, 80.0e3],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            0.0,
            false,
            nothing
        )
        @test isapprox(rhos_piece[1], ρ_piece_low; atol=0.0, rtol=1e-12)
        @test isapprox(rhos_piece[2], ρ_piece_high; atol=0.0, rtol=1e-12)
        @test isapprox(rhos_piece[3], ρ_mid * exp(-(80.0e3 - 20.0e3) / 12.0e3); atol=0.0, rtol=1e-12)
        @test all(==(210.0), Ts_piece)
        @test all(==(SVector{3, Float64}(0.0, 0.0, 0.0)), winds_piece)

        @test_throws ArgumentError PiecewiseExponentialAtmosphereModel([0.0], [1.0], [8.5e3])
        @test_throws ArgumentError PiecewiseExponentialAtmosphereModel([0.0, 20.0e3, 10.0e3], [1.0, 0.1], [8.5e3, 7.0e3])

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=-1.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SpiceEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=false,
            topo_degree=8,
            topo_order=8,
            wind=false
        )

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SpiceEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=-1,
            topo_order=8,
            wind=false
        )

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=120.0,
            density_model=NoAtmosphereModel(),
            ephemerides_model=SpiceEphemeridesModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=8,
            topo_order=-1,
            wind=false
        )
    end

    ic = InitialCondition()
    @test ic isa InitialCondition
    @test ic.a == 0.0
    @test ic.e == 0.0

    link = Link()
    @test link isa Link{0}

    joint = Joint()
    @test joint isa Joint

    sc = SpacecraftModel()
    @test sc isa SpacecraftModel
    @test sc.root.root

    custom_root = Link{0}(root=true, m=100.0)
    sc_custom = SpacecraftModel(root=custom_root, id=42)
    @test sc_custom.id == 42
    @test sc_custom.root === custom_root
    @test any(link -> link === custom_root, sc_custom.links)

    nbody = NBodyGravityModel(["Sun"], "Earth", SPICE_PATH)
    @test nbody isa NBodyGravityModel
    @test nbody.primary_body_name == "Earth"
    @test nbody.body_names == ("Sun",)

    nbody_mars = NBodyGravityModel(["Sun"], "Mars", SPICE_PATH)
    nbody_venus = NBodyGravityModel(["Sun"], "Venus", SPICE_PATH)
    nbody_titan = NBodyGravityModel(["Sun"], "Titan", SPICE_PATH)
    @test lowercase(nbody_mars.planet.name) == "mars"
    @test lowercase(nbody_venus.planet.name) == "venus"
    @test lowercase(nbody_titan.planet.name) == "titan"
    @test_throws ArgumentError NBodyGravityModel(["Sun"], "Pluto", SPICE_PATH)

    nbody_jupiter = NBodyGravityModel(["Jupiter"], "Earth", SPICE_PATH)
    nbody_state = [EARTH.Rp_e + 500e3, 0.0, 0.0, 0.0, 0.0, 0.0, 500.0]
    sc_nbody = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=175.0)
    args_nbody = build_config(
        spacecraft=sc_nbody,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    force_nbody, torque_nbody = calcForceTorque(nbody_jupiter, nbody_state, ODEParams(n_sats=1, args=args_nbody), 1)
    @test all(isfinite, force_nbody)
    @test torque_nbody == SVector{3, Float64}(0.0, 0.0, 0.0)

    harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
    harmonics_l20 = GravitationalHarmonicsModel(20, 20, harmonics_file, EARTH)
    @test size(harmonics_l20.C) == (21, 21)
    @test size(harmonics_l20.S) == (21, 21)
    @test harmonics_l20.coefficient_normalization == :full
    @test harmonics_l20.include_central == true
    @test GravitationalHarmonicsModel(2, 0, harmonics_file, EARTH; include_central=false).include_central == false
    @test_throws ArgumentError GravitationalHarmonicsModel(10, 11, harmonics_file, EARTH)

    child_link = Link(root=false, q=MVector{4, Float64}(sin(pi / 4), 0.0, 0.0, cos(pi / 4)))
    rot_child = rotate_to_body(child_link)
    @test size(rot_child) == (3, 3)
    @test isapprox(det(Matrix(rot_child)), 1.0; atol=1e-12)
    @test norm(Matrix(rot_child) - Matrix{Float64}(I, 3, 3)) > 0.1

    @testset "Quaternion DCM Conversion Negative-Trace Branch" begin
        dcm_180_x = SMatrix{3, 3, Float64}(1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0)
        q_neg_trace = SimulationModel.dcm_to_quaternion(dcm_180_x)
        @test isapprox(norm(q_neg_trace), 1.0; atol=1e-12, rtol=0.0)

        dcm_roundtrip = SimulationModel.rot(q_neg_trace)
        @test isapprox(Matrix(dcm_roundtrip), Matrix(dcm_180_x); atol=1e-12, rtol=0.0)
    end

    @testset "Effector Rate Validation" begin
        @test GuidanceModel((), Float64[]) isa GuidanceModel
        @test NavigationModel((), Float64[]) isa NavigationModel
        @test ControlModel((), Float64[]) isa ControlModel

        @test_throws ArgumentError GuidanceModel((:g1,), Float64[])
        @test_throws ArgumentError NavigationModel((:n1,), Float64[])
        @test_throws ArgumentError ControlModel((:c1,), Float64[])

        @test_throws ArgumentError GuidanceModel((:g1,), [0.0])
        @test_throws ArgumentError NavigationModel((:n1,), [-1.0])
        @test_throws ArgumentError ControlModel((:c1,), [Inf])

        sc1 = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
        sc2 = make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        args_bad_slots = build_config_multi(
            spacecraft=[sc1, sc2],
            density_model=NoAtmosphereModel(),
            orientation_sim=false,
            mission_time=60.0,
            EI_km=120.0,
            dynamic_effectors=(InverseSquaredGravityModel(),),
            control_effectors=(make_base_thruster_model(thrust=1.0, Δv=1.0, start_burn_time=-1.0, stop_burn_time=-1.0),),
            control_rates=[1.0],
            keplerian=true
        )
        @test_throws ArgumentError run_case_silent(args_bad_slots)
    end
end
