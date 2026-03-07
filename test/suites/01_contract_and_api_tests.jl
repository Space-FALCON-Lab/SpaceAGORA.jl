@testset "SimulationModel Export Contract" begin
    sandbox = EXPORT_IMPORT_SANDBOX
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
    @test_nowarn Core.eval(sandbox, :(using .SimulationModel))

    required_public_names = [
        :SimulationConfiguration,
        :InitialCondition,
        :InitialTime,
        :MissionConfiguration,
        :EnvironmentModel,
        :IntegrationTolerances,
        :ControlModel,
        :GuidanceModel,
        :NavigationModel,
        :DynamicsModel
    ]

    for sym in required_public_names
        @test Core.eval(sandbox, :(isdefined(@__MODULE__, $(QuoteNode(sym)))))
    end
end

@testset "Include-Order + Name Ambiguity Smoke" begin
    sandbox = INCLUDE_ORDER_SANDBOX

    Base.include_string(sandbox, """
    module ConflictingExports
        export SimulationConfiguration
        struct SimulationConfiguration end
    end
    using .ConflictingExports
    """)

    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "core", "simulation_model.jl"))
    Core.eval(sandbox, :(const quat_mult = SimulationModel.quat_mult))
    @test_nowarn Base.include(sandbox, joinpath(REPO_ROOT, "src", "simulation", "engine", "simulation_engine.jl"))
    @test isdefined(sandbox, :SimulationEngine)
    @test Core.eval(sandbox, :(isdefined(SimulationEngine, :run_simulation)))
end

@testset "Simulation Filename Canonical Contract" begin
    engine_path = joinpath(REPO_ROOT, "src", "simulation", "engine", "execution.jl")
    legacy_exec_dir = joinpath(REPO_ROOT, "src", "simulation", "execution")
    legacy_run_path = joinpath(legacy_exec_dir, "run_simulation.jl")
    legacy_elements_path = joinpath(legacy_exec_dir, "simulation_elements.jl")
    legacy_execution_path = joinpath(legacy_exec_dir, "simulation_execution.jl")
    legacy_aerobraking_path = joinpath(REPO_ROOT, "src", "simulation", "Aerobraking.jl")
    legacy_complete_path = joinpath(REPO_ROOT, "src", "simulation", "Complete_passage.jl")

    @test isfile(engine_path)
    @test !isfile(legacy_run_path)
    @test !isfile(legacy_elements_path)
    @test !isfile(legacy_execution_path)
    @test !isfile(legacy_aerobraking_path)
    @test !isfile(legacy_complete_path)
end











function ensure_guidance_sandbox_loaded!()
    if GUIDANCE_SANDBOX_LOADED[]
        return
    end

    Core.eval(GUIDANCE_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "gnc", "guidance", "thruster_guidance", "thruster_guidance_models.jl"))))
    Core.eval(GUIDANCE_SANDBOX, :(include(joinpath(Main.REPO_ROOT, "src", "gnc", "guidance", "thruster_guidance", "thruster_guidance_functions.jl"))))
    GUIDANCE_SANDBOX_LOADED[] = true
end









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
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=EARTH),
            topography=true,
            topo_degree=8,
            topo_order=8,
            wind=false
        )
        @test env_ok.EI == 120.0

        @test_throws ArgumentError EnvironmentModel(
            planet=EARTH,
            EI=-1.0,
            density_model=NoAtmosphereModel(),
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
    force_nbody, torque_nbody = calcForceTorque(nbody_jupiter, nbody_state, ODEParams{1}(args=args_nbody), 1)
    @test all(isfinite, force_nbody)
    @test torque_nbody == SVector{3, Float64}(0.0, 0.0, 0.0)

    harmonics_file = joinpath(REPO_ROOT, "data/Gravity_harmonics_data", "EarthGGM05C.csv")
    harmonics_l20 = GravitationalHarmonicsModel(20, 20, harmonics_file, EARTH)
    @test size(harmonics_l20.C) == (21, 21)
    @test size(harmonics_l20.S) == (21, 21)
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

@testset "Run Simulation Metadata Return" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=20.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        keplerian=true
    )

    metadata = mktempdir() do tmp
        cd(tmp) do
            run_simulation(args; return_solution=true, return_solver_metadata=true)
        end
    end

    @test metadata isa NamedTuple
    @test hasproperty(metadata, :solution)
    @test hasproperty(metadata, :solver_mode)
    @test hasproperty(metadata, :solver_trace)
    @test hasproperty(metadata, :parallel_policy)
    @test hasproperty(metadata, :spice_counters)
end

@testset "RHS Completeness: Mass Derivative" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=120.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )

    u = build_initial_conditions(args)
    du = copy(u)
    du.sc[1].mass = 789.0
    p = ODEParams{1}(args=args)
    spacecraft_dynamics!(du, u, p, 0.0)
    @test du.sc[1].mass == 0.0

    du_inactive = copy(u)
    du_inactive.sc[1].mass = 123.0
    p_inactive = ODEParams{1}(args=args, is_active=[false])
    spacecraft_dynamics!(du_inactive, u, p_inactive, 0.0)
    @test du_inactive.sc[1].mass == 0.0
end







@testset "Guidance Thruster Campaign Helpers" begin
    ensure_guidance_sandbox_loaded!()
    sandbox = GUIDANCE_SANDBOX

    guidance_model = sandbox.AerobrakingCampaignPropulsiveManeuverGuidanceModel(
        maneuver_orbit_number=[2, 3],
        maneuver_Δv=[-5.0, 20.0]
    )

    thruster = BaseThrusterModel(
        thrust=[500.0],
        direction=[0.0],
        Δv=[0.0],
        start_burn_time=[-1.0],
        stop_burn_time=[-1.0],
        Isp=[300.0]
    )
    args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thruster,),
        control_rates=[1.0],
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p = ODEParams{1}(args=args)
    u = ComponentVector(pos=[0.0, 0.0, 0.0], vel=[0.0, 0.0, 0.0], mass=1.0, heat_loads=[0.0])

    p.orbit_counter[1] = 3
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test thruster.Δv[1] == 20.0
    @test isapprox(thruster.direction[1], π; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 2
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test thruster.Δv[1] == 5.0
    @test isapprox(thruster.direction[1], 0.0; atol=1e-12, rtol=0.0)

    p.orbit_counter[1] = 1
    sandbox.calcGuidanceEffect!(guidance_model, u, p, 0.0, 1)
    @test thruster.Δv[1] == 0.0
end

