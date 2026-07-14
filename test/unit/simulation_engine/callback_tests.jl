@testset "Callback Internal Branch Coverage" begin
    mutable struct MockCallbackOpts
        dtmax::Float64
        reltol::Float64
        abstol::Float64
    end

    mutable struct MockCallbackIntegrator{P, U, O}
        p::P
        u::U
        t::Float64
        opts::O
        tdir::Int
        tstop_max::Float64
    end
    DiffEqBase.get_tstops_max(integrator::MockCallbackIntegrator) = integrator.tstop_max
    DiffEqBase.add_tstop!(integrator::MockCallbackIntegrator, tstop) = begin
        integrator.tstop_max = max(integrator.tstop_max, Float64(tstop))
        nothing
    end

    args_base = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=true, generate_plots=false, normalize=false)
    )
    args_orient = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=true, generate_plots=false, normalize=false)
    )

    mission_orbits = MissionConfiguration(
        mission_type=MissionOrbits,
        keplerian=true,
        number_of_orbits=1,
        mission_time=120.0,
        orientation_sim=false,
        num_steps_to_save=100
    )
    args_orbits = SimulationConfiguration(
        file_paths=args_base.file_paths,
        simulation_settings=SimulationSettings(results=false, verbose=true, generate_plots=false, normalize=false),
        mission_configuration=mission_orbits,
        environment_model=args_base.environment_model,
        dynamics_model=args_base.dynamics_model,
        guidance_model=args_base.guidance_model,
        navigation_model=args_base.navigation_model,
        control_model=args_base.control_model,
        initial_time=args_base.initial_time,
        integration_tolerances=args_base.integration_tolerances
    )

    _ = SimulationModel.SimulationCallbacks.get_callbacks(
        1,
        args_orbits.dynamics_model.dynamic_effectors,
        args_orbits
    )

    thruster_control = make_base_thruster_model(
        thrust=2.0,
        Δv=20.0,
        start_burn_time=-1.0,
        stop_burn_time=-1.0,
        Isp=300.0
    )
    args_control = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(thruster_control,),
        control_rates=[1.0],
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    control_cbs = SimulationModel.SimulationCallbacks.get_control_callbacks(1, args_control)
    p_control = ODEParams(n_sats=1, args=args_control)
    u_control = build_initial_conditions(args_control)
    integrator_control = MockCallbackIntegrator(
        p_control,
        u_control,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    control_cbs[1].affect!(integrator_control)
    @test isfinite(thruster_control.start_burn_time[1])
    @test isfinite(thruster_control.stop_burn_time[1])
    @test integrator_control.tstop_max >= thruster_control.start_burn_time[1]

    orbit_cb = SimulationModel.SimulationCallbacks.get_orbit_end_callback(1)
    p_orbit = ODEParams(n_sats=1, args=args_orbits)
    u_orbit = build_initial_conditions(args_orbits)
    integrator_orbit = MockCallbackIntegrator(
        p_orbit,
        u_orbit,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    out_orbit = zeros(1)
    orbit_cb.condition(out_orbit, u_orbit, 0.0, integrator_orbit)
    @test isfinite(out_orbit[1])
    orbit_count_before = p_orbit.orbit_counter[1]
    orbit_cb.affect!(integrator_orbit, 1)
    @test p_orbit.orbit_counter[1] == orbit_count_before + 1

    args_impact = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    impact_cb = SimulationModel.SimulationCallbacks.get_impact_callback(2)
    p_impact = ODEParams(n_sats=2, args=args_impact)
    u_impact = build_initial_conditions(args_impact)
    integrator_impact = MockCallbackIntegrator(
        p_impact,
        u_impact,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    impact_out = zeros(2)
    impact_cb.condition(impact_out, u_impact, 0.0, integrator_impact)
    @test all(impact_out .> 0.0)
    @test impact_cb.affect! === nothing
    impact_cb.affect_neg!(integrator_impact, 1)
    @test p_impact.is_active[1] == false
    @test p_impact.is_active[2] == true

    args_drag = build_config(
        spacecraft=make_spacecraft(ra_alt_m=220e3, rp_alt_m=100e3, ν_deg=180.0),
        density_model=ConstantDensityModel(1e-6, 240.0),
        orientation_sim=false,
        mission_time=180.0,
        EI_km=140.0,
        dynamic_effectors=(InverseSquaredGravityModel(), AerodynamicCoefficientfM()),
        keplerian=false,
        simulation_settings=SimulationSettings(results=false, verbose=true, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_atmosphere=1e-8,
            abstol_atmosphere=1e-8,
            dt_max_orbit=5.0,
            dt_max_atmosphere=0.2
        )
    )
    drag_cb = SimulationModel.SimulationCallbacks.get_drag_state_callback(1)
    p_drag = ODEParams(n_sats=1, args=args_drag)
    u_drag = build_initial_conditions(args_drag)
    integrator_drag = MockCallbackIntegrator(
        p_drag,
        u_drag,
        12.0,
        MockCallbackOpts(
            args_drag.integration_tolerances.dt_max_atmosphere,
            args_drag.integration_tolerances.reltol_atmosphere,
            args_drag.integration_tolerances.abstol_atmosphere
        ),
        1,
        Inf
    )
    drag_output = ""
    mktemp() do path, io
        redirect_stdout(io) do
            drag_cb.affect!(integrator_drag, 1)
        end
        flush(io)
        seekstart(io)
        drag_output = read(io, String)
    end
    @test occursin("Switching to space integration", drag_output)
    @test integrator_drag.opts.dtmax == args_drag.integration_tolerances.dt_max_orbit
    @test integrator_drag.opts.reltol == args_drag.integration_tolerances.reltol_orbit
    @test integrator_drag.opts.abstol == args_drag.integration_tolerances.abstol_orbit

    quat_proj_cb = SimulationModel.SimulationCallbacks.get_quaternion_projection_callback(1, args_orient)
    p_orient = ODEParams(n_sats=1, args=args_orient)
    u_orient = build_initial_conditions(args_orient)
    u_orient.sc[1].q .= [0.0, 0.0, 0.0, 2.0]
    integrator_orient = MockCallbackIntegrator(
        p_orient,
        u_orient,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    @test quat_proj_cb.condition(u_orient, 0.0, integrator_orient) == true
    quat_proj_cb.affect!(integrator_orient)
    @test isapprox(norm(integrator_orient.u.sc[1].q), 1.0; atol=1e-12, rtol=0.0)

    u_orient_unit = build_initial_conditions(args_orient)
    integrator_orient_unit = MockCallbackIntegrator(
        p_orient,
        u_orient_unit,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    @test quat_proj_cb.condition(u_orient_unit, 0.0, integrator_orient_unit) == true

    counting_navigation = CountingNavigationModel([0])
    args_navigation = SimulationConfiguration(
        file_paths=args_base.file_paths,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=args_base.mission_configuration,
        environment_model=args_base.environment_model,
        dynamics_model=args_base.dynamics_model,
        guidance_model=args_base.guidance_model,
        navigation_model=NavigationModel(navigation_effectors=(counting_navigation,), navigation_rates=[1.0]),
        control_model=args_base.control_model,
        initial_time=args_base.initial_time,
        integration_tolerances=args_base.integration_tolerances
    )
    navigation_cbs = SimulationModel.SimulationCallbacks.get_navigation_callbacks(1, args_navigation)
    p_navigation = ODEParams(n_sats=1, args=args_navigation)
    u_navigation = build_initial_conditions(args_navigation)
    integrator_navigation = MockCallbackIntegrator(
        p_navigation,
        u_navigation,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    navigation_cbs[1].affect!.affect!(integrator_navigation)
    @test counting_navigation.hits == [1]
    withenv("SPACEAGORA_DEV_HOT_RELOAD" => "1") do
        navigation_cbs_hot = SimulationModel.SimulationCallbacks.get_navigation_callbacks(1, args_navigation)
        navigation_cbs_hot[1].affect!.affect!(integrator_navigation)
    end
    @test counting_navigation.hits == [2]

    counting_guidance = CountingGuidanceModel([0, 0])
    args_guidance_base = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    args_guidance = SimulationConfiguration(
        file_paths=args_guidance_base.file_paths,
        simulation_settings=args_guidance_base.simulation_settings,
        mission_configuration=args_guidance_base.mission_configuration,
        environment_model=args_guidance_base.environment_model,
        dynamics_model=args_guidance_base.dynamics_model,
        guidance_model=GuidanceModel(guidance_effectors=(counting_guidance,), guidance_rates=[1.0]),
        navigation_model=args_guidance_base.navigation_model,
        control_model=args_guidance_base.control_model,
        initial_time=args_guidance_base.initial_time,
        integration_tolerances=args_guidance_base.integration_tolerances
    )
    guidance_cbs = SimulationModel.SimulationCallbacks.get_guidance_callbacks(2, args_guidance)
    p_guidance = ODEParams(n_sats=2, args=args_guidance)
    u_guidance = build_initial_conditions(args_guidance)
    integrator_guidance = MockCallbackIntegrator(
        p_guidance,
        u_guidance,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    guidance_cbs[1].affect!.affect!(integrator_guidance)
    @test counting_guidance.hits == [1, 1]
    withenv("SPACEAGORA_DEV_HOT_RELOAD" => "1") do
        guidance_cbs_hot = SimulationModel.SimulationCallbacks.get_guidance_callbacks(2, args_guidance)
        guidance_cbs_hot[1].affect!.affect!(integrator_guidance)
    end
    @test counting_guidance.hits == [2, 2]

    control_model = CountingControlModel([0, 0])
    args_control = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(control_model,),
        control_rates=[1.0],
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_control = ODEParams(n_sats=2, args=args_control)
    u_control = build_initial_conditions(args_control)
    integrator_control = MockCallbackIntegrator(
        p_control,
        u_control,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    withenv("SPACEAGORA_DEV_HOT_RELOAD" => "1") do
        control_cbs = SimulationModel.SimulationCallbacks.get_control_callbacks(2, args_control)
        control_cbs[1].affect!.affect!(integrator_control)
    end
    @test control_model.hits == [1, 1]

    requires_density = SimulationModel.SimulationCallbacks._requires_density_callback
    requires_orbit_end = SimulationModel.SimulationCallbacks._requires_orbit_end_callback
    requires_drag_state = SimulationModel.SimulationCallbacks._requires_drag_state_callback
    requires_quat_projection = SimulationModel.SimulationCallbacks._requires_quaternion_projection_callback
    density_use_threads = SimulationModel.SimulationCallbacks._density_callback_use_threads
    control_use_threads = SimulationModel.SimulationCallbacks._control_callback_use_threads

    @test requires_density((InverseSquaredGravityModel(),), args_base) == false
    @test requires_density((InverseSquaredGravityModel(), AerodynamicCoefficientfM()), args_base) == true
    @test requires_density((InverseSquaredGravityModel(),), args_drag) == true
    @test requires_orbit_end(args_base) == false
    @test requires_orbit_end(args_orbits) == true
    @test requires_drag_state((InverseSquaredGravityModel(),), args_base) == false
    @test requires_drag_state((InverseSquaredGravityModel(), AerodynamicCoefficientfM()), args_drag) == true
    @test requires_quat_projection(args_base) == false
    @test requires_quat_projection(args_orient) == true

    has_worker_threads = Threads.nthreads() > 1
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1"
    ) do
        @test density_use_threads(args_control, 8) == false
    end
    # Pin the auto-mode budget floor so this branch check is independent of the
    # default (16), which exceeds the CI thread count.
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET" => "2"
    ) do
        @test density_use_threads(args_control, 8) == has_worker_threads
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "64"
    ) do
        @test density_use_threads(args_control, 2) == has_worker_threads
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "0"
    ) do
        @test density_use_threads(args_control, 8) == false
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_AUTO_THREAD_MIN_BUDGET" => "2",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1"
    ) do
        @test density_use_threads(args_control, 8) == has_worker_threads
    end

    control_thruster = BaseThrusterModel(
        thrust=[0.5, 0.6],
        direction=[0.0, π],
        Δv=[0.0, 0.0],
        start_burn_time=[0.0, 0.0],
        stop_burn_time=[10.0, 10.0],
        Isp=[300.0, 300.0]
    )
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "1"
    ) do
        @test control_use_threads(control_thruster, 8, false) == false
    end
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "1"
    ) do
        @test control_use_threads(control_thruster, 8, false) == has_worker_threads
        @test control_use_threads(control_thruster, 8, true) == false
    end
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "0"
    ) do
        @test control_use_threads(control_thruster, 8, false) == false
    end
    withenv(
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1"
    ) do
        @test control_use_threads(control_thruster, 8, false) == has_worker_threads
    end

    n_parallel_sats = 4
    threaded_spacecraft = [
        make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0 - 5.0 * (k - 1))
        for k in 1:n_parallel_sats
    ]
    threaded_control = CountingControlModel(zeros(Int, n_parallel_sats))
    args_control_parallel = build_config_multi(
        spacecraft=threaded_spacecraft,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        control_effectors=(threaded_control,),
        control_rates=[1.0],
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_control_parallel = ODEParams(n_sats=n_parallel_sats, args=args_control_parallel)
    u_control_parallel = build_initial_conditions(args_control_parallel)
    integrator_control_parallel = MockCallbackIntegrator(
        p_control_parallel,
        u_control_parallel,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    withenv(
        "SPACEAGORA_DEV_HOT_RELOAD" => "0",
        "SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_CONTROL_CALLBACK_ASSUME_THREADSAFE" => "1"
    ) do
        control_cbs_parallel = SimulationModel.SimulationCallbacks.get_control_callbacks(
            n_parallel_sats,
            args_control_parallel
        )
        control_cbs_parallel[1].affect!.affect!(integrator_control_parallel)
    end
    @test threaded_control.hits == ones(Int, n_parallel_sats)

    args_density_parallel = build_config_multi(
        spacecraft=threaded_spacecraft,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_density_parallel = ODEParams(n_sats=n_parallel_sats, args=args_density_parallel)
    u_density_parallel = build_initial_conditions(args_density_parallel)
    integrator_density_parallel = MockCallbackIntegrator(
        p_density_parallel,
        u_density_parallel,
        0.0,
        MockCallbackOpts(1.0, 1e-8, 1e-8),
        1,
        Inf
    )
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1"
    ) do
        density_cb_parallel = SimulationModel.SimulationCallbacks.get_density_callback(
            n_parallel_sats,
            args_density_parallel
        )
        density_cb_parallel.affect!(integrator_density_parallel)
    end
    @test all(isfinite, p_density_parallel.shared_buffers.densities)
    @test all(ρ -> ρ >= 0.0, p_density_parallel.shared_buffers.densities)
end

@testset "Callback Env Helper Branch Coverage" begin
    callbacks = SimulationModel.SimulationCallbacks

    @test callbacks._parse_bool_env("SPACEAGORA_CB_TEST_BOOL", false) == false
    withenv("SPACEAGORA_CB_TEST_BOOL" => "yes") do
        @test callbacks._parse_bool_env("SPACEAGORA_CB_TEST_BOOL", false) == true
    end
    withenv("SPACEAGORA_CB_TEST_BOOL" => "off") do
        @test callbacks._parse_bool_env("SPACEAGORA_CB_TEST_BOOL", true) == false
    end
    withenv("SPACEAGORA_CB_TEST_BOOL" => "invalid") do
        @test_throws ArgumentError callbacks._parse_bool_env("SPACEAGORA_CB_TEST_BOOL", false)
    end

    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off") do
        @test callbacks._density_callback_parallel_mode() == :off
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on") do
        @test callbacks._density_callback_parallel_mode() == :on
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto") do
        @test callbacks._density_callback_parallel_mode() == :auto
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "invalid") do
        @test_throws ArgumentError callbacks._density_callback_parallel_mode()
    end

    withenv("SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "4") do
        @test callbacks._density_callback_thread_threshold() == 4
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "0") do
        @test callbacks._density_callback_thread_threshold() == 1
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "oops") do
        @test_throws ArgumentError callbacks._density_callback_thread_threshold()
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1") do
        @test callbacks._density_callback_allow_with_outer() == true
    end
    withenv("SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "invalid") do
        @test_throws ArgumentError callbacks._density_callback_allow_with_outer()
    end

    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "off") do
        @test callbacks._control_callback_parallel_mode() == :off
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "on") do
        @test callbacks._control_callback_parallel_mode() == :on
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "auto") do
        @test callbacks._control_callback_parallel_mode() == :auto
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => "invalid") do
        @test_throws ArgumentError callbacks._control_callback_parallel_mode()
    end

    withenv("SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "4") do
        @test callbacks._control_callback_thread_threshold() == 4
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "0") do
        @test callbacks._control_callback_thread_threshold() == 1
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "oops") do
        @test_throws ArgumentError callbacks._control_callback_thread_threshold()
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1") do
        @test callbacks._control_callback_allow_with_outer() == true
    end
    withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "invalid") do
        @test_throws ArgumentError callbacks._control_callback_allow_with_outer()
    end

    @test callbacks.density_model_threadsafe(NoAtmosphereModel()) == true
    @test callbacks.density_model_threadsafe(ConstantDensityModel(1e-6, 220.0)) == false

    thruster_model = make_base_thruster_model(thrust=0.1, Δv=0.0, start_burn_time=0.0, stop_burn_time=1.0)
    @test callbacks.control_model_threadsafe(thruster_model) == true
    @test callbacks.control_model_threadsafe(CountingControlModel([0])) == false

    callbacks._gram_runtime_stats_reset!()
    @test callbacks._gram_runtime_stats_update!(s -> begin
        s.density_calls += 1
        s.direct_calls += 1
    end) === nothing
    stats_snap = callbacks._gram_runtime_stats_snapshot()
    @test stats_snap.density_calls == 1
    @test stats_snap.direct_calls == 1

    withenv("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW" => "off") do
        @test callbacks._gram_track_cache_ignore_time_window() == false
    end
    withenv("SPACEAGORA_GRAM_TRACK_CACHE_TARGET_USE_J2" => "on") do
        @test callbacks._gram_track_cache_target_use_j2() == true
    end

    withenv("SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => "on") do
        @test callbacks._thermal_callback_parallel_mode() == :on
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL" => nothing
    ) do
        @test callbacks._thermal_callback_parallel_mode() == :auto
    end
    withenv("SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD" => "5") do
        @test callbacks._thermal_callback_thread_threshold() == 5
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "6",
        "SPACEAGORA_THERMAL_CALLBACK_THREAD_THRESHOLD" => nothing
    ) do
        @test callbacks._thermal_callback_thread_threshold() == 6
    end
    withenv("SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1") do
        @test callbacks._thermal_callback_allow_with_outer() == true
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "1",
        "SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => nothing
    ) do
        @test callbacks._thermal_callback_allow_with_outer() == true
    end
    withenv("SPACEAGORA_THERMAL_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "invalid") do
        @test_throws ArgumentError callbacks._thermal_callback_allow_with_outer()
    end

    @test callbacks._is_gram_density_model(NoAtmosphereModel()) == false

    args_density_lookup = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    p_density_lookup = ODEParams(n_sats=1, args=args_density_lookup)
    gram_model_lookup = SimulationModel.GRAMAtmosphereModel(nothing)
    push!(p_density_lookup.shared_buffers.density_models, gram_model_lookup)
    @test callbacks._density_model_for_sat(p_density_lookup, 1) === gram_model_lookup
    @test callbacks._density_model_for_sat(p_density_lookup, 2) isa NoAtmosphereModel

    withenv("SPACEAGORA_CB_TEST_FLOAT" => "oops") do
        @test_throws ArgumentError callbacks._parse_float_env("SPACEAGORA_CB_TEST_FLOAT", 1.0)
    end
    withenv("SPACEAGORA_CB_TEST_FLOAT_OPT" => "oops") do
        @test_throws ArgumentError callbacks._parse_float_env_optional("SPACEAGORA_CB_TEST_FLOAT_OPT")
    end
    withenv("SPACEAGORA_CB_TEST_INT_OPT" => "oops") do
        @test_throws ArgumentError callbacks._parse_int_env_optional("SPACEAGORA_CB_TEST_INT_OPT")
    end

    withenv("SPACEAGORA_GRAM_TRACK_CACHE" => "on") do
        @test callbacks._gram_track_cache_mode() == :on
    end
    withenv("SPACEAGORA_GRAM_TRACK_CACHE" => "auto") do
        @test callbacks._gram_track_cache_mode() == :auto
    end
    withenv("SPACEAGORA_GRAM_TRACK_CACHE" => "invalid") do
        @test_throws ArgumentError callbacks._gram_track_cache_mode()
    end

    cfg = callbacks.GramTrackCacheConfig(
        :on,
        2.0,
        100.0,
        deg2rad(0.5),
        8,
        10.0,
        1000.0,
        deg2rad(4.0),
        32,
        20_000.0
    )
    @test callbacks._gram_track_cache_enabled(cfg, NoAtmosphereModel()) == false
    entry_profile = callbacks._gram_track_cache_profile(cfg, p_density_lookup, 130e3)
    orbit_profile = callbacks._gram_track_cache_profile(cfg, p_density_lookup, 200e3)
    @test entry_profile[1] == cfg.entry_horizon_s
    @test orbit_profile[1] == cfg.orbit_horizon_s

    @test callbacks._lerp(0.0, 10.0, 0.25) == 2.5
    @test isapprox(callbacks._angdiff_rad(0.0, 2pi - 0.1), 0.1; atol=1e-12, rtol=0.0)
    @test isapprox(callbacks._lerp_angle_rad(0.0, Float64(pi), 0.5), Float64(pi) / 2; atol=1e-12, rtol=0.0)

    cache = callbacks.GramTrackCache()
    cache.valid = true
    cache.t0 = 0.0
    cache.t1 = 10.0
    cache.index_hint = 1
    cache.times = [0.0, 5.0, 10.0]
    cache.alts = [1000.0, 2000.0, 3000.0]
    cache.lats = [0.0, 0.1, 0.2]
    cache.lons = [0.0, 0.1, 0.2]
    cache.rhos = [1.0, 2.0, 3.0]
    cache.Ts = [200.0, 220.0, 240.0]
    cache.winds = [SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(1.0, 0.0, 0.0), SVector{3, Float64}(2.0, 0.0, 0.0)]

    withenv("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW" => "0") do
        @test callbacks._gram_track_cache_segment(cache, -1.0) === nothing
    end
    withenv("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW" => "1") do
        seg_clamped = callbacks._gram_track_cache_segment(cache, -1.0)
        @test seg_clamped == (1, 0.0)
    end

    seg_hit = callbacks._gram_track_cache_ready(cache, 2.5, 1500.0, 0.05, 0.05, 1e-9, 1e-9)
    @test seg_hit == (1, 0.5)
    @test callbacks._gram_track_cache_ready(cache, 2.5, 2000.0, 0.05, 0.05, 1e-9, 1e-9) === nothing
    ρ_eval, T_eval, wind_eval = callbacks._gram_track_cache_eval(cache, 1, 0.25)
    @test isapprox(ρ_eval, 1.25; atol=1e-12, rtol=0.0)
    @test isapprox(T_eval, 205.0; atol=1e-12, rtol=0.0)
    @test isapprox(wind_eval[1], 0.25; atol=1e-12, rtol=0.0)

    cache_flat = callbacks.GramTrackCache()
    cache_flat.valid = true
    cache_flat.t0 = 2.0
    cache_flat.t1 = 2.0
    cache_flat.index_hint = 1
    cache_flat.times = [2.0, 2.0]
    cache_flat.alts = [1000.0, 1000.0]
    cache_flat.lats = [0.0, 0.0]
    cache_flat.lons = [0.0, 0.0]
    cache_flat.rhos = [1.0, 1.0]
    cache_flat.Ts = [200.0, 200.0]
    cache_flat.winds = [SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0)]
    withenv("SPACEAGORA_GRAM_TRACK_CACHE_IGNORE_TIME_WINDOW" => "1") do
        @test callbacks._gram_track_cache_segment(cache_flat, 2.0) == (1, 0.0)
    end

    @test callbacks._uses_j2_gravity_effector((InverseSquaredGravityModel(), InverseSquaredJ2GravityModel())) == true

    tol_args = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=true,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        tolerances=IntegrationTolerances(
            reltol_orbit=1e-8,
            abstol_orbit=1e-8,
            reltol_quaternion=1e-9,
            abstol_quaternion=1e-9,
            reltol_angular_rate=4e-7,
            abstol_angular_rate=5e-8
        )
    )
    u_tol = build_initial_conditions(tol_args)
    template_reltol, template_abstol = _build_solver_tolerances(u_tol, tol_args)
    reltol_phase, abstol_phase = callbacks._callback_tolerances_for_phase(template_reltol, template_abstol, tol_args, false)
    @test all(isapprox.(reltol_phase.sc[1].q, tol_args.integration_tolerances.reltol_quaternion; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_phase.sc[1].q, tol_args.integration_tolerances.abstol_quaternion; atol=0.0, rtol=0.0))
    @test all(isapprox.(reltol_phase.sc[1].ω, tol_args.integration_tolerances.reltol_angular_rate; atol=0.0, rtol=0.0))
    @test all(isapprox.(abstol_phase.sc[1].ω, tol_args.integration_tolerances.abstol_angular_rate; atol=0.0, rtol=0.0))

    q_id = SVector{4, Float64}(0.0, 0.0, 0.0, 1.0)
    cache_lpi = PlanetFrameEphemerisCache(
        [0.0, 5.0, 10.0],
        [q_id, q_id, q_id]
    )
    @test callbacks._planet_lpi_from_cache(cache_lpi, -1.0) === nothing
    @test callbacks._planet_lpi_from_cache(cache_lpi, NaN) isa SMatrix{3, 3, Float64}
    @test callbacks._planet_lpi_from_cache(cache_lpi, 10.0) isa SMatrix{3, 3, Float64}
    cache_lpi_flip = PlanetFrameEphemerisCache(
        [0.0, 10.0],
        [q_id, -q_id]
    )
    @test callbacks._planet_lpi_from_cache(cache_lpi_flip, 5.0) isa SMatrix{3, 3, Float64}

    args_density_stats = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    p_density_stats = ODEParams(n_sats=1, args=args_density_stats)
    u_density_stats = build_initial_conditions(args_density_stats)
    withenv(
        "SPACEAGORA_GRAM_PROFILE" => "1",
        "SPACEAGORA_GRAM_TRACK_CACHE" => "off"
    ) do
        callbacks._gram_runtime_stats_reset!()
        density_cb_stats = callbacks.get_density_callback(1, args_density_stats)
        density_cb_stats.affect!((p=p_density_stats, u=u_density_stats, t=0.0))
        stats_after_density = callbacks._gram_runtime_stats_snapshot()
        @test stats_after_density.density_calls >= 1
        @test stats_after_density.direct_calls >= 1
    end

    args_thermal_branches = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=60.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    p_thermal_branches = ODEParams(n_sats=1, args=args_thermal_branches)
    u_thermal_branches = build_initial_conditions(args_thermal_branches)
    thermal_cb_branches = callbacks.get_thermal_callback(1, args_thermal_branches)

    p_thermal_branches.shared_buffers.heat_rates[1] = [1.0, 2.0, 3.0]
    p_thermal_branches.shared_buffers.densities[1] = 0.0
    p_thermal_branches.shared_buffers.temperatures[1] = 250.0
    thermal_cb_branches.affect!((p=p_thermal_branches, u=u_thermal_branches, t=0.0))
    @test length(p_thermal_branches.shared_buffers.heat_rates[1]) == length(args_thermal_branches.dynamics_model.spacecraft[1].links)

    p_thermal_branches.shared_buffers.densities[1] = 1e-6
    p_thermal_branches.shared_buffers.temperatures[1] = 250.0
    p_thermal_branches.shared_buffers.winds[1] = SVector{3, Float64}(NaN, NaN, NaN)
    thermal_cb_branches.affect!((p=p_thermal_branches, u=u_thermal_branches, t=0.0))

    p_thermal_branches.shared_buffers.winds[1] = SVector{3, Float64}(0.0, 0.0, 0.0)
    p_thermal_branches.shared_buffers.densities[1] = 1e-6
    p_thermal_branches.shared_buffers.temperatures[1] = 250.0
    args_thermal_branches.dynamics_model.spacecraft[1].links[1].α = NaN
    thermal_cb_branches.affect!((p=p_thermal_branches, u=u_thermal_branches, t=0.0))
    @test p_thermal_branches.shared_buffers.heat_rates[1][1] == 0.0

    mission_orbits = MissionConfiguration(
        mission_type=MissionOrbits,
        keplerian=true,
        number_of_orbits=1,
        mission_time=120.0,
        orientation_sim=false,
        num_steps_to_save=100
    )
    args_orbit_multi_base = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0),
            make_spacecraft(ra_alt_m=550e3, rp_alt_m=500e3, ν_deg=160.0)
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=120.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true
    )
    args_orbit_multi = SimulationConfiguration(
        file_paths=args_orbit_multi_base.file_paths,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=mission_orbits,
        environment_model=args_orbit_multi_base.environment_model,
        dynamics_model=args_orbit_multi_base.dynamics_model,
        guidance_model=args_orbit_multi_base.guidance_model,
        navigation_model=args_orbit_multi_base.navigation_model,
        control_model=args_orbit_multi_base.control_model,
        initial_time=args_orbit_multi_base.initial_time,
        integration_tolerances=args_orbit_multi_base.integration_tolerances
    )
    p_orbit_multi = ODEParams(n_sats=2, args=args_orbit_multi)
    p_orbit_multi.orbit_counter .= [2, 1]
    p_orbit_multi.is_active .= [true, true]
    orbit_cb_multi = callbacks.get_orbit_end_callback(2)
    mutable struct OrbitStopIntegrator{P, U}
        p::P
        u::U
        t::Float64
        terminated::Bool
    end
    DiffEqBase.terminate!(integrator::OrbitStopIntegrator) = begin
        integrator.terminated = true
        nothing
    end
    orbit_integrator = OrbitStopIntegrator(
        p_orbit_multi,
        build_initial_conditions(args_orbit_multi),
        0.0,
        false
    )
    orbit_cb_multi.affect!(orbit_integrator, 1)
    @test orbit_integrator.terminated == false
    p_orbit_multi.orbit_counter .= [2, 2]
    orbit_cb_multi.affect!(orbit_integrator, 1)
    @test orbit_integrator.terminated == true
end

