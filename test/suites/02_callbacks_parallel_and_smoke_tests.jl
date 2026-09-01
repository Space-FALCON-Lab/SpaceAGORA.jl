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
    # default (16), which exceeds the CI thread count. args_control uses
    # NoAtmosphereModel, which routes through the lock-free density-callback
    # source (see the source= selection in
    # src/simulation/callbacks/density_callbacks/config.jl), so the
    # lock-free-specific knob is what actually gates this decision.
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_LOCKFREE_AUTO_THREAD_MIN_BUDGET" => "2"
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
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL_ALLOW_WITH_OUTER" => "0"
    ) do
        @test density_use_threads(args_control, 8) == false
    end
    withenv(
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "auto",
        "SPACEAGORA_DENSITY_CALLBACK_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_DENSITY_CALLBACK_LOCKFREE_AUTO_THREAD_MIN_BUDGET" => "2",
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

@testset "Coverage Threaded Probe Driver" begin
    if Base.JLOptions().code_coverage == 0
        @test true
    else
        probe_script = joinpath(REPO_ROOT, "test", "probes", "coverage_threaded_probes.jl")
        cmd = `$(Base.julia_cmd()) --startup-file=no --depwarn=error --project=$(REPO_ROOT) --code-coverage=user --threads=2 $(probe_script)`
        cmd = addenv(
            cmd,
            "SPACEAGORA_WARN_DEPRECATED_CONFIG" => "0",
            "SPACEAGORA_WARN_NORMALIZE" => "0"
        )

        output = IOBuffer()
        proc = run(pipeline(ignorestatus(cmd), stdout=output, stderr=output))
        text = String(take!(output))
        if !success(proc)
            println(text)
        end

        @test success(proc)
        @test occursin("coverage_threaded_probes_ok", text)
    end
end

@testset "Multibody Parallel Policy Gates" begin
    use_threads = SimulationModel.DynamicEffectors._multibody_use_threads
    has_worker_threads = Threads.nthreads() > 1

    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
        "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1"
    ) do
        @test use_threads(64) == false
    end

    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
        "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0"
    ) do
        @test use_threads(64) == has_worker_threads
    end

    withenv(
        "SPACEAGORA_MULTIBODY_PARALLEL" => "on",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1"
    ) do
        @test use_threads(64) == has_worker_threads
    end
end

@testset "Parallel Policy Adaptive Controller" begin
    policy = SimulationModel.ParallelPolicy
    policy.reset_policy_telemetry!()

    # Forced region: budget <= 1 or num_items <= 1 means the *next* decision
    # cannot read the controller (thread_policy_decision returns use_threads
    # false and allotment 1 before consulting desire), so record_policy_observation!
    # must not move it. This block used to drive the AIMD arithmetic through
    # exactly that region -- budget 1, num_items 1/1/0 -- and assert the
    # classifications it produced. Those classifications were unreadable by
    # construction, and worse, the arithmetic there is degenerate: at budget 1
    # `allotment` collapses to 1 and `useful` to min(num_items, 1), so any
    # num_items >= 1 scores utilization 1.0 and inflates desire on every call,
    # which then leaks into a later higher-budget solve in the same process as
    # a desire that was never earned. The controller is exercised for real in
    # the live-budget block below; here the contract is that the forced region
    # is inert.
    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
        "SPACEAGORA_PARALLEL_POLICY_WINDOW" => "1",
        "SPACEAGORA_PARALLEL_POLICY_TRIM_QUANTA" => "1",
        "SPACEAGORA_PARALLEL_POLICY_DELTA" => "0.8",
        "SPACEAGORA_PARALLEL_POLICY_RHO" => "1.5",
        "SPACEAGORA_INNER_THREAD_BUDGET" => "1"
    ) do
        @test policy.use_threads_policy(
            8;
            mode=:auto,
            threshold=1,
            source=:density_callback
        ) == false

        for (items, ns) in ((1, 10), (1, 11), (0, 12))
            policy.record_policy_observation!(
                :density_callback;
                mode=:auto,
                num_items=items,
                use_threads=false,
                elapsed_ns=ns
            )
        end
        snap = policy.policy_telemetry_snapshot()
        # Raw accounting still advances -- only the controller update is skipped.
        @test snap.observations_total >= 3
        @test snap.serial_elapsed_ns_total >= 33
        @test snap.adaptation_updates_total == 0
        @test snap.quantums_total == 0
        @test snap.last_classification == "none"
        @test snap.accounted_fraction_proxy >= 0.0
        @test snap.trimmed_accounted_fraction_proxy >= 0.0
    end

    # Live budget: the controller can actually be read, so the three AIMD
    # classifications are reachable and are asserted here instead.
    #
    # `inefficient` needs utilization < delta, i.e. num_items < delta * allotment
    # with num_items >= 2 (below that the forced-region guard applies), so it
    # needs an allotment of at least 3 and therefore a budget of at least 4.
    # `efficient_deprived` needs desire > budget, which the rho ramp reaches
    # after a few satisfied windows since the desire cap is ceil(rho * budget).
    if Threads.nthreads() >= 4
        policy.reset_policy_telemetry!()
        budget = Threads.nthreads()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_WINDOW" => "1",
            "SPACEAGORA_PARALLEL_POLICY_TRIM_QUANTA" => "1",
            "SPACEAGORA_PARALLEL_POLICY_DELTA" => "0.8",
            "SPACEAGORA_PARALLEL_POLICY_RHO" => "1.5",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget)
        ) do
            saturating_items = 4 * budget
            classifications = String[]
            elapsed_total = 0
            # Saturating windows: utilization pins at 1.0, so desire ramps by rho
            # until it exceeds the budget and the allotment starts being clamped.
            for i in 1:8
                decision = policy.thread_policy_decision(
                    saturating_items;
                    mode=:auto,
                    threshold=1,
                    source=:other_source
                )
                elapsed_total += 10 + i
                policy.record_policy_observation!(
                    :other_source;
                    mode=:auto,
                    num_items=saturating_items,
                    use_threads=decision.use_threads,
                    elapsed_ns=10 + i
                )
                push!(classifications, policy.policy_telemetry_snapshot().last_classification)
            end
            @test "efficient_satisfied" in classifications
            @test "efficient_deprived" in classifications
            @test policy.policy_telemetry_snapshot().last_desire >= 2

            # Starve the same controller: two items against an allotment that has
            # ramped to at least 3 scores below delta.
            starved = policy.thread_policy_decision(
                2;
                mode=:auto,
                threshold=1,
                source=:other_source
            )
            policy.record_policy_observation!(
                :other_source;
                mode=:auto,
                num_items=2,
                use_threads=starved.use_threads,
                elapsed_ns=13
            )
            elapsed_total += 13
            snap = policy.policy_telemetry_snapshot()
            @test snap.last_classification == "inefficient"
            @test snap.last_utilization < 0.8
            @test snap.adaptation_updates_total >= 3
            @test snap.quantum_length == 1
            @test snap.trim_quanta_budget == 1
            @test snap.quantums_total >= 3
            @test snap.quantums_inefficient >= 1
            @test snap.quantums_efficient_satisfied >= 1
            @test snap.quantums_efficient_deprived >= 1
            @test snap.elapsed_ns_total >= elapsed_total
        end
    end

    if Threads.nthreads() > 1
        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_WINDOW" => "1",
            "SPACEAGORA_PARALLEL_POLICY_DELTA" => "0.8",
            "SPACEAGORA_PARALLEL_POLICY_RHO" => "1.5",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            use_history = Bool[]
            for _ in 1:8
                decision = policy.thread_policy_decision(
                    6;
                    mode=:auto,
                    threshold=1,
                    source=:other_source
                )
                push!(use_history, decision.use_threads)
                policy.record_policy_observation!(
                    :other_source;
                    mode=:auto,
                    num_items=6,
                    use_threads=decision.use_threads,
                    elapsed_ns=1
                )
            end
            @test any(use_history)
            @test use_history[end] == true
        end

        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS" => "1",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            decision = policy.thread_policy_decision(
                8;
                mode=:auto,
                threshold=2,
                source=:control_callback
            )
            @test decision.use_threads == true
            @test decision.allotment >= 2
        end

        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS" => "0",
            "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD" => "0",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            decision = policy.thread_policy_decision(
                8;
                mode=:auto,
                threshold=2,
                source=:control_callback
            )
            @test decision.use_threads == false
            @test decision.allotment == 1
        end

        policy.reset_policy_telemetry!()
        withenv(
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "1",
            "SPACEAGORA_PARALLEL_POLICY_BOOTSTRAP_THREADS" => "0",
            "SPACEAGORA_PARALLEL_POLICY_CONTROL_TAIL_GUARD" => "1",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads())
        ) do
            decision = policy.thread_policy_decision(
                8;
                mode=:auto,
                threshold=2,
                source=:control_callback
            )
            @test decision.use_threads == true
            @test decision.allotment == min(Threads.nthreads(), 8)
        end
    end

    withenv(
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_INNER_THREAD_BUDGET" => "1"
    ) do
        @test policy.use_threads_policy(4; mode=:on, threshold=1, source=:other_source) == false
        policy.record_policy_observation!(
            :other_source;
            mode=:on,
            num_items=4,
            use_threads=true,
            elapsed_ns=UInt(5)
        )
        snap = policy.policy_telemetry_snapshot()
        @test snap.threaded_elapsed_ns_total >= 5
        @test snap.other_decisions >= 1
    end

    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => "oops") do
        @test_throws ArgumentError policy.effective_inner_thread_budget()
    end
end

@testset "Parallel Policy threaded_reduce" begin
    policy = SimulationModel.ParallelPolicy
    budget = max(1, Threads.nthreads())
    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => string(budget)) do
        reduced = policy.threaded_reduce(
            16,
            budget,
            () -> MVector{2, Int}(0, 0),
            (local_acc, idx) -> begin
                local_acc[1] += idx
                local_acc[2] += 1
                return nothing
            end,
            (dest, src) -> begin
                dest[1] += src[1]
                dest[2] += src[2]
                return nothing
            end
        )
        @test reduced[1] == sum(1:16)
        @test reduced[2] == 16

        reduced_empty = policy.threaded_reduce(
            0,
            budget,
            () -> Ref(0),
            (local_acc, idx) -> begin
                local_acc[] += idx
                return nothing
            end,
            (dest, src) -> begin
                dest[] += src[]
                return nothing
            end
        )
        @test reduced_empty[] == 0
    end
end

@testset "Parallel Policy threaded_foreach_persistent" begin
    policy = SimulationModel.ParallelPolicy
    budget = max(1, Threads.nthreads())

    withenv(
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
        "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "1"
    ) do
        acc = Base.Threads.Atomic{Int}(0)
        policy.with_policy_context() do
            policy.threaded_foreach_persistent(:density_callback, 16, budget) do idx
                Base.Threads.atomic_add!(acc, idx)
            end
            policy.threaded_foreach_persistent(:density_callback, 16, budget) do idx
                Base.Threads.atomic_add!(acc, idx)
            end
        end
        @test acc[] == 2 * sum(1:16)
    end

    withenv(
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
        "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "0"
    ) do
        acc = Base.Threads.Atomic{Int}(0)
        policy.threaded_foreach_persistent(:density_callback, 8, budget) do idx
            Base.Threads.atomic_add!(acc, idx)
        end
        @test acc[] == sum(1:8)
    end

    withenv(
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(budget),
        "SPACEAGORA_CALLBACK_PERSISTENT_WORKERS" => "1"
    ) do
        err = try
            policy.with_policy_context() do
                policy.threaded_foreach_persistent(:thermal_callback, 8, budget) do idx
                    if idx == 3
                        error("threaded_foreach_persistent_probe")
                    end
                end
            end
            nothing
        catch e
            e
        end
        @test err !== nothing
    end
end

@testset "Aerodynamic Helper Branch Coverage" begin
    dynamic_effectors = SimulationModel.DynamicEffectors

    @test dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", false) == false
    withenv("SPACEAGORA_TEST_BOOL_PARSE" => "yes") do
        @test dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", false) == true
    end
    withenv("SPACEAGORA_TEST_BOOL_PARSE" => "off") do
        @test dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", true) == false
    end
    withenv("SPACEAGORA_TEST_BOOL_PARSE" => "invalid") do
        @test_throws ArgumentError dynamic_effectors._parse_bool_env("SPACEAGORA_TEST_BOOL_PARSE", false)
    end

    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "off") do
        @test dynamic_effectors._multibody_parallel_mode() == :off
    end
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "on") do
        @test dynamic_effectors._multibody_parallel_mode() == :on
    end
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "auto") do
        @test dynamic_effectors._multibody_parallel_mode() == :auto
    end
    withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "invalid") do
        @test_throws ArgumentError dynamic_effectors._multibody_parallel_mode()
    end

    withenv("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "4") do
        @test dynamic_effectors._multibody_thread_threshold() == 4
    end
    withenv("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "0") do
        @test dynamic_effectors._multibody_thread_threshold() == 1
    end
    withenv("SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "oops") do
        @test_throws ArgumentError dynamic_effectors._multibody_thread_threshold()
    end
    withenv("SPACEAGORA_MULTIBODY_MAX_THREADS" => "3") do
        @test dynamic_effectors._multibody_max_threads() == 3
    end
    withenv("SPACEAGORA_MULTIBODY_MAX_THREADS" => "oops") do
        @test_throws ArgumentError dynamic_effectors._multibody_max_threads()
    end

    withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1") do
        @test dynamic_effectors._multibody_outer_parallel_hint() == true
    end
    withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0") do
        @test dynamic_effectors._multibody_outer_parallel_hint() == false
    end
    withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "oops") do
        @test_throws ArgumentError dynamic_effectors._multibody_outer_parallel_hint()
    end

    @test dynamic_effectors._multibody_use_threads(1) == false
    if Threads.nthreads() > 1
        withenv("SPACEAGORA_MULTIBODY_PARALLEL" => "on") do
            @test dynamic_effectors._multibody_use_threads(64) == true
        end
        withenv(
            "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
            "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1",
            "SPACEAGORA_MULTIBODY_PARALLEL_ALLOW_WITH_OUTER" => "0"
        ) do
            @test dynamic_effectors._multibody_use_threads(64) == false
        end
        withenv(
            "SPACEAGORA_MULTIBODY_PARALLEL" => "auto",
            "SPACEAGORA_MULTIBODY_THREAD_THRESHOLD" => "2",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
            "SPACEAGORA_MULTIBODY_PARALLEL_HEAVY_ONLY" => "1"
        ) do
            @test dynamic_effectors._multibody_use_threads(64; heavy_work=false) == false
        end
    end

    @test dynamic_effectors._threadid_capacity() >= Threads.maxthreadid()

    body_a = Link{0}(root=true)
    body_b = Link{0}(root=false)
    body_a.net_force .= SVector{3, Float64}(1.0, 2.0, 3.0)
    body_a.net_torque .= SVector{3, Float64}(4.0, 5.0, 6.0)
    body_b.net_force .= SVector{3, Float64}(-0.5, 0.0, 0.5)
    body_b.net_torque .= SVector{3, Float64}(1.0, -1.0, 0.0)
    force_sum, torque_sum = dynamic_effectors.collect_and_reset_link_wrenches!([body_a, body_b])

    @test force_sum == SVector{3, Float64}(0.5, 2.0, 3.5)
    @test torque_sum == SVector{3, Float64}(5.0, 4.0, 6.0)
    @test body_a.net_force == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test body_a.net_torque == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test body_b.net_force == SVector{3, Float64}(0.0, 0.0, 0.0)
    @test body_b.net_torque == SVector{3, Float64}(0.0, 0.0, 0.0)
end




@testset "Typed Planet Constructors + Topography Workspace" begin
    mars = Mars("", SPICE_PATH)
    venus = Venus("", SPICE_PATH)
    titan = Titan("", SPICE_PATH)

    @test mars.name == "Mars"
    @test venus.name == "Venus"
    @test titan.name == "Titan"
    @test mars.μ > 0.0
    @test venus.μ > 0.0
    @test titan.μ > 0.0

    earth = Earth("", SPICE_PATH)
    moon = Moon("", SPICE_PATH)
    earth_radii_km = lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
        bodvrd("EARTH", "RADII")
    end
    moon_radii_km = lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
        bodvrd("MOON", "RADII")
    end
    @test isapprox(earth.Rp_e, earth_radii_km[1] * 1e3; atol=1e-6, rtol=0.0)
    @test isapprox(earth.Rp_p, earth_radii_km[3] * 1e3; atol=1e-6, rtol=0.0)
    @test isapprox(moon.Rp_e, moon_radii_km[1] * 1e3; atol=1e-6, rtol=0.0)
    @test isapprox(moon.Rp_p, moon_radii_km[3] * 1e3; atol=1e-6, rtol=0.0)

    mktempdir() do tmp
        topo_file = joinpath(tmp, "topo_harmonics.csv")
        write(topo_file, "degree,order,C,S\n0,0,1.0,0.0\n1,0,0.1,0.0\n1,1,0.05,0.02\n")
        SimulationModel.Planets.TopographyHarmonicsWorkspace!(topo_file, earth)
        @test size(earth.topography_workspace.Clm) == (2, 2)
        @test size(earth.topography_workspace.Slm) == (2, 2)
        @test isapprox(earth.topography_workspace.Clm[2, 2], 0.05; atol=0.0, rtol=0.0)
        @test isapprox(earth.topography_workspace.Slm[2, 2], 0.02; atol=0.0, rtol=0.0)
    end
end

@testset "Run Simulation State Isolation" begin
    sc = make_spacecraft(ra_alt_m=500e3, rp_alt_m=450e3, ν_deg=170.0)
    args = build_config(
        spacecraft=sc,
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=45.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredGravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=10.0)
    )

    l_pi_before = args.environment_model.planet.L_PI
    @test_nowarn run_simulation(args)
    @test args.environment_model.planet.L_PI == l_pi_before

    args_no_isolation = deepcopy(args)
    @test_nowarn run_simulation(args_no_isolation; isolate_state=false)

    if Threads.nthreads() >= 2
        args_parallel = deepcopy(args)
        t1 = Threads.@spawn run_simulation(args_parallel)
        t2 = Threads.@spawn run_simulation(args_parallel)
        @test_nowarn fetch(t1)
        @test_nowarn fetch(t2)
    end
end

# Test doubles for the cost probe. Defined at file scope because the probe
# dispatches on their concrete types.
#
# _ProbeWorkEffector stands in for a user-defined effector with no
# effector_cost_terms method: the whole point of the probe is that such an
# effector needs to implement nothing to be routable. _ProbeThrowingEffector
# stands in for one that fails when evaluated, which must degrade the model to
# "cannot predict" rather than take the solve down with it.
struct _ProbeWorkEffector <: SimulationModel.AbstractForceTorqueModel
    iterations::Int
end

function SimulationModel.calcForceTorque(
    model::_ProbeWorkEffector, x::AbstractVector{Float64}, param::ODEParams, i::Int64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    acc = 0.0
    @inbounds for j in 1:model.iterations
        acc = muladd(x[1], 1.0000001, acc)
    end
    return (SVector{3, Float64}(acc * 0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
end

struct _ProbeThrowingEffector <: SimulationModel.AbstractForceTorqueModel end

function SimulationModel.calcForceTorque(
    ::_ProbeThrowingEffector, ::AbstractVector{Float64}, ::ODEParams, ::Int64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    error("deliberate probe failure")
end

# Survives the warm-up RHS, then throws once the probe starts hammering it.
# This is the failure mode the per-effector try/catch exists for: an effector
# that works in the RHS but not under isolated re-evaluation. The
# always-throwing double above cannot exercise it, because it takes the warm-up
# down first and nothing is ever probed.
const _PROBE_FLAKY_CALLS = Ref(0)

struct _ProbeFlakyEffector <: SimulationModel.AbstractForceTorqueModel
    succeed_for::Int
end

function SimulationModel.calcForceTorque(
    model::_ProbeFlakyEffector, ::AbstractVector{Float64}, ::ODEParams, ::Int64
)::Tuple{SVector{3, Float64}, SVector{3, Float64}}
    _PROBE_FLAKY_CALLS[] += 1
    _PROBE_FLAKY_CALLS[] > model.succeed_for && error("deliberate probe failure")
    return (SVector{3, Float64}(0.0, 0.0, 0.0), SVector{3, Float64}(0.0, 0.0, 0.0))
end

@testset "Effector Cost Probe" begin
    PC = SimulationModel.ParallelCost
    harmonics_file = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    harmonics = GravitationalHarmonicsModel(4, 4, harmonics_file, EARTH)
    cheap = _ProbeWorkEffector(64)
    dear  = _ProbeWorkEffector(16_384)

    args_probe = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3 + 10e3 * i, rp_alt_m=480e3 + 5e3 * i, ν_deg=120.0 + 5.0 * i)
            for i in 1:4
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(harmonics, cheap, dear),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_probe = ODEParams(n_sats=4, args=args_probe)
    _initialize_heat_rate_buffers!(p_probe)
    _initialize_harmonics_workspace_buffers!(p_probe)
    SimulationEngine._initialize_density_model_instances!(p_probe)
    SimulationEngine._initialize_density_cache_buffers!(p_probe)
    SimulationEngine._initialize_gram_isolated_pool_buffers!(p_probe)
    _initialize_aero_workspace_buffers!(p_probe)
    _initialize_nbody_workspace_buffers!(p_probe)
    u_probe = build_initial_conditions(args_probe)

    probes = SimulationEngine.probe_effector_costs!(p_probe, u_probe, args_probe)

    @test length(probes) == 3
    @test all(haskey(probes, i) for i in 1:3)

    # Harmonics declares analytic terms, so it is probed only to validate.
    @test probes[1].declared
    @test !probes[2].declared
    @test !probes[3].declared

    # A user-defined effector with no trait method is measured, and the
    # measurement is real work rather than a folded-away constant.
    floor_ns = PC.timing_noise_floor()
    for i in 2:3
        @test isfinite(probes[i].ns_per_sat)
        @test probes[i].ns_per_sat > floor_ns
        @test isfinite(probes[i].reference_ratio)
        @test probes[i].reference_ratio > 0.0
    end

    # The probe must order effectors by cost -- that is the only property the
    # model actually consumes from it. 256x the inner iterations must read as
    # materially dearer; the bound is loose because SIMD width and loop overhead
    # make the true factor machine-specific.
    @test probes[3].ns_per_sat > 4 * probes[2].ns_per_sat

    # probe_ns_map feeds constellation_work_counts and must omit declared
    # effectors, whose cost comes from their counts -- including them would
    # double-count.
    ns_map = PC.probe_ns_map(probes)
    @test !haskey(ns_map, 1)
    @test haskey(ns_map, 2) && haskey(ns_map, 3)

    # Second call, post-compilation: the probe runs at every solve setup, so
    # its own cost has to stay negligible against the sweep it replaces (which
    # evaluated the whole RHS 15 times for each of up to 10 candidate plans).
    probe_elapsed_s = @elapsed SimulationEngine.probe_effector_costs!(p_probe, u_probe, args_probe)
    @info "Effector cost probe overhead" effectors=3 seconds=probe_elapsed_s
    @test probe_elapsed_s < 0.5

    counts = PC.constellation_work_counts(args_probe, 4; probe=ns_map)
    @test counts.unknown_effectors == 2
    @test counts.probe_ns ≈ ns_map[2] + ns_map[3]
    @test counts.simd_terms > 0.0                     # harmonics contributed analytically
    # Harmonics is pre-pass handled; the two doubles are not, so only they
    # produce queue nodes.
    @test counts.queue_nodes == 4 * 2
    @test counts.in_domain
end

@testset "Effector Cost Probe: failures do not break the run" begin
    PC = SimulationModel.ParallelCost
    args_bad = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3 + 10e3 * i, rp_alt_m=480e3 + 5e3 * i, ν_deg=120.0 + 5.0 * i)
            for i in 1:2
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(InverseSquaredJ2GravityModel(), _ProbeThrowingEffector()),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_bad = ODEParams(n_sats=2, args=args_bad)
    _initialize_heat_rate_buffers!(p_bad)
    _initialize_harmonics_workspace_buffers!(p_bad)
    SimulationEngine._initialize_density_model_instances!(p_bad)
    SimulationEngine._initialize_density_cache_buffers!(p_bad)
    SimulationEngine._initialize_gram_isolated_pool_buffers!(p_bad)
    _initialize_aero_workspace_buffers!(p_bad)
    _initialize_nbody_workspace_buffers!(p_bad)
    u_bad = build_initial_conditions(args_bad)

    # An effector that cannot be evaluated at all takes the warm-up RHS down,
    # so nothing is probed. That must return empty rather than throw: a solve
    # whose RHS fails once was not going to run regardless, and the probe is
    # not the place to surface it.
    probes = SimulationEngine.probe_effector_costs!(p_bad, u_bad, args_bad)
    @test isempty(probes)

    # And the model must decline rather than treat the unmeasured effector as
    # free -- the whole point of routing an unknown effector to the probe is
    # lost if a failed probe silently becomes zero cost.
    counts = PC.constellation_work_counts(args_bad, 2; probe=PC.probe_ns_map(probes))
    @test counts.unknown_effectors == 1
    @test !counts.in_domain
end

@testset "Effector Cost Probe: isolated-evaluation failure is contained" begin
    PC = SimulationModel.ParallelCost
    _PROBE_FLAKY_CALLS[] = 0
    args_flaky = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3 + 10e3 * i, rp_alt_m=480e3 + 5e3 * i, ν_deg=120.0 + 5.0 * i)
            for i in 1:2
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        # succeed_for covers the warm-up RHS (one call per active satellite),
        # then the probe's own repeated calls trip the failure.
        dynamic_effectors=(InverseSquaredJ2GravityModel(), _ProbeFlakyEffector(2)),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_flaky = ODEParams(n_sats=2, args=args_flaky)
    _initialize_heat_rate_buffers!(p_flaky)
    _initialize_harmonics_workspace_buffers!(p_flaky)
    SimulationEngine._initialize_density_model_instances!(p_flaky)
    SimulationEngine._initialize_density_cache_buffers!(p_flaky)
    SimulationEngine._initialize_gram_isolated_pool_buffers!(p_flaky)
    _initialize_aero_workspace_buffers!(p_flaky)
    _initialize_nbody_workspace_buffers!(p_flaky)
    u_flaky = build_initial_conditions(args_flaky)

    probes = SimulationEngine.probe_effector_costs!(p_flaky, u_flaky, args_flaky)
    @test haskey(probes, 2)
    @test isnan(probes[2].ns_per_sat)
    @test isnan(probes[2].reference_ratio)

    # NaN must not leak into the counts as a number.
    counts = PC.constellation_work_counts(args_flaky, 2; probe=PC.probe_ns_map(probes))
    @test !counts.in_domain
    @test counts.probe_ns == 0.0
    @test !isnan(counts.probe_ns)
end

@testset "Declaration mismatch is warned, not enforced" begin
    PC = SimulationModel.ParallelCost
    # Agreement inside the factor, disagreement outside it.
    @test PC.validate_declaration(100.0, 150.0)
    @test PC.validate_declaration(150.0, 100.0)
    @test !PC.validate_declaration(100.0, 10.0)
    @test !PC.validate_declaration(10.0, 100.0)
    # An unmeasurable probe cannot invalidate a declaration.
    @test PC.validate_declaration(NaN, 100.0)

    probes = Dict(
        1 => PC.EffectorProbe(100.0, 5.0, true),    # declared, agrees
        2 => PC.EffectorProbe(100.0, 5.0, true),    # declared, disagrees 10x
        3 => PC.EffectorProbe(100.0, 5.0, false),   # undeclared, never checked
    )
    predicted = Dict(1 => 120.0, 2 => 10.0, 3 => 1.0)
    mismatches = SimulationEngine.warn_on_declaration_mismatch(probes, predicted)
    @test mismatches == 1
end

@testset "Successive-halving sweep" begin
    # The sweep now eliminates candidates in rounds instead of giving every one
    # the full timed block. What must hold is that it still returns a plan drawn
    # from the candidate set, that it survives a budget of 1, and that a
    # candidate set it cannot evaluate degrades to "no override" rather than to
    # an exception -- the calibration is an optimisation and must never be able
    # to take a solve down.
    Threads.nthreads() > 1 || return

    harmonics_file = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    args_sweep = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3 + 10e3 * i, rp_alt_m=480e3 + 5e3 * i,
                            ν_deg=120.0 + 5.0 * i)
            for i in 1:16
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(GravitationalHarmonicsModel(8, 8, harmonics_file, EARTH),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_sweep = ODEParams(n_sats=16, args=args_sweep)
    _initialize_heat_rate_buffers!(p_sweep)
    _initialize_harmonics_workspace_buffers!(p_sweep)
    SimulationEngine._initialize_density_model_instances!(p_sweep)
    SimulationEngine._initialize_density_cache_buffers!(p_sweep)
    SimulationEngine._initialize_gram_isolated_pool_buffers!(p_sweep)
    _initialize_aero_workspace_buffers!(p_sweep)
    _initialize_nbody_workspace_buffers!(p_sweep)
    SimulationEngine._initialize_runtime_env_config!(p_sweep)
    u_sweep = build_initial_conditions(args_sweep)

    effs = args_sweep.dynamics_model.dynamic_effectors
    withenv("SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER" => "1") do
        SimulationEngine._initialize_runtime_env_config!(p_sweep)
        candidates = SimulationEngine._rhs_plan_candidates(p_sweep, effs)
        @test length(candidates) > 1

        plan, elapsed = SimulationEngine._run_rhs_sweep!(p_sweep, u_sweep, effs, false)
        @test plan !== nothing
        @test elapsed > 0.0
        @test any(c -> c.mode === plan.mode && c.allotment == plan.allotment &&
                       c.scheduler === plan.scheduler, candidates)

        # A one-call budget must still terminate and still choose. Successive
        # halving doubles its repetition count each round, so a budget that is
        # already reached on the first round has to break immediately rather
        # than loop.
        withenv("SPACEAGORA_RHS_CALIBRATE_N_TIMED" => "1",
                "SPACEAGORA_RHS_CALIBRATE_N_WARMUP" => "1") do
            plan_min, elapsed_min = SimulationEngine._run_rhs_sweep!(p_sweep, u_sweep, effs, false)
            @test plan_min !== nothing
            @test elapsed_min > 0.0
        end
    end
    # The sweep must leave no override behind: it restores the Ref on every path,
    # including the error path, or a candidate would leak into the real solve.
    @test p_sweep.shared_buffers.rhs_plan_override[] === nothing
end

@testset "Flat-queue reduction: parallel path is bit-identical to serial" begin
    # The cross-worker partials reduction was serial O(N*W) and now splits over
    # the worker pool. Each worker owns a disjoint satellite range and still
    # accumulates worker_id 1..W in order for every satellite it owns, so the
    # floating-point summation order per satellite is UNCHANGED and the result
    # must be bit-identical -- not merely close. Anything else would mean the
    # optimisation silently altered trajectories.
    Threads.nthreads() > 1 || return

    harmonics_file = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    args_red = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3 + 10e3 * (i % 7), rp_alt_m=480e3 + 5e3 * (i % 5),
                            ν_deg=(120.0 + 7.0 * i) % 360.0)
            for i in 1:160
        ],
        # Harmonics is pre-pass handled; the aerodynamic effector is flat-queue
        # only, which is what makes the partials reduction actually run.
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(GravitationalHarmonicsModel(8, 8, harmonics_file, EARTH),
                           AerodynamicCoefficientConstant()),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_red = ODEParams(n_sats=160, args=args_red)
    _initialize_heat_rate_buffers!(p_red)
    _initialize_harmonics_workspace_buffers!(p_red)
    SimulationEngine._initialize_density_model_instances!(p_red)
    SimulationEngine._initialize_density_cache_buffers!(p_red)
    SimulationEngine._initialize_gram_isolated_pool_buffers!(p_red)
    _initialize_aero_workspace_buffers!(p_red)
    _initialize_nbody_workspace_buffers!(p_red)
    SimulationEngine._initialize_runtime_env_config!(p_red)
    u_red = build_initial_conditions(args_red)
    du_red = zero(u_red)

    function derivative_at(allot::Int, min_sats::String)
        withenv("SPACEAGORA_RHS_REDUCTION_MIN_SATS_PER_WORKER" => min_sats) do
            SimulationEngine._initialize_runtime_env_config!(p_red)
            p_red.shared_buffers.rhs_plan_override[] =
                SimulationEngine._make_calib_flat_plan(allot, :static)
            try
                SimulationEngine.spacecraft_dynamics!(du_red, u_red, p_red, 0.0)
                return [du_red.sc[i][j] for i in 1:160 for j in 4:6]
            finally
                p_red.shared_buffers.rhs_plan_override[] = nothing
            end
        end
    end

    # "999999" forces the serial path regardless of satellite count.
    reference = derivative_at(1, "999999")
    @test any(!iszero, reference)

    for allot in (2, 4, Threads.nthreads())
        serial_result = derivative_at(allot, "999999")
        parallel_result = derivative_at(allot, "1")
        @test serial_result == reference
        @test parallel_result == reference
    end
end

@testset "Light-work effector short-circuit keeps the cost model fed" begin
    # _dynamic_effector_thread_decision returns before consulting the policy when
    # heavy-only rejects light work. That return MUST keep policy_applied true:
    # at the call site it gates both the elapsed-time measurement and
    # _update_effector_cost_model!, and that cost model is what computes
    # heavy_work in the first place. Returning false freezes the estimate, so
    # work classified light once can never be reclassified, and the guard latches
    # permanently on its first observation.
    args_light = build_config_multi(
        spacecraft=[make_spacecraft(ra_alt_m=500e3 + 10e3 * i, rp_alt_m=480e3 + 5e3 * i,
                                    ν_deg=120.0 + 5.0 * i) for i in 1:2],
        density_model=ExponentialAtmosphereModel(EARTH),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        # Two cheap effectors, both thread-safe so the decision actually reaches
        # the heavy-work guard. AerodynamicCoefficientConstant would not: it is
        # not declared thread-safe, so _dynamic_effectors_parallel_supported
        # rejects the set earlier and the guard under test never runs.
        dynamic_effectors=(InverseSquaredGravityModel(), InverseSquaredJ2GravityModel()),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_light = ODEParams(n_sats=2, args=args_light)
    _initialize_heat_rate_buffers!(p_light)
    _initialize_harmonics_workspace_buffers!(p_light)
    SimulationEngine._initialize_density_model_instances!(p_light)
    SimulationEngine._initialize_density_cache_buffers!(p_light)
    SimulationEngine._initialize_gram_isolated_pool_buffers!(p_light)
    _initialize_aero_workspace_buffers!(p_light)
    _initialize_nbody_workspace_buffers!(p_light)
    SimulationEngine._initialize_runtime_env_config!(p_light)
    u_light = build_initial_conditions(args_light)
    du_light = zero(u_light)

    withenv("SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
            "SPACEAGORA_EFFECTOR_PARALLEL_HEAVY_ONLY" => "1") do
        SimulationEngine._initialize_runtime_env_config!(p_light)
        decision = SimulationEngine._dynamic_effector_thread_decision(
            args_light, p_light, args_light.dynamics_model.dynamic_effectors, 2
        )
        # Light work must not thread, but must still be timed.
        @test decision.use_threads == false
        @test decision.allotment == 1
        @test decision.policy_applied == true

        # And the cost model must actually accumulate samples across RHS calls,
        # which is the observable consequence.
        before = p_light.shared_buffers.effector_cost_samples[]
        for _ in 1:5
            SimulationEngine.spacecraft_dynamics!(du_light, u_light, p_light, 0.0)
        end
        @test p_light.shared_buffers.effector_cost_samples[] > before
    end
end

@testset "RHS Plan Calibration" begin
    # ── Env var parsing ──────────────────────────────────────────────────────
    withenv("SPACEAGORA_RHS_CALIBRATE" => "off") do
        @test SimulationEngine._rhs_calibration_mode() == :off
    end
    withenv("SPACEAGORA_RHS_CALIBRATE" => "auto") do
        @test SimulationEngine._rhs_calibration_mode() == :auto
    end
    withenv("SPACEAGORA_RHS_CALIBRATE" => "force") do
        @test SimulationEngine._rhs_calibration_mode() == :force
    end
    withenv("SPACEAGORA_RHS_CALIBRATE" => "invalid_xyz") do
        @test_throws ArgumentError SimulationEngine._rhs_calibration_mode()
    end

    withenv("SPACEAGORA_RHS_CALIBRATE_N_WARMUP" => "3") do
        @test SimulationEngine._rhs_calibrate_n_warmup() == 3
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_N_WARMUP" => "0") do
        @test SimulationEngine._rhs_calibrate_n_warmup() == 1
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_N_TIMED" => "7") do
        @test SimulationEngine._rhs_calibrate_n_timed() == 7
    end
    withenv("SPACEAGORA_RHS_CALIBRATE_N_TIMED" => "0") do
        @test SimulationEngine._rhs_calibrate_n_timed() == 1
    end

    # ── Machine label stability ───────────────────────────────────────────────
    # The label is cached in a module-level Ref after the first call; test idempotency only.
    label1 = SimulationEngine._calib_machine_label()
    label2 = SimulationEngine._calib_machine_label()
    @test label1 == label2
    @test !isempty(label1)

    # ── Plan constructors return expected field layout ────────────────────────
    sat_batch_plan = SimulationEngine._make_calib_satellite_batch_plan()
    @test sat_batch_plan.mode == :satellite_batch
    @test sat_batch_plan.allotment == 1
    @test sat_batch_plan.dominant_axis == :satellite
    @test sat_batch_plan.policy_applied == true

    flat_plan = SimulationEngine._make_calib_flat_plan(4)
    @test flat_plan.mode == :flat_constellation_effector_queue
    @test flat_plan.allotment == 4
    @test flat_plan.dominant_axis == :flat_effector
    @test flat_plan.policy_applied == true

    flat_plan_floor = SimulationEngine._make_calib_flat_plan(0)
    @test flat_plan_floor.allotment == 1

    # ── In-memory store / lookup round-trip ──────────────────────────────────
    test_sig = "v3|machine=test_gate|budget=8|sats=2_4|effs=1|harm=1"
    SimulationEngine._rhs_calib_store!(test_sig, sat_batch_plan, 1.5e6)
    retrieved = SimulationEngine._rhs_calib_lookup(test_sig)
    @test retrieved !== nothing
    @test retrieved.mode == :satellite_batch

    SimulationEngine._rhs_calib_store!(test_sig, flat_plan, 0.9e6)
    retrieved_flat = SimulationEngine._rhs_calib_lookup(test_sig)
    @test retrieved_flat !== nothing
    @test retrieved_flat.mode == :flat_constellation_effector_queue
    @test retrieved_flat.allotment == 4

    # Scheduler survives store/lookup: it is a swept axis now, so a cache entry
    # that dropped it would silently hand the solve back to the env var.
    for sched in (:static, :dynamic)
        sig_sched = "v3|machine=test_gate_$(sched)|budget=8|sats=2_4|effs=1|harm=1"
        SimulationEngine._rhs_calib_store!(
            sig_sched, SimulationEngine._make_calib_flat_plan(4, sched), 0.8e6
        )
        round_tripped = SimulationEngine._rhs_calib_lookup(sig_sched)
        @test round_tripped !== nothing
        @test round_tripped.scheduler == sched
    end

    # ── TOML persistence round-trip (temp directory) ─────────────────────────
    mktempdir() do tmp
        calib_path = joinpath(tmp, "calib_test.toml")
        withenv("SPACEAGORA_RHS_CALIBRATION_PATH" => calib_path) do
            # Store and save to disk
            sig_disk = "v3|machine=disktest|budget=4|sats=2_4|effs=1|harm=1"
            SimulationEngine._rhs_calib_store!(sig_disk, flat_plan, 2.0e6)
            SimulationEngine._rhs_calib_save!()
            @test isfile(calib_path)

            parsed = TOML.parsefile(calib_path)
            @test haskey(parsed, "calibrations")
            @test parsed["schema_version"] == 1
            rows = parsed["calibrations"]
            @test any(r -> get(r, "signature", "") == sig_disk, rows)
            disk_row = only(filter(r -> get(r, "signature", "") == sig_disk, rows))
            @test get(disk_row, "scheduler", "") == String(flat_plan.scheduler)
        end
    end

    # ── Calibration guard: SPACEAGORA_RHS_CALIBRATE=off skips entirely ───────
    harmonics_file = joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv")
    harmonics_model = GravitationalHarmonicsModel(4, 4, harmonics_file, EARTH)

    args_calib = build_config_multi(
        spacecraft=[
            make_spacecraft(ra_alt_m=500e3 + 10e3 * i, rp_alt_m=480e3 + 5e3 * i, ν_deg=120.0 + 5.0 * i)
            for i in 1:4
        ],
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(harmonics_model,),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_calib = ODEParams(n_sats=4, args=args_calib)
    _initialize_heat_rate_buffers!(p_calib)
    _initialize_harmonics_workspace_buffers!(p_calib)
    SimulationEngine._initialize_density_model_instances!(p_calib)
    SimulationEngine._initialize_density_cache_buffers!(p_calib)
    SimulationEngine._initialize_gram_isolated_pool_buffers!(p_calib)
    _initialize_aero_workspace_buffers!(p_calib)
    _initialize_nbody_workspace_buffers!(p_calib)
    u_calib = build_initial_conditions(args_calib)

    withenv("SPACEAGORA_RHS_CALIBRATE" => "off") do
        p_calib.shared_buffers.rhs_plan_override[] = nothing
        SimulationEngine._calibrate_rhs_plan_if_needed!(p_calib, u_calib, args_calib)
        @test p_calib.shared_buffers.rhs_plan_override[] === nothing
    end

    # ── Calibration guard: single satellite skips calibration ────────────────
    args_single_sat = build_config(
        spacecraft=make_spacecraft(ra_alt_m=500e3, rp_alt_m=480e3, ν_deg=120.0),
        density_model=NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=30.0,
        EI_km=120.0,
        dynamic_effectors=(harmonics_model,),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false)
    )
    p_single = ODEParams(n_sats=1, args=args_single_sat)
    _initialize_heat_rate_buffers!(p_single)
    _initialize_harmonics_workspace_buffers!(p_single)
    SimulationEngine._initialize_density_model_instances!(p_single)
    SimulationEngine._initialize_density_cache_buffers!(p_single)
    SimulationEngine._initialize_gram_isolated_pool_buffers!(p_single)
    _initialize_aero_workspace_buffers!(p_single)
    _initialize_nbody_workspace_buffers!(p_single)
    u_single = build_initial_conditions(args_single_sat)
    withenv("SPACEAGORA_RHS_CALIBRATE" => "auto") do
        p_single.shared_buffers.rhs_plan_override[] = nothing
        SimulationEngine._calibrate_rhs_plan_if_needed!(p_single, u_single, args_single_sat)
        @test p_single.shared_buffers.rhs_plan_override[] === nothing
    end

    # ── Signature stability ───────────────────────────────────────────────────
    sig_a = SimulationEngine._rhs_calib_signature(
        p_calib, args_calib.dynamics_model.dynamic_effectors,
        args_calib.environment_model.density_model
    )
    sig_b = SimulationEngine._rhs_calib_signature(
        p_calib, args_calib.dynamics_model.dynamic_effectors,
        args_calib.environment_model.density_model
    )
    @test sig_a == sig_b
    # v4: bumped again for the effector-identity and outer-split terms. v3
    # entries key on effector COUNT and are blind to the outer split, so
    # replaying one is exactly the collision those terms exist to prevent --
    # measured: gravity_4096sat_l50 calibrates to "retain the heuristic" 16/16
    # with no outer split and to a pinned flat plan 16/16 with one, and under v3
    # those two shared a single cache key.
    #
    # (v3 itself was the no-regret floor bump: a cached entry short-circuits the
    # sweep, and the floor only runs inside a sweep, so pre-floor entries would
    # keep pinning plans the floor exists to reject, indefinitely.)
    #
    # v5: bumped for the density-model term -- v4 keys were blind to the density
    # model, so a cached entry could replay a GRAM plan the sweep had refused to
    # rank.
    @test startswith(sig_a, "v5|machine=")
    @test occursin("|sats=", sig_a)
    @test occursin("|effs=1|", sig_a)
    @test occursin("|harm=1", sig_a)
    # Effector identity, not just count: two shapes with the same effector count
    # must not share a cache key.
    @test occursin("|eff=", sig_a)
    @test occursin("|outer=", sig_a)

    # ── Multi-sat + harmonics calibration (requires worker threads) ───────────
    if Threads.nthreads() > 1
        # Pin the SIMD batch floor so the 4-sat fixture yields >= 2 viable
        # workers and the sweep produces an override regardless of the
        # default SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER (4).
        # Point the calibration path at a throwaway file so the force-mode
        # save cannot persist this suite's synthetic signatures into the
        # machine-default calibration state (which poisons later suite runs).
        withenv(
            "SPACEAGORA_RHS_CALIBRATE" => "force",
            "SPACEAGORA_RHS_CALIBRATE_N_WARMUP" => "1",
            "SPACEAGORA_RHS_CALIBRATE_N_TIMED" => "2",
            "SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER" => "1",
            "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(mktempdir(), "calib_force_gate.toml"),
            "SPACEAGORA_RHS_CALIBRATE_OVERRIDE_MARGIN" => "0.0"
        ) do
            # End to end, EITHER outcome is correct and which one occurs is a
            # timing race: at margin zero the heuristic is still retained
            # whenever it is faster by any amount, which on a four-satellite
            # fixture it often is. An earlier version of this asserted an
            # override was installed and passed only by luck of the draw. So
            # assert the outcome is well formed rather than which way it went.
            p_calib.shared_buffers.rhs_plan_override[] = nothing
            SimulationEngine._calibrate_rhs_plan_if_needed!(p_calib, u_calib, args_calib)
            override = p_calib.shared_buffers.rhs_plan_override[]
            if override !== nothing
                @test override.mode ∈ (:satellite_batch, :flat_constellation_effector_queue)
                @test override.policy_applied == true
            else
                @test override === nothing
            end
        end

        # That the sweep CAN produce a plan is tested directly, without the
        # floor, so it does not depend on which side wins a timing race:
        # passing no `args` means no heuristic candidate and therefore no floor.
        withenv(
            "SPACEAGORA_RHS_CALIBRATE_N_WARMUP" => "1",
            "SPACEAGORA_RHS_CALIBRATE_N_TIMED" => "2",
            "SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER" => "1"
        ) do
            plan, elapsed = SimulationEngine._run_rhs_sweep!(
                p_calib, u_calib, args_calib.dynamics_model.dynamic_effectors, false
            )
            @test plan !== nothing
            @test plan.mode ∈ (:satellite_batch, :flat_constellation_effector_queue)
            @test elapsed > 0.0
            @test p_calib.shared_buffers.rhs_plan_override[] === nothing
        end

        # ── No-regret floor ──────────────────────────────────────────────────
        # The sweep may only displace the runtime heuristic when it beats it by
        # more than the sweep's own resolution. Without this the sweep could pin
        # a plan worse than doing nothing and no measurement could detect it,
        # because the heuristic's answer was never one of its candidates:
        # measured on gravity_4096sat_l50_vacuum_1hr and heavy_1024sat_l50_6hr
        # at eight threads, where policy_decisions_total is zero for every
        # profile and the whole 10-12% gap was a pinned plan losing to the
        # heuristic it was never compared against.
        withenv(
            "SPACEAGORA_RHS_CALIBRATE" => "force",
            "SPACEAGORA_RHS_CALIBRATE_N_WARMUP" => "1",
            "SPACEAGORA_RHS_CALIBRATE_N_TIMED" => "2",
            "SPACEAGORA_HARMONICS_BATCH_MIN_SATS_PER_WORKER" => "1",
            "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(mktempdir(), "calib_floor_gate.toml"),
            # An unreachable margin: nothing can clear it, so the heuristic must
            # always be retained regardless of what the sweep measured.
            "SPACEAGORA_RHS_CALIBRATE_OVERRIDE_MARGIN" => "0.89"
        ) do
            p_calib.shared_buffers.rhs_plan_override[] = nothing
            SimulationEngine._calibrate_rhs_plan_if_needed!(p_calib, u_calib, args_calib)
            @test p_calib.shared_buffers.rhs_plan_override[] === nothing
        end

        @test SimulationEngine._rhs_calibrate_override_margin() ≈ 0.10
        withenv("SPACEAGORA_RHS_CALIBRATE_OVERRIDE_MARGIN" => "0.25") do
            @test SimulationEngine._rhs_calibrate_override_margin() ≈ 0.25
        end
        # Clamped into [0, 0.9] so a mistyped value cannot disable calibration
        # entirely or make every sweep override.
        withenv("SPACEAGORA_RHS_CALIBRATE_OVERRIDE_MARGIN" => "5.0") do
            @test SimulationEngine._rhs_calibrate_override_margin() ≈ 0.9
        end
    end
end

function _no_gram_ring_sat(planet_ng, id::Int, raan_deg::Float64)
    return make_single_link_spacecraft(
        ra_alt_m=700e3, rp_alt_m=650e3, i_deg=45.0, ω_deg=0.0, Ω_deg=raan_deg, ν_deg=10.0,
        planet=planet_ng, id=id, m=100.0, ref_area=2.0
    )
end

@testset "Constellation Ensemble Campaign" begin
    planet_ng = SimulationModel.make_no_gram_planet(:earth)

    sats = [
        _no_gram_ring_sat(planet_ng, 11, 0.0),
        _no_gram_ring_sat(planet_ng, 22, 40.0),
        _no_gram_ring_sat(planet_ng, 33, 80.0),
    ]
    results_dir = mktempdir()
    cfg = build_config_multi(
        spacecraft=sats,
        density_model=SimulationModel.NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=600.0,
        EI_km=120.0,
        dynamic_effectors=(SimulationModel.InverseSquaredJ2GravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(
            results=true,
            verbose=false,
            generate_plots=false,
            normalize=false,
            results_directory=results_dir
        ),
        ephemerides_model=SimulationModel.SimpleEphemeridesModel(),
        planet=planet_ng
    )

    ensemble_threads = min(2, Threads.nthreads())
    res = SimulationCampaigns.run_constellation_ensemble(cfg; threads=ensemble_threads, return_solution=true)
    @test res isa SimulationCampaigns.MonteCarloResult
    @test length(res.samples) == 3
    @test isempty(res.failed)
    @test [s.seed for s in res.samples] == [1, 2, 3]
    @test all(s -> s.value !== nothing, res.samples)

    # Each ensemble member writes into its own per-satellite results directory.
    @test isdir(joinpath(results_dir, "sat_1_id_11"))
    @test isdir(joinpath(results_dir, "sat_2_id_22"))
    @test isdir(joinpath(results_dir, "sat_3_id_33"))

    # An ensemble member must reproduce a direct single-satellite propagation.
    single_cfg = SimulationCampaigns._ensemble_member_configuration(cfg, sats[3], "solo_direct")
    sol_direct = run_simulation(single_cfg; return_solution=true)
    u_direct = collect(sol_direct.u[end])
    u_member = collect(res.samples[3].value.u[end])
    @test maximum(abs.(u_direct .- u_member)) <= 1e-9

    # Uncoupled guard: GNC effectors are rejected unless explicitly allowed.
    cfg_ctrl = build_config_multi(
        spacecraft=sats,
        density_model=SimulationModel.NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=600.0,
        EI_km=120.0,
        dynamic_effectors=(SimulationModel.InverseSquaredJ2GravityModel(),),
        control_effectors=(make_base_thruster_model(thrust=0.0),),
        control_rates=[1.0],
        keplerian=true,
        ephemerides_model=SimulationModel.SimpleEphemeridesModel(),
        planet=planet_ng
    )
    @test_throws ArgumentError SimulationCampaigns.run_constellation_ensemble(cfg_ctrl)
    @test SimulationCampaigns._validate_ensemble_uncoupled(cfg_ctrl, true) === nothing

    # Empty constellation is rejected.
    cfg_empty = build_config_multi(
        spacecraft=SpacecraftModel[],
        density_model=SimulationModel.NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=600.0,
        EI_km=120.0,
        dynamic_effectors=(SimulationModel.InverseSquaredJ2GravityModel(),),
        keplerian=true,
        ephemerides_model=SimulationModel.SimpleEphemeridesModel(),
        planet=planet_ng
    )
    @test_throws ArgumentError SimulationCampaigns.run_constellation_ensemble(cfg_empty)

    # Member splits alias the parent environment (cheap split); parallel isolation
    # therefore relies on the per-worker deepcopy, which must sever mutable model
    # state such as the planet object.
    member_split = SimulationCampaigns._ensemble_member_configuration(cfg, sats[1], "alias_probe")
    @test member_split.environment_model === cfg.environment_model
    member_isolated = deepcopy(member_split)
    @test member_isolated.environment_model !== cfg.environment_model
    @test member_isolated.environment_model.planet !== cfg.environment_model.planet

    # An explicit isolate_state=false must not break the parallel path: the runner
    # overrides it with its own per-worker copy.
    res_noiso = SimulationCampaigns.run_constellation_ensemble(
        cfg; threads=ensemble_threads, return_solution=true, isolate_state=false
    )
    @test isempty(res_noiso.failed)
    u_noiso = collect(res_noiso.samples[3].value.u[end])
    @test maximum(abs.(u_direct .- u_noiso)) <= 1e-9

    # Checkpoint path splitting: any checkpoint interaction gets a per-member path.
    settings_ck_default = SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false,
        results_directory="outdir", checkpoint_enabled=true
    )
    split_ck_default = SimulationCampaigns._ensemble_member_settings(settings_ck_default, "sat_1_id_11")
    @test split_ck_default.results_directory == joinpath("outdir", "sat_1_id_11")

    settings_resume_explicit = SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false,
        results_directory="outdir", checkpoint_directory="ckdir", resume_from_checkpoint=true
    )
    split_resume = SimulationCampaigns._ensemble_member_settings(settings_resume_explicit, "sat_2_id_22")
    @test split_resume.checkpoint_directory == joinpath("ckdir", "sat_2_id_22")
    @test split_resume.results_directory == "outdir"

    settings_plain = SimulationSettings(
        results=false, verbose=false, generate_plots=false, normalize=false
    )
    @test SimulationCampaigns._ensemble_member_settings(settings_plain, "sat_3_id_33") === settings_plain
end

@testset "GRAM Lock Scope" begin
    EMt = SimulationModel.EnvironmentModels

    withenv("SPACEAGORA_GRAM_LOCK_SCOPE" => nothing) do
        @test EMt._gram_lock_scope() === :global
    end
    withenv("SPACEAGORA_GRAM_LOCK_SCOPE" => "global") do
        @test EMt._gram_lock_scope() === :global
    end
    for token in ("model", "per_model", "per-model", "instance", "MODEL")
        withenv("SPACEAGORA_GRAM_LOCK_SCOPE" => token) do
            @test EMt._gram_lock_scope() === :model
        end
    end
    withenv("SPACEAGORA_GRAM_LOCK_SCOPE" => "bogus") do
        @test_throws ArgumentError EMt._gram_lock_scope()
    end

    # Every wrapper construction carries its own instance lock so distinct
    # models never contend under SPACEAGORA_GRAM_LOCK_SCOPE=model.
    m1 = EMt.GRAMAtmosphereModel(:dummy_core_a)
    m2 = EMt.GRAMAtmosphereModel(:dummy_core_b)
    @test m1.instance_lock isa ReentrantLock
    @test m1.instance_lock !== m2.instance_lock
    @test m1.core === :dummy_core_a
    @test :instance_lock in propertynames(m1)
end

@testset "Adaptive Campaign Routing" begin
    # ── Features from an explicit campaign shape ──────────────────────────────
    feats = SimulationCampaigns.campaign_route_features(
        samples=12, n_sats=2, density_family="exponential", mission_time_s=1200.0
    )
    @test feats isa ParallelProfiles.OuterRouteFeatures
    @test feats.category == "montecarlo"
    @test feats.montecarlo_samples == 12
    @test feats.n_sats == 2
    sig = ParallelProfiles.outer_route_signature(feats)
    @test occursin("cat=montecarlo", sig)
    @test occursin("dens=exp", sig)
    @test occursin("mission=short", sig)
    @test_throws ArgumentError SimulationCampaigns.campaign_route_features(samples=-1)
    @test_throws ArgumentError SimulationCampaigns.campaign_route_features(samples=4, n_sats=0)

    # ── Features derived from a SimulationConfiguration ───────────────────────
    planet_ng = SimulationModel.make_no_gram_planet(:earth)

    adaptive_sats = [_no_gram_ring_sat(planet_ng, 1, 0.0), _no_gram_ring_sat(planet_ng, 2, 120.0)]
    cfg = build_config_multi(
        spacecraft=adaptive_sats,
        density_model=SimulationModel.NoAtmosphereModel(),
        orientation_sim=false,
        mission_time=300.0,
        EI_km=120.0,
        dynamic_effectors=(SimulationModel.InverseSquaredJ2GravityModel(),),
        keplerian=true,
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        ephemerides_model=SimulationModel.SimpleEphemeridesModel(),
        planet=planet_ng
    )
    cfg_feats = SimulationCampaigns.campaign_route_features(cfg; samples=2, n_sats=1)
    @test cfg_feats.density_family == "none"
    @test cfg_feats.mission_time_s == 300.0
    @test cfg_feats.n_sats == 1
    @test cfg_feats.montecarlo_samples == 2
    @test cfg_feats.dynamic_effector_count == 1
    @test cfg_feats.has_control == false
    @test cfg_feats.orientation_on == false
    full_feats = SimulationCampaigns.campaign_route_features(cfg; samples=2)
    @test full_feats.n_sats == 2

    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => nothing, "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => nothing) do
        # ── threads=:auto Monte Carlo: route selection, env split, feedback ───
        state = ParallelProfiles.OuterRouteState()
        n_seeds = 6
        inner_budget_seen = Vector{String}(undef, n_seeds)
        outer_active_seen = Vector{String}(undef, n_seeds)
        runner = seed -> begin
            inner_budget_seen[seed] = get(ENV, "SPACEAGORA_INNER_THREAD_BUDGET", "")
            outer_active_seen[seed] = get(ENV, "SPACEAGORA_OUTER_PARALLEL_ACTIVE", "")
            seed * 10
        end
        res1 = SimulationCampaigns.run_monte_carlo(runner, 1:n_seeds; threads=:auto, route_state=state)
        @test res1 isa SimulationCampaigns.MonteCarloResult
        @test length(res1.successful) == n_seeds
        @test [s.value for s in res1.samples] == [10, 20, 30, 40, 50, 60]
        # The adaptive env overrides are scoped to the campaign run.
        @test isempty(get(ENV, "SPACEAGORA_INNER_THREAD_BUDGET", ""))

        expected_workers = Threads.nthreads() > 1 ? min(n_seeds, Threads.nthreads()) : 1
        @test res1.threads == expected_workers
        if Threads.nthreads() > 1
            expected_budget = string(max(1, fld(Threads.nthreads(), expected_workers)))
            @test all(==(expected_budget), inner_budget_seen)
            @test all(==("1"), outer_active_seen)
        end

        auto_sig = ParallelProfiles.outer_route_signature(
            SimulationCampaigns.campaign_route_features(samples=n_seeds)
        )
        route1 = Threads.nthreads() > 1 ? :threads : :none
        snap1 = ParallelProfiles.outer_route_stats_snapshot(state, auto_sig)
        @test haskey(snap1, route1)
        @test snap1[route1].samples == n_seeds
        @test snap1[route1].success_rate == 1.0
        # Feedback stores amortized campaign wall time per sample.
        @test isapprox(snap1[route1].mean_s, res1.elapsed_s / n_seeds; rtol=1e-6)

        # A second identical campaign explores the under-sampled serial route,
        # so repeated campaigns accumulate stats for every feasible allocation.
        if Threads.nthreads() > 1
            res2 = SimulationCampaigns.run_monte_carlo(runner, 1:n_seeds; threads=:auto, route_state=state)
            @test res2.threads == 1
            snap2 = ParallelProfiles.outer_route_stats_snapshot(state, auto_sig)
            @test snap2[:none].samples == n_seeds
            @test snap2[:threads].samples == n_seeds
        end

        # Caller-provided route features are patched with the actual sample count.
        if Threads.nthreads() > 1
            override_state = ParallelProfiles.OuterRouteState()
            stale = SimulationCampaigns.campaign_route_features(samples=0, density_family="exponential")
            res_o = SimulationCampaigns.run_monte_carlo(
                x -> x, 1:n_seeds;
                threads=:auto, route_features=stale, route_state=override_state
            )
            @test res_o.threads == expected_workers
            snap_o = ParallelProfiles.outer_route_stats_snapshot(
                override_state,
                ParallelProfiles.outer_route_signature(
                    SimulationCampaigns.campaign_route_features(samples=n_seeds, density_family="exponential")
                )
            )
            @test sum(info.samples for info in values(snap_o)) == n_seeds
        end

        # Failures are recorded and lower the stored success rate.
        fail_state = ParallelProfiles.OuterRouteState()
        res_fail = SimulationCampaigns.run_monte_carlo(1:4; threads=:auto, route_state=fail_state) do seed
            seed == 2 && error("seed two failed")
            seed
        end
        @test length(res_fail.failed) == 1
        snap_fail = ParallelProfiles.outer_route_stats_snapshot(
            fail_state,
            ParallelProfiles.outer_route_signature(SimulationCampaigns.campaign_route_features(samples=4))
        )
        @test sum(info.samples for info in values(snap_fail)) == 4
        @test any(info.success_rate == 0.75 for info in values(snap_fail))

        # Nested adaptive campaigns yield to an enclosing outer split: serial
        # execution, no route feedback (contended timings must not poison the
        # shared statistics).
        nested_state = ParallelProfiles.OuterRouteState()
        res_nested = withenv("SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "1") do
            SimulationCampaigns.run_monte_carlo(x -> x, 1:n_seeds; threads=:auto, route_state=nested_state)
        end
        @test res_nested.threads == 1
        @test length(res_nested.successful) == n_seeds
        @test isempty(nested_state.history)

        # Caller-provided features with a non-montecarlo category are normalized
        # before routing and recording: other categories' default-route rules can
        # answer :process for shapes whose candidates are [:none, :threads],
        # which the selector would clamp to a serial default on larger machines.
        if Threads.nthreads() > 1
            cat_state = ParallelProfiles.OuterRouteState()
            odd_features = ParallelProfiles.OuterRouteFeatures(
                category="scaling",
                n_sats=2,
                mission_time_s=3600.0,
                has_control=true,
                montecarlo_samples=0
            )
            res_cat = SimulationCampaigns.run_monte_carlo(
                x -> x, 1:n_seeds;
                threads=:auto, route_features=odd_features, route_state=cat_state
            )
            @test res_cat.threads == expected_workers
            normalized_sig = ParallelProfiles.outer_route_signature(
                ParallelProfiles.OuterRouteFeatures(
                    category="montecarlo",
                    n_sats=2,
                    mission_time_s=3600.0,
                    has_control=true,
                    montecarlo_samples=n_seeds
                )
            )
            snap_cat = ParallelProfiles.outer_route_stats_snapshot(cat_state, normalized_sig)
            @test sum(info.samples for info in values(snap_cat)) == n_seeds
        end

        # ── threads=:auto constellation ensemble ──────────────────────────────
        ens_state = ParallelProfiles.OuterRouteState()
        res_ens = SimulationCampaigns.run_constellation_ensemble(
            cfg; threads=:auto, route_state=ens_state, return_solution=true
        )
        @test res_ens isa SimulationCampaigns.MonteCarloResult
        @test length(res_ens.samples) == 2
        @test isempty(res_ens.failed)
        @test all(s -> s.value !== nothing, res_ens.samples)
        ens_sig = ParallelProfiles.outer_route_signature(
            SimulationCampaigns.campaign_route_features(cfg; samples=2, n_sats=1)
        )
        ens_snap = ParallelProfiles.outer_route_stats_snapshot(ens_state, ens_sig)
        @test sum(info.samples for info in values(ens_snap)) == 2
        @test all(info.success_rate == 1.0 for info in values(ens_snap))
    end

    # ── Argument validation ────────────────────────────────────────────────────
    guard_state = ParallelProfiles.OuterRouteState()
    @test_throws ArgumentError SimulationCampaigns.run_monte_carlo(identity, 1:2; threads=:fast)
    @test_throws ArgumentError SimulationCampaigns.run_monte_carlo(identity, 1:2; threads=1, route_state=guard_state)
    @test_throws ArgumentError SimulationCampaigns.run_constellation_ensemble(cfg; threads=:fast)
    @test_throws ArgumentError SimulationCampaigns.run_constellation_ensemble(cfg; threads=1, route_state=guard_state)

    # Empty seed sets return an empty result without consulting the bandit.
    empty_state = ParallelProfiles.OuterRouteState()
    res_empty = SimulationCampaigns.run_monte_carlo(identity, Int[]; threads=:auto, route_state=empty_state)
    @test isempty(res_empty.samples)
    @test isempty(empty_state.history)

    # The process-global state is a stable, inspectable OuterRouteState.
    @test SimulationCampaigns.campaign_outer_route_state() isa ParallelProfiles.OuterRouteState
    @test SimulationCampaigns.campaign_outer_route_state() === SimulationCampaigns.campaign_outer_route_state()
end
