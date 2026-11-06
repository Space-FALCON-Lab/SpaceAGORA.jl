module EventOutputTests

using Test, StaticArrays, LinearAlgebra, OrdinaryDiffEq, DiffEqCallbacks, SpaceAGORA
const SM = SpaceAGORA.SimulationModel
const CB = SM.SimulationCallbacks

function battery_run(; rate=1.0)
    planet = make_no_gram_planet(:earth)
    radius = planet.Rp_e + 700e3
    root = Link(root=true, m=100.0)
    initial = CartesianInitialCondition(SVector(radius, 0.0, 0.0),
        SVector(0.0, sqrt(planet.μ / radius), 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SpacecraftModel(root=root, links=[root], inertia_tensor=root.inertia,
        initial_condition=initial, id=1)
    add_power!(sc, root, Battery(capacity_wh=1.0, initial_soc=0.5))
    add_power!(sc, root, PowerLoad(power_w=0.5 * 3600 / 6.1))
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, generate_plots=false, verbose=false),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, mission_time=12.0,
            orientation_sim=false, data_rate=rate),
        environment_model=make_no_gram_environment(planet=planet),
        dynamics_model=DynamicsModel([sc], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]), navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]), initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-8, dt_max_orbit=0.7),
        solver_config=SM.SolverConfig(solver_mode=:tsit5))
end

expected_battery(t) = max(0.5 * (1.0 - t / 6.1), 0.0)

@testset "Output before an event uses the incoming trajectory" begin
    # Results-only runs do not retain dense solver history. The same correction
    # must work there, as well as when a caller also requests the solution.
    for full_solution in (false, true)
        result = run_simulation(battery_run(); return_results=true, return_solution=full_solution)
        @test result.table.time == collect(0.0:1.0:12.0)
        @test result.table.sc1_battery_wh ≈ expected_battery.(result.table.time) atol=1e-10
        @test result.table.sc1_battery_wh[7] ≈ expected_battery(6.0) atol=1e-10
    end
    # An event after several reinitializations must drain the CURRENT inner
    # table saver. Noncommensurate checkpoints leave nominal5s/6s samples
    # inside the4.6..6.9s segment that contains depletion at6.1s.
    checkpoint_args = SM.SimConfig._with_configuration(battery_run();
        simulation_settings=SimulationSettings(results=false, generate_plots=false,
            verbose=false, normalize=false, results_directory=mktempdir(),
            checkpoint_enabled=true, checkpoint_interval_s=2.3))
    checkpoint_result = run_simulation(checkpoint_args; return_results=true, return_solution=false)
    @test all(t -> t in checkpoint_result.table.time, 0.0:1.0:12.0)
    @test checkpoint_result.table.sc1_battery_wh ≈
        expected_battery.(checkpoint_result.table.time) atol=1e-10
    six = findfirst(==(6.0), checkpoint_result.table.time)
    @test checkpoint_result.table.sc1_battery_wh[six] ≈ expected_battery(6.0) atol=1e-10

    # The preallocated recorder shares the correction even when it is the
    # only saver. Its custom fields still use the same SaveField getters.
    args = battery_run()
    fields = [CB.SaveField(:held_energy, (u, t, i) -> Float64(u.sc[1].battery_wh))]
    recorder = CB.TrajectoryRecorder(args; save_fields=fields)
    run_simulation(args; extra_callbacks=(CB.get_trajectory_recorder_callback(recorder),))
    @test collect(CB.trajectory_times(recorder)) == collect(0.0:1.0:12.0)
    @test collect(CB.trajectory_field(recorder, :held_energy)) ≈
        expected_battery.(CB.trajectory_times(recorder)) atol=1e-10

    # The event itself keeps its post-effect sample. Initializers/finalizers
    # must receive the original callback, whose affect object they own.
    for vector_event in (false, true)
        initialized, finalized = Ref(false), Ref(false)
        jump!(i) = (i.u[1] += 10.0)
        jump_vector!(i, events) = (events[1] == 0 || jump!(i); nothing)
        effect = vector_event ? jump_vector! : jump!
        init(cb, u, t, i) = (initialized[] = cb.affect! === effect; nothing)
        finish(cb, u, t, i) = (finalized[] = cb.affect! === effect; nothing)
        condition!(out, u, t, i) = (out[1] = t - 0.5; nothing)
        event = vector_event ? VectorContinuousCallback(condition!, effect, 1;
            initialize=init, finalize=finish, save_positions=(false, false)) :
            ContinuousCallback((u, t, i) -> t - 0.5, effect, nothing;
                initialize=init, finalize=finish, save_positions=(false, false))
        saved = SavedValues(Float64, Float64)
        saver = SavingCallback((u, t, i) -> u[1], saved; saveat=0.25, save_everystep=false)
        callbacks = CB._save_before_continuous_events(CallbackSet(event, saver), (saver.affect!,))
        rhs!(du, u, p, t) = (du[1] = 1.0)
        sol = solve(ODEProblem(rhs!, [0.0], (0.0, 1.0)), Tsit5(); callback=callbacks,
            dt=0.7, abstol=1e-12, reltol=1e-12)
        @test initialized[] && finalized[]
        @test saved.t == collect(0.0:0.25:1.0)
        @test saved.saveval ≈ [0.0, 0.25, 10.5, 10.75, 11.0] atol=1e-10
        @test sol.u[end][1] ≈ 11.0 atol=1e-10
    end
end

end
