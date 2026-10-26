module CoupledMultirateTests

using Test, LinearAlgebra, StaticArrays, SpaceAGORA, OrdinaryDiffEq, DiffEqCallbacks, Serialization
using Arrow, DataFrames
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PLANET = make_no_gram_planet(:earth)

# Equal-mass rigid spacecraft, with a fixed data rate or a strong momentum
# link. The data case has an analytic transfer and transmitter battery drain.
function fixture(mode; duration=10.0, cap=0.5, initial_bits=1e6, kind=:data,
    orientation=false, sensor=false,
    settings=SimulationSettings(results=false, generate_plots=false))
    radius = PLANET.Rp_e + 700e3
    speed = sqrt(PLANET.μ / radius)
    spacecraft = map(1:2) do index
        bus = Link(root=true, m=100.0, dims=MVector(1.0, 1.0, 1.0))
        theta = (index - 1) * 1000 / radius
        initial = CartesianInitialCondition(
            SVector(radius * cos(theta), radius * sin(theta), 0.0),
            SVector(-speed * sin(theta), speed * cos(theta), 0.0);
            q=SVector(0.0, 0.0, 0.0, 1.0), ang_vel=SVector(0.0, 0.0, 0.01))
        sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=initial, id=index, n_terminal=kind === :laser ? 1 : 0)
        if kind === :data
            add_actuator!(sc, bus, RadioTerminal(boresight=[0.0, 1.0, 0.0],
                max_power_w=5.0, gain_dbi=6.0, efficiency=0.5))
        end
        if kind === :mounted
            add_actuator!(sc, bus, LaserTerminal(boresight=[0.0, index == 1 ? 1.0 : -1.0, 0.0],
                location=[1.0, 0.0, 0.0], max_power_w=50.0, efficiency=0.5))
        end
        if kind === :power
            terminal = index == 1 ?
                LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=100.0, efficiency=0.5) :
                PowerReceiver(normal=[0.0, -1.0, 0.0], radius_m=0.25, efficiency=0.4)
            add_actuator!(sc, bus, terminal)
        end
        sensor && add_sensor!(sc, bus, RateGyro(period_s=0.5))
        add_power!(sc, bus, Battery(capacity_wh=10.0, initial_soc=0.5))
        add_data!(sc, bus, DataStore(capacity_bits=1e7,
            initial_bits=index == 1 ? initial_bits : 0.0))
        sc
    end
    model = InterlinkModel()
    link = kind === :data ?
        DataLink(transmit_power_w=1.0, budget=FixedRateBudget(rate_bps=1000.0), range_m=2000e3) :
        kind === :power ? PowerBeamLink(power_w=80.0, range_m=50e3) :
        kind === :mounted ? LaserBeamLink(power_w=40.0, magnification=100.0,
            range_m=2000e3, recoil=:mounted) :
        LaserMomentumLink(power_w=1e8, magnification=100.0, range_m=2000e3)
    register_link!(model, spacecraft, (1, 1), (2, 1), link)
    return SimulationConfiguration(
        simulation_settings=settings,
        mission_configuration=MissionConfiguration(mission_type=MissionTime,
            mission_time=duration, orientation_sim=orientation, data_rate=1.0),
        environment_model=make_no_gram_environment(planet=PLANET),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]),
        navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]), initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10,
            abstol_orbit=1e-8, dt_max_orbit=cap),
        solver_config=SolverConfig(solver_mode=mode, multirate_slow_dt_s=cap),
        interlink_model=model,
        scheduling_policy_model=kind === :data ? DataThroughputPolicy() :
            kind === :power ? PowerDeliveryPolicy() : SemimajorAxisRatePolicy(targets=2))
end

struct Interrupted <: Exception end

withenv("SPACEAGORA_RHS_CALIBRATE"=>"off", "SPACEAGORA_RHS_IDENTIFY"=>"0",
    "SPACEAGORA_RHS_EXECUTION_MODE"=>"serial", "SPACEAGORA_VISUALIZATION"=>"0") do
@testset "Coupled MRI integration" begin
    @test SE._parse_solver_mode_sym("multirate_mri") === :multirate_mri

    @testset "Complete partition and reference trajectory: $kind" for kind in (:data, :laser, :power, :mounted)
        reference=run_simulation(fixture(:dp8; kind, orientation=true); return_results=true,return_solution=true)
        p=reference.solution.prob.p
        u=copy(reference.solution.u[end])
        for inactive in (false,true)
            p.is_active .= true
            inactive && (p.is_active[2]=false)
            whole=zero(u); fast=fill!(similar(u),NaN); slow=fill!(similar(u),NaN)
            SE.spacecraft_dynamics!(whole,u,p,10.0)
            SE.spacecraft_dynamics_mri_fast!(fast,u,p,10.0)
            SE.spacecraft_dynamics_mri_slow!(slow,u,p,10.0)
            @test all(isfinite,fast) && all(isfinite,slow)
            @test whole ≈ fast+slow rtol=1e-12 atol=1e-12
            if kind === :mounted && !inactive
                @test norm(fast.sc[1].ω) > 1e-8
                @test fast.sc[1].ω ≈ whole.sc[1].ω rtol=1e-12 atol=1e-12
            end
        end
        result=run_simulation(fixture(:multirate_mri; kind, orientation=true);return_results=true,return_solution=true)
        @test result.solution.t[end]==10.0
        @test result.solution.u[end] ≈ reference.solution.u[end] rtol=1e-9 atol=1e-6
        history=result.configuration.interlink_model.history
        @test all(diff([x.time for x in history]).>0)
        @test length(history)==result.solution.stats.naccept+1
        if kind===:data
            @test result.solution.u[end].sc[1].data_sent_bits ≈ 10000.0 atol=1e-6
            @test result.solution.u[end].sc[2].data_received_bits ≈ 10000.0 atol=1e-6
            @test result.solution.u[end].sc[1].battery_wh ≈ 5-2*10/3600 atol=1e-10
        elseif kind === :mounted
            @test norm(result.solution.u[end].sc[1].ω - SVector(0.0, 0.0, 0.01)) > 1e-8
            @test result.solution.u[end].sc[1].ω ≈ reference.solution.u[end].sc[1].ω rtol=1e-9 atol=1e-12
            @test norm(result.solution.u[end].sc[1].interlink_dv) > 0
        elseif kind === :power
            @test result.solution.u[end].sc[1].battery_wh ≈ 5 - 160 * 10 / 3600 atol=1e-10
            @test result.solution.u[end].sc[2].battery_wh ≈ 5 + 32 * 10 / 3600 atol=1e-10
        else
            @test norm(result.solution.u[end].sc[1].interlink_dv)>0
            @test result.solution.u[end].sc[2].interlink_delta_sma>0
        end
    end

    @testset "An empty producer stays connected under backlog scheduling" begin
        finals = []
        for mode in (:dp8, :multirate_mri)
            args = fixture(mode; initial_bits=0.0)
            sender = args.dynamics_model.spacecraft[1]
            add_data!(sender, sender.root, DataSource(rate_bps=100.0))
            args = SM.SimConfig._with_configuration(args;
                scheduling_policy_model=DataBacklogPolicy(horizon_s=5.0))
            result = run_simulation(args; return_results=true, return_solution=true)
            final = result.solution.u[end]
            history = result.configuration.interlink_model.history
            active = [((1, 1), (2, 1))]
            @test result.solution.t[end] == 10.0
            @test length(history) > 1
            @test first(history).active == active
            @test all(sample -> sample.active == active, history)
            @test final.sc[1].data_stored_bits ≈ 0.0 atol=1e-7
            @test final.sc[1].data_produced_bits ≈ 1000.0 atol=1e-7
            @test final.sc[1].data_sent_bits ≈ 1000.0 atol=1e-7
            @test final.sc[2].data_received_bits ≈ 1000.0 atol=1e-7
            @test final.sc[2].data_stored_bits ≈ 1000.0 atol=1e-7
            @test final.sc[1].data_lost_bits ≈ 0.0 atol=1e-7
            @test final.sc[2].data_lost_bits ≈ 0.0 atol=1e-7
            push!(finals, final)
        end
        @test finals[1] ≈ finals[2] rtol=1e-9 atol=1e-6
    end

    @testset "One physical event and periodic callback timeline" begin
        ticks=Float64[]; events=Float64[]
        fast!(du,u,p,t)=(du[1]=-10u[1]+10u[2];du[2]=0.0;nothing)
        slow!(du,u,p,t)=(du[1]=0.0;du[2]=-u[2];nothing)
        exact(t)=[10/9*exp(-t)-1/9*exp(-10t),exp(-t)]
        periodic=PeriodicCallback(i->push!(ticks,i.t),0.2;initial_affect=false,final_affect=false,
            save_positions=(false,false))
        event=ContinuousCallback((u,t,i)->t-0.53,i->(push!(events,i.t);i.u .*=0.5),nothing;
            save_positions=(false,false))
        problem=SplitODEProblem(fast!,slow!,ones(2),(0.0,1.0);callback=CallbackSet(event,periodic))
        config=SolverConfig(solver_mode=:multirate_mri,multirate_slow_dt_s=0.02)
        args=(integration_tolerances=IntegrationTolerances(dt_max_orbit=0.02),)
        sol,meta=SE._solve_with_solver_policy(problem,config,args,1e-10,1e-12)
        @test ticks ≈ [0.2,0.4,0.6,0.8] atol=1e-12
        @test length(events)==1
        @test only(events) ≈ 0.53 atol=1e-12
        @test sol.u[end] ≈ 0.5exact(1.0) atol=1e-8
        @test all(diff(sol.t).>0)
        @test meta.solver=="Multirate(MRIGARKERK45a)" && !meta.fallback_used
    end

    @testset "An empty-store event preserves incoming output" begin
        for full_solution in (false,true)
            result=run_simulation(fixture(:multirate_mri;initial_bits=6100.0,cap=0.7);
                return_results=true,return_solution=full_solution)
            @test result.table.time==collect(0.0:1.0:10.0)
            @test result.table.sc1_data_stored_bits ≈ max.(6100.0 .- 1000result.table.time,0.0) atol=1e-5
            @test result.table.sc1_data_sent_bits ≈ min.(1000result.table.time,6100.0) atol=1e-5
            @test result.table.sc2_data_received_bits ≈ result.table.sc1_data_sent_bits atol=1e-5
        end
    end

    @testset "Ground-station counters belong to the fast component" begin
        args = fixture(:multirate_mri)
        ephemerides = args.environment_model.ephemerides_model
        epoch = SM.EphemeridesModels.ephemerides_epoch_seconds(ephemerides, args.initial_time)
        frame = SM.EphemeridesModels.ephemerides_frame_lpi(ephemerides, PLANET, epoch)
        initial = SE.build_initial_conditions(args)
        fixed = frame * initial.sc[1].pos
        station = GroundStation(name="overhead", lat_deg=asind(fixed[3] / norm(fixed)),
            lon_deg=rad2deg(atan(fixed[2], fixed[1])), min_elevation_deg=0.0,
            antennas=[RadioTerminal(boresight=[0.0, 0.0, 1.0], max_power_w=5.0,
                gain_dbi=6.0, efficiency=0.5)])
        model = InterlinkModel()
        add_ground_station!(model, station)
        register_link!(model, args.dynamics_model.spacecraft, (1, 1), (-1, 1),
            DataLink(transmit_power_w=1.0, budget=FixedRateBudget(rate_bps=1000.0), range_m=2000e3))
        args = SM.SimConfig._with_configuration(args; interlink_model=model)
        result = run_simulation(args; return_results=true, return_solution=true)
        final = result.solution.u[end]
        @test final.gs[1].received_bits ≈ 10000.0 atol=1e-6
        @test final.gs[1].received_bits ≈ final.sc[1].data_sent_bits atol=1e-6
        p = result.solution.prob.p
        whole, fast, slow = zero(final), zero(final), zero(final)
        SE.spacecraft_dynamics!(whole, final, p, 10.0)
        SE.spacecraft_dynamics_mri_fast!(fast, final, p, 10.0)
        SE.spacecraft_dynamics_mri_slow!(slow, final, p, 10.0)
        @test whole ≈ fast + slow rtol=1e-12 atol=1e-12
        @test slow.gs[1].received_bits == slow.gs[1].sent_bits == 0.0
    end

    @testset "Mounted sensor sampling" begin
        result=run_simulation(fixture(:multirate_mri;orientation=true,sensor=true);
            return_results=true,return_solution=true)
        for sc in result.configuration.dynamics_model.spacecraft
            gyro=only(SM.mounted_sensors(sc)).sensor
            @test sensor_measurement(gyro) ≈ SVector(0.0,0.0,0.01) atol=1e-9
        end
    end

    @testset "Environment-selected mode preserves its cap" begin
        withenv("SPACEAGORA_SOLVER_MODE"=>"multirate_mri",
            "SPACEAGORA_MULTIRATE_SLOW_DT_S"=>"0.2") do
            args = SM.SimConfig._with_configuration(fixture(:tsit5); solver_config=nothing)
            result = run_simulation(args; return_results=true, return_solution=true)
            @test result.configuration.solver_config.solver_mode === :multirate_mri
            @test result.configuration.solver_config.multirate_slow_dt_s == 0.2
            @test maximum(diff(result.solution.t)) <= 0.2 + 1e-12
            @test result.solution.u[end].sc[1].data_sent_bits ≈ 10000.0 atol=1e-6
        end
    end

    @testset "Phase changes retain the configured macro-step cap" begin
        args = fixture(:multirate_mri; cap=0.3)
        p = SM.ODEParams(n_sats=2, args=args)
        reltol, abstol = SE._build_solver_tolerances(SE.build_initial_conditions(args), args)
        for inside in (false, true)
            p.shared_buffers.in_atmosphere .= inside
            cap, _, _ = SM.SimulationCallbacks._active_phase_solver_settings(p, reltol, abstol)
            @test cap == 0.3
        end
        tighter = SM.SimConfig._with_configuration(args;
            integration_tolerances=IntegrationTolerances(dt_max_orbit=1.0, dt_max_atmosphere=0.1))
        p = SM.ODEParams(n_sats=2, args=tighter)
        p.shared_buffers.in_atmosphere .= true
        cap, _, _ = SM.SimulationCallbacks._active_phase_solver_settings(p, reltol, abstol)
        @test cap == 0.1
    end

    @testset "MRI streamed outputs survive interruption" begin
        mktempdir() do directory
            settings(path; resume=false) = SimulationSettings(results=true, generate_plots=false,
                results_layout=:long, save_csv=false, results_directory=path,
                checkpoint_enabled=true, checkpoint_interval_s=3.0,
                checkpoint_directory=path, resume_from_checkpoint=resume)
            reference = run_simulation(fixture(:multirate_mri;
                settings=settings(joinpath(directory, "reference"))); return_results=true)
            resumed_dir = joinpath(directory, "resumed")
            stop = DiscreteCallback((u,t,i)->t>6.2, i->throw(Interrupted()); save_positions=(false,false))
            @test_throws Interrupted run_simulation(fixture(:multirate_mri;
                settings=settings(resumed_dir)); extra_callbacks=(stop,))
            resumed = run_simulation(fixture(:multirate_mri;
                settings=settings(resumed_dir; resume=true)); return_results=true)
            @test resumed.table.time == reference.table.time
            @test resumed.table.spacecraft == reference.table.spacecraft
            for column in (:data_sent_bits, :data_received_bits, :data_stored_bits)
                @test resumed.table[!, column] ≈ reference.table[!, column] rtol=1e-10 atol=1e-7
            end
            sender = resumed.table[resumed.table.spacecraft .== 1, :]
            @test last(sender.time) == 10.0
            @test last(sender.data_sent_bits) ≈ 10000.0 atol=1e-6
            prefix = joinpath(resumed_dir, "simulation_results")
            @test isequal(DataFrame(Arrow.Table(prefix * "_spacecraft.feather")), resumed.table)
            @test isequal(DataFrame(Arrow.Table(prefix * "_link_schedule.feather")), link_schedule(resumed))
        end
    end

    @testset "Checkpoint interruption uses the same mission horizon" begin
        mktempdir() do directory
            settings(resume)=SimulationSettings(results=false,generate_plots=false,
                checkpoint_enabled=true,checkpoint_interval_s=3.0,checkpoint_directory=directory,
                results_directory=directory,resume_from_checkpoint=resume)
            reference=run_simulation(fixture(:multirate_mri;settings=settings(false));return_solution=true)
            stop=DiscreteCallback((u,t,i)->t>6.2,i->throw(Interrupted());save_positions=(false,false))
            @test_throws Interrupted run_simulation(fixture(:multirate_mri;settings=settings(false));
                extra_callbacks=(stop,))
            resumed_args = fixture(:multirate_mri; settings=settings(true))
            payload = open(deserialize, SM.IOConfig._checkpoint_paths(resumed_args).data)
            @test payload.t == 6.0
            resumed=run_simulation(resumed_args;return_solution=true)
            # The API returns only the final checkpoint segment, 9 to 10 s.
            @test first(resumed.t)==9.0
            @test resumed.t[end]==10.0
            @test resumed.u[end] ≈ reference.u[end] rtol=1e-10 atol=1e-8
            @test resumed.u[end].sc[1].data_sent_bits ≈ 10000.0 atol=1e-6
        end
    end
end
end
end
