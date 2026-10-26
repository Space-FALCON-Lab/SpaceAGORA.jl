module InterlinkMatchingResumeTests
using Test, SpaceAGORA, StaticArrays, LinearAlgebra
using OrdinaryDiffEq: DiscreteCallback
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PLANET = make_no_gram_planet(:earth)
const UPLINK = ((-1, 1), (1, 1))
struct MatchingInterrupted <: Exception end

function configuration(directory; resume=false)
    radius = PLANET.Rp_e + 700e3
    env = make_no_gram_environment(planet=PLANET)
    epoch = SM.EphemeridesModels.ephemerides_epoch_seconds(env.ephemerides_model, InitialTime())
    frame = SMatrix{3,3,Float64}(SM.EphemeridesModels.ephemerides_frame_lpi(env.ephemerides_model, PLANET, epoch))
    up = transpose(frame) * SVector(1.0,0.0,0.0)
    pole = transpose(frame) * SVector(0.0,0.0,1.0)
    radial = cosd(10.0)*up - sind(10.0)*cross(pole,up)
    root = Link(root=true, m=100.0)
    craft = SpacecraftModel(root=root, links=[root], inertia_tensor=root.inertia, id=1,
        initial_condition=CartesianInitialCondition(radius*radial, sqrt(PLANET.μ/radius)*cross(pole,radial)))
    radio() = RadioTerminal(boresight=[0.0,0.0,1.0], max_power_w=5.0, gain_dbi=0.0, efficiency=0.5)
    for _ in 1:2
        add_actuator!(craft, root, radio())
    end
    add_data!(craft, root, DataStore(capacity_bits=1e7, initial_bits=1e7))
    model = InterlinkModel()
    add_ground_station!(model, GroundStation(name="relay", lat_deg=0.0, lon_deg=0.0,
        min_elevation_deg=10.0, antennas=[radio(), radio()]))
    link(rate) = DataLink(transmit_power_w=0.0, budget=FixedRateBudget(rate_bps=rate), range_m=3e6)
    register_link!(model, [craft], (-1,1), (1,1), link(2e4))
    register_link!(model, [craft], (1,2), (-1,2), link(1e4))
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, generate_plots=false, verbose=false,
            results_directory=directory, checkpoint_directory=directory,
            checkpoint_enabled=true, checkpoint_interval_s=60.0, resume_from_checkpoint=resume),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, mission_time=160.0,
            orientation_sim=false, data_rate=1.0),
        environment_model=env, dynamics_model=DynamicsModel([craft], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((),Float64[]), navigation_model=NavigationModel((),Float64[]),
        control_model=ControlModel((),Float64[]), initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-8, dt_max_orbit=5.0),
        solver_config=SolverConfig(solver_mode=:tsit5), interlink_model=model,
        scheduling_policy_model=DataBacklogPolicy(horizon_s=60.0))
end
function observer(records)
    DiscreteCallback((u,t,i)->false, i->nothing; save_positions=(false,false),
        initialize=(cb,u,t,i)->begin
            if t == 60.0
                args = i.p.args
                push!(records, (active=copy(last(args.interlink_model.history).active),
                    score=args.scheduling_policy_model.matching.score[1],
                    stored=Float64(u.sc[1].data_stored_bits)))
            end
            nothing
        end)
end
@testset "Resume preserves the previous matching for receiver outflow" begin
    baseline = Any[]
    whole = run_simulation(configuration(mktempdir()); return_solution=true,
        extra_callbacks=(observer(baseline),))
    @test length(baseline) == 1
    @test only(baseline).stored ≈ 1e7 atol=1e-6
    @test UPLINK in only(baseline).active
    @test only(baseline).score ≈ 1e4 atol=1e-6
    for keep_history in (true, false)
        mktempdir() do directory
            crash = DiscreteCallback((u,t,i)->t>65.0, i->throw(MatchingInterrupted()); save_positions=(false,false))
            @test_throws MatchingInterrupted run_simulation(configuration(directory); extra_callbacks=(crash,))
            checkpoint = SE._load_checkpoint(configuration(directory))
            @test checkpoint.t == 60.0
            snapshot = checkpoint.runtime_state.interlinks
            @test Set(snapshot.active_links) == Set(only(baseline).active)
            # Reject unknown/duplicate candidates before mutating the destination.
            for invalid in ([UPLINK, UPLINK], [((99,1),(100,1))])
                fresh = configuration(directory)
                flags = trues(1)
                damaged = merge(checkpoint, (runtime_state=merge(checkpoint.runtime_state,
                    (interlinks=merge(snapshot, (active_links=invalid,)),)),))
                @test_throws ArgumentError SE._restore_interlink_checkpoint_state!(fresh, damaged, flags)
                @test all(flags) && !any(fresh.interlink_model.graph.active)
            end
            # A distinct registered edge can still conflict at a terminal.
            fresh = configuration(directory)
            register_link!(fresh.interlink_model, fresh.dynamics_model.spacecraft,
                (-1,1), (1,2), DataLink(transmit_power_w=0.0,
                    budget=FixedRateBudget(rate_bps=1e4), range_m=3e6))
            bad_runtime = SE._interlink_checkpoint_runtime_state(fresh, [false])
            conflicting = merge(checkpoint, (runtime_state=(interlinks=merge(bad_runtime.interlinks,
                (active_links=[UPLINK, ((-1,1),(1,2))],)),),))
            flags = trues(1)
            @test_throws r"must not share terminals" SE._restore_interlink_checkpoint_state!(fresh, conflicting, flags)
            @test all(flags) && !any(fresh.interlink_model.graph.active)
            legacy = merge(checkpoint, (runtime_state=merge(checkpoint.runtime_state,
                (interlinks=Base.structdiff(snapshot, NamedTuple{(:active_links,)}),)),))
            fresh = configuration(directory)
            @test_logs (:warn, r"no previous active interlink matching") SE._restore_interlink_checkpoint_state!(fresh, legacy, trues(1))
            @test !any(fresh.interlink_model.graph.active)
            if !keep_history
                history_path = joinpath(directory, "simulation_checkpoint.interlinks.bin")
                @test isfile(history_path)
                rm(history_path)
                @test !isfile(history_path)
            end
            records = Any[]
            resumed = run_simulation(configuration(directory; resume=true); return_solution=true,
                extra_callbacks=(observer(records),))
            @test records == baseline
            history = resumed.prob.p.args.interlink_model.history
            @test all(sample -> UPLINK in sample.active, filter(sample -> sample.time >= 60.0, history))
            @test collect(resumed.u[end]) == collect(whole.u[end])
            @test iszero(resumed.u[end].sc[1].data_lost_bits)
        end
    end
end
end
