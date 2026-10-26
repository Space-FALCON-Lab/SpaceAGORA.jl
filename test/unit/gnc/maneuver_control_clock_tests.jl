module ManeuverControlClockTests
using Test, SpaceAGORA, CSV, DataFrames
using OrdinaryDiffEq: DiscreteCallback
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const L = SpaceAGORA.SimulationLifecycle
const CH = SM.ControlHooks
const PLANET = make_no_gram_planet(:earth)

struct TimedRequest <: AbstractGuidanceModel
    after_s::Float64
    ticks::Vector{Float64}
end
function SpaceAGORA.update_guidance!(m::TimedRequest, u, p, t::Float64, i::Int)
    push!(m.ticks, t)
    if t >= m.after_s
        p.shared_buffers.maneuver_commands[i] = SM.PropulsiveManeuverCommand(valid=true,
            delta_v_mps=1.0, direction_rad=0.0, source_orbit=1)
    end
    return nothing
end
L.initialize_guidance!(m::TimedRequest, u, p, t) = empty!(m.ticks)
L.checkpoint_state(m::TimedRequest) = (ticks=copy(m.ticks),)
L.restore_checkpoint_state!(m::TimedRequest, data) = copy!(m.ticks, data.ticks)

function configuration(; duration=400.0, guidance_period=3.0, request_after=5.0,
        directory="", interval=120.0, resume=false, mode=:dp8, configured_dv=0.0)
    root = Link(root=true, m=100.0)
    sc = SpacecraftModel(root=root, links=[root], inertia_tensor=root.inertia, prop_mass=50.0,
        initial_condition=InitialCondition(ra=PLANET.Rp_e+2000e3, rp=PLANET.Rp_e+400e3,
            i=30.0, Ω=20.0, ω=0.0, ν=170.0), id=1)
    ctl = BaseThrusterModel(thrust=[10.0], direction=[0.0], Δv=[configured_dv],
        start_burn_time=[-1.0], stop_burn_time=[-1.0], Isp=[300.0])
    env = EnvironmentModel(planet=PLANET, EI=120.0, density_model=NoAtmosphereModel(),
        ephemerides_model=SimpleEphemeridesModel(), thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=PLANET),
        topography=false, wind=false)
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false,
            checkpoint_enabled=!isempty(directory), checkpoint_directory=directory,
            checkpoint_interval_s=interval, resume_from_checkpoint=resume),
        mission_configuration=MissionConfiguration(mission_type=MissionTime, keplerian=true, mission_time=duration,
            orientation_sim=false, num_steps_to_save=100),
        environment_model=env, dynamics_model=DynamicsModel([sc], (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((TimedRequest(request_after, Float64[]),), [guidance_period]),
        navigation_model=NavigationModel((), Float64[]), control_model=ControlModel((ctl,), [10.0]),
        solver_config=SolverConfig(solver_mode=mode, multirate_slow_dt_s=1.0),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end
function traced(args; path=tempname()*".csv")
    result = withenv("SPACEAGORA_MANEUVER_TRACE_CSV" => path) do
        run_simulation(args; return_solution=true)
    end
    return result, isfile(path) ? CSV.read(path, DataFrame) : DataFrame(), path
end
scheduler(sol) = only(sol.prob.p.args.control_model.control_effectors)
plans(trace) = trace[trace.event .== "schedule_set", :]

@testset "Maneuver requests follow the control clock" begin
    for guidance_period in (3.0, 10.0)
        sol, trace, _ = traced(configuration(; guidance_period))
        accepted = plans(trace)
        @test nrow(accepted) == 1 && accepted.t_s == [10.0]
        @test only(sol.prob.p.args.guidance_model.guidance_effectors).ticks == collect(guidance_period:guidance_period:399.0)
        @test count(==("burn_start"), trace.event) == 1
        @test count(==("burn_end"), trace.event) == 1
        @test !any(==("schedule_update"), trace.event)
        @test scheduler(sol).consumed_command == [(1.0, 0.0, 1)]
        expected = 150.0 * (1-exp(-1/(300*9.80665)))
        @test 150.0-sol.u[end].sc[1].mass ≈ expected rtol=1e-8
    end
    @testset "checkpoint before acceptance, before ignition, during burn and at a control tick" begin
        for mode in (:dp8, :multirate_mri)
            baseline, events, _ = traced(configuration(; mode))
            accepted = only(eachrow(plans(events)))
            # All windows must remain unchanged once accepted, even when every
            # checkpoint introduces a new solver segment.
            for split in (8.0, 10.0, 120.0, (accepted.start_burn_s+accepted.stop_burn_s)/2)
                mktempdir() do directory
                    _, _, tracepath = traced(configuration(; mode, directory, interval=split, duration=split))
                    resumed, trace, _ = traced(configuration(; mode, directory, interval=split, resume=true); path=tracepath)
                    @test nrow(plans(trace)) == 1 && plans(trace).t_s == [10.0]
                    @test only(plans(trace).start_burn_s) ≈ accepted.start_burn_s atol=1e-6
                    @test only(plans(trace).stop_burn_s) ≈ accepted.stop_burn_s atol=1e-6
                    @test count(==("burn_start"), trace.event) == 1 && count(==("burn_end"), trace.event) == 1
                    @test scheduler(resumed).consumed_command == scheduler(baseline).consumed_command
                    @test only(resumed.prob.p.args.guidance_model.guidance_effectors).ticks ==
                        only(baseline.prob.p.args.guidance_model.guidance_effectors).ticks
                    @test collect(resumed.u[end]) ≈ collect(baseline.u[end]) rtol=1e-8 atol=1e-5
                end
            end
        end
    end
    @testset "interruption preserves the original mission horizon" begin
        mktempdir() do directory
            whole, _, _ = traced(configuration(directory=mktempdir(), interval=120.0))
            path = tempname()*".csv"
            stop = DiscreteCallback((u,t,integrator) -> t > 125.0,
                integrator -> error("injected interruption"); save_positions=(false,false))
            @test_throws Exception withenv("SPACEAGORA_MANEUVER_TRACE_CSV" => path) do
                run_simulation(configuration(; directory); extra_callbacks=(stop,))
            end
            resumed, trace, _ = traced(configuration(; directory, resume=true); path)
            @test nrow(plans(trace)) == 1
            @test collect(resumed.u[end]) == collect(whole.u[end])
            @test only(resumed.prob.p.args.guidance_model.guidance_effectors).ticks ==
                only(whole.prob.p.args.guidance_model.guidance_effectors).ticks
        end
    end
    @testset "accepted requests, new configured requests and explicit replacement plans" begin
        args = configuration(configured_dv=1.0)
        p = SM.ODEParams(n_sats=1, args=args)
        u = SE.build_initial_conditions(args)
        ctl = only(args.control_model.control_effectors)
        L.initialize_control!(ctl, u, p, 0.0)
        update_control!(ctl, u, p, 10.0, 1)
        first_window = (ctl.start_burn_time[1], ctl.stop_burn_time[1])
        ctl.Δv[1] = 2.0
        update_control!(ctl, u, p, 20.0, 1)
        @test (ctl.start_burn_time[1], ctl.stop_burn_time[1]) == first_window
        @test ctl.consumed_configured[1] == (1.0, 0.0)
        update_control!(ctl, u, p, first_window[2]+1, 1)
        @test ctl.consumed_configured[1] == (2.0, 0.0)
        @test p.shared_buffers.maneuver_burn_plans[1].delta_v_mps == 2.0
        for after_stop in (false, true)
            oldstart, oldstop = CH._effective_burn_window(ctl, p, 1)
            # Give the scheduler a stale own-window entry for both controls.
            ctl.start_burn_time[1], ctl.stop_burn_time[1] = oldstart, oldstop
            t = after_stop ? oldstop+1 : (oldstart+oldstop)/2
            replacement = PropulsiveBurnPlan(valid=true, delta_v_mps=3.0, direction_rad=0.0,
                source_orbit=2, start_burn_s=t+50, stop_burn_s=t+60, thrust_n=10.0, isp_s=300.0)
            p.shared_buffers.maneuver_burn_plans[1] = replacement
            update_control!(ctl, u, p, t, 1)
            @test p.shared_buffers.maneuver_burn_plans[1] === replacement
            @test CH._effective_burn_window(ctl, p, 1) == (t+50, t+60)
        end
        # An unchanged guidance request cannot repeat; the same Δv on a new
        # source orbit is a new request.
        end_s = CH._effective_burn_window(ctl, p, 1)[2]
        p.shared_buffers.maneuver_commands[1] = SM.PropulsiveManeuverCommand(valid=true,
            delta_v_mps=3.0, direction_rad=0.0, source_orbit=2)
        update_control!(ctl, u, p, end_s+1, 1)
        @test !p.shared_buffers.maneuver_burn_plans[1].valid
        p.shared_buffers.maneuver_commands[1] = SM.PropulsiveManeuverCommand(valid=true,
            delta_v_mps=3.0, direction_rad=0.0, source_orbit=3)
        update_control!(ctl, u, p, end_s+2, 1)
        @test p.shared_buffers.maneuver_burn_plans[1].valid
        @test ctl.consumed_command[1] == (3.0, 0.0, 3)
        snapshot = L.checkpoint_state(ctl)
        restored = deepcopy(ctl)
        L.initialize_control!(restored, u, p, 0.0)
        L.restore_checkpoint_state!(restored, snapshot)
        @test restored.consumed_command == ctl.consumed_command
        @test restored.consumed_configured == ctl.consumed_configured
        legacy = Base.structdiff(snapshot, NamedTuple{(:version, :consumed_command, :consumed_configured)})
        @test_throws r"BaseThrusterModel.*predates control-clock.*even when no burn was configured" L.restore_checkpoint_state!(restored, legacy)
    end
end
# A real burn-plan binding installs consecutive plans between control ticks.
# Direct stale-model-window replacement is covered in the unit controls above.
struct OffGridPlanner <: AbstractBurnPlanner end
function SpaceAGORA.GNCAlgorithms.propose!(state, ::OffGridPlanner, r::BurnPlanRequest, rng)
    sized = CH._burn_plan_sizing(r.mass_kg, 0.1, 0.0, r.thrust_n, r.isp_s, nothing)
    start = r.time_s == 0.0 ? 5.0 : r.time_s+1
    plan = PropulsiveBurnPlan(valid=true, delta_v_mps=0.1, direction_rad=0.0, source_orbit=2,
        start_burn_s=start, stop_burn_s=start+sized.burn_duration_s, thrust_n=r.thrust_n, isp_s=r.isp_s,
        commanded_impulse_n_s=sized.commanded_impulse_n_s, propellant_required_kg=sized.propellant_required_kg)
    return AlgorithmResult(request_id=r.request_id, status=:candidate, termination=:planned, product=plan)
end
@testset "Explicit plans installed off the control grid keep their ignition" begin
    args = configuration(duration=40.0)
    binding = AlgorithmGuidance(algorithm=OffGridPlanner(), spacecraft_index=1,
        events=[AlgorithmEvent(0.0), AlgorithmEvent(15.0; token=2)])
    args = SM.SimConfig._with_configuration(args; guidance_model=GuidanceModel((binding,), [3.0]))
    sol, trace, _ = traced(args)
    plan = algorithm_run_report(sol).bindings[1].product
    starts = trace[trace.event .== "burn_start", :].t_s
    stops = trace[trace.event .== "burn_end", :].t_s
    @test starts == [5.0, 16.0]
    first_propellant = 150.0 * (1-exp(-0.1/(300*9.80665)))
    @test stops ≈ [5.0 + first_propellant*(300*9.80665)/10.0, plan.stop_burn_s] atol=1e-10
    @test 16.0 in sol.t && plan.stop_burn_s in sol.t
    @test 150.0-sol.u[end].sc[1].mass ≈ first_propellant+plan.propellant_required_kg rtol=1e-8
end

end
