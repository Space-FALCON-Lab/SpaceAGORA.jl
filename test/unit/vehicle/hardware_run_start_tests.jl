# The run state of mounted hinges, thruster gimbals and thrusters with firing
# dynamics at a run's start (G-8): runs of one configuration with
# `isolate_state=false` start equal, a user's edit between two runs is kept and
# becomes the start of later runs, and a run resumed from a checkpoint continues
# the checkpoint's state.
module HardwareRunStartTests

using Test, Logging, StaticArrays, SpaceAGORA
using OrdinaryDiffEq: DiscreteCallback
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const ZERO3 = SVector(0.0, 0.0, 0.0)

# A control effector that commands mounted actuators at constant levels and
# keeps no state of its own.
struct HoldLevels <: AbstractControlEffectorModel
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::HoldLevels) = :mounted
SpaceAGORA.actuator_levels(c::HoldLevels, i::Int) = c.levels
SpaceAGORA.update_control!(::HoldLevels, u, p, t::Float64, i::Int64) = nothing

# A spacecraft far from any planet with a panel link turned by a hinge (at most
# 0.1 rad/s), a 1 N jet along body x with firing dynamics (a ramp slow enough
# that its ramp factor carries over several periods) turned by a gimbal (at
# most 0.05 rad/s), and a plain 1 N jet along body y. A controller holds the
# jets at 0.6 N and 0.3 N and the hinge and gimbal at 0.5 rad, ticking every
# second for 3 s. `hardware` picks the parts mounted (all by default); `hinge`,
# `gimbal` and `κ` are their configured starts.
function configuration(; hardware=(:hinge, :gimbal, :firing), hinge=0.0, gimbal=0.0, κ=0.0, mission_time=3.0, rate=1.0)
    root = SM.Link(root=true, m=100.0, dims=MVector(1.0, 1.0, 1.0))
    panel = SM.Link(m=5.0, r=MVector(0.0, 1.5, 0.0), dims=MVector(0.02, 2.0, 1.0))
    ic = SM.CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), ZERO3; q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SM.SpacecraftModel(links=[root, panel], root=root, prop_mass=1.0, inertia_tensor=root.inertia, initial_condition=ic, id=1)
    jet = add_actuator!(sc, root, Thruster(max_thrust=1.0, direction=MVector(1.0, 0.0, 0.0), Isp=200.0, cutoff_frequency=2.0,
        firing_dynamics=:firing in hardware, κ=κ))
    plain = add_actuator!(sc, root, Thruster(max_thrust=1.0, direction=MVector(0.0, 1.0, 0.0), Isp=200.0))
    levels = [0.6, 0.3]
    if :hinge in hardware
        add_actuator!(sc, panel, PanelHinge(axis=SVector(0.0, 1.0, 0.0), max_rate_rad_s=0.1, angle_rad=hinge))
        push!(levels, 0.5)
    end
    if :gimbal in hardware
        add_actuator!(sc, root, ThrusterGimbal(jet; axis=SVector(0.0, 1.0, 0.0), max_rate_rad_s=0.05, angle_rad=gimbal))
        push!(levels, 0.5)
    end
    base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=mission_time,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
        keplerian=false, EI_km=1.0, verbose=false, results=false)
    return SM.SimConfig._with_configuration(base;
        control_model=SM.ControlModel(control_effectors=(HoldLevels(levels),), control_rates=[rate]))
end

spacecraft(args) = only(args.dynamics_model.spacecraft)
hinge_of(args) = only(spacecraft(args).links[2].hinges)
gimbal_of(args) = only(spacecraft(args).root.gimbals)
jets(args) = spacecraft(args).root.thrusters
panel_q(args) = collect(spacecraft(args).links[2].q)

# Everything a run leaves in the parts and their geometry.
hardware_state(args) = (
    hinges=[(h.angle_rad, h.last_update_s) for l in spacecraft(args).links for h in l.hinges],
    gimbals=[(g.angle_rad, g.last_update_s, g.previous_angle_rad) for g in spacecraft(args).root.gimbals],
    firing=[(t.κ, t.thrust) for t in jets(args)],
    panel=panel_q(args),
    turn=collect(jets(args)[1].gimbal_rotation),
)

# The parts' state after every accepted step of a run.
function traced(args; kw...)
    trace = Any[]
    record = DiscreteCallback((u, t, integrator) -> true, integrator -> push!(trace, (integrator.t, hardware_state(args)));
        save_positions=(false, false))
    sol = reused(args; extra_callbacks=(record,), kw...)
    return sol, trace
end

quiet(f) = with_logger(f, ConsoleLogger(stderr, Logging.Error))
reused(args; kw...) = quiet(() -> run_simulation(args; isolate_state=false, return_solution=true, kw...))
bits(sol) = [reinterpret(UInt64, x) for u in sol.u for x in collect(Float64, u)]
same_run(a, b) = a.t == b.t && bits(a) == bits(b)

# `args` checkpointed every 2 s into `directory`, or resumed from there.
function checkpointed(args, directory; resume=false)
    s = args.simulation_settings
    names = fieldnames(typeof(s))
    fields = merge(NamedTuple{names}(map(name -> getfield(s, name), names)), (checkpoint_enabled=true,
        checkpoint_interval_s=2.0, checkpoint_directory=String(directory), resume_from_checkpoint=resume))
    return SM.SimConfig._with_configuration(args; simulation_settings=SM.SimulationSettings(; fields...))
end

# Thrown at the first step past the checkpoint at 2 s, so the run stops there as
# an interrupted run would; unarmed, the callback leaves a run unchanged.
struct Interrupted <: Exception end
interrupter(armed) = DiscreteCallback((u, t, integrator) -> armed[] && t > 2.0, integrator -> throw(Interrupted());
    save_positions=(false, false))
unwrapped(err) = err isa TaskFailedException ? unwrapped(err.task.exception) : err
outcome(f) = try
    f()
    nothing
catch err
    unwrapped(err)
end

@testset "Hardware run state at a run's start (G-8)" begin
    @testset "three runs of one configuration start equal: $(join(parts, ", "))" for parts in
            ((:hinge,), (:gimbal,), (:firing,), (:hinge, :gimbal, :firing))
        args = configuration(hardware=parts)
        start = hardware_state(args)
        first_run, first_trace = traced(args)
        ended = hardware_state(args)
        @test !isequal(ended, start)                                           # the run moved them
        runs = [traced(args) for _ in 1:2]
        @test all(r -> same_run(r[1], first_run) && isequal(r[2], first_trace), runs)
        @test isequal(hardware_state(args), ended)
        # Each run ends as a fresh configuration's run does.
        fresh = configuration(hardware=parts)
        @test same_run(reused(fresh), first_run) && isequal(hardware_state(fresh), ended)
    end

    @testset "the hinge of the example: 0.1 rad after each 3 s run" begin
        args = configuration(hardware=(:hinge,))
        angles = Float64[]
        for _ in 1:3
            reused(args)
            push!(angles, hinge_of(args).angle_rad)
        end
        @test angles ≈ [0.1, 0.1, 0.1] atol=1e-12
        @test panel_q(args) ≈ [0.0, sin(0.05), 0.0, cos(0.05)]
    end

    @testset "an edit between two runs is kept and becomes the start: $(part)" for (part, edit!, built) in (
            (:hinge, args -> (hinge_of(args).angle_rad = 0.25), () -> configuration(hinge=0.25)),
            (:gimbal, args -> (gimbal_of(args).angle_rad = 0.02), () -> configuration(gimbal=0.02)),
            (:firing, args -> (jets(args)[1].κ = 0.5), () -> configuration(κ=0.5)))
        # (Without an atmosphere or radiation pressure here, a hinge moves no
        # state of the run: its angle shows in the parts' state at the end.)
        args = configuration()
        first_run = reused(args)
        first_end = hardware_state(args)
        edit!(args)
        edited_run = reused(args)
        @test !isequal(hardware_state(args), first_end)
        # As a fresh configuration with the edited value, and so is the run after it.
        fresh = built()
        fresh_run = reused(fresh)
        @test same_run(edited_run, fresh_run)
        edited_end = hardware_state(args)
        @test isequal(edited_end, hardware_state(fresh))
        @test same_run(reused(args), edited_run) && isequal(hardware_state(args), edited_end)
    end

    @testset "an edited hinge or gimbal angle turns its link or thrust from the run's start" begin
        args = configuration(hardware=(:hinge, :gimbal))
        reused(args)
        hinge_of(args).angle_rad, gimbal_of(args).angle_rad = 0.25, 0.02
        observed = Ref{Any}(nothing)
        watch = DiscreteCallback((u, t, integrator) -> true, integrator -> nothing;
            initialize=(cb, u, t, integrator) -> (observed[] === nothing && (observed[] = (panel_q(args), collect(jets(args)[1].gimbal_rotation)))),
            save_positions=(false, false))
        reused(args; extra_callbacks=(watch,))
        fresh = configuration(hardware=(:hinge, :gimbal), hinge=0.25, gimbal=0.02)
        @test observed[][1] == panel_q(fresh)
        @test observed[][2] == collect(jets(fresh)[1].gimbal_rotation)
    end

    @testset "a resumed run continues the checkpoint's state" begin
        armed = Ref(false)
        callback = interrupter(armed)
        straight = checkpointed(configuration(mission_time=4.0, rate=0.5), mktempdir())
        straight_sol = reused(straight; extra_callbacks=(callback,))
        directory = mktempdir()
        interrupted = checkpointed(configuration(mission_time=4.0, rate=0.5), directory)
        armed[] = true
        stopped = try
            outcome(() -> reused(interrupted; extra_callbacks=(callback,)))
        finally
            armed[] = false
        end
        @test stopped isa Interrupted
        saved = only(SE._load_checkpoint(interrupted).runtime_state.callbacks.spacecraft)
        @test 0.0 < only(saved.hinges).angle_rad && 0.0 < only(saved.thrusters).κ < 1.0   # a history to continue
        # The interrupted configuration itself resumes: a resume is not a run
        # start, so the checkpoint's state is not reset.
        resumed = checkpointed(interrupted, directory; resume=true)
        resumed_sol = reused(resumed; extra_callbacks=(callback,))
        @test same_run(resumed_sol, straight_sol)
        @test isequal(hardware_state(resumed), hardware_state(straight))
        # A run from the start after the resume starts from the configured state.
        plain = SM.SimConfig._with_configuration(resumed; simulation_settings=configuration().simulation_settings)
        after = reused(plain)
        fresh = configuration(mission_time=4.0, rate=0.5)
        @test same_run(after, reused(fresh))
        @test isequal(hardware_state(plain), hardware_state(fresh))
    end

    @testset "negative controls" begin
        @testset "a thruster without firing dynamics keeps no run start" begin
            args = configuration()
            plain = jets(args)[2]
            plain.κ, plain.thrust = 0.7, 0.4                                   # fields no tick writes
            reused(args)
            @test plain.run_written === nothing && (plain.κ, plain.thrust) == (0.7, 0.4)
        end
        @testset "parts no controller drives keep their state" begin
            args = SM.SimConfig._with_configuration(configuration();
                control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]))
            hinge_of(args).angle_rad = 0.3
            reused(args)
            reused(args)
            @test hinge_of(args).angle_rad == 0.3 && hinge_of(args).run_written === nothing
        end
        @testset "an edit equal to the last write is not seen" begin
            args = configuration(hardware=(:hinge,))
            reused(args)
            hinge_of(args).angle_rad = hinge_of(args).angle_rad                # set to the value the run wrote
            reused(args)
            @test hinge_of(args).run_start.angle_rad == 0.0                    # the start stays the configured one
        end
        @testset "an isolated run leaves the configuration untouched" begin
            args = configuration()
            before = hardware_state(args)
            quiet(() -> run_simulation(args))
            @test isequal(hardware_state(args), before) && hinge_of(args).run_written === nothing
        end
    end
end

end # module HardwareRunStartTests
