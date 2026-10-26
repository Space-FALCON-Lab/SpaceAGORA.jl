# A mounted thruster with firing dynamics in checkpoint files: its ramp factor
# `κ` and held average thrust `thrust` are saved with the engine's callback
# snapshot and restored on a resume, so a resumed run ends as an uninterrupted
# run with the same checkpoints, bit for bit; a snapshot for another thruster
# configuration, or one written before the firing state was saved, is refused.
module ThrusterFiringCheckpointTests

using Test, Logging, Serialization, StaticArrays, SpaceAGORA
using OrdinaryDiffEq: DiscreteCallback
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const C = SpaceAGORA.Conformance
const ZERO3 = SVector(0.0, 0.0, 0.0)

# A control effector that commands mounted actuators at constant levels and
# keeps no state of its own.
struct HoldLevels <: AbstractControlEffectorModel
    levels::Vector{Float64}
end
SpaceAGORA.actuation_mode(::HoldLevels) = :mounted
SpaceAGORA.actuator_levels(c::HoldLevels, i::Int) = c.levels
SpaceAGORA.update_control!(::HoldLevels, u, p, t::Float64, i::Int64) = nothing

# A spacecraft far from any planet with two 1 N thrusters along body x,
# commanded at a constant 0.6 N and 0.3 N every 0.5 s for 4 s. The first (or,
# with `swap`, the second) has firing dynamics with a ramp of `cutoff` rad/s,
# slow enough that its ramp factor carries over several control periods.
# `fault` adds a failure of the second thruster after the run's end. With
# `on=:root` or `on=:panel` the spacecraft also has a rigid panel link 1.5 m off
# its root, and both thrusters are mounted on the link named: their positions in
# the mounting order are 1 and 2 either way.
function configuration(; firing=true, swap=false, cutoff=2.0, fault=false, on=nothing)
    root = SM.Link(root=true, m=100.0, dims=MVector(1.0, 1.0, 1.0))
    ic = SM.CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), ZERO3; q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=1.0, inertia_tensor=root.inertia, initial_condition=ic, id=1)
    panel = SM.Link(m=2.0, r=MVector(0.0, 1.5, 0.0), dims=MVector(0.02, 2.0, 1.0))
    on === nothing || push!(sc.links, panel)
    link = on === :panel ? panel : root
    jets = [add_actuator!(sc, link, Thruster(max_thrust=1.0, direction=MVector(1.0, 0.0, 0.0), Isp=200.0,
                cutoff_frequency=cutoff, firing_dynamics=firing && k == (swap ? 2 : 1))) for k in 1:2]
    fault && add_fault!(sc, ActuatorFault(jets[2]; kind=:failed, start_s=100.0, stop_s=200.0))
    base = make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=4.0,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
        keplerian=false, EI_km=1.0, verbose=false, results=false)
    return SM.SimConfig._with_configuration(base;
        control_model=SM.ControlModel(control_effectors=(HoldLevels([0.6, 0.3]),), control_rates=[0.5]))
end

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

# A checkpointed run warns that it returns its last segment's solution.
quiet(f) = with_logger(f, ConsoleLogger(stderr, Logging.Error))
run!(args, callback) = quiet(() -> run_simulation(args; isolate_state=false, return_solution=true, extra_callbacks=(callback,)))
outcome(f) = try
    f()
    nothing
catch err
    unwrapped(err)
end

bits(sol) = [reinterpret(UInt64, x) for u in sol.u for x in collect(Float64, u)]
jets(args) = only(args.dynamics_model.spacecraft).root.thrusters
firing_state(args) = [(t.κ, t.thrust) for t in jets(args)]

# The uninterrupted run, and a run interrupted after its checkpoint at 2 s (its
# directory then holds the checkpoint), each made by `make()`.
function interrupted_pair(make)
    armed = Ref(false)
    callback = interrupter(armed)
    straight = checkpointed(make(), mktempdir())
    straight_sol = run!(straight, callback)
    directory = mktempdir()
    interrupted = checkpointed(make(), directory)
    armed[] = true
    stopped = try
        outcome(() -> run!(interrupted, callback))
    finally
        armed[] = false
    end
    return (; straight, straight_sol, interrupted, stopped, directory, callback)
end

resume!(pair, args) = run!(checkpointed(args, pair.directory; resume=true), pair.callback)

# The checkpoint file of `args` with its callbacks' snapshot replaced by
# `f(snapshot)` (removed when `f` returns nothing).
function rewrite_snapshot!(args, f)
    path = SE._checkpoint_paths(args).data
    payload = deserialize(path)
    snapshot = f(get(payload.runtime_state, :callbacks, nothing))
    runtime = snapshot === nothing ? Base.structdiff(payload.runtime_state, NamedTuple{(:callbacks,)}) :
        merge(payload.runtime_state, (callbacks=snapshot,))
    open(io -> serialize(io, merge(payload, (runtime_state=runtime,))), path, "w")
    return nothing
end
# A snapshot as SpaceAGORA 0.4 wrote it: version 1, no thrusters, and a faulted
# actuator placed by its position in the mounting order alone.
as_version_1(snapshot) = (version=1,
    configuration=[merge(Base.structdiff(c, NamedTuple{(:thrusters,)}), (faults=[(f[1], first(f[2]), f[3:end]...) for f in c.faults],))
                   for c in snapshot.configuration],
    spacecraft=[Base.structdiff(s, NamedTuple{(:thrusters,)}) for s in snapshot.spacecraft])
refusal(f) = (err = outcome(f); err isa ArgumentError ? sprint(showerror, err) : err)

@testset "Thruster firing state in checkpoint files" begin
    @testset "a resumed run continues the firing history" begin
        pair = interrupted_pair(configuration)
        @test pair.stopped isa Interrupted
        checkpoint = SE._load_checkpoint(pair.interrupted)
        @test checkpoint.t == 2.0
        snapshot = checkpoint.runtime_state.callbacks
        @test snapshot.version == 2
        # Mounting position 1, on link 1 as its first actuator, and its firing parameters.
        @test only(snapshot.configuration).thrusters == [((1, 1, 1), 1.0, 2.0, 0.0, 0.75, 0.25)]
        saved = only(only(snapshot.spacecraft).thrusters)
        @test 0.0 < saved.κ < 1.0 && saved.thrust > 0.0                                    # a firing history to continue
        @test (saved.κ, saved.thrust) == first(firing_state(pair.interrupted))
        resumed = checkpointed(configuration(), pair.directory; resume=true)
        @test firing_state(resumed) == [(0.0, 0.0), (0.0, 0.0)]                            # a fresh configuration
        resumed_sol = run!(resumed, pair.callback)
        @test resumed_sol.t == pair.straight_sol.t
        @test bits(resumed_sol) == bits(pair.straight_sol)
        @test firing_state(resumed) == firing_state(pair.straight)
    end

    @testset "without firing dynamics, as before" begin
        pair = interrupted_pair(() -> configuration(firing=false))
        @test pair.stopped isa Interrupted
        @test !hasproperty(SE._load_checkpoint(pair.interrupted).runtime_state, :callbacks)   # nothing to save
        resumed_sol = resume!(pair, configuration(firing=false))
        @test bits(resumed_sol) == bits(pair.straight_sol)
    end

    @testset "the conformance round trip passes" begin
        for firing in (true, false)
            args = configuration(; firing)
            check = quiet(() -> C._gnc_checkpoint_round_trip(only(args.control_model.control_effectors), :control_model, args))
            @test check.status === :pass
        end
    end

    @testset "another thruster configuration is refused" begin
        pair = interrupted_pair(configuration)
        for (name, args) in ("another ramp" => configuration(cutoff=3.0),
                             "firing dynamics on the other thruster" => configuration(swap=true),
                             "no firing dynamics" => configuration(firing=false))
            @testset "$name" begin
                msg = refusal(() -> resume!(pair, args))
                @test msg isa String && occursin(r"thrusters with firing dynamics.*differ from the checkpoint's", msg)
            end
        end
        @test refusal(() -> resume!(pair, configuration())) === nothing                    # the same configuration resumes
    end

    @testset "an actuator moved to another link is refused" begin
        # Both thrusters move from the root to the panel: their positions in the
        # mounting order stay 1 and 2, their links change.
        on_root, on_panel = configuration(on=:root), configuration(on=:panel)
        root_jets = only(on_root.dynamics_model.spacecraft).root.thrusters
        panel_jets = only(on_panel.dynamics_model.spacecraft).links[2].thrusters
        @test SE._actuator_place(only(on_root.dynamics_model.spacecraft), root_jets[1]) == (1, 1, 1)
        @test SE._actuator_place(only(on_panel.dynamics_model.spacecraft), panel_jets[1]) == (1, 2, 1)
        pair = interrupted_pair(() -> configuration(on=:root))
        msg = refusal(() -> resume!(pair, configuration(on=:panel)))
        @test msg isa String && occursin(r"thrusters with firing dynamics.*differ from the checkpoint's", msg)
        @test bits(resume!(pair, configuration(on=:root))) == bits(pair.straight_sol)
        # Restored onto a thruster on the panel, a run continues as well.
        panel_pair = interrupted_pair(() -> configuration(on=:panel))
        @test bits(resume!(panel_pair, configuration(on=:panel))) == bits(panel_pair.straight_sol)
        # A faulted actuator, without firing dynamics, likewise.
        fault_pair = interrupted_pair(() -> configuration(firing=false, fault=true, on=:root))
        @test hasproperty(SE._load_checkpoint(fault_pair.interrupted).runtime_state, :callbacks)
        msg = refusal(() -> resume!(fault_pair, configuration(firing=false, fault=true, on=:panel)))
        @test msg isa String && occursin("differ from the checkpoint's", msg)
        @test refusal(() -> resume!(fault_pair, configuration(firing=false, fault=true, on=:root))) === nothing
    end

    @testset "a checkpoint without the firing state is refused" begin
        pair = interrupted_pair(configuration)
        rewrite_snapshot!(pair.interrupted, as_version_1)
        msg = refusal(() -> resume!(pair, configuration()))
        @test msg isa String && occursin("thrusters with firing dynamics, whose firing state this checkpoint does not hold", msg) &&
            occursin("SpaceAGORA 0.4 or earlier", msg)
        rewrite_snapshot!(pair.interrupted, snapshot -> nothing)
        msg = refusal(() -> resume!(pair, configuration()))
        @test msg isa String && occursin(r"thrusters with firing dynamics.*whose run state this checkpoint does not hold", msg)
        rewrite_snapshot!(pair.interrupted, snapshot -> (version=3,))
        msg = refusal(() -> resume!(pair, configuration()))
        @test msg isa String && occursin("version 1 or 2 is required", msg)
    end

    @testset "an earlier snapshot still resumes a run without firing dynamics" begin
        args = configuration(firing=false, fault=true)
        snapshot = SE._callback_checkpoint_state(args)
        @test only(snapshot.configuration).thrusters == [] && only(snapshot.spacecraft).thrusters == []
        fault = only(only(args.dynamics_model.spacecraft).root.faults.faults)
        saved = as_version_1(merge(snapshot, (spacecraft=[merge(only(snapshot.spacecraft),
            (faults=[(active=true, held_level=0.25)],))],)))
        p = (shared_buffers=(callback_resume=Ref(false),),)
        @test SE._restore_callback_checkpoint_state!(args, (t=2.0, runtime_state=(callbacks=saved,)), p) === nothing
        @test fault.active && fault.held_level == 0.25 && p.shared_buffers.callback_resume[]
        # With firing dynamics, the same earlier snapshot restores nothing.
        firing = configuration(fault=true)
        firing_fault = only(only(firing.dynamics_model.spacecraft).root.faults.faults)
        q = (shared_buffers=(callback_resume=Ref(false),),)
        @test_throws ArgumentError SE._restore_callback_checkpoint_state!(firing, (t=2.0, runtime_state=(callbacks=saved,)), q)
        @test !firing_fault.active && !q.shared_buffers.callback_resume[]
    end
end

end # module ThrusterFiringCheckpointTests
