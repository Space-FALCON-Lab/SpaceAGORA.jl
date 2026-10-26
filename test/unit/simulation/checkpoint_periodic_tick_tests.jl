module CheckpointPeriodicTickTests
using Test, Serialization, StaticArrays, SpaceAGORA
import SpaceAGORA.TelemetryVerification: make_example_config
const SM = SpaceAGORA.SimulationModel
const SC = SM.SimulationCallbacks
const SE = SpaceAGORA.SimulationEngine
const L = SpaceAGORA.SimulationLifecycle

# A mounted controller that commands its wheel 0.001 * t N m at each tick and
# records the tick, so a lost or shifted tick changes the commanded level and
# the wheel momentum. It keeps the held level in checkpoint files, which a
# resumed run needs until its first tick.
mutable struct RampWheel <: AbstractControlEffectorModel
    levels::Vector{Float64}
    ticks::Vector{Float64}
end
RampWheel() = RampWheel([0.0], Float64[])
SpaceAGORA.actuation_mode(::RampWheel) = :mounted
SpaceAGORA.actuator_levels(c::RampWheel, i::Int) = i == 1 ? c.levels : nothing
function SpaceAGORA.update_control!(c::RampWheel, u, p, t::Float64, i::Int64)
    c.levels[1] = 0.001 * t
    push!(c.ticks, t)
    return nothing
end
L.checkpoint_state(c::RampWheel) = (levels=copy(c.levels),)
L.restore_checkpoint_state!(c::RampWheel, data) = (c.levels .= data.levels; nothing)

# Guidance, navigation and control effectors that only record their ticks (and,
# for control, the clock the engine reports), and a navigation effector that
# fails at its first tick after a given time, to interrupt a run.
struct TickRecorder
    ticks::Vector{Float64}
end
TickRecorder() = TickRecorder(Float64[])
SpaceAGORA.update_guidance!(r::TickRecorder, u, p, t::Float64, i::Int) = (push!(r.ticks, t); nothing)
SpaceAGORA.update_navigation!(r::TickRecorder, u, p, t::Float64, i::Int) = (push!(r.ticks, t); nothing)
struct ControlRecorder <: AbstractControlEffectorModel
    ticks::Vector{Float64}
    clock::Vector{Tuple{Float64, Float64}}
end
ControlRecorder() = ControlRecorder(Float64[], Tuple{Float64, Float64}[])
SpaceAGORA.update_control!(r::ControlRecorder, u, p, t::Float64, i::Int64) = (push!(r.ticks, t); nothing)
SpaceAGORA.control_wrench(::ControlRecorder, x, p, t) = (SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0))
SpaceAGORA.gnc_clock_started!(r::ControlRecorder, period_s, origin_s) = (push!(r.clock, (period_s, origin_s)); nothing)
struct FailAfter
    t_s::Float64
end
SpaceAGORA.update_navigation!(f::FailAfter, u, p, t::Float64, i::Int) = t > f.t_s ? error("interrupted at $t s") : nothing

# A 100 kg spacecraft far from Earth with one wheel and no forces.
function wheel_config(mission_time)
    root = SM.Link(root=true, m=100.0, dims=MVector(1.0, 1.0, 1.0))
    ic = SM.CartesianInitialCondition(SVector(4.0e9, 0.0, 0.0), SVector(0.0, 0.0, 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SM.SpacecraftModel(links=[root], root=root, prop_mass=1.0, inertia_tensor=root.inertia, initial_condition=ic, id=1)
    add_actuator!(sc, root, ReactionWheel(spin_axis=SVector(1.0, 0.0, 0.0), max_torque_nm=1.0, max_momentum_nms=10.0))
    return make_example_config(planet=SM.Earth(), spacecraft=sc, mission_time=mission_time,
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        dynamic_effectors=(), ephemerides_model=SM.SimpleEphemeridesModel(), orientation_sim=true,
        keplerian=false, EI_km=1.0, verbose=false, results=false)
end

# The ramp controller at `period`, a second control effector, guidance and
# navigation at their own periods; `fail_after` adds a navigation effector that
# interrupts the run.
function ramp_setup(; interval=nothing, period, data_rate=period, guidance_period=0.25, navigation_period=0.3,
                    second_control_period=0.3, mission_time=5.0, dir=mktempdir(), resume=false, fail_after=nothing,
                    results=false, layout=:wide)
    control, second, guidance, navigation = RampWheel(), ControlRecorder(), TickRecorder(), TickRecorder()
    navigation_effectors = fail_after === nothing ? (navigation,) : (navigation, FailAfter(fail_after))
    args = SM.SimConfig._with_configuration(wheel_config(mission_time);
        simulation_settings=SM.SimulationSettings(results=results, verbose=false, generate_plots=false,
            normalize=false, save_csv=results, results_directory=dir, results_layout=layout,
            checkpoint_enabled=interval !== nothing, checkpoint_interval_s=something(interval, 300.0),
            checkpoint_directory=joinpath(dir, "checkpoints"), resume_from_checkpoint=resume),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=false,
            mission_time=mission_time, orientation_sim=true, data_rate=data_rate),
        control_model=SM.ControlModel(control_effectors=(control, second), control_rates=[period, second_control_period]),
        guidance_model=SM.GuidanceModel(guidance_effectors=(guidance,), guidance_rates=[guidance_period]),
        navigation_model=SM.NavigationModel(navigation_effectors=navigation_effectors,
            navigation_rates=fill(navigation_period, length(navigation_effectors))))
    return (args=args, control=control, second=second, guidance=guidance, navigation=navigation)
end

# The final state of a checkpointed run is its last checkpoint. `factories` adds
# callbacks to the run (`extra_callback_factories`).
function ramp_run(; interval=nothing, factories=(), kw...)
    s = ramp_setup(; interval=interval, kw...)
    out = run_simulation(s.args; return_results=true, return_solution=interval === nothing, isolate_state=false,
        extra_callback_factories=factories)
    final = interval === nothing ? out.solution.u[end].sc[1] : SE._load_checkpoint(s.args).u.sc[1]
    return (table=out.table, final=final, control=s.control.ticks, second=s.second.ticks, guidance=s.guidance.ticks,
            navigation=s.navigation.ticks, clock=s.second.clock, args=s.args)
end

const TICKS = (:control, :second, :guidance, :navigation)

# The segment boundaries of a checkpointed run from `start`, its end included.
function checkpoint_boundaries(interval, mission_time; start=0.0)
    boundaries = Float64[]
    t = start
    while t < mission_time
        t = SE._checkpoint_segment_end(t, interval, mission_time)
        push!(boundaries, t)
    end
    return boundaries
end

# The ticks of a run in one piece as a checkpointed run takes them: a tick that
# differs from a segment boundary only by rounding (`_same_run_time`) at the boundary.
function at_boundaries(ticks, boundaries)
    return map(ticks) do t
        i = findfirst(b -> SC._same_run_time(t, b), boundaries)
        i === nothing ? t : boundaries[i]
    end
end

function same_ticks_and_state(split, unsplit)
    settings = split.args.simulation_settings
    boundaries = settings.checkpoint_enabled ?
        checkpoint_boundaries(settings.checkpoint_interval_s, split.args.mission_configuration.mission_time) : Float64[]
    @test all(field -> !isempty(getfield(unsplit, field)), TICKS)
    for field in TICKS
        @test getfield(split, field) == at_boundaries(getfield(unsplit, field), boundaries)
    end
    @test split.final.h_wheels[1] ≈ unsplit.final.h_wheels[1] rtol=1e-12
    @test split.final.ω ≈ unsplit.final.ω rtol=1e-12
end

# The run state of a checkpoint (`nothing` when it has none).
runtime_state(args) = get(SE._load_checkpoint(args), :runtime_state, nothing)

# The periods of the clocks a checkpoint lists as not having taken the tick due at
# its time (version 2 entry, at that time).
function pending_periods(args)
    entry = runtime_state(args).periodic_ticks
    @test entry.version == 2 && entry.time_s == SE._load_checkpoint(args).t
    return sort([clock.period_s for clock in entry.pending])
end

@testset "checkpointed runs tick as runs in one piece" begin
    @testset "segment boundaries on the tick grid keep their ticks" begin
        # Exact binary times: every boundary and saved row lands on a tick.
        unsplit = ramp_run(period=0.125)
        split = ramp_run(period=0.125, interval=1.0)
        @test unsplit.control == 0.125:0.125:4.875
        same_ticks_and_state(split, unsplit)
        @test count(t -> t in (1.0, 2.0, 3.0, 4.0), split.control) == 4
        @test count(t -> t in (1.0, 2.0, 3.0, 4.0), split.guidance) == 4
        @test split.table.time == unsplit.table.time
        @test split.table.sc1_actuator_level_1 == unsplit.table.sc1_actuator_level_1
        # The level saved at a boundary is the one commanded there.
        @test split.table.sc1_actuator_level_1[findfirst(==(1.0), split.table.time)] == 0.001
        # Each segment reports the run's clock origin, as the run's start does.
        @test unsplit.clock == [(0.3, 0.0), (0.3, 0.0)]
        @test split.clock == fill((0.3, 0.0), 1 + 5)
    end

    @testset "the reported 0.1 s control period with 1 s segments" begin
        unsplit = ramp_run(period=0.1)
        split = ramp_run(period=0.1, interval=1.0)
        # 10, 20, 30 and 40 times 0.1 are exactly 1.0 to 4.0, so each boundary has
        # its tick (boundaries off by an ulp are tested below).
        @test count(t -> t in (1.0, 2.0, 3.0, 4.0), split.control) == 4
        same_ticks_and_state(split, unsplit)
        # Within a segment tick times stay exact: 3 * 0.1 is 0.30000000000000004,
        # an ulp after the second controller's tick at 0.3, and keeps its own time.
        @test 3 * 0.1 in split.control && 0.3 in split.second
    end

    @testset "segment boundaries off the tick grid do not shift it" begin
        unsplit = ramp_run(period=0.1)
        split = ramp_run(period=0.1, interval=0.75)
        same_ticks_and_state(split, unsplit)
        # Several rates, some ticking on the 1.5 s boundaries and some not.
        rates = (period=0.125, guidance_period=0.4, navigation_period=0.75, second_control_period=0.2)
        unsplit = ramp_run(; rates...)
        split = ramp_run(; rates..., interval=1.5)
        same_ticks_and_state(split, unsplit)
        @test 1.5 in split.navigation && 3.0 in split.navigation && 1.5 in split.control
        @test !any(t -> t in (1.5, 3.0, 4.5), split.guidance)
    end

    @testset "a run resumed after an interruption ticks as one in one piece" begin
        unsplit = ramp_run(period=0.125)
        dir = mktempdir()
        # Interrupted after its checkpoint at 2 s, which it wrote after that
        # time's ticks: a checkpoint at a segment boundary lists no pending clock.
        interrupted = ramp_setup(period=0.125, interval=1.0, dir=dir, fail_after=2.5)
        err = try
            run_simulation(interrupted.args; isolate_state=false)
            nothing
        catch e
            e
        end
        @test err !== nothing && occursin("interrupted at", sprint(showerror, err))
        @test 2.0 in interrupted.control.ticks && SE._load_checkpoint(interrupted.args).t == 2.0
        @test isempty(pending_periods(interrupted.args))
        resumed = ramp_run(period=0.125, interval=1.0, dir=dir, resume=true, fail_after=Inf)
        @test resumed.control[1] == 2.125
        for field in TICKS
            @test getfield(resumed, field) == filter(>(2.0), getfield(unsplit, field))
        end
        @test all(c -> c[2] == 0.0, resumed.clock)
        @test resumed.final.h_wheels[1] ≈ unsplit.final.h_wheels[1] rtol=1e-12
        @test resumed.table.sc1_actuator_level_1[end] == unsplit.table.sc1_actuator_level_1[end]
    end

    @testset "a run stopped by an event at a tick takes the tick before its checkpoint" begin
        # An event at 2.5 s, on the control and guidance grids, ends the checkpointed
        # solve with terminate!, as an impact or an orbit count does. The solver applies
        # the step's discrete callbacks after the event's, so the ticks due then run and
        # the checkpoint written there lists no pending clock.
        unsplit = ramp_run(period=0.125)
        dir = mktempdir()
        stop = configuration -> SC.ContinuousCallback((u, t, integrator) -> t - 2.5, integrator -> SC.terminate!(integrator))
        stopped = ramp_run(period=0.125, interval=1.0, dir=dir, factories=(stop,))
        @test last(stopped.control) == 2.5 && last(stopped.guidance) == 2.5 && last(stopped.navigation) < 2.5
        checkpoint = SE._load_checkpoint(stopped.args)
        @test checkpoint.t == 2.5 && isempty(pending_periods(stopped.args))
        resumed = ramp_run(period=0.125, interval=1.0, dir=dir, resume=true)
        @test resumed.control[1] == 2.625
        for field in TICKS
            @test vcat(getfield(stopped, field), getfield(resumed, field)) == getfield(unsplit, field)
        end
        @test resumed.final.h_wheels[1] ≈ unsplit.final.h_wheels[1] rtol=1e-12
    end

    @testset "a streamed long-layout run extended from its end saves the row due there again" begin
        # The streams hold every row but the last at each checkpoint; the resumed
        # run takes the ticks due at 2 s and saves the row there again, as the
        # run never interrupted does.
        whole_dir, dir = mktempdir(), mktempdir()
        whole = ramp_run(period=0.125, interval=1.0, dir=whole_dir, results=true, layout=:long)
        first_leg = ramp_run(period=0.125, interval=1.0, dir=dir, mission_time=2.0, results=true, layout=:long)
        @test SE._load_checkpoint(first_leg.args).runtime_state.results_rows.version == 2
        @test first_leg.table.time[end] == 2.0 && first_leg.table.actuator_level_1[end] == 0.001875
        resumed = ramp_run(period=0.125, interval=1.0, dir=dir, resume=true, results=true, layout=:long)
        @test isequal(resumed.table, whole.table)
        row = findfirst(==(2.0), resumed.table.time)
        @test resumed.table.time[1] == 0.0 && resumed.table.actuator_level_1[row] == 0.002
        @test read(joinpath(dir, "simulation_results_spacecraft.csv")) ==
            read(joinpath(whole_dir, "simulation_results_spacecraft.csv"))
    end

    @testset "a run extended from the end of a finished run takes the ticks due there" begin
        unsplit = ramp_run(period=0.125)
        dir, again_dir = mktempdir(), mktempdir()
        first_leg = ramp_run(period=0.125, interval=1.0, dir=dir, mission_time=2.0)
        @test last(first_leg.control) == 1.875                               # no tick at the run's end
        @test pending_periods(first_leg.args) == [0.125, 0.25]                # the grids that meet 2 s
        cp(joinpath(dir, "checkpoints"), joinpath(again_dir, "checkpoints"))
        resumed = ramp_run(period=0.125, interval=1.0, dir=dir, resume=true)
        # Only the rates whose grid meets 2 s tick there.
        @test resumed.control[1] == 2.0 && resumed.guidance[1] == 2.0
        @test resumed.second[1] > 2.0 && resumed.navigation[1] > 2.0
        for field in TICKS
            @test vcat(getfield(first_leg, field), getfield(resumed, field)) == getfield(unsplit, field)
        end
        @test resumed.final.h_wheels[1] ≈ unsplit.final.h_wheels[1] rtol=1e-12
        # The resumed run's table holds the whole run, as a checkpointed run never
        # interrupted does. The first leg's last row, saved before the ticks due at
        # 2 s, is saved again after them and holds the level commanded there.
        whole = ramp_run(period=0.125, interval=1.0)
        @test isequal(resumed.table, whole.table)
        @test first_leg.table.time[end] == 2.0 && first_leg.table.sc1_actuator_level_1[end] == 0.001875
        row = findfirst(==(2.0), resumed.table.time)
        @test resumed.table.time[1] == 0.0 && resumed.table.sc1_actuator_level_1[row] == 0.002
        # The rows were written anew without the first leg's last row: a resume of
        # the finished run solves nothing and returns the same table.
        @test isequal(ramp_run(period=0.125, interval=1.0, dir=dir, resume=true).table, whole.table)
        # Resumed with checkpoint writing off (a resume is still checkpointed, here in
        # one segment): the same ticks.
        again = ramp_run(period=0.125, dir=again_dir, resume=true)
        @test again.control == resumed.control && again.guidance == resumed.guidance
        @test again.final.h_wheels[1] ≈ unsplit.final.h_wheels[1] rtol=1e-12
        # The ticks a resumed run takes at its start, by entry: a version 2 entry
        # lists the clocks; 0.4.0's version 1 entry (at the run's end) and its
        # checkpoints without an entry (at a boundary or an early stop) leave pending,
        # for every clock, what 0.4.0 did not take: the ticks at the run's end, and
        # those just after the time. An entry that does not belong to the checkpoint's
        # time, or of another shape, is refused.
        written = SE._load_checkpoint(first_leg.args)                        # the extended run's end
        listed = SE._pending_periodic_ticks(written)
        @test listed.at_s == written.t && sort(last.(listed.clocks)) == [0.125, 0.25]
        with_entry(entry) = merge(written, (t=2.0, runtime_state=(periodic_ticks=entry,),))
        @test SE._pending_periodic_ticks(with_entry((version=2, time_s=2.0, pending=[(clock=2, period_s=0.125)]))) ==
            (at_s=2.0, from_s=-Inf, clocks=[(2, 0.125)])
        @test SE._pending_periodic_ticks(with_entry((version=2, time_s=2.0, pending=[]))).clocks == Tuple{Int, Float64}[]
        @test SE._pending_periodic_ticks(with_entry((version=1, pending_s=2.0))) == (at_s=2.0, from_s=2.0, clocks=nothing)
        @test SE._pending_periodic_ticks((t=2.0, u=nothing)) == (at_s=2.0, from_s=nextfloat(2.0), clocks=nothing)
        for bad in ((version=2, time_s=1.0, pending=[]), (version=2, time_s=2.0, pending=true),
                    (version=2, time_s=2.0, pending=[(clock=1.0, period_s=0.1)]),
                    (version=2, time_s=2.0, pending=[(position=1, period_s=0.1)]),
                    (version=1, pending_s=1.0), (version=3, time_s=2.0, pending=[]), (pending_s=2.0,))
            @test_throws ArgumentError SE._pending_periodic_ticks(with_entry(bad))
        end
    end

    @testset "an extension that fails keeps the finished run's rows and results" begin
        whole = ramp_run(period=0.125, interval=1.0, fail_after=Inf)
        # Its first checkpoint cannot be written. The rows file is only appended
        # to, so the finished run's checkpoint still finds its own rows, and a
        # later resume holds the whole run.
        dir = mktempdir()
        ramp_run(period=0.125, interval=1.0, dir=dir, mission_time=2.0, fail_after=Inf)
        checkpoints = joinpath(dir, "checkpoints")
        chmod(checkpoints, 0o555)
        read_only = try
            touch(joinpath(checkpoints, "probe"))
            false
        catch
            true
        end
        try
            if read_only
                @test_throws Exception ramp_run(period=0.125, interval=1.0, dir=dir, resume=true, fail_after=Inf)
            else
                @test_skip "the checkpoint directory stays writable here"
            end
        finally
            chmod(checkpoints, 0o755)
            rm(joinpath(checkpoints, "probe"); force=true)
        end
        @test isequal(ramp_run(period=0.125, interval=1.0, dir=dir, resume=true, fail_after=Inf).table, whole.table)
        # Its first segment fails: the partial results it saves are the finished
        # run's, its row at 2 s included, which the extension had not yet saved again.
        dir = mktempdir()
        ramp_run(period=0.125, interval=1.0, dir=dir, mission_time=2.0, fail_after=Inf, results=true)
        csv = joinpath(dir, "simulation_results.csv")
        finished = read(csv)
        @test_throws Exception ramp_run(period=0.125, interval=1.0, dir=dir, resume=true, fail_after=2.3, results=true)
        @test read(csv) == finished
        @test isequal(ramp_run(period=0.125, interval=1.0, dir=dir, resume=true, fail_after=Inf).table, whole.table)
    end

    @testset "decimal boundaries take the tick they round away from, once" begin
        # 3 * 0.1 is 0.30000000000000004, and a 0.3 s interval ends the first segment
        # at 0.3. The next segment started an ulp before that tick, and the solver
        # aborted the step to it as unstable after two navigation ticks.
        rates = (period=0.125, navigation_period=0.1, mission_time=1.5)
        unsplit = ramp_run(; rates...)
        split = ramp_run(; rates..., interval=0.3)
        @test checkpoint_boundaries(0.3, 1.5) == [0.3, 0.6, 3 * 0.3, 1.2, 1.5]
        @test length(split.navigation) == length(unsplit.navigation) == 14
        @test all(b -> count(==(b), split.navigation) == 1, (0.3, 0.6, 3 * 0.3, 1.2))
        @test 3 * 0.1 in unsplit.navigation && 3 * 0.1 ∉ split.navigation
        same_ticks_and_state(split, unsplit)
        @test pending_periods(split.args) == [0.1, 0.125, 0.25, 0.3]          # every grid meets 1.5 s

        # Boundaries an ulp before ticks of the 0.25 s grid; the third falls an ulp
        # before the run's end, which ends that segment instead.
        rates = (period=0.125, navigation_period=0.25, mission_time=1.5)
        unsplit = ramp_run(; rates...)
        split = ramp_run(; rates..., interval=prevfloat(0.5))
        @test checkpoint_boundaries(prevfloat(0.5), 1.5) == [prevfloat(0.5), prevfloat(1.0), 1.5]
        @test split.navigation == [0.25, prevfloat(0.5), 0.75, prevfloat(1.0), 1.25]
        same_ticks_and_state(split, unsplit)

        # Ten 0.1 s intervals reach 0.9999999999999999, an ulp before the run's end;
        # a last segment that short aborted as unstable as well.
        rates = (period=0.125, navigation_period=0.1, mission_time=1.0)
        unsplit = ramp_run(; rates...)
        split = ramp_run(; rates..., interval=0.1)
        @test length(checkpoint_boundaries(0.1, 1.0)) == 10 && last(checkpoint_boundaries(0.1, 1.0)) == 1.0
        same_ticks_and_state(split, unsplit)
        @test SE._load_checkpoint(split.args).t == 1.0

        # A run resumed after an interruption at the 0.6 s boundary: the tick at
        # 6 * 0.1 came before that checkpoint, at 0.6, and does not come again.
        rates = (period=0.125, navigation_period=0.1, mission_time=1.5)
        split = ramp_run(; rates..., interval=0.3)
        dir = mktempdir()
        interrupted = ramp_setup(; rates..., interval=0.3, dir=dir, fail_after=0.65)
        err = try
            run_simulation(interrupted.args; isolate_state=false)
            nothing
        catch e
            e
        end
        @test err !== nothing && occursin("interrupted at", sprint(showerror, err))
        @test SE._load_checkpoint(interrupted.args).t == 0.6 && count(==(0.6), interrupted.navigation.ticks) == 1
        resumed = ramp_run(; rates..., interval=0.3, dir=dir, resume=true, fail_after=Inf)
        @test resumed.navigation[1] == 7 * 0.1
        for field in TICKS
            @test getfield(resumed, field) == filter(>(0.6), getfield(split, field))
        end
        @test resumed.final.h_wheels[1] ≈ split.final.h_wheels[1] rtol=1e-12

        # A run ending an ulp before a tick leaves it pending; the run extended from
        # its end takes it there, once.
        unsplit = ramp_run(; rates..., mission_time=0.6)
        dir = mktempdir()
        first_leg = ramp_run(; rates..., interval=0.3, dir=dir, mission_time=0.3)
        @test first_leg.navigation == [0.1, 0.2]
        @test pending_periods(first_leg.args) == [0.1, 0.3]                    # 3 * 0.1 and 0.3
        extended = ramp_run(; rates..., interval=0.3, dir=dir, resume=true, mission_time=0.6)
        @test extended.navigation == [0.3, 0.4, 0.5]
        @test vcat(first_leg.navigation, extended.navigation) == at_boundaries(unsplit.navigation, [0.3, 0.6])
        @test extended.final.h_wheels[1] ≈ unsplit.final.h_wheels[1] rtol=1e-12
    end

    @testset "a checkpoint from before the entry leaves pending the tick 0.4.0 did not take" begin
        # SpaceAGORA 0.4.0 wrote no entry at a segment boundary and compared tick
        # times exactly, so at a boundary an ulp before 3 * 0.1 it had not taken that
        # tick. Such a checkpoint is made here from a run ending at 0.3 s, where no
        # other model ticks, with its entry removed. Resumed, it takes the tick at its
        # start, once, and the two runs tick as one run in one piece; so does the same
        # checkpoint with its entry.
        rates = (period=0.125, navigation_period=0.1, second_control_period=0.2, mission_time=1.5)
        unsplit = ramp_run(; rates...)
        boundaries = checkpoint_boundaries(0.3, 1.5)
        for written_by_040 in (true, false)
            dir = mktempdir()
            first_leg = ramp_run(; rates..., interval=0.3, dir=dir, mission_time=0.3)
            @test !any(field -> any(t -> SC._same_run_time(t, 0.3), getfield(first_leg, field)), TICKS)
            if written_by_040
                paths = SE._checkpoint_paths(first_leg.args)
                payload = deserialize(paths.data)
                payload = merge(payload, (runtime_state=Base.structdiff(payload.runtime_state, NamedTuple{(:periodic_ticks,)}),))
                open(io -> serialize(io, payload), paths.data, "w")
                @test !hasproperty(runtime_state(first_leg.args), :periodic_ticks)
            end
            resumed = ramp_run(; rates..., interval=0.3, dir=dir, resume=true)
            @test resumed.navigation[1] == 0.3 && count(t -> SC._same_run_time(t, 0.3), resumed.navigation) == 1
            for field in TICKS
                @test vcat(getfield(first_leg, field), getfield(resumed, field)) ==
                    at_boundaries(getfield(unsplit, field), boundaries)
            end
            @test resumed.final.h_wheels[1] ≈ unsplit.final.h_wheels[1] rtol=1e-12
        end

        # A run that went on past that boundary took the tick there and wrote it as
        # taken: resumed from that checkpoint, it does not take it again.
        dir = mktempdir()
        interrupted = ramp_setup(; rates..., interval=0.3, dir=dir, fail_after=0.35)
        err = try
            run_simulation(interrupted.args; isolate_state=false)
            nothing
        catch e
            e
        end
        @test err !== nothing && occursin("interrupted at", sprint(showerror, err))
        @test isempty(pending_periods(interrupted.args))
        before = filter(<=(0.3), interrupted.navigation.ticks)
        @test count(==(0.3), before) == 1
        resumed = ramp_run(; rates..., interval=0.3, dir=dir, resume=true, fail_after=Inf)
        @test resumed.navigation[1] == 4 * 0.1
        @test vcat(before, resumed.navigation) == at_boundaries(unsplit.navigation, boundaries)
        @test resumed.final.h_wheels[1] ≈ unsplit.final.h_wheels[1] rtol=1e-12
    end

    @testset "an early stop records which clocks had not taken the tick due then" begin
        # An extra callback, which runs after the periodic ones, stops the solve.
        stop_at(stop_s) = configuration -> SC.DiscreteCallback((u, t, integrator) -> t == stop_s,
            integrator -> SC.terminate!(integrator); initialize=(c, u, t, integrator) -> SC.add_tstop!(integrator, stop_s))
        unsplit = ramp_run(period=0.1, mission_time=1.5)

        # Stopped at 0.3 s: the 0.3 s navigation and second control ticks ran there,
        # before the stop, but the 0.1 s control tick at 3 * 0.1, an ulp later, did
        # not. The checkpoint lists that clock alone. The resumed run takes its tick
        # at its start, once, repeats no other clock's, and the two runs tick and end
        # as the run in one piece.
        dir = mktempdir()
        stopped = ramp_run(period=0.1, interval=1.0, mission_time=1.5, dir=dir, factories=(stop_at(0.3),))
        @test SE._load_checkpoint(stopped.args).t == 0.3 && pending_periods(stopped.args) == [0.1]
        @test last(stopped.navigation) == 0.3 && last(stopped.second) == 0.3 && last(stopped.control) == 0.2
        resumed = ramp_run(period=0.1, interval=1.0, mission_time=1.5, dir=dir, resume=true)
        @test resumed.control[1] == 0.3 && resumed.navigation[1] == 0.6 && resumed.second[1] == 0.6
        for field in TICKS
            @test vcat(getfield(stopped, field), getfield(resumed, field)) ==
                at_boundaries(getfield(unsplit, field), vcat(0.3, checkpoint_boundaries(1.0, 1.5; start=0.3)))
        end
        @test resumed.final.h_wheels[1] ≈ unsplit.final.h_wheels[1] rtol=1e-12

        # Stopped an ulp before the run's 1 s end: the ticks at 1 s belong to the run's
        # end, which the run does not take. The checkpoint lists the clocks whose grids
        # meet 1 s, and the run extended to 1.5 s takes their ticks at its start.
        dir = mktempdir()
        stopped = ramp_run(period=0.1, interval=1.0, mission_time=1.0, dir=dir, factories=(stop_at(prevfloat(1.0)),))
        @test SE._load_checkpoint(stopped.args).t == prevfloat(1.0) && pending_periods(stopped.args) == [0.1, 0.25]
        @test last(stopped.control) == 0.9 && last(stopped.guidance) == 0.75
        resumed = ramp_run(period=0.1, interval=1.0, mission_time=1.5, dir=dir, resume=true)
        @test resumed.control[1] == prevfloat(1.0) && resumed.guidance[1] == prevfloat(1.0)
        @test resumed.navigation[1] == 4 * 0.3 && resumed.second[1] == 4 * 0.3
        for field in TICKS
            @test vcat(getfield(stopped, field), getfield(resumed, field)) ==
                at_boundaries(getfield(unsplit, field), [prevfloat(1.0)])
        end
        @test resumed.final.h_wheels[1] ≈ unsplit.final.h_wheels[1] rtol=1e-12
    end
end

@testset "periodic GNC callback outside a checkpointed run is PeriodicCallback" begin
    function direct_ticks(make_callback, tspan; origin=NaN, horizon=NaN, pending=NaN, pending_from=-Inf,
                          pending_clocks=nothing, before=())
        p = SM.ODEParams(n_sats=1, args=wheel_config(5.0))
        p.shared_buffers.periodic_clock_origin_s[] = origin
        p.shared_buffers.periodic_clock_horizon_s[] = horizon
        p.shared_buffers.periodic_clock_pending_s[] = pending
        p.shared_buffers.periodic_clock_pending_from_s[] = pending_from
        p.shared_buffers.periodic_clock_pending_clocks[] = pending_clocks
        ticks = Float64[]
        cb = make_callback(integrator -> push!(ticks, integrator.t))
        prob = SC.ODEProblem((du, u, p, t) -> fill!(du, 0.0), [1.0], tspan, p)
        sol = SC.solve(prob, SC.Tsit5(); callback=SC.CallbackSet(before..., cb))
        return ticks, sol.t
    end
    for tspan in ((0.0, 1.0), (12.345, 13.0), (0.0, 0.95))
        reference = direct_ticks(f -> SC.PeriodicCallback(f, 0.1), tspan)
        @test direct_ticks(f -> SC._run_clock_periodic_callback(f, 0.1), tspan) == reference
        @test !isempty(first(reference)) && last(tspan) ∉ first(reference)
    end
    # A later segment of a checkpointed run: the run's grid, and its tick at
    # the segment's end, which is before the run's end.
    origins = Float64[]
    recorded(f) = SC._run_clock_periodic_callback(f, 0.125; on_clock=(origin, integrator) -> push!(origins, origin))
    ticks, _ = direct_ticks(recorded, (1.0, 2.0); origin=0.0, horizon=5.0)
    @test ticks == 1.125:0.125:2.0
    @test origins == [0.0]
    ticks, _ = direct_ticks(recorded, (4.0, 5.0); origin=0.0, horizon=5.0)
    @test ticks == 4.125:0.125:4.875
    _ = direct_ticks(recorded, (12.345, 13.0))
    @test origins == [0.0, 0.0, 12.345]
    # The ticks pending at the start of a resumed run's first solve: only on the
    # grid, and only at the pending time.
    ticks, _ = direct_ticks(recorded, (2.0, 3.0); origin=0.0, horizon=5.0, pending=2.0)
    @test ticks == 2.0:0.125:3.0
    ticks, _ = direct_ticks(recorded, (2.0, 3.0); origin=0.0, horizon=5.0, pending=1.0)
    @test ticks == 2.125:0.125:3.0
    ticks, _ = direct_ticks(f -> SC._run_clock_periodic_callback(f, 0.3), (2.0, 3.0); origin=0.0, horizon=5.0, pending=2.0)
    @test ticks == [k * 0.3 for k in 7:10]
    # A solve terminated at a tick, by an event or by an earlier discrete callback,
    # takes that tick (PeriodicCallback threw there, scheduling a stop after
    # terminate! had emptied the stops).
    event = SC.ContinuousCallback((u, t, integrator) -> t - 2.5, integrator -> SC.terminate!(integrator))
    at_stop = SC.DiscreteCallback((u, t, integrator) -> t == 2.5, integrator -> SC.terminate!(integrator);
        initialize=(c, u, t, integrator) -> SC.add_tstop!(integrator, 2.5))
    for terminator in (event, at_stop), clock in ((origin=0.0, horizon=5.0), (origin=NaN, horizon=NaN))
        ticks, times = direct_ticks(recorded, (2.0, 3.0); before=(terminator,), clock...)
        @test ticks == 2.125:0.125:2.5 && last(times) == 2.5
    end

    # Run times that differ only by rounding are one time to the run clock: within
    # 100 ulps of the larger, OrdinaryDiffEqCore's tolerance for reaching a stop.
    @test SC._same_run_time(0.3, 3 * 0.1) && SC._same_run_time(3 * 0.1, 0.3)
    @test SC._same_run_time(1.0e7, 1.0e7 + 100 * eps(1.0e7)) && !SC._same_run_time(1.0e7, 1.0e7 + 101 * eps(1.0e7))
    @test SE._checkpoint_segment_end(0.0, 0.3, 1.5) == 0.3 && SE._checkpoint_segment_end(1.2, 0.3, 1.5) == 1.5
    @test SE._checkpoint_segment_end(0.8999999999999999, 0.1, 1.0) == 1.0
    @test SE._checkpoint_segment_end(0.0, 1.0 - 1e-9, 1.0) == 1.0 - 1e-9
    # Segments ending near 3 * 0.1, which is 0.30000000000000004. Within that
    # tolerance of the boundary, before or after it, the tick fires at the boundary,
    # at the first segment's end, and the second segment starts after it; farther,
    # it fires at its own time. It fires once either way, and no solve has a step
    # shorter than the tolerance.
    tick = 3 * 0.1
    tolerance = 100 * eps(tick)
    clock = (origin=0.0, horizon=5.0)
    every_tenth(f) = SC._run_clock_periodic_callback(f, 0.1)
    for (boundary, at_boundary) in ((tick - 1e-9, false), (tick - 101 * eps(tick), false), (tick - tolerance, true),
                                    (prevfloat(tick), true), (tick, true), (nextfloat(tick), true),
                                    (tick + tolerance, true), (tick + 101 * eps(tick), false), (tick + 1e-9, false))
        before, before_times = direct_ticks(every_tenth, (0.0, boundary); clock...)
        after, after_times = direct_ticks(every_tenth, (boundary, 0.55); clock...)
        @test vcat(before, after) == [0.1, 0.2, at_boundary ? boundary : tick, 0.4, 0.5]
        @test minimum(diff(unique(before_times))) > tolerance && minimum(diff(unique(after_times))) > tolerance
    end
    # A run ending that close to a tick, before or after it, leaves the tick pending,
    # and a run resumed there takes it at its start, once.
    for run_end in (prevfloat(tick), tick + tolerance)
        ticks, _ = direct_ticks(every_tenth, (0.0, run_end); origin=0.0, horizon=run_end)
        @test ticks == [0.1, 0.2]
        ticks, _ = direct_ticks(every_tenth, (run_end, 0.55); clock..., pending=run_end)
        @test ticks == [run_end, 0.4, 0.5]
    end
    # A checkpoint written by 0.4.0 near that tick, read as it left it: without an
    # entry (a boundary or an early stop) the tick is pending only when it falls just
    # after the checkpoint's time, since 0.4.0 took one at or just before it; with
    # 0.4.0's entry at the run's end, also when it falls at that time.
    legacy(t) = SE._pending_periodic_ticks((t=t, u=nothing))
    version1(t) = SE._pending_periodic_ticks((t=t, u=nothing, runtime_state=(periodic_ticks=(version=1, pending_s=t),)))
    for (read, start, pending) in ((legacy, prevfloat(tick), true), (legacy, tick, false), (legacy, nextfloat(tick), false),
                                   (version1, prevfloat(tick), true), (version1, tick, true), (version1, nextfloat(tick), false))
        due = read(start)
        ticks, _ = direct_ticks(every_tenth, (start, 0.55); clock..., pending=due.at_s, pending_from=due.from_s)
        @test ticks == (pending ? [start, 0.4, 0.5] : [0.4, 0.5])
    end
    # A checkpoint that lists its clocks: only a listed clock, at its position and
    # with its period, takes the tick at the start; the callback here is the run's
    # first clock.
    for (listed, pending) in (([(1, 0.1)], true), ([(1, 0.2)], false), ([(2, 0.1)], false),
                              (Tuple{Int, Float64}[], false))
        ticks, _ = direct_ticks(every_tenth, (prevfloat(tick), 0.55); clock..., pending=prevfloat(tick),
            pending_clocks=listed)
        @test ticks == (pending ? [prevfloat(tick), 0.4, 0.5] : [0.4, 0.5])
    end
    # A solve stopped by `terminate!` that close to a tick takes it where it stopped,
    # so a run resumed there starts after it. Without the run clock it does not, as
    # PeriodicCallback does not.
    stop_s = tick - 50 * eps(tick)
    near_stop = SC.DiscreteCallback((u, t, integrator) -> t == stop_s, integrator -> SC.terminate!(integrator);
        initialize=(c, u, t, integrator) -> SC.add_tstop!(integrator, stop_s))
    ticks, times = direct_ticks(every_tenth, (0.0, 1.0); before=(near_stop,), clock...)
    @test ticks == [0.1, 0.2, stop_s] && last(times) == stop_s
    @test first(direct_ticks(every_tenth, (stop_s, 0.55); clock...)) == [0.4, 0.5]
    reference = first(direct_ticks(f -> SC.PeriodicCallback(f, 0.1), (0.0, 1.0); before=(near_stop,)))
    @test first(direct_ticks(every_tenth, (0.0, 1.0); before=(near_stop,))) == reference == [0.1, 0.2]
end

end
