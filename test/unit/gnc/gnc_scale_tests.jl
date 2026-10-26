# How the periodic GNC ticks visit spacecraft: each tick reads a model's scope once
# (a linked consumer's twice, once before its producers prepare) and ticks only the
# positions it names, in ascending order, so the work of a tick follows the scope
# rather than the number of spacecraft. Counted, not timed: the ticks and releases
# per spacecraft, the scope reads per tick, the order of every tick and stop
# request, and the linked producers' preparation.
module GNCScaleTests
using Test, StaticArrays, SpaceAGORA
using DiffEqCallbacks: PresetTimeCallback
const SM = SpaceAGORA.SimulationModel
const CB = SM.SimulationCallbacks
const ZERO3 = SVector(0.0, 0.0, 0.0)

# Every tick, release and stop request, in order: (time, kind, model id, spacecraft).
const LOG = Tuple{Float64, Symbol, Int, Int}[]
record!(t, kind, id, i) = push!(LOG, (Float64(t), kind, id, i))

# The bound copies a template makes, so a test can read them after the run.
mutable struct Copies
    models::Vector{Any}
end
Copies() = Copies(Any[])

# A model of each kind that counts its scope reads, ticks and releases. Index 0 is
# a template, bound to a spacecraft by `bind_spacecraft`; `scope` is the scope a
# model declares when it is given (configuration level).
mutable struct Counter
    reads::Int
    ticks::Vector{Int}
    released::Vector{Int}
end
Counter(n) = Counter(0, zeros(Int, n), zeros(Int, n))

for (T, Abstract, tick, kind) in ((:CountNav, :AbstractNavigationModel, :update_navigation!, :nav),
                                  (:CountGd, :AbstractGuidanceModel, :update_guidance!, :gd),
                                  (:CountCtl, :(SM.AbstractControlEffectorModel), :update_control!, :ctl))
    @eval begin
        mutable struct $T <: $Abstract
            id::Int
            spacecraft_index::Int
            scope::Any
            skip::Bool
            n::Int
            counter::Counter
            copies::Copies
        end
        $T(id, n; scope=nothing, skip=true, copies=Copies()) = $T(id, 0, scope, skip, n, Counter(n), copies)
        function SpaceAGORA.$tick(m::$T, u, p, t::Float64, i::Int)
            m.counter.ticks[i] += 1
            record!(t, $(QuoteNode(kind)), m.id, i)
            return nothing
        end
        function SpaceAGORA.spacecraft_scope(m::$T)
            m.counter.reads += 1
            m.scope === nothing || return m.scope
            return m.spacecraft_index == 0 ? AllSpacecraft() : BoundSpacecraft((m.spacecraft_index,))
        end
        SpaceAGORA.skip_inactive(m::$T) = m.skip
        function SpaceAGORA.release_command!(m::$T, t, i)
            m.counter.released[i] += 1
            record!(t, :release, m.id, i)
            return nothing
        end
        function SpaceAGORA.bind_spacecraft(m::$T, k::Int)
            bound = $T(m.id + k, k, nothing, m.skip, m.n, Counter(m.n), m.copies)
            push!(m.copies.models, bound)
            return bound
        end
    end
end
SpaceAGORA.control_wrench(::CountCtl, x::StateSample, p, t::Float64) = (ZERO3, ZERO3)
function SpaceAGORA.register_tstops!(integrator, m::CountCtl, i)
    record!(integrator.t, :stop, m.id, i)
    return nothing
end

# A controller whose initializer widens its scope past the run, after the setup
# check: the ticks skip the positions the run does not have.
mutable struct WidenedCtl <: SM.AbstractControlEffectorModel
    id::Int
    indices::Tuple{Vararg{Int}}
    widened::Tuple{Vararg{Int}}
end
SpaceAGORA.update_control!(m::WidenedCtl, u, p, t::Float64, i::Int) = (record!(t, :ctl, m.id, i); nothing)
SpaceAGORA.control_wrench(::WidenedCtl, x::StateSample, p, t::Float64) = (ZERO3, ZERO3)
SpaceAGORA.spacecraft_scope(m::WidenedCtl) = BoundSpacecraft(m.indices)
SpaceAGORA.initialize_control!(m::WidenedCtl, u, p, t) = (m.indices = m.widened; nothing)
SpaceAGORA.register_tstops!(integrator, m::WidenedCtl, i) = (record!(integrator.t, :stop, m.id, i); nothing)

# A linked pair: a guidance model that prepares the plan a bound controller reads.
struct Plan
    prepared::Base.RefValue{Int}
end
struct PlanGd <: AbstractGuidanceModel
    plan::Plan
end
SpaceAGORA.update_guidance!(::PlanGd, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.gnc_links(g::PlanGd) = (Produces(g.plan),)
SpaceAGORA.before_reference_control!(g::PlanGd, c, u, p, t) = (g.plan.prepared[] += 1; nothing)
struct PlanCtl <: SM.AbstractControlEffectorModel
    spacecraft_index::Int
    skip::Bool
    plan::Plan
    ticks::Vector{Int}
end
SpaceAGORA.update_control!(m::PlanCtl, u, p, t::Float64, i::Int) = (m.ticks[i] += 1; nothing)
SpaceAGORA.control_wrench(::PlanCtl, x::StateSample, p, t::Float64) = (ZERO3, ZERO3)
SpaceAGORA.spacecraft_scope(m::PlanCtl) = BoundSpacecraft((m.spacecraft_index,))
SpaceAGORA.skip_inactive(m::PlanCtl) = m.skip
SpaceAGORA.gnc_links(m::PlanCtl) = (Consumes(m.plan),)

# A controller whose producer retargets it before each tick, which changes its
# scope during the run: the tick reads the scope after the producers prepared it.
mutable struct RetargetCtl <: SM.AbstractControlEffectorModel
    spacecraft_index::Int
    plan::Plan
    ticks::Vector{Int}
end
SpaceAGORA.update_control!(m::RetargetCtl, u, p, t::Float64, i::Int) = (m.ticks[i] += 1; nothing)
SpaceAGORA.control_wrench(::RetargetCtl, x::StateSample, p, t::Float64) = (ZERO3, ZERO3)
SpaceAGORA.spacecraft_scope(m::RetargetCtl) = BoundSpacecraft((m.spacecraft_index,))
SpaceAGORA.gnc_links(m::RetargetCtl) = (Consumes(m.plan),)
SpaceAGORA.before_reference_control!(g::PlanGd, c::RetargetCtl, u, p, t) = (c.spacecraft_index = 3; nothing)

# A controller that counts its scope reads, linked to a producer or not.
mutable struct ReadCountCtl{S} <: SM.AbstractControlEffectorModel
    scope::S
    plan::Union{Nothing, Plan}
    reads::Int
end
SpaceAGORA.update_control!(m::ReadCountCtl, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.control_wrench(::ReadCountCtl, x::StateSample, p, t::Float64) = (ZERO3, ZERO3)
SpaceAGORA.spacecraft_scope(m::ReadCountCtl) = (m.reads += 1; m.scope)
SpaceAGORA.gnc_links(m::ReadCountCtl) = m.plan === nothing ? () : (Consumes(m.plan),)

# The bytes one call of `_ticks_any` allocates, once compiled.
ticks_any_bytes(m, p, n) = (CB._ticks_any(m, p, n); @allocated CB._ticks_any(m, p, n))

# A thread-safe controller, for every spacecraft or for the positions given.
struct ThreadedCtl <: SM.AbstractControlEffectorModel
    ticks::Vector{Int}
    scope::Any
end
ThreadedCtl(n; scope=AllSpacecraft()) = ThreadedCtl(zeros(Int, n), scope)
SpaceAGORA.update_control!(m::ThreadedCtl, u, p, t::Float64, i::Int) = (m.ticks[i] += 1; nothing)
SpaceAGORA.control_wrench(::ThreadedCtl, x::StateSample, p, t::Float64) = (ZERO3, ZERO3)
SpaceAGORA.spacecraft_scope(m::ThreadedCtl) = m.scope
CB.control_model_threadsafe(::ThreadedCtl) = true

const PLANET = SM.make_no_gram_planet(:earth)

# `n` spacecraft, each declaring `per_spacecraft` (a guidance, navigation and
# control tuple of templates, bound to it at the run's start), plus the
# configuration-level models. Every GNC rate is 1 s.
function scale_args(n; per_spacecraft=((), (), ()), control=(), guidance=(), navigation=(), mission_time=4.5)
    g, nav, c = per_spacecraft
    sats = map(1:n) do id
        root = SM.Link(root=true, m=100.0, ref_area=2.0)
        ic = SM.InitialCondition(ra=PLANET.Rp_e + 700e3, rp=PLANET.Rp_e + 650e3, i=45.0, ω=0.0, Ω=30.0id, ν=10.0)
        SM.SpacecraftModel(joints=SM.Joint[], links=[root], root=root, prop_mass=0.0, inertia_tensor=root.inertia,
            initial_condition=ic, id=id,
            guidance=SM.GuidanceModel(Tuple(g), fill(1.0, length(g))),
            navigation=SM.NavigationModel(Tuple(nav), fill(1.0, length(nav))),
            control=SM.ControlModel(Tuple(c), fill(1.0, length(c))))
    end
    env = SM.EnvironmentModel(planet=PLANET, EI=120.0, density_model=SM.NoAtmosphereModel(),
        ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=PLANET), topography=false, wind=false)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false, generate_plots=false, normalize=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=false, num_steps_to_save=100),
        environment_model=env,
        dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[sats...], (SM.InverseSquaredJ2GravityModel(),)),
        guidance_model=SM.GuidanceModel(Tuple(guidance), fill(1.0, length(guidance))),
        navigation_model=SM.NavigationModel(Tuple(navigation), fill(1.0, length(navigation))),
        control_model=SM.ControlModel(control_effectors=Tuple(control), control_rates=fill(1.0, length(control))),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

deactivate(i, t_s) = PresetTimeCallback([t_s], integrator -> (integrator.p.is_active[i] = false); save_positions=(false, false))

# A run of `n` spacecraft, each with a bound navigation, guidance and control model,
# beside one configuration-level controller for every spacecraft. Returns the
# templates' bound copies by kind and the shared controller.
function bound_run(n; mission_time=4.5)
    templates = (CountGd(200, n), CountNav(100, n), CountCtl(300, n))
    shared = CountCtl(1, n)
    run_simulation(scale_args(n; per_spacecraft=((templates[1],), (templates[2],), (templates[3],)), control=(shared,),
        mission_time=mission_time); isolate_state=false)
    copies = map(t -> sort(t.copies.models; by=m -> m.spacecraft_index), templates)
    return (guidance=copies[1], navigation=copies[2], control=copies[3], shared=shared)
end

# The order a run ticks in, from the declarations alone: at each tick time the
# navigation models, then the control models (each followed by its stop
# requests), then the guidance models, each in configuration order (configuration
# level first, then each spacecraft's own), and each over its scope's positions in
# ascending order.
function expected_log(times, navigation, control, guidance)
    log = Tuple{Float64, Symbol, Int, Int}[]
    for t in times
        for (id, positions) in navigation, i in positions
            push!(log, (t, :nav, id, i))
        end
        for (id, positions) in control
            for i in positions
                push!(log, (t, :ctl, id, i))
            end
            for i in positions
                push!(log, (t, :stop, id, i))
            end
        end
        for (id, positions) in guidance, i in positions
            push!(log, (t, :gd, id, i))
        end
    end
    return log
end

@testset "GNC ticks follow the scope" begin
    @testset "the positions a scope ticks" begin
        tp = CB._tick_positions
        @test collect(tp(AllSpacecraft(), 4)) == [1, 2, 3, 4]
        @test tp(BoundSpacecraft((2, 4)), 4) === (2, 4)                 # increasing and in range: as given
        @test tp(BoundSpacecraft((3, 1, 3)), 4) == [1, 3]               # sorted, each once
        @test tp(BoundSpacecraft((5, 2, 0)), 4) == [2]                  # only positions the run has
        @test isempty(tp(BoundSpacecraft(()), 4))
        @test tp(CoupledSpacecraft(3, (1,)), 4) === (3,)                # the primary only
        @test isempty(tp(CoupledSpacecraft(7, (1,)), 4))
        # The same positions `_ticks_for` accepts, for every scope of a small run.
        for scope in (AllSpacecraft(), BoundSpacecraft((2,)), BoundSpacecraft((4, 2, 2, 9)), CoupledSpacecraft(1, (3,)))
            @test collect(tp(scope, 4)) == [i for i in 1:4 if CB._ticks_for(scope, i)]
        end
    end

    @testset "a bound model ticks its own spacecraft, and reads its scope once per tick" begin
        short, long = bound_run(3; mission_time=4.5), bound_run(3; mission_time=8.5)
        k = short.shared.counter.ticks[1]
        @test k == 4 && short.shared.counter.ticks == fill(k, 3)          # ticks at 1, 2, 3 and 4 s
        for kind in (:guidance, :navigation, :control), (j, m) in enumerate(getproperty(short, kind))
            @test m.spacecraft_index == j
            @test m.counter.ticks == [i == j ? k : 0 for i in 1:3]
        end
        # Four more ticks, four more scope reads: one per tick, not one per spacecraft.
        for kind in (:guidance, :navigation, :control), j in 1:3
            @test getproperty(long, kind)[j].counter.reads - getproperty(short, kind)[j].counter.reads == 4
        end
        @test long.shared.counter.reads - short.shared.counter.reads == 4
        # And the reads per model do not grow with the constellation.
        wide = bound_run(9; mission_time=4.5)
        for kind in (:guidance, :navigation, :control)
            @test all(m -> m.counter.reads == getproperty(short, kind)[1].counter.reads, getproperty(wide, kind))
        end
        @test wide.shared.counter.reads == short.shared.counter.reads
        # Negative control: the step that visits one spacecraft at a time reads the
        # scope at each call, so the counter sees per-spacecraft reads.
        probe = CountCtl(9, 3)
        p = (is_active=trues(3),)
        before = probe.counter.reads
        for i in 1:3
            CB._gnc_tick!((m, u, p, t, i) -> nothing, probe, nothing, p, 0.0, i)
        end
        @test probe.counter.reads - before == 3
    end

    @testset "ticks and stop requests come in the declared order" begin
        empty!(LOG)
        n = 5
        all_c = CountCtl(1, n)
        odd = CountCtl(2, n; scope=BoundSpacecraft((3, 1, 3)))
        coupled = CountCtl(3, n; scope=CoupledSpacecraft(4, (2,)))
        widened = WidenedCtl(4, (2,), (5, 2, 99))
        shared_nav = CountNav(5, n; scope=BoundSpacecraft((5, 2)))
        shared_gd = CountGd(6, n)
        templates = (CountGd(200, n), CountNav(100, n), CountCtl(300, n))
        run_simulation(scale_args(n; per_spacecraft=((templates[1],), (templates[2],), (templates[3],)),
            control=(all_c, odd, coupled, widened), navigation=(shared_nav,), guidance=(shared_gd,));
            isolate_state=false)
        times = [1.0, 2.0, 3.0, 4.0]
        expected = expected_log(times,
            [(5, [2, 5]); [(100 + j, [j]) for j in 1:n]],
            [(1, 1:n); (2, [1, 3]); (3, [4]); (4, [2, 5]); [(300 + j, [j]) for j in 1:n]],
            [(6, 1:n); [(200 + j, [j]) for j in 1:n]])
        @test LOG == expected
        @test !any(e -> e[4] == 99, LOG)
    end

    @testset "register_tstops! is called only for the spacecraft in the scope" begin
        empty!(LOG)
        n = 4
        models = (CountCtl(1, n; scope=BoundSpacecraft((2, 4))), CountCtl(2, n; scope=CoupledSpacecraft(3, (1,))),
                  CountCtl(3, n))
        in_scope = Dict(1 => [2, 4], 2 => [3], 3 => [1, 2, 3, 4])
        run_simulation(scale_args(n; control=models); isolate_state=false)
        stops = [(e[3], e[4]) for e in LOG if e[2] === :stop]
        @test !isempty(stops)
        @test all(((id, i),) -> i in in_scope[id], stops)                    # none outside the scope
        for (id, positions) in in_scope, i in 1:n
            @test count(==((id, i)), stops) == (i in positions ? 4 : 0)      # one per tick in the scope
        end
    end

    @testset "an inactive spacecraft is released, and linked producers prepare only an acting tick" begin
        empty!(LOG)
        n = 3
        skipping = CountCtl(300, n)
        ticking = CountCtl(400, n; skip=false)
        run_simulation(scale_args(n; per_spacecraft=((), (), (skipping, ticking))); isolate_state=false,
            extra_callbacks=(deactivate(2, 2.5),))
        bound_skip = sort(skipping.copies.models; by=m -> m.spacecraft_index)
        bound_tick = sort(ticking.copies.models; by=m -> m.spacecraft_index)
        @test bound_skip[2].counter.ticks == [0, 2, 0] && bound_skip[2].counter.released == [0, 2, 0]
        @test bound_skip[1].counter.ticks == [4, 0, 0] && iszero(bound_skip[1].counter.released)
        @test bound_tick[2].counter.ticks == [0, 4, 0] && iszero(bound_tick[2].counter.released)
        # Releases come where the ticks would, and are followed by their stop requests.
        @test filter(e -> e[3] == 302, LOG) == [(1.0, :ctl, 302, 2), (1.0, :stop, 302, 2), (2.0, :ctl, 302, 2),
            (2.0, :stop, 302, 2), (3.0, :release, 302, 2), (3.0, :stop, 302, 2), (4.0, :release, 302, 2),
            (4.0, :stop, 302, 2)]
        # The producer prepares the controller only at ticks that act on its spacecraft.
        for (skip, prepared) in ((true, 2), (false, 4))
            plan = Plan(Ref(0))
            consumer = PlanCtl(2, skip, plan, zeros(Int, n))
            run_simulation(scale_args(n; control=(consumer,), guidance=(PlanGd(plan),)); isolate_state=false,
                extra_callbacks=(deactivate(2, 2.5),))
            @test plan.prepared[] == prepared
            @test consumer.ticks == [0, skip ? 2 : 4, 0]
        end
        # A producer that retargets its consumer before the tick changes the
        # consumer's scope during the run, which stops it (SpaceAGORA 0.6 followed the
        # change): the tick reads the scope after the producers prepared and compares
        # it with the one the run took when it started.
        plan = Plan(Ref(0))
        retargeted = RetargetCtl(1, plan, zeros(Int, n))
        @test_throws "control model 1 (RetargetCtl)'s spacecraft_scope changed during the run" run_simulation(
            scale_args(n; control=(retargeted,), guidance=(PlanGd(plan),)); isolate_state=false)
        @test retargeted.ticks == [0, 0, 0]
    end

    @testset "a linked consumer reads its scope twice per tick, once before its producers prepare" begin
        function reads(n, mission_time; linked)
            plan = Plan(Ref(0))
            consumer = ReadCountCtl(BoundSpacecraft((3, 1, 3)), linked ? plan : nothing, 0)
            run_simulation(scale_args(n; control=(consumer,), guidance=linked ? (PlanGd(plan),) : (), mission_time);
                isolate_state=false)
            linked && @test plan.prepared[] == round(Int, mission_time - 0.5)   # it prepared at every tick
            return consumer.reads
        end
        # Four more ticks: eight more reads for a linked consumer, four for another.
        @test reads(3, 8.5; linked=true) - reads(3, 4.5; linked=true) == 8
        @test reads(3, 8.5; linked=false) - reads(3, 4.5; linked=false) == 4
        # Neither grows with the constellation.
        @test reads(9, 4.5; linked=true) == reads(3, 4.5; linked=true)
        @test reads(9, 4.5; linked=false) == reads(3, 4.5; linked=false)
        # The check before the producers prepare answers what the tick's positions
        # would, for unsorted, repeated and out-of-run indices, and builds no copy.
        for scope in (AllSpacecraft(), BoundSpacecraft((3, 1, 3)), BoundSpacecraft((9, 2)), BoundSpacecraft((9,)),
                      CoupledSpacecraft(2, (1,)), CoupledSpacecraft(7, (1,))),
            active in (trues(3), BitVector([false, true, false]), falses(3))
            probe, p = ReadCountCtl(scope, nothing, 0), (is_active=active,)
            @test CB._ticks_any(probe, p, 3) == any(i -> active[i], CB._tick_positions(scope, 3))
            @test ticks_any_bytes(probe, p, 3) == 0
        end
    end

    @testset "threaded ticks visit the same spacecraft" begin
        counts = map(("off", "on")) do mode
            law = ThreadedCtl(6)
            some = ThreadedCtl(6; scope=BoundSpacecraft((2, 4, 6)))
            bound = CountCtl(300, 6)
            withenv("SPACEAGORA_CONTROL_CALLBACK_PARALLEL" => mode, "SPACEAGORA_CONTROL_CALLBACK_THREAD_THRESHOLD" => "2") do
                run_simulation(scale_args(6; per_spacecraft=((), (), (bound,)), control=(law, some)); isolate_state=false)
            end
            (law.ticks, some.ticks, [m.counter.ticks for m in sort(bound.copies.models; by=m -> m.spacecraft_index)])
        end
        @test counts[1] == counts[2]
        @test counts[1][1] == fill(4, 6)
        @test counts[1][2] == [0, 4, 0, 4, 0, 4]
    end
end

end
