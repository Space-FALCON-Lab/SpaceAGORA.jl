# A GNC model's per-spacecraft hooks follow its scope: the output hooks, the levels
# and `atmosphere_transition!` are called only for the spacecraft the model's scope
# names. The run takes each model's scope when it starts, after the run-start
# initializers, and stops when a later read differs; a SpacecraftSelection scopes a
# GNC model and is resolved once. Counted, not timed: hook calls per spacecraft and
# per evaluation, against the same runs with the member lists switched off, which
# also give the same bits.
module GNCScopeReadersTests
using Test, StaticArrays, LinearAlgebra, SpaceAGORA
using DiffEqCallbacks: PresetTimeCallback
import DiffEqBase
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const CB = SM.SimulationCallbacks
const Z3 = SVector(0.0, 0.0, 0.0)
const PLANET = SM.make_no_gram_planet(:earth)

# The right-hand side on one thread, so the hook counters need no locks.
serial(f) = withenv(f, "SPACEAGORA_RHS_EXECUTION_MODE" => "serial")

# A run with the member lists switched off: every model's hooks are called for
# every spacecraft, as before 0.7.0.
function without_lists(f)
    SE._GNC_MEMBER_LISTS[] = false
    try
        return f()
    finally
        SE._GNC_MEMBER_LISTS[] = true
    end
end

same_solution(a, b) = a.t == b.t && length(a.u) == length(b.u) && all(isequal(collect(x), collect(y)) for (x, y) in zip(a.u, b.u))

# A controller that holds a small push along its spacecraft's velocity and counts
# its control_wrench calls per spacecraft. Index 0 acts on every spacecraft (a
# template, or a configuration-level model); `scope` overrides the scope. Its
# output is zero for a spacecraft it never ticks.
mutable struct PushCtl <: SM.AbstractControlEffectorModel
    spacecraft_index::Int
    scope::Any
    held::Vector{SVector{3, Float64}}
    calls::Vector{Int}
    copies::Vector{Any}
end
PushCtl(n; index=0, scope=nothing) = PushCtl(index, scope, fill(Z3, n), zeros(Int, n), Any[])
SpaceAGORA.spacecraft_scope(m::PushCtl) = m.scope !== nothing ? m.scope :
    m.spacecraft_index == 0 ? AllSpacecraft() : BoundSpacecraft((m.spacecraft_index,))
SpaceAGORA.update_control!(m::PushCtl, u, p, t::Float64, i::Int) =
    (m.held[i] = 1.0e-3 * normalize(SVector{3, Float64}(u.sc[i].vel)); nothing)
function SpaceAGORA.control_wrench(m::PushCtl, x::StateSample, p, t::Float64)
    m.calls[x.spacecraft_index] += 1
    return m.held[x.spacecraft_index], Z3
end
SpaceAGORA.initialize_control!(m::PushCtl, u, p, t) = (fill!(m.held, Z3); fill!(m.calls, 0); nothing)
function SpaceAGORA.bind_spacecraft(m::PushCtl, k::Int)
    bound = PushCtl(length(m.held); index=k)
    push!(m.copies, bound)
    return bound
end

# Every control_wrench call, in order: (spacecraft, model id).
const CALLS = Tuple{Int, Int}[]

# A controller that logs its control_wrench calls. `id` names it; `scope` is its
# scope, and a template's bound copies are bound to their spacecraft.
struct OrderCtl <: SM.AbstractControlEffectorModel
    id::Int
    scope::Any
end
SpaceAGORA.spacecraft_scope(m::OrderCtl) = m.scope
SpaceAGORA.update_control!(m::OrderCtl, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.control_wrench(m::OrderCtl, x::StateSample, p, t::Float64) = (push!(CALLS, (x.spacecraft_index, m.id)); (Z3, Z3))
SpaceAGORA.bind_spacecraft(m::OrderCtl, k::Int) = OrderCtl(m.id + k, BoundSpacecraft((k,)))

# `n` spacecraft with J2, each declaring `declared` and `declared_guidance` (control
# and guidance templates bound to it when the run starts), beside the
# configuration-level `control` and `guidance`. With `drag`, an exponential
# atmosphere and orbits that dip into it shortly after the start.
function scope_args(n; declared=(), declared_guidance=(), control=(), guidance=(), mission_time=6.0, drag=false,
                    spacecraft=nothing, settings=(;))
    sats = spacecraft !== nothing ? spacecraft : map(1:n) do id
        root = SM.Link(root=true, m=100.0, ref_area=2.0)
        ic = SM.InitialCondition(ra=PLANET.Rp_e + (drag ? 2000e3 : 700e3), rp=PLANET.Rp_e + (drag ? 140e3 : 650e3),
            i=45.0, ω=0.0, Ω=30.0id, ν=drag ? 320.0 + 2.0id : 10.0)
        SM.SpacecraftModel(joints=SM.Joint[], links=[root], root=root, prop_mass=10.0, inertia_tensor=root.inertia,
            initial_condition=ic, id=id, control=SM.ControlModel(Tuple(declared), fill(1.0, length(declared))),
            guidance=SM.GuidanceModel(Tuple(declared_guidance), fill(1.0, length(declared_guidance))))
    end
    density = drag ? SM.ExponentialAtmosphereModel(1.0e-9, 1.4e5, 7.0e3; temperature_k=800.0,
        valid_min_altitude_m=0.0, valid_max_altitude_m=2.0e6) : SM.NoAtmosphereModel()
    effectors = drag ? (SM.InverseSquaredJ2GravityModel(), SM.AerodynamicCoefficientConstant()) :
        (SM.InverseSquaredJ2GravityModel(),)
    env = SM.EnvironmentModel(planet=PLANET, EI=200.0, density_model=density, ephemerides_model=SM.SimpleEphemeridesModel(),
        thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=PLANET), topography=false, wind=false)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(; results=false, verbose=false, generate_plots=false, normalize=false,
            settings...),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=mission_time, orientation_sim=false, num_steps_to_save=100),
        environment_model=env,
        dynamics_model=SM.DynamicsModel(SM.SpacecraftModel[sats...], effectors),
        guidance_model=SM.GuidanceModel(Tuple(guidance), fill(1.0, length(guidance))),
        control_model=SM.ControlModel(control_effectors=Tuple(control), control_rates=fill(1.0, length(control))),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0))
end

# The run's GNC scopes, from a solution's parameters.
run_scopes(solution) = solution.prob.p.shared_buffers.gnc_scopes[]

# A controller with no spacecraft until its initializer picks the one with id
# `id`, as the RPO pair finds its chaser.
mutable struct FindsCtl <: SM.AbstractControlEffectorModel
    id::Int
    spacecraft_index::Int
    calls::Vector{Int}
end
SpaceAGORA.spacecraft_scope(m::FindsCtl) = m.spacecraft_index == 0 ? AllSpacecraft() : BoundSpacecraft((m.spacecraft_index,))
SpaceAGORA.update_control!(m::FindsCtl, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.control_wrench(m::FindsCtl, x::StateSample, p, t::Float64) = (m.calls[x.spacecraft_index] += 1; (Z3, Z3))
function SpaceAGORA.initialize_control!(m::FindsCtl, u, p, t)
    m.spacecraft_index = findfirst(sc -> sc.id == m.id, p.args.dynamics_model.spacecraft)
    fill!(m.calls, 0)
    return nothing
end

# A controller that sets `next` as its scope in its tick at `at` s or later.
mutable struct TurnsCtl <: SM.AbstractControlEffectorModel
    scope::Any
    next::Any
    at::Float64
end
SpaceAGORA.spacecraft_scope(m::TurnsCtl) = m.scope
SpaceAGORA.update_control!(m::TurnsCtl, u, p, t::Float64, i::Int) = (t >= m.at && (m.scope = m.next); nothing)
SpaceAGORA.control_wrench(m::TurnsCtl, x::StateSample, p, t::Float64) = (Z3, Z3)

# A controller whose selection's positions its tick changes in place at `at` s.
mutable struct MutatesCtl <: SM.AbstractControlEffectorModel
    scope::SpacecraftSelection
    at::Float64
end
SpaceAGORA.spacecraft_scope(m::MutatesCtl) = m.scope
SpaceAGORA.update_control!(m::MutatesCtl, u, p, t::Float64, i::Int) =
    (t >= m.at && (m.scope.positions[1] = 3); nothing)
SpaceAGORA.control_wrench(m::MutatesCtl, x::StateSample, p, t::Float64) = (Z3, Z3)

# An event-driven controller whose own callback changes its scope at `at` s.
mutable struct EventCtl <: SM.AbstractControlEffectorModel
    scope::Any
    at::Float64
end
SpaceAGORA.tick_schedule(::EventCtl) = :event_driven
SpaceAGORA.spacecraft_scope(m::EventCtl) = m.scope
SpaceAGORA.update_control!(m::EventCtl, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.control_wrench(m::EventCtl, x::StateSample, p, t::Float64) = (Z3, Z3)
SpaceAGORA.gnc_callbacks(m::EventCtl, num_sats, args, period_s) =
    (PresetTimeCallback([m.at], integrator -> (m.scope = BoundSpacecraft((2,))); save_positions=(false, false)),)

# Calls of a selection's `where`; a global, which the run's copy of the
# configuration does not copy.
const ASKED = Ref(0)

# A controller that counts its ticks per spacecraft, with the scope given.
mutable struct TickCtl <: SM.AbstractControlEffectorModel
    scope::Any
    ticks::Vector{Int}
    calls::Vector{Int}
end
TickCtl(n, scope) = TickCtl(scope, zeros(Int, n), zeros(Int, n))
SpaceAGORA.spacecraft_scope(m::TickCtl) = m.scope
SpaceAGORA.update_control!(m::TickCtl, u, p, t::Float64, i::Int) = (m.ticks[i] += 1; nothing)
SpaceAGORA.control_wrench(m::TickCtl, x::StateSample, p, t::Float64) = (m.calls[x.spacecraft_index] += 1; (Z3, Z3))
SpaceAGORA.initialize_control!(m::TickCtl, u, p, t) = (fill!(m.ticks, 0); fill!(m.calls, 0); nothing)

# A guidance model that records the atmosphere crossings it is told about.
mutable struct CrossingGd <: AbstractGuidanceModel
    spacecraft_index::Int
    crossings::Vector{Tuple{Int, Bool}}
    copies::Vector{Any}
end
CrossingGd() = CrossingGd(0, Tuple{Int, Bool}[], Any[])
SpaceAGORA.spacecraft_scope(m::CrossingGd) = m.spacecraft_index == 0 ? AllSpacecraft() : BoundSpacecraft((m.spacecraft_index,))
SpaceAGORA.update_guidance!(m::CrossingGd, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.SimulationLifecycle.requires_atmosphere_events(::CrossingGd) = true
SpaceAGORA.SimulationLifecycle.atmosphere_transition!(m::CrossingGd, u, p, t, i, inside) = (push!(m.crossings, (i, inside)); nothing)
function SpaceAGORA.bind_spacecraft(m::CrossingGd, k::Int)
    bound = CrossingGd(k, Tuple{Int, Bool}[], Any[])
    push!(m.copies, bound)
    return bound
end

# Mounted actuation: spacecraft with a reaction wheel and a battery, so the
# dynamics, the power snapshot and the result columns read the levels (with a
# thruster and a fault too, the thruster-level column and the fault reader). Each
# bound copy holds a wheel command for its own spacecraft, the last of its
# `actuators` levels (thrusters mount first), and counts every actuator_levels
# call per spacecraft.
mutable struct WheelCtl <: SM.AbstractControlEffectorModel
    spacecraft_index::Int
    levels::Vector{Float64}
    calls::Vector{Int}
    copies::Vector{Any}
end
WheelCtl(n; index=0, actuators=1) = WheelCtl(index, zeros(actuators), zeros(Int, n), Any[])
SpaceAGORA.actuation_mode(::WheelCtl) = :mounted
SpaceAGORA.spacecraft_scope(m::WheelCtl) = m.spacecraft_index == 0 ? AllSpacecraft() : BoundSpacecraft((m.spacecraft_index,))
SpaceAGORA.update_control!(m::WheelCtl, u, p, t::Float64, i::Int) = (m.levels[end] = 1.0e-3 * i; nothing)
function SpaceAGORA.actuator_levels(m::WheelCtl, i::Int)
    m.calls[i] += 1
    return i == m.spacecraft_index ? m.levels : nothing
end
SpaceAGORA.initialize_control!(m::WheelCtl, u, p, t) = (fill!(m.levels, 0.0); fill!(m.calls, 0); nothing)
function SpaceAGORA.bind_spacecraft(m::WheelCtl, k::Int)
    bound = WheelCtl(length(m.calls); index=k, actuators=length(m.levels))
    push!(m.copies, bound)
    return bound
end

function wheel_spacecraft(id, template; thruster=false, fault=false)
    r = PLANET.Rp_e + 500e3 + 10e3 * id
    ic = SM.CartesianInitialCondition(SVector(r, 0.0, 0.0), SVector(0.0, sqrt(PLANET.μ / r), 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    bus = SM.Link(root=true, m=10.0, dims=MVector(0.1, 0.2, 0.3))
    sc = SM.SpacecraftModel(links=[bus], root=bus, inertia_tensor=bus.inertia, initial_condition=ic, id=id,
        control=SM.ControlModel((template,), [1.0]))
    wheel = add_actuator!(sc, bus, ReactionWheel(spin_axis=SVector(0.0, 0.0, 1.0)))
    thruster && add_actuator!(sc, bus, Thruster(max_thrust=1.0, direction=MVector(1.0, 0.0, 0.0), Isp=300.0))
    fault && add_fault!(sc, ActuatorFault(wheel; kind=:degraded, value=0.5, start_s=1.0, stop_s=3.0))
    add_power!(sc, bus, Battery(capacity_wh=5.0, initial_soc=0.5))
    add_power!(sc, bus, PowerLoad(power_w=1.0, name="avionics"))
    return sc
end

function wheel_args(n, template; mission_time=5.0, thruster=false, fault=false, settings=(;))
    args = scope_args(n; spacecraft=[wheel_spacecraft(id, template; thruster, fault) for id in 1:n],
        mission_time=mission_time, settings=settings)
    mission = args.mission_configuration
    return SM.SimConfig._with_configuration(args; mission_configuration=SM.MissionConfiguration(
        mission_type=SM.MissionTime, keplerian=true, number_of_orbits=1, mission_time=mission_time,
        orientation_sim=true, num_steps_to_save=mission.num_steps_to_save))
end

# An event-driven controller that counts, per spacecraft, the reschedules an
# atmosphere exit asks of it.
mutable struct ExitCtl <: SM.AbstractControlEffectorModel
    spacecraft_index::Int
    exits::Vector{Int}
    copies::Vector{Any}
end
ExitCtl(n; index=0) = ExitCtl(index, zeros(Int, n), Any[])
SpaceAGORA.tick_schedule(::ExitCtl) = :event_driven
SpaceAGORA.spacecraft_scope(m::ExitCtl) = m.spacecraft_index == 0 ? AllSpacecraft() : BoundSpacecraft((m.spacecraft_index,))
SpaceAGORA.update_control!(m::ExitCtl, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.control_wrench(m::ExitCtl, x::StateSample, p, t::Float64) = (Z3, Z3)
SpaceAGORA.gnc_callbacks(m::ExitCtl, num_sats, args, period_s) = ()
CB._reschedule_event_driven!(integrator, m::ExitCtl, i::Int) = (m.exits[i] += 1; nothing)
function SpaceAGORA.bind_spacecraft(m::ExitCtl, k::Int)
    bound = ExitCtl(length(m.exits); index=k)
    push!(m.copies, bound)
    return bound
end

# A controller whose scope method builds its selection, with a new closure, at
# every call.
struct RebuildsCtl <: SM.AbstractControlEffectorModel
    ids::Vector{Int}
end
SpaceAGORA.spacecraft_scope(m::RebuildsCtl) = (ids = copy(m.ids); SpacecraftSelection(where=(i, sc) -> sc.id in ids))
SpaceAGORA.update_control!(m::RebuildsCtl, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.control_wrench(m::RebuildsCtl, x::StateSample, p, t::Float64) = (Z3, Z3)

# A controller bound to spacecraft 1 until its initializer moves it to 2, and that
# pushes every spacecraft.
mutable struct MovesCtl <: SM.AbstractControlEffectorModel
    spacecraft_index::Int
end
SpaceAGORA.spacecraft_scope(m::MovesCtl) = BoundSpacecraft((m.spacecraft_index,))
SpaceAGORA.update_control!(m::MovesCtl, u, p, t::Float64, i::Int) = nothing
SpaceAGORA.control_wrench(m::MovesCtl, x::StateSample, p, t::Float64) = (SVector(1.0e-6, 0.0, 0.0), Z3)
SpaceAGORA.initialize_control!(m::MovesCtl, u, p, t) = (m.spacecraft_index = 2; nothing)

@testset "GNC hooks follow the scope" begin
    @testset "the output hooks are called only for the spacecraft in the scope" begin
        n = 4
        template, shared = PushCtl(n), PushCtl(n)
        coupled = PushCtl(n; scope=CoupledSpacecraft(2, (3,)))
        serial() do
            run_simulation(scope_args(n; declared=(template,), control=(shared, coupled)); isolate_state=false)
        end
        evaluations = shared.calls[1]
        @test evaluations > 0 && shared.calls == fill(evaluations, n)        # every spacecraft, once per evaluation
        @test length(template.copies) == n
        for bound in template.copies
            k = bound.spacecraft_index
            @test bound.calls == [i == k ? evaluations : 0 for i in 1:n]
        end
        @test coupled.calls == [i == 2 ? evaluations : 0 for i in 1:n]   # the primary only
        # Negative control: without the member lists every copy is called everywhere.
        template, shared = PushCtl(n), PushCtl(n)
        serial() do
            without_lists(() -> run_simulation(scope_args(n; declared=(template,), control=(shared,)); isolate_state=false))
        end
        @test all(bound -> bound.calls == fill(shared.calls[1], n), template.copies)
    end

    @testset "and give the bits of a run that calls every model everywhere" begin
        for (label, args) in (("bound copies, a coupled and a shared model",
                () -> scope_args(4; declared=(PushCtl(4),), control=(PushCtl(4), PushCtl(4; scope=CoupledSpacecraft(2, (3,)))))),
                ("mounted wheels with batteries", () -> wheel_args(2, WheelCtl(2))))
            with_lists = run_simulation(args(); return_solution=true)
            everywhere = without_lists(() -> run_simulation(args(); return_solution=true))
            @test same_solution(with_lists, everywhere)
            @test length(with_lists.t) > 2
        end
    end

    @testset "each spacecraft's models are called in configuration order" begin
        # The control tuple is (1, 2, 11, 12): spacecraft 1 has models 1 and 11 and
        # spacecraft 2 has 1, 2 and 12, in that order at every evaluation.
        empty!(CALLS)
        serial(() -> run_simulation(scope_args(2; declared=(OrderCtl(10, AllSpacecraft()),),
            control=(OrderCtl(1, AllSpacecraft()), OrderCtl(2, BoundSpacecraft((2,)))))))
        for (k, expected) in ((1, [1, 11]), (2, [1, 2, 12]))
            ids = [id for (i, id) in CALLS if i == k]
            @test !isempty(ids) && length(ids) % length(expected) == 0
            @test ids == repeat(expected, length(ids) ÷ length(expected))
        end
    end

    @testset "hook calls per evaluation follow the scopes, not the number of spacecraft" begin
        per_evaluation = map((3, 9)) do n
            template, shared = PushCtl(n), PushCtl(n)
            serial(() -> run_simulation(scope_args(n; declared=(template,), control=(shared,)); isolate_state=false))
            sum(sum(bound.calls) for bound in template.copies) / shared.calls[1]
        end
        @test per_evaluation == (3.0, 9.0)                               # one per bound copy: n, not n^2
    end

    @testset "the mounted levels are read only for the spacecraft in the scope" begin
        n = 2
        template = WheelCtl(n)
        serial(() -> run_simulation(wheel_args(n, template); isolate_state=false))
        for bound in template.copies
            k = bound.spacecraft_index
            @test bound.calls[k] > 0 && bound.calls[3 - k] == 0
        end
        # Negative control: without the member lists the dynamics, the power snapshot
        # and the result columns read every copy's levels for both spacecraft.
        template = WheelCtl(n)
        serial(() -> without_lists(() -> run_simulation(wheel_args(n, template); isolate_state=false)))
        @test all(bound -> all(>(0), bound.calls), template.copies)
        # With a thruster and a wheel fault, the thruster-level column (saved with the
        # visualization scene's fields) and the fault reader read the levels too, and
        # only for the spacecraft in the scope.
        template = WheelCtl(n; actuators=2)
        table = serial(() -> run_simulation(wheel_args(n, template; thruster=true, fault=true,
            settings=(save_visualization_scene=true,)); isolate_state=false, return_results=true)).table
        @test "sc1_thruster_level_1" in names(table) && "sc2_thruster_level_1" in names(table)
        for bound in template.copies
            k = bound.spacecraft_index
            @test bound.calls[k] > 0 && bound.calls[3 - k] == 0
        end
    end

    @testset "atmosphere_transition! is called for the models whose scope names the spacecraft" begin
        template, shared = CrossingGd(), CrossingGd()
        run_simulation(scope_args(2; declared_guidance=(template,), guidance=(shared,), drag=true, mission_time=900.0);
            isolate_state=false)
        @test !isempty(shared.crossings) && Set(first.(shared.crossings)) == Set([1, 2])
        for bound in template.copies
            @test !isempty(bound.crossings) && all(c -> first(c) == bound.spacecraft_index, bound.crossings)
        end
    end

    @testset "an atmosphere exit reschedules the event-driven models whose scope names the spacecraft" begin
        template, shared = ExitCtl(2), ExitCtl(2)
        run_simulation(scope_args(2; declared=(template,), control=(shared,), drag=true, mission_time=900.0);
            isolate_state=false)
        @test all(>(0), shared.exits)
        @test length(template.copies) == 2
        for bound in template.copies
            k = bound.spacecraft_index
            @test bound.exits[k] > 0 && bound.exits[3 - k] == 0
        end
    end

    @testset "the run takes each scope when it starts" begin
        n = 3
        # Every model acts on every spacecraft: no member lists, today's path.
        all_only = run_simulation(scope_args(n; control=(PushCtl(n),)); return_solution=true)
        scopes = run_scopes(all_only)
        @test scopes.control_members === nothing && scopes.guidance_members === nothing
        @test scopes.control[1].positions == 1:n && scopes.control[1].owner == "control model 1 (PushCtl)"
        # Bound copies: each spacecraft lists the shared model, then its own copy.
        bound = run_simulation(scope_args(n; declared=(PushCtl(n),), control=(PushCtl(n),)); return_solution=true)
        @test run_scopes(bound).control_members == [[1, 2], [1, 3], [1, 4]]
        # A scope an initializer sets, as the RPO pair's, is the one the run takes.
        finds = FindsCtl(3, 0, zeros(Int, n))
        solution = serial(() -> run_simulation(scope_args(n; control=(finds,)); return_solution=true, isolate_state=false))
        @test run_scopes(solution).control[1].positions == [3]
        @test finds.calls[1] == 0 && finds.calls[2] == 0 && finds.calls[3] > 0
        # The checkpoint's record of them is the gnc_scopes entry #476 writes.
        args = SE._flatten_spacecraft_gnc(scope_args(n; declared=(PushCtl(n),),
            control=(PushCtl(n), PushCtl(n; scope=CoupledSpacecraft(2, (3,))), PushCtl(n; scope=BoundSpacecraft((3, 1))))))
        @test SE._gnc_scope_record(args, SE._gnc_run_scopes(args)) == SE._gnc_scope_baseline(args)
    end

    @testset "a scope that changes during the run stops it" begin
        n = 3
        message = "spacecraft_scope changed during the run"
        # A periodic model that changes its own scope in its tick.
        @test_throws message run_simulation(scope_args(n;
            control=(TurnsCtl(BoundSpacecraft((1,)), BoundSpacecraft((2,)), 3.0),)))
        @test_throws "control model 1 (TurnsCtl)" run_simulation(scope_args(n;
            control=(TurnsCtl(BoundSpacecraft((1,)), BoundSpacecraft((2,)), 3.0),)))
        # The scopes as users write them, without their type parameters.
        @test_throws "it was BoundSpacecraft((1,)) when the run started and is BoundSpacecraft((2,)) now" run_simulation(
            scope_args(n; control=(TurnsCtl(BoundSpacecraft((1,)), BoundSpacecraft((2,)), 3.0),)))
        @test CB._scope_text(CoupledSpacecraft(1, (2, 3))) == "CoupledSpacecraft(1, (2, 3))"
        @test CB._scope_text(SpacecraftSelection(ids=[2], where=(i, sc) -> true)) ==
            "SpacecraftSelection(ids=[2], where=<function>)"
        # A model that acts on every spacecraft and then on one: the run took its scope too.
        @test_throws message run_simulation(scope_args(n; control=(TurnsCtl(AllSpacecraft(), BoundSpacecraft((1,)), 3.0),)))
        # An event-driven model whose own callback changes its scope.
        @test_throws message run_simulation(scope_args(n; control=(EventCtl(BoundSpacecraft((1,)), 2.5),)))
        # A selection changed in place: the run compares its own copy.
        @test_throws message run_simulation(scope_args(n;
            control=(MutatesCtl(SpacecraftSelection(positions=[1, 2]), 3.0),)))
        # A selection built anew at each call, with a new where closure: the message
        # says to return one object.
        @test_throws "another where function" run_simulation(scope_args(n; control=(RebuildsCtl([2, 3]),)))
        # The same spacecraft in another order is the same scope.
        @test run_simulation(scope_args(n; control=(TurnsCtl(BoundSpacecraft((2, 1)), BoundSpacecraft((1, 2)), 3.0),))) ===
            nothing
        # An event-driven model that keeps its scope runs to the end.
        @test run_simulation(scope_args(n; control=(EventCtl(BoundSpacecraft((2,)), 2.5),))) === nothing
    end

    @testset "a SpacecraftSelection scopes a GNC model" begin
        factory = (plane, slot, index) -> begin
            root = SM.Link(root=true, m=100.0, ref_area=2.0)
            SM.SpacecraftModel(joints=SM.Joint[], links=[root], root=root, prop_mass=10.0, inertia_tensor=root.inertia,
                initial_condition=SM.InitialCondition(a=PLANET.Rp_e + 5.5e5, e=0.0, i=53.0, Ω=0.0, ω=0.0, ν=0.0),
                id=100 + index)
        end
        walker = walker_constellation(factory; total=6, planes=3, phasing=1, altitude=5.5e5, planet=PLANET, i=53)
        plane_two = [i for (i, sc) in enumerate(walker) if sc.constellation_slot.plane == 2]
        @test length(plane_two) == 2
        model = TickCtl(6, SpacecraftSelection(planes=2))
        serial(() -> run_simulation(scope_args(6; spacecraft=walker, control=(model,)); isolate_state=false))
        @test [model.ticks[i] > 0 for i in 1:6] == [i in plane_two for i in 1:6]
        @test [model.calls[i] > 0 for i in 1:6] == [i in plane_two for i in 1:6]
        # The setup check accepts a selection; the run resolves it once, `where`
        # called once per spacecraft, in a checkpointed run too.
        @test CB.check_spacecraft_scopes(scope_args(6; spacecraft=walker, control=(model,))) === nothing
        ASKED[] = 0
        choose = SpacecraftSelection(where=(i, sc) -> (ASKED[] += 1; i in (2, 5)))
        run_simulation(scope_args(6; spacecraft=walker, control=(TickCtl(6, choose),)))
        @test ASKED[] == 6
        ASKED[] = 0
        mktempdir() do dir
            settings = (checkpoint_enabled=true, checkpoint_interval_s=2.0, checkpoint_directory=dir)
            run_simulation(scope_args(6; spacecraft=walker, control=(TickCtl(6, choose),), settings=settings))
        end
        @test ASKED[] == 6
        # A selection that cannot be resolved is refused when the run starts.
        @test_throws "selects spacecraft id 99" run_simulation(scope_args(6; spacecraft=walker,
            control=(TickCtl(6, SpacecraftSelection(ids=[99])),)))
        # Outside a run nothing can say which spacecraft it names.
        @test_throws "resolved when the run starts" CB._tick_positions(SpacecraftSelection(planes=2), 6)
    end

    @testset "the conformance kit checks the outputs outside the scope" begin
        C = SpaceAGORA.Conformance
        status(report) = only(check.status for check in report.checks if check.id === :outputs_in_scope)
        bound = PushCtl(2; index=1)
        @test status(C.check_control_model(bound; fixture=C.conformance_fixture(control=(bound,), spacecraft_count=2))) ===
            :pass
        shared = PushCtl(2)
        @test status(C.check_control_model(shared; fixture=C.conformance_fixture(control=(shared,), spacecraft_count=2))) ===
            :skip
        # A one-spacecraft fixture is checked with a copy of its spacecraft as the second.
        bound = PushCtl(2; index=1)
        @test status(C.check_control_model(bound)) === :pass
        shared = PushCtl(2)
        @test status(C.check_control_model(shared)) === :skip
        # The spacecraft checked is outside the scope the run took, after the
        # initializer moved the model to spacecraft 2.
        moves = MovesCtl(1)
        report = C.check_control_model(moves; fixture=C.conformance_fixture(control=(moves,), spacecraft_count=2))
        check = only(c for c in report.checks if c.id === :outputs_in_scope)
        @test check.status === :advisory_fail && occursin("for spacecraft 1", check.detail)
        # A selection the fixture cannot resolve: a report, not an exception.
        planes = TickCtl(1, SpacecraftSelection(planes=2))
        check = only(c for c in C.check_control_model(planes).checks if c.id === :outputs_in_scope)
        @test check.status === :skip && occursin("constellation", check.detail)
    end
end

end
