module InterlinkPolicyTests

# Scheduling policies as types: the built-in policies' constructors, a policy of
# your own through `candidate_score`, the setup check, the refusal of the
# symbol forms removed in 0.6.0, and the checkpoint entry of a policy.

using Test
using LinearAlgebra
using Serialization
using StaticArrays
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const SP = SM.SchedulingPolicies
const PLANET = make_no_gram_planet(:earth)
const RADIUS = PLANET.Rp_e + 700e3
const SPEED = sqrt(PLANET.μ / RADIUS)

# A policy of your own: prefer short links, `1 / separation` (1/m). It records
# the times it is called at, so the tests can count its calls.
struct ShortLinkPolicy <: AbstractSchedulingPolicy
    calls::Vector{Float64}
end
ShortLinkPolicy() = ShortLinkPolicy(Float64[])
function SpaceAGORA.candidate_score(policy::ShortLinkPolicy, candidate::LinkCandidate, ctx::SchedulingContext, t::Float64)
    push!(policy.calls, t)
    @assert isequal(t, ctx.link_context.time_s)
    return 1 / norm(candidate.a.pos_ii - candidate.b.pos_ii)
end

struct NoScorePolicy <: AbstractSchedulingPolicy end

struct ConstantPolicy <: AbstractSchedulingPolicy
    value::Float64
end
SpaceAGORA.candidate_score(policy::ConstantPolicy, candidate, ctx, t) = policy.value

# Same-named policies in two modules, and a policy with nested and vector fields.
module Left
import SpaceAGORA
struct Same <: SpaceAGORA.AbstractSchedulingPolicy end
SpaceAGORA.candidate_score(::Same, candidate, ctx, t) = 1.0
end
module Right
import SpaceAGORA
struct Same <: SpaceAGORA.AbstractSchedulingPolicy end
SpaceAGORA.candidate_score(::Same, candidate, ctx, t) = 1.0
end
struct NestedPolicy <: AbstractSchedulingPolicy
    inner::AbstractSchedulingPolicy
    weights::Vector{Float64}
    label::Symbol
end
SpaceAGORA.candidate_score(policy::NestedPolicy, candidate, ctx, t) = 1.0

# Policies with a function field: a closure from one factory expression, which
# captures its gain or a vector; a closure that captures a reassigned variable
# (a `Core.Box`); and a list that can contain itself.
struct CallablePolicy{F} <: AbstractSchedulingPolicy
    f::F
end
# Counts every call in the process, through the run's own copy of the policy.
const CALLABLE_SCORES = Ref(0)
SpaceAGORA.candidate_score(policy::CallablePolicy, candidate, ctx, t) = (CALLABLE_SCORES[] += 1; policy.f(1.0))
gain_policy(gain) = CallablePolicy(x -> gain * x)
vector_policy(v) = CallablePolicy(x -> v[1] * x)
function boxed_policy()
    gain = 1.0
    f = x -> gain * x
    gain = 2.0
    return CallablePolicy(f)
end
struct ListPolicy <: AbstractSchedulingPolicy
    items::Vector{Any}
end
# Values recorded field by field: a callable struct that is not a function, a
# closure that captures a `Ref`, and a table of arrays held by a policy.
struct Scale
    v::Vector{Float64}
end
(scale::Scale)(x) = scale.v[1] * x
ref_policy(x) = (held = Ref(x); CallablePolicy(y -> held[] * y))
struct Table
    xs::Vector{Float64}
    ys::Matrix{Float64}
end
struct HeldPolicy <: AbstractSchedulingPolicy
    held::Any
end
mutable struct Node
    next::Any
end
struct Handle
    pointer::Ptr{Cvoid}
end
SpaceAGORA.candidate_score(policy::HeldPolicy, candidate, ctx, t) = 1.0
SpaceAGORA.candidate_score(policy::ListPolicy, candidate, ctx, t) = 1.0

# A 100 kg bus `arc_m` along the orbit, with a laser transmitter (`sender`) or a
# power receiver, and a battery of `capacity_wh` at `soc`.
function beam_sat(index, arc_m; sender, capacity_wh=10.0, soc=0.5, battery=true, momentum=false)
    bus = Link(root=true, m=100.0)
    θ = arc_m / RADIUS
    ic = CartesianInitialCondition(SVector(RADIUS * cos(θ), RADIUS * sin(θ), 0.0),
        SVector(-SPEED * sin(θ), SPEED * cos(θ), 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, initial_condition=ic, id=index, n_terminal=0)
    if sender || momentum
        add_actuator!(sc, bus, LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=100.0, efficiency=0.5))
    else
        add_actuator!(sc, bus, PowerReceiver(normal=[0.0, -1.0, 0.0], radius_m=0.25, efficiency=0.4))
    end
    battery && add_power!(sc, bus, Battery(capacity_wh=capacity_wh, initial_soc=soc))
    return sc
end
const BEAM = PowerBeamLink(power_w=80.0, range_m=50e3)
# A directed link kind of your own whose receiver gets `power_w`, zero or negative.
struct FixedPowerLink <: AbstractInterlink
    power_w::Float64
end
SpaceAGORA.link_directed(::FixedPowerLink) = true
SpaceAGORA.link_received_power_w(link::FixedPowerLink, receiver, a, b, ctx) = link.power_w

# The state a scheduler sees: each spacecraft's position, velocity, mass and, with
# a battery, its stored energy.
function beam_state(spacecraft)
    return (sc=[begin
        ic = sc.initial_condition
        b = SM.PowerSystems.battery(sc)
        base = (pos=collect(ic.pos), vel=collect(ic.vel), mass=100.0)
        b === nothing ? base : merge(base, (battery_wh=b.capacity_wh * b.initial_soc,))
    end for sc in spacecraft],)
end

# One sender and two receivers, a beam to each from the sender's one terminal.
function energy_fixture(; sender=(;), first=(;), second=(;))
    spacecraft = [beam_sat(1, 0.0; sender=true, soc=0.9, sender...), beam_sat(2, 5000.0; sender=false, soc=0.2, first...),
                  beam_sat(3, -5000.0; sender=false, soc=0.6, second...)]
    model = InterlinkModel()
    register_link!(model, spacecraft, (1, 1), (2, 1), BEAM)
    register_link!(model, spacecraft, (1, 1), (3, 1), BEAM)
    state = beam_state(spacecraft)
    update_availability!(model, spacecraft, state)
    return spacecraft, model, state
end

# The power a beam on `edge` delivers, from the hook, at `state`.
function delivered_w(model, spacecraft, state, edge)
    (tx, txt), (rx, rxt) = map(v -> model.graph.vertices[v], model.graph.edges[edge])
    a = SM.InterLinkModels.link_endpoint(state.sc[tx], spacecraft[tx], tx)
    b = SM.InterLinkModels.link_endpoint(state.sc[rx], spacecraft[rx], rx)
    receiver = SM.InterLinkModels._terminal_device(spacecraft[rx], rxt)
    return link_received_power_w(model.graph.edge_properties[edge], receiver, a, b, LinkContext())
end
scores_of(model, policy, state, spacecraft=nothing) = (scheduler = SchedulingPolicyModel(policy);
    score_candidates!(model, scheduler, state, PLANET.μ; spacecraft); copy(SP.matching_workspace!(model, scheduler).score))
selected_by(model, policy, state, spacecraft=nothing) = (scheduler = SchedulingPolicyModel(policy);
    score_candidates!(model, scheduler, state, PLANET.μ; spacecraft); select_interlinks!(model, scheduler))

function sat(index, arc_m)
    bus = Link(root=true, m=100.0)
    θ = arc_m / RADIUS
    ic = CartesianInitialCondition(SVector(RADIUS * cos(θ), RADIUS * sin(θ), 0.0),
        SVector(-SPEED * sin(θ), SPEED * cos(θ), 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    return SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, initial_condition=ic, id=index,
        n_terminal=1)
end

function run_args(policy; duration=5.0, settings=SimulationSettings(results=false, generate_plots=false))
    spacecraft = [sat(1, 0.0), sat(2, 2000.0), sat(3, 9000.0)]
    model = InterlinkModel()
    link = LaserMomentumLink(power_w=40.0, magnification=1.0, range_m=1e4)
    register_link!(model, spacecraft, (1, 1), (2, 1), link)
    register_link!(model, spacecraft, (1, 1), (3, 1), link)
    return SimulationConfiguration(
        simulation_settings=settings,
        mission_configuration=MissionConfiguration(mission_type=MissionTime, mission_time=duration,
            orientation_sim=false, data_rate=1.0),
        environment_model=make_no_gram_environment(planet=PLANET),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]), navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]), initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-8, dt_max_orbit=1.0),
        solver_config=SM.SolverConfig(solver_mode=:tsit5), interlink_model=model,
        scheduling_policy_model=policy)
end

# A state of two spacecraft 1 km apart and a third out of range, with one link to each.
function scoring_fixture()
    spacecraft = [(n_terminal=1,), (n_terminal=1,), (n_terminal=1,)]
    model = InterlinkModel()
    link = LaserMomentumLink(power_w=40.0, magnification=1.0, range_m=1e4)
    register_link!(model, spacecraft, (1, 1), (2, 1), link)
    register_link!(model, spacecraft, (1, 1), (3, 1), link)
    state = (sc=[(pos=[RADIUS, 0.0, 0.0], vel=[0.0, SPEED, 0.0], mass=100.0),
                 (pos=[RADIUS, 1000.0, 0.0], vel=[0.0, SPEED, 0.0], mass=100.0),
                 (pos=[RADIUS, 5e4, 0.0], vel=[0.0, SPEED, 0.0], mass=100.0)],)
    update_availability!(model, spacecraft, state)
    return spacecraft, model, state
end

# The removed symbol forms and the typed call each message names.
const REMOVED = [
    ((:gve_sma, (;)), "SemimajorAxisRatePolicy()"),
    ((:gve_sma, (target_index=2,)), "SemimajorAxisRatePolicy(targets=2)"),
    ((:gve_sma, (target_idx=2,)), "SemimajorAxisRatePolicy(targets=2)"),
    ((:gve_sma, (target_index=nothing,)), "SemimajorAxisRatePolicy(targets=:both)"),
    ((:gve_sma, (target_idx=nothing,)), "SemimajorAxisRatePolicy(targets=:both)"),
    ((:gve_sma, (target_indices=[31, 32],)), "SemimajorAxisRatePolicy(targets=[31, 32])"),
    ((:gve_sma, (target_idxs=[31, 32],)), "SemimajorAxisRatePolicy(targets=[31, 32])"),
    ((:gve_sma, (target_indices=31:40,)), "SemimajorAxisRatePolicy(targets=31:40)"),
    ((:gve_sma, (target_index=99, target_indices=[2, 1])), "SemimajorAxisRatePolicy(targets=[2, 1])"),
    ((:gve_sma_first_index, (;)), "SemimajorAxisRatePolicy(targets=:lower_index)"),
    ((:gve_sma_first_index, (target_index=3,)), "SemimajorAxisRatePolicy(targets=:lower_index)"),
    ((:gve_eccentricity, (;)), "EccentricityRatePolicy()"),
    ((:gve_eccentricity, (target_index=3,)), "EccentricityRatePolicy(targets=3)"),
    ((:gve_eccentricity, (target_idx=nothing,)), "EccentricityRatePolicy(targets=:both)"),
    ((:power_delivery, (;)), "PowerDeliveryPolicy()"),
    ((:power_delivery, (target_index=2,)), "PowerDeliveryPolicy()"),
    ((:data_throughput, (;)), "DataThroughputPolicy()"),
]

refusal(f) = try f(); nothing catch err err end

@testset "Interlink scheduling policies" begin
    @testset "built-in constructors" begin
        @test SemimajorAxisRatePolicy().targets == 1
        @test SemimajorAxisRatePolicy(targets=3).targets === 3
        @test SemimajorAxisRatePolicy(targets=[3, 1]).targets == [1, 3]
        @test SemimajorAxisRatePolicy(targets=31:33).targets == [31, 32, 33]
        @test SemimajorAxisRatePolicy(targets=:both).targets === :both
        @test SemimajorAxisRatePolicy(targets=:lower_index).targets === :lower_index
        for bad in (0, -1, Int[], [0], [1, 1], :all, true, [true], 1.5, "2", nothing)
            @test_throws ArgumentError SemimajorAxisRatePolicy(targets=bad)
        end
        @test EccentricityRatePolicy().targets === 1
        @test EccentricityRatePolicy(targets=:both).targets === :both
        for bad in (0, [1], :lower_index, true, nothing)
            @test_throws ArgumentError EccentricityRatePolicy(targets=bad)
        end
        @test PowerDeliveryPolicy() isa AbstractSchedulingPolicy && DataThroughputPolicy() isa AbstractSchedulingPolicy
    end

    @testset "the scheduler and the configuration" begin
        scheduler = SchedulingPolicyModel(PowerDeliveryPolicy())
        @test scheduler.objective isa PowerDeliveryPolicy && scheduler.matching === nothing
        @test fieldnames(SchedulingPolicyModel) == (:objective, :matching)
        args = run_args(DataThroughputPolicy())
        @test args.scheduling_policy_model isa SchedulingPolicyModel && args.scheduling_policy_model.objective isa DataThroughputPolicy
        args = run_args(scheduler)
        @test args.scheduling_policy_model === scheduler
        default = SimulationConfiguration(environment_model=args.environment_model, dynamics_model=args.dynamics_model,
            initial_time=args.initial_time).scheduling_policy_model
        @test default.objective isa SemimajorAxisRatePolicy && default.objective.targets === 1
    end

    @testset "a policy of your own" begin
        spacecraft, model, state = scoring_fixture()
        policy = ShortLinkPolicy()
        scheduler = SchedulingPolicyModel(policy)
        score_candidates!(model, scheduler, state, PLANET.μ; context=LinkContext(time_s=12.0, planet=PLANET))
        # Only the available candidate is scored; the other scores 0 without a call.
        @test policy.calls == [12.0]
        @test SP.matching_workspace!(model, scheduler).score == [1 / 1000.0, 0.0]
        @test select_interlinks!(model, scheduler) == [((1, 1), (2, 1))]
        score_candidates!(model, scheduler, state, PLANET.μ)
        @test isnan(last(policy.calls))
        # A nonfinite score is refused, naming the policy and the candidate.
        err = refusal(() -> score_candidates!(model, SchedulingPolicyModel(ConstantPolicy(NaN)), state, PLANET.μ))
        @test err isa ArgumentError && occursin("ConstantPolicy", err.msg) && occursin("((1, 1), (2, 1))", err.msg)
        # The positive-score rule holds for any policy.
        negative = SchedulingPolicyModel(ConstantPolicy(-1.0))
        score_candidates!(model, negative, state, PLANET.μ)
        @test isempty(select_interlinks!(model, negative))

        # In a run: selected at every accepted step, at the step's time.
        policy = ShortLinkPolicy()
        solution = run_simulation(run_args(policy); return_solution=true)
        history = solution.prob.p.args.interlink_model.history
        @test !isempty(history) && all(sample -> sample.active == [((1, 1), (2, 1))], history)
        @test solution.prob.p.args.scheduling_policy_model.objective isa ShortLinkPolicy
        # The run's policy is a copy (state isolation), whose calls cover the run's steps.
        calls = solution.prob.p.args.scheduling_policy_model.objective.calls
        @test !isempty(calls) && all(t -> 0.0 <= t <= 5.0, calls)
    end

    @testset "a policy without candidate_score is refused at setup" begin
        err = refusal(() -> run_simulation(run_args(NoScorePolicy())))
        @test err isa ArgumentError && occursin("NoScorePolicy has no candidate_score method", err.msg)
        @test_throws MethodError score_candidates!(scoring_fixture()[2], SchedulingPolicyModel(NoScorePolicy()),
            scoring_fixture()[3], PLANET.μ)
    end

    @testset "power and data policies need the spacecraft models" begin
        spacecraft, model, state = scoring_fixture()
        for policy in (PowerDeliveryPolicy(), DataThroughputPolicy())
            scheduler = SchedulingPolicyModel(policy)
            score_candidates!(model, scheduler, state, PLANET.μ)
            @test all(iszero, SP.matching_workspace!(model, scheduler).score)
        end
    end

    @testset "a candidate carries the hardware at its terminal numbers" begin
        # Two virtual terminals and two mounted ones on the first spacecraft, one
        # of each on the second: each end's hardware is the terminal its number names.
        function hardware_sat(index, arc_m, n_terminal; laser)
            bus = Link(root=true, m=100.0)
            θ = arc_m / RADIUS
            ic = CartesianInitialCondition(SVector(RADIUS * cos(θ), RADIUS * sin(θ), 0.0),
                SVector(-SPEED * sin(θ), SPEED * cos(θ), 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
            sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, initial_condition=ic, id=index,
                n_terminal=n_terminal)
            laser && add_actuator!(sc, bus, LaserTerminal(boresight=[0.0, 1.0, 0.0], max_power_w=100.0, efficiency=0.5))
            add_actuator!(sc, bus, PowerReceiver(normal=[0.0, -1.0, 0.0], radius_m=0.25, efficiency=0.4))
            return sc
        end
        spacecraft = [hardware_sat(1, 0.0, 2; laser=true), hardware_sat(2, 5000.0, 1; laser=false)]
        model = InterlinkModel()
        register_link!(model, spacecraft, (1, 1), (2, 1), LaserMomentumLink(power_w=40.0, magnification=1.0, range_m=1e4))
        register_link!(model, spacecraft, (1, 3), (2, 2), PowerBeamLink(power_w=80.0, range_m=1e4))
        state = (sc=[(pos=collect(sc.initial_condition.pos), vel=collect(sc.initial_condition.vel), mass=100.0)
            for sc in spacecraft],)
        update_availability!(model, spacecraft, state)
        ctx = SchedulingContext(PLANET.μ, LinkContext(), spacecraft, state, nothing)
        for edge in eachindex(model.graph.edges)
            candidate = SP._link_candidate(model, edge, ctx)
            for (endpoint, terminal) in zip(candidate.key, (candidate.terminal_a, candidate.terminal_b))
                @test terminal === SM.InterLinkModels._terminal_device(spacecraft[endpoint[1]], endpoint[2])
            end
        end
        @test SP._link_candidate(model, 2, ctx).terminal_a isa LaserTerminal
        @test SP._link_candidate(model, 1, ctx).terminal_a === nothing
    end

    @testset "the symbol forms are refused with their replacement" begin
        for ((symbol, keywords), replacement) in REMOVED
            err = refusal(() -> SchedulingPolicyModel(symbol; keywords...))
            @test err isa ArgumentError
            @test occursin("was removed in SpaceAGORA 0.6.0", err.msg)
            @test occursin("Use $(replacement), for example scheduling_policy_model=$(replacement).", err.msg)
            # The named call builds the policy the symbol form stood for.
            policy = Core.eval(@__MODULE__, Meta.parse(replacement))
            expected = SP._symbol_scheduling_policy(symbol, get(keywords, :target_index, get(keywords, :target_idx, 1)),
                get(keywords, :target_indices, get(keywords, :target_idxs, nothing)))
            @test typeof(policy) === typeof(expected) && isequal(getfield.(Ref(policy), fieldnames(typeof(policy))),
                getfield.(Ref(expected), fieldnames(typeof(expected))))
        end
        err = refusal(() -> SchedulingPolicyModel(:invalid))
        @test err isa ArgumentError && occursin(":invalid is not a scheduling policy", err.msg) &&
            all(name -> occursin(name, err.msg), ("SemimajorAxisRatePolicy", "EccentricityRatePolicy",
                "PowerDeliveryPolicy", "DataThroughputPolicy", "EnergyBalancingPolicy", "WeightedPolicy",
                "AbstractSchedulingPolicy"))
        err = refusal(() -> SchedulingPolicyModel(:gve_sma; target_index=2, mode=:fast))
        @test err isa ArgumentError && occursin("targets=2", err.msg) && occursin("The keyword mode has no replacement.", err.msg)
        # A combination the old constructor refused still names a typed call, whose own check refuses it.
        err = refusal(() -> SchedulingPolicyModel(:gve_eccentricity; target_indices=[1]))
        @test err isa ArgumentError && occursin("EccentricityRatePolicy(targets=[1])", err.msg)
        @test_throws ArgumentError EccentricityRatePolicy(targets=[1])
        err = refusal(() -> SchedulingPolicyModel(:gve_sma; target_indices=Int[]))
        @test err isa ArgumentError && occursin("SemimajorAxisRatePolicy(targets=Int64[])", err.msg)
        # The removed properties refuse, naming their replacement; other names are fields.
        scheduler = SchedulingPolicyModel(SemimajorAxisRatePolicy())
        for name in (:policy, :target_index, :target_indices, :target_idx, :target_idxs)
            err = refusal(() -> getproperty(scheduler, name))
            @test err isa ArgumentError && occursin("SchedulingPolicyModel.$(name) was removed in SpaceAGORA 0.6.0", err.msg) &&
                occursin(name === :policy ? "`.objective`" : "`.objective.targets`", err.msg)
            err = refusal(() -> setproperty!(scheduler, name, 2))
            @test err isa ArgumentError && occursin("removed in SpaceAGORA 0.6.0", err.msg)
        end
        @test scheduler.objective.targets == 1
        @test_throws FieldError scheduler.unknown
        scheduler.matching = nothing
        @test scheduler.matching === nothing
        # The same refusal under --depwarn=error: it is no deprecation.
        script = """
            using SpaceAGORA
            @assert Base.JLOptions().depwarn == 2
            err = try SchedulingPolicyModel(:gve_sma; target_idx=2); nothing catch e e end
            println("REFUSAL ", nameof(typeof(err)), " ", occursin("SemimajorAxisRatePolicy(targets=2)", err.msg))
            """
        output = read(pipeline(`$(Base.julia_cmd()) --startup-file=no --depwarn=error
            --project=$(Base.active_project()) -e $script`; stderr=devnull), String)
        @test occursin("REFUSAL ArgumentError true", output)
    end

    @testset "checkpoint entry" begin
        spacecraft = [(n_terminal=2,) for _ in 1:3]
        model = InterlinkModel()
        register_link!(model, spacecraft, (1, 1), (2, 1), LaserMomentumLink(power_w=12.0, magnification=3.0, range_m=100.0))
        args(policy) = (interlink_model=model, scheduling_policy_model=SchedulingPolicyModel(policy),
            dynamics_model=(spacecraft=spacecraft,))
        entry(policy) = (c = SE._interlink_checkpoint_configuration(args(policy)); (c.policy, c.target_idx, c.target_idxs))
        # The policies 0.5 chose by a symbol are recorded as 0.5 recorded them.
        @test entry(SemimajorAxisRatePolicy()) == (:gve_sma, 1, nothing)
        @test entry(SemimajorAxisRatePolicy(targets=2)) == (:gve_sma, 2, nothing)
        @test entry(SemimajorAxisRatePolicy(targets=:both)) == (:gve_sma, nothing, nothing)
        @test entry(SemimajorAxisRatePolicy(targets=[3, 1])) == (:gve_sma, 1, [1, 3])
        @test entry(SemimajorAxisRatePolicy(targets=:lower_index)) == (:gve_sma_first_index, 1, nothing)
        @test entry(EccentricityRatePolicy(targets=3)) == (:gve_eccentricity, 3, nothing)
        @test entry(EccentricityRatePolicy(targets=:both)) == (:gve_eccentricity, nothing, nothing)
        @test entry(PowerDeliveryPolicy()) == (:power_delivery, 1, nothing)
        @test entry(DataThroughputPolicy()) == (:data_throughput, 1, nothing)
        @test keys(SE._interlink_checkpoint_configuration(args(PowerDeliveryPolicy()))) ==
            (:vertices, :candidates, :forbidden_pairs, :policy, :target_idx, :target_idxs, :terminal_counts)
        # Another policy by its kind and every field: floats with their bits, a policy
        # field by its own identity, a vector element by element.
        nested = NestedPolicy(ConstantPolicy(2.0), [0.5, -0.0], :x)
        policy, target_idx, target_idxs = entry(nested)
        @test target_idx === nothing && target_idxs === nothing
        @test policy.kind.name === :NestedPolicy && policy.kind.module_path == fullname(@__MODULE__)
        @test policy.fields.inner.kind.name === :ConstantPolicy
        @test policy.fields.inner.fields.value.bits == reinterpret(UInt64, 2.0)
        @test [element.bits for element in policy.fields.weights.elements] == reinterpret.(UInt64, [0.5, -0.0])

        # Resume: the snapshot written with a policy, and with the targets 0.5 wrote.
        flags = Bool[true, true, true]
        state = (sc=[(pos=[Float64(i), 0.0, 0.0],) for i in 1:3],)
        checkpoint(policy; configuration=identity) = (u=state, runtime_state=(interlinks=merge(
            SE._interlink_checkpoint_runtime_state(args(policy), flags).interlinks,
            (configuration=configuration(SE._interlink_checkpoint_configuration(args(policy))),)),))
        resumes(policy, snapshot) = try
            SE._restore_interlink_checkpoint_state!(args(policy), snapshot, trues(3)); true
        catch err
            err isa ArgumentError || rethrow(); false
        end
        recorded(symbol, index, indices) = c -> merge(c, (policy=symbol, target_idx=index, target_idxs=indices))
        base = SemimajorAxisRatePolicy()
        for (snapshot_entry, policy) in (
                ((:gve_sma, 1, nothing), SemimajorAxisRatePolicy()),
                ((:gve_sma, nothing, [2, 3]), SemimajorAxisRatePolicy(targets=[2, 3])),
                ((:gve_sma, 99, [1, 2]), SemimajorAxisRatePolicy(targets=[1, 2])),
                ((:gve_sma_first_index, 3, nothing), SemimajorAxisRatePolicy(targets=:lower_index)),
                ((:gve_eccentricity, nothing, nothing), EccentricityRatePolicy(targets=:both)),
                ((:power_delivery, 7, nothing), PowerDeliveryPolicy()),
                ((:data_throughput, nothing, nothing), DataThroughputPolicy()))
            @test resumes(policy, checkpoint(base; configuration=recorded(snapshot_entry...)))
        end
        for (snapshot_entry, policy) in (
                ((:gve_sma, 1, nothing), SemimajorAxisRatePolicy(targets=2)),
                ((:gve_sma, 1, [1, 2]), SemimajorAxisRatePolicy(targets=[1, 3])),
                ((:gve_sma, 1, nothing), EccentricityRatePolicy()),
                ((:gve_sma_first_index, 1, nothing), SemimajorAxisRatePolicy(targets=:both)),
                ((:power_delivery, 1, nothing), DataThroughputPolicy()),
                ((:radio, 1, nothing), SemimajorAxisRatePolicy()),
                ((:gve_eccentricity, 1, [1, 2]), EccentricityRatePolicy()))
            @test !resumes(policy, checkpoint(base; configuration=recorded(snapshot_entry...)))
        end
        # A policy of your own resumes itself and refuses any other field, kind or module.
        @test resumes(nested, checkpoint(nested))
        @test resumes(ConstantPolicy(2.0), checkpoint(ConstantPolicy(2.0)))
        @test !resumes(ConstantPolicy(2.5), checkpoint(ConstantPolicy(2.0)))
        @test !resumes(ConstantPolicy(0.0), checkpoint(ConstantPolicy(-0.0)))
        @test !resumes(NestedPolicy(ConstantPolicy(2.0), [0.5, 0.0], :x), checkpoint(nested))
        @test !resumes(NestedPolicy(ConstantPolicy(3.0), [0.5, -0.0], :x), checkpoint(nested))
        @test !resumes(Right.Same(), checkpoint(Left.Same()))
        @test !resumes(SemimajorAxisRatePolicy(), checkpoint(ConstantPolicy(1.0)))
        @test !resumes(ConstantPolicy(1.0), checkpoint(SemimajorAxisRatePolicy()))

        # A function field: a closure is recorded with every value it captures, so
        # two closures of one factory expression differ by their captures.
        @test typeof(gain_policy(1.0)) === typeof(gain_policy(2.0))
        @test resumes(gain_policy(2.0), checkpoint(gain_policy(2.0)))
        @test !resumes(gain_policy(1.0), checkpoint(gain_policy(2.0)))
        @test resumes(vector_policy([2.0]), checkpoint(vector_policy([2.0])))
        @test !resumes(vector_policy([1.0]), checkpoint(vector_policy([2.0])))
        captured = entry(gain_policy(2.0))[1].fields.f
        @test captured.captures.gain.bits == reinterpret(UInt64, 2.0)
        # A named function is its type: another function is refused.
        @test resumes(CallablePolicy(sin), checkpoint(CallablePolicy(sin)))
        @test !resumes(CallablePolicy(cos), checkpoint(CallablePolicy(sin)))
        @test entry(CallablePolicy(sin))[1].fields.f == (callable=SE._interlink_kind_identity(typeof(sin)),)
        # A reassigned capture and a payload that contains itself cannot be recorded.
        err = refusal(() -> entry(boxed_policy()))
        @test err isa ArgumentError && occursin("reassigned after the closure was made", err.msg)
        cycle = Any[]
        push!(cycle, cycle)
        err = refusal(() -> entry(ListPolicy(cycle)))
        @test err isa ArgumentError && occursin("contain themselves", err.msg)
        shared = [1.0]
        @test entry(ListPolicy(Any[shared, shared]))[2] === nothing   # shared, not cyclic
        # Any other value is recorded by its type and fields, so an equal value
        # resumes, also once written and read back, and another value is refused.
        roundtrip(x) = (io = IOBuffer(); serialize(io, x); seekstart(io); deserialize(io))
        table(y) = Table([1.0, 2.0], [0.5 y; -0.0 4.0])
        for (same, other) in ((CallablePolicy(Scale([2.0])), CallablePolicy(Scale([1.0]))),
                (ref_policy(2.0), ref_policy(1.0)), (HeldPolicy(table(3.0)), HeldPolicy(table(-3.0))),
                (HeldPolicy((gain=[2.0], name="a")), HeldPolicy((gain=[2.0], name="b"))),
                (HeldPolicy(big(2)^70), HeldPolicy(big(2)^71)))
            @test resumes(same, checkpoint(same))
            @test !resumes(other, checkpoint(same))
            @test isequal(roundtrip(entry(same)), entry(same))
        end
        @test entry(HeldPolicy(table(3.0)))[1].fields.held.fields.ys.size == (2, 2)
        # A box inside a named tuple, a dictionary, a set, a stream and an
        # unassigned value cannot be recorded.
        err = refusal(() -> entry(HeldPolicy((gain=Core.Box(1.0),))))
        @test err isa ArgumentError && occursin("reassigned after the closure was made", err.msg)
        for held in (Dict(:gain => 2.0), Set([1.0]), IOBuffer())
            err = refusal(() -> entry(HeldPolicy(held)))
            @test err isa ArgumentError && occursin("holds a value of type `$(repr(typeof(held)))`, which a checkpoint " *
                "cannot record by value", err.msg)
        end
        # An unassigned element or field is named as such.
        err = refusal(() -> entry(HeldPolicy(Vector{Any}(undef, 2))))
        @test err isa ArgumentError && occursin("a value of type `Vector{Any}` with an unassigned element", err.msg)
        err = refusal(() -> entry(HeldPolicy(Ref{Any}())))
        @test err isa ArgumentError && occursin("whose field `x` is not assigned", err.msg)
        # A pointer's address differs from one process to the next: a pointer, a bits
        # value holding one, a `Regex` and a closure capturing a pointer are refused.
        for held in (C_NULL, Handle(C_NULL), r"abc", (let p = C_NULL; x -> p; end))
            err = refusal(() -> entry(HeldPolicy(held)))
            @test err isa ArgumentError && occursin("which holds a pointer", err.msg)
        end
        node = Node(nothing)
        node.next = node
        err = refusal(() -> entry(HeldPolicy(node)))
        @test err isa ArgumentError && occursin("contain themselves", err.msg)
        # A checkpointed run refuses such a policy at setup, before it starts.
        mktempdir() do directory
            settings = SimulationSettings(results=false, generate_plots=false, checkpoint_enabled=true,
                checkpoint_interval_s=1.0, checkpoint_directory=directory)
            CALLABLE_SCORES[] = 0
            err = refusal(() -> run_simulation(run_args(boxed_policy(); settings)))
            @test err isa ArgumentError && occursin("reassigned after the closure was made", err.msg)
            @test isempty(readdir(directory))
            @test CALLABLE_SCORES[] == 0   # refused before the scheduler's first call
            err = refusal(() -> run_simulation(run_args(HeldPolicy(Dict(:gain => 2.0)); settings)))
            @test err isa ArgumentError && occursin("cannot record by value", err.msg)
            @test isempty(readdir(directory))
        end
    end

    @testset "state of charge" begin
        spacecraft, model, state = energy_fixture(second=(battery=false,))
        scheduler = SchedulingPolicyModel(EnergyBalancingPolicy())
        ctx = SchedulingContext(PLANET.μ, LinkContext(), spacecraft, state, nothing)
        @test state_of_charge(ctx, 1) === state.sc[1].battery_wh / 10.0
        @test state_of_charge(ctx, 2) === state.sc[2].battery_wh / 10.0
        @test isnan(state_of_charge(ctx, 3))                     # no battery
        @test isnan(state_of_charge(SchedulingContext(PLANET.μ, LinkContext(), nothing, state, nothing), 1))
    end

    @testset "energy balancing" begin
        spacecraft, model, state = energy_fixture()
        s = [state.sc[i].battery_wh / 10.0 for i in 1:3]
        scores = scores_of(model, EnergyBalancingPolicy(), state, spacecraft)
        @test scores == [delivered_w(model, spacecraft, state, 1) * (s[1] - s[2]),
                         delivered_w(model, spacecraft, state, 2) * (s[1] - s[3])]
        @test all(>(0), scores)
        # The emptier receiver wins the sender's one terminal.
        @test selected_by(model, EnergyBalancingPolicy(), state, spacecraft) == [((1, 1), (2, 1))]
        # A full receiver, and one fuller than its sender, score negative; no battery at either end scores 0.
        spacecraft, model, state = energy_fixture(first=(soc=1.0,), second=(soc=0.95,))
        @test all(<(0), scores_of(model, EnergyBalancingPolicy(), state, spacecraft))
        @test isempty(selected_by(model, EnergyBalancingPolicy(), state, spacecraft))
        # Negative control: PowerDeliveryPolicy beams into the full battery.
        @test selected_by(model, PowerDeliveryPolicy(), state, spacecraft) != InterLinkKey[]
        # A battery within 1e-9 of its capacity, where the power system calls it full,
        # counts as full: a full sender and full receivers tie at 0, also when the
        # sender is a hair fuller; just outside the tolerance the beams score.
        for (tx, rx) in ((1.0, 1 - 1e-12), (1.0, 1 - 1e-10), (1 - 1e-12, 1 - 5e-10))
            spacecraft, model, state = energy_fixture(sender=(soc=tx,), first=(soc=rx,), second=(soc=rx,))
            @test all(iszero, scores_of(model, EnergyBalancingPolicy(), state, spacecraft))
            @test isempty(selected_by(model, EnergyBalancingPolicy(), state, spacecraft))
            @test !isempty(selected_by(model, PowerDeliveryPolicy(), state, spacecraft))
        end
        spacecraft, model, state = energy_fixture(sender=(soc=1.0,), first=(soc=1 - 1e-8,), second=(soc=1 - 1e-8,))
        @test all(>(0), scores_of(model, EnergyBalancingPolicy(), state, spacecraft))
        # Around the threshold, a beam into a receiver scores exactly when the power
        # system's boundary event does not flag the receiver full.
        for capacity in (10.0, 0.033450261609188466, 7300.0)
            spacecraft, model, state = energy_fixture(sender=(soc=1.0,), first=(capacity_wh=capacity,))
            receiver = SM.PowerSystems.battery(spacecraft[2])
            threshold = capacity * (1 - 1e-9)
            for shift in -4:4
                energy = nextfloat(threshold, shift)
                receiver.full = receiver.empty = false
                SM.PowerSystems.battery_boundary!(spacecraft[2], energy)
                at = (sc=[i == 2 ? merge(sc, (battery_wh=energy,)) : sc for (i, sc) in enumerate(state.sc)],)
                @test (scores_of(model, EnergyBalancingPolicy(), at, spacecraft)[1] > 0) == !receiver.full
            end
        end
        # A link that would deliver no power, or a negative power, scores 0, also from
        # an emptier battery to a fuller one.
        spacecraft, model, state = energy_fixture(first=(soc=0.95,))
        ctx = SchedulingContext(PLANET.μ, LinkContext(), spacecraft, state, nothing)
        beam = SP._link_candidate(model, 1, ctx)
        for power_w in (0.0, -5.0)
            candidate = LinkCandidate(beam.key, FixedPowerLink(power_w), beam.a, beam.b, beam.terminal_a, beam.terminal_b)
            @test candidate_score(EnergyBalancingPolicy(), candidate, ctx, 0.0) === 0.0
        end
        spacecraft, model, state = energy_fixture(sender=(battery=false,))
        @test all(iszero, scores_of(model, EnergyBalancingPolicy(), state, spacecraft))
        spacecraft, model, state = energy_fixture(first=(battery=false,))
        @test scores_of(model, EnergyBalancingPolicy(), state, spacecraft)[1] == 0.0
        @test all(iszero, scores_of(model, EnergyBalancingPolicy(), state))   # no spacecraft models
        # A momentum link scores 0 under energy balancing.
        spacecraft = [beam_sat(1, 0.0; sender=true), beam_sat(2, -2000.0; sender=false, momentum=true)]
        model = InterlinkModel()
        register_link!(model, spacecraft, (1, 1), (2, 1), LaserMomentumLink(only(mounted_link_terminals(spacecraft[1])).terminal))
        state = beam_state(spacecraft)
        update_availability!(model, spacecraft, state)
        @test scores_of(model, EnergyBalancingPolicy(), state, spacecraft) == [0.0]

        # In a run: the emptier receiver is charged first; once it passes the other, the other is.
        spacecraft = [beam_sat(1, 0.0; sender=true, soc=0.9), beam_sat(2, 5000.0; sender=false, capacity_wh=0.05, soc=0.30),
                      beam_sat(3, -5000.0; sender=false, soc=0.35)]
        model = InterlinkModel()
        register_link!(model, spacecraft, (1, 1), (2, 1), BEAM)
        register_link!(model, spacecraft, (1, 1), (3, 1), BEAM)
        args = run_args(EnergyBalancingPolicy(); duration=10.0)
        args = SM.SimConfig._with_configuration(args; interlink_model=model,
            dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)))
        solution = run_simulation(args; return_solution=true)
        history = solution.prob.p.args.interlink_model.history
        @test first(history).active == [((1, 1), (2, 1))]
        @test any(sample -> sample.active == [((1, 1), (3, 1))], history)

        # In a run, a sender that its solar array keeps full beams to a receiver until
        # the receiver is full, and never again; PowerDeliveryPolicy keeps beaming.
        ephemerides = SM.SimpleEphemeridesModel()
        sun = SM.ephemerides_sun_ii(ephemerides, PLANET, SM.ephemerides_epoch_seconds(ephemerides, InitialTime()))
        function full_receiver_run(policy; capacity_wh=0.05, soc=0.9)
            sender = beam_sat(1, 0.0; sender=true, soc=1.0)
            add_power!(sender, sender.root, SolarArray(normal=normalize(sun), area_m2=3.0, efficiency=0.3))
            spacecraft = [sender, beam_sat(2, 5000.0; sender=false, capacity_wh, soc)]
            model = InterlinkModel()
            register_link!(model, spacecraft, (1, 1), (2, 1), BEAM)
            args = SM.SimConfig._with_configuration(run_args(policy; duration=20.0); interlink_model=model,
                dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)))
            solution = run_simulation(args; return_solution=true)
            full_at = findfirst(u -> u.sc[2].battery_wh >= capacity_wh * (1 - 1e-9), solution.u)
            full_at === nothing && return (full=false, beamed_after=0)
            history = solution.prob.p.args.interlink_model.history
            return (full=SM.PowerSystems.battery(solution.prob.p.args.dynamics_model.spacecraft[2]).full,
                beamed_after=count(sample -> sample.time > solution.t[full_at] && !isempty(sample.active), history))
        end
        balancing = full_receiver_run(EnergyBalancingPolicy())
        @test balancing.full && balancing.beamed_after == 0
        @test full_receiver_run(PowerDeliveryPolicy()).beamed_after > 0
        # A receiver that starts exactly at the threshold is full from the start, and
        # is never beamed into.
        at_threshold = (capacity_wh=0.033450261609188466, soc=1 - 1e-9)
        balancing = full_receiver_run(EnergyBalancingPolicy(); at_threshold...)
        @test balancing.full && balancing.beamed_after == 0
        @test full_receiver_run(PowerDeliveryPolicy(); at_threshold...).beamed_after > 0
    end

    @testset "weighted mixes" begin
        spacecraft, model, state = scoring_fixture()
        base = SemimajorAxisRatePolicy(targets=:both)
        # One term of weight 1 and scale 1 scores, and selects, as its policy does.
        for policy in (base, SemimajorAxisRatePolicy(targets=:lower_index), ShortLinkPolicy())
            alone, mixed = scores_of(model, policy, state), scores_of(model, WeightedPolicy(PolicyTerm(policy; scale=1.0)), state)
            @test mixed == alone
            @test selected_by(model, WeightedPolicy(PolicyTerm(policy; scale=1.0)), state) == selected_by(model, policy, state)
        end
        # A scale and a weight scaled by one power of two leave the scores exactly as they were.
        @test scores_of(model, WeightedPolicy(PolicyTerm(base; scale=2.0^-3, weight=2.0^-3)), state) == scores_of(model, base, state)
        # The terms are summed in order: 1e16 + 1 + 1 rounds to 1e16, 1 + 1 + 1e16 does not.
        ordered = WeightedPolicy(PolicyTerm(ConstantPolicy(1e16); scale=1.0), PolicyTerm(ConstantPolicy(1.0); scale=1.0),
                                 PolicyTerm(ConstantPolicy(1.0); scale=1.0))
        @test scores_of(model, ordered, state)[1] === 1e16
        # Nested mixes.
        nested = WeightedPolicy(PolicyTerm(WeightedPolicy(PolicyTerm(base; scale=1.0)); scale=1.0))
        @test scores_of(model, nested, state) == scores_of(model, base, state)
        # A nonfinite term is refused, naming its policy.
        err = refusal(() -> scores_of(model, WeightedPolicy(PolicyTerm(ConstantPolicy(Inf); scale=1.0)), state))
        @test err isa ArgumentError && occursin("ConstantPolicy", err.msg)
        for bad in ((scale=0.0,), (scale=-1.0,), (scale=Inf,), (scale=NaN,), (scale=1.0, weight=Inf), (scale=1.0, weight=NaN))
            @test_throws ArgumentError PolicyTerm(base; bad...)
        end
        @test_throws ArgumentError WeightedPolicy()

        # Momentum against power on one terminal: the weights decide.
        spacecraft = [beam_sat(1, 0.0; sender=true), beam_sat(2, -2000.0; sender=false, momentum=true),
                      beam_sat(3, 5000.0; sender=false)]
        model = InterlinkModel()
        register_link!(model, spacecraft, (1, 1), (2, 1), LaserMomentumLink(only(mounted_link_terminals(spacecraft[1])).terminal))
        register_link!(model, spacecraft, (1, 1), (3, 1), BEAM)
        state = beam_state(spacecraft)
        update_availability!(model, spacecraft, state)
        sma, power = scores_of(model, SemimajorAxisRatePolicy(), state, spacecraft), scores_of(model, PowerDeliveryPolicy(), state, spacecraft)
        @test sma[1] > 0 && sma[2] == 0 && power[1] == 0 && power[2] > 0
        mix(w_sma, w_power) = WeightedPolicy(PolicyTerm(SemimajorAxisRatePolicy(); scale=sma[1], weight=w_sma),
                                             PolicyTerm(PowerDeliveryPolicy(); scale=power[2], weight=w_power))
        @test selected_by(model, mix(2.0, 1.0), state, spacecraft) == [((1, 1), (2, 1))]
        @test selected_by(model, mix(1.0, 2.0), state, spacecraft) == [((1, 1), (3, 1))]
        @test scores_of(model, mix(2.0, 1.0), state, spacecraft) == [2.0 * (sma[1] / sma[1]), 1.0 * (power[2] / power[2])]

        # An eccentricity term keeps its near-circular refusal over the mix's candidates.
        circular = [(n_terminal=2,), (n_terminal=1,), (n_terminal=1,)]
        model = InterlinkModel()
        register_link!(model, circular, (1, 1), (2, 1), LaserMomentumLink(power_w=40.0, magnification=1.0, range_m=1e4))
        register_link!(model, circular, (1, 2), (3, 1), LaserMomentumLink(power_w=40.0, magnification=1.0, range_m=1e4))
        state = (sc=[(pos=[RADIUS, 0.0, 0.0], vel=[0.0, SPEED, 0.0], mass=100.0),
                     (pos=[RADIUS, 1000.0, 0.0], vel=[0.0, SPEED, 0.0], mass=100.0),
                     (pos=[RADIUS, -1000.0, 0.0], vel=[0.0, SPEED, 0.0], mass=100.0)],)
        update_availability!(model, circular, state)
        @test_throws ArgumentError scores_of(model, EccentricityRatePolicy(), state)
        err = refusal(() -> scores_of(model, WeightedPolicy(PolicyTerm(EccentricityRatePolicy(); scale=1e-9),
            PolicyTerm(ConstantPolicy(0.0); scale=1.0)), state))
        @test err isa ArgumentError && occursin("near-circular", err.msg)

        # The setup check reaches every term.
        err = refusal(() -> run_simulation(run_args(WeightedPolicy(PolicyTerm(ShortLinkPolicy(); scale=1.0),
            PolicyTerm(NoScorePolicy(); scale=1.0)))))
        @test err isa ArgumentError && occursin("NoScorePolicy has no candidate_score method", err.msg)
        # A mix in a run.
        solution = run_simulation(run_args(WeightedPolicy(PolicyTerm(ShortLinkPolicy(); scale=1e-3))); return_solution=true)
        @test all(sample -> sample.active == [((1, 1), (2, 1))], solution.prob.p.args.interlink_model.history)
    end

    @testset "checkpoint entry of energy balancing and mixes" begin
        spacecraft = [(n_terminal=2,) for _ in 1:3]
        model = InterlinkModel()
        register_link!(model, spacecraft, (1, 1), (2, 1), LaserMomentumLink(power_w=12.0, magnification=3.0, range_m=100.0))
        args(policy) = (interlink_model=model, scheduling_policy_model=SchedulingPolicyModel(policy),
            dynamics_model=(spacecraft=spacecraft,))
        state = (sc=[(pos=[Float64(i), 0.0, 0.0],) for i in 1:3],)
        checkpoint(policy) = (u=state, runtime_state=SE._interlink_checkpoint_runtime_state(args(policy), trues(3)))
        resumes(policy, snapshot) = try
            SE._restore_interlink_checkpoint_state!(args(policy), snapshot, trues(3)); true
        catch err
            err isa ArgumentError || rethrow(); false
        end
        entry = SE._interlink_checkpoint_configuration(args(EnergyBalancingPolicy()))
        @test entry.policy.kind.name === :EnergyBalancingPolicy && entry.target_idx === nothing
        @test resumes(EnergyBalancingPolicy(), checkpoint(EnergyBalancingPolicy()))
        @test !resumes(PowerDeliveryPolicy(), checkpoint(EnergyBalancingPolicy()))
        mixture(; weight=0.5, scale=100.0, targets=[1, 2]) = WeightedPolicy(PolicyTerm(SemimajorAxisRatePolicy(; targets); scale=1.0),
            PolicyTerm(EnergyBalancingPolicy(); scale, weight))
        entry = SE._interlink_checkpoint_configuration(args(mixture()))
        @test entry.policy.kind.name === :WeightedPolicy
        @test entry.policy.fields.terms[2].fields.weight.bits == reinterpret(UInt64, 0.5)
        @test entry.policy.fields.terms[1].fields.policy.kind.name === :SemimajorAxisRatePolicy
        @test entry.policy.fields.terms[1].kind.name === :PolicyTerm   # a term is recorded as a policy is
        @test resumes(mixture(), checkpoint(mixture()))
        @test !resumes(mixture(weight=0.25), checkpoint(mixture()))
        @test !resumes(mixture(scale=50.0), checkpoint(mixture()))
        @test !resumes(mixture(targets=[1, 3]), checkpoint(mixture()))
        reordered = WeightedPolicy(PolicyTerm(EnergyBalancingPolicy(); scale=100.0, weight=0.5),
            PolicyTerm(SemimajorAxisRatePolicy(targets=[1, 2]); scale=1.0))
        @test !resumes(reordered, checkpoint(mixture()))
    end
end

end # module InterlinkPolicyTests
