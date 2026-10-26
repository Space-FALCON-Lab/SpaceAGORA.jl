module InterlinkStorageTests

# Onboard data storage: the store and the sources, the setup checks, the flow
# control at a store (an empty store passes what comes in on, a full one takes
# in no more than it sends), relays, production lost into a full store, a source
# of your own, the columns, the checkpoint state, scheduling by backlog, and a
# negative control.

using Test
using LinearAlgebra
using StaticArrays
using TOML
using DataFrames
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PLANET = make_no_gram_planet(:earth)
const RADIUS = PLANET.Rp_e + 700e3

# A source of your own: a rate that rises and falls with time.
struct WaveSource <: AbstractDataSource
    rate_bps::Float64
    period_s::Float64
end
SpaceAGORA.data_production_bps(source::WaveSource, x, ctx, t) = source.rate_bps * (1 + sin(2π * t / source.period_s)) / 2
# A source kind without a production rate, and one with a negative rate.
struct SilentSource <: AbstractDataSource end
struct DrainSource <: AbstractDataSource end
SpaceAGORA.data_production_bps(::DrainSource, x, ctx, t) = -1.0

# Record which directed endpoint and current context a backlog policy samples.
struct BacklogProbeSource <: AbstractDataSource
    rate_bps::Float64
    calls::Vector{Any}
end
function SpaceAGORA.data_production_bps(source::BacklogProbeSource, x, ctx, t)
    push!(source.calls, (x, ctx, t))
    return t < 5.0 ? source.rate_bps * x.mass_kg / 100.0 : 0.0
end

# A 100 kg bus `arc_m` along a circular orbit, with `n_terminal` virtual
# terminals, and a data store and sources when given.
function data_sat(index, arc_m; n_terminal=1, store=nothing, sources=())
    bus = Link(root=true, m=100.0)
    θ = arc_m / RADIUS
    speed = sqrt(PLANET.μ / RADIUS)
    ic = CartesianInitialCondition(SVector(RADIUS * cos(θ), RADIUS * sin(θ), 0.0),
        SVector(-speed * sin(θ), speed * cos(θ), 0.0); q=SVector(0.0, 0.0, 0.0, 1.0))
    sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia, initial_condition=ic, id=index,
        n_terminal=n_terminal)
    store === nothing || add_data!(sc, bus, store)
    foreach(source -> add_data!(sc, bus, source), sources)
    return sc
end

fixed(rate) = DataLink(transmit_power_w=0.0, budget=FixedRateBudget(rate_bps=rate), range_m=1e7)

# Links `(from, to, rate)` between terminals: each spacecraft's first terminal,
# or `(spacecraft, terminal)` pairs.
function data_model(spacecraft, links)
    model = InterlinkModel()
    endpoint(e) = e isa Integer ? (e, 1) : e
    for (from, to, rate) in links
        register_link!(model, spacecraft, endpoint(from), endpoint(to), fixed(rate))
    end
    return model
end

function run_args(spacecraft, model; duration=60.0, solver_mode=:tsit5,
                  settings=SimulationSettings(results=false, generate_plots=false))
    return SimulationConfiguration(
        simulation_settings=settings,
        mission_configuration=MissionConfiguration(mission_type=MissionTime, mission_time=duration,
            orientation_sim=false, data_rate=1.0),
        environment_model=make_no_gram_environment(planet=PLANET),
        dynamics_model=DynamicsModel(spacecraft, (InverseSquaredGravityModel(),)),
        guidance_model=GuidanceModel((), Float64[]), navigation_model=NavigationModel((), Float64[]),
        control_model=ControlModel((), Float64[]), initial_time=InitialTime(),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-10, abstol_orbit=1e-8, dt_max_orbit=1.0),
        solver_config=SM.SolverConfig(solver_mode=solver_mode), interlink_model=model,
        scheduling_policy_model=DataThroughputPolicy())
end

run(spacecraft, model; kwargs...) = run_simulation(run_args(spacecraft, model; kwargs...);
    return_results=true, return_solution=true)
refusal(f) = try f(); nothing catch err err end
run_store(solution, i) = SM.data_store(solution.prob.p.args.dynamics_model.spacecraft[i])

# The store's balance at every saved row: stored = initial + produced + received − sent − lost.
function balanced(table, i, initial; atol)
    column(name) = getproperty(table, Symbol("sc$(i)_", name))
    received = hasproperty(table, Symbol("sc$(i)_data_received_bits")) ? column(:data_received_bits) : zero(table.time)
    sent = hasproperty(table, Symbol("sc$(i)_data_sent_bits")) ? column(:data_sent_bits) : zero(table.time)
    residual = column(:data_stored_bits) .- (initial .+ column(:data_produced_bits) .+ received .- sent .- column(:data_lost_bits))
    return maximum(abs, residual) <= atol
end

# One sender holding data and producing, one receiver.
function one_link(; c=2e6, production=1e5, initial=1e7)
    spacecraft = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e8, initial_bits=initial),
                           sources=(DataSource(rate_bps=production),)),
                  data_sat(2, 100e3; store=DataStore(capacity_bits=5e7))]
    return spacecraft, data_model(spacecraft, [(1, 2, c)])
end

@testset "Data storage" begin
    @testset "parts" begin
        @test_throws ArgumentError DataStore(capacity_bits=0.0)
        @test_throws ArgumentError DataStore(capacity_bits=Inf)
        @test_throws ArgumentError DataStore(capacity_bits=1e6, initial_bits=2e6)
        @test_throws ArgumentError DataStore(capacity_bits=1e6, initial_bits=-1.0)
        @test_throws ArgumentError DataSource(rate_bps=-1.0)
        @test_throws ArgumentError DataSource(rate_bps=NaN)
        sc = data_sat(1, 0.0)
        @test !SM.has_data_store(sc) && SM.data_store(sc) === nothing
        store = add_data!(sc, sc.root, DataStore(capacity_bits=1e6, initial_bits=1e5))
        @test SM.has_data_store(sc) && SM.data_store(sc) === store && !store.full && !store.empty
        @test_throws ArgumentError add_data!(sc, sc.root, DataStore(capacity_bits=1e6))   # at most one
        @test_throws ArgumentError add_data!(sc, Link(root=true, m=1.0), DataSource(rate_bps=1.0))   # not its link
        source = add_data!(sc, sc.root, DataSource(rate_bps=5.0, name="camera"))
        @test only(sc.root.data_sources) === source && source.name == "camera"
        x = SM.InterLinkModels.link_endpoint((pos=[RADIUS, 0.0, 0.0],), nothing, 1)
        @test data_production_bps(source, x, LinkContext(), 0.0) == 5.0
        @test data_production_bps(WaveSource(4.0, 8.0), x, LinkContext(), 2.0) ≈ 4.0
        # The spacecraft with data parts are claimed by the constellation builder's sharing check.
        @test SM.DataStorage._has_data_sources(sc)
    end

    @testset "setup refusals" begin
        no_store = [data_sat(1, 0.0; sources=(DataSource(rate_bps=1.0),)), data_sat(2, 100e3)]
        err = refusal(() -> run_simulation(run_args(no_store, data_model(no_store, [(1, 2, 1e6)]))))
        @test err isa ArgumentError && occursin("Spacecraft 1 has a data source but no data store", err.msg)
        silent = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e6), sources=(SilentSource(),))]
        err = refusal(() -> run_simulation(run_args(silent, nothing)))
        @test err isa ArgumentError && occursin("SilentSource has no data_production_bps method", err.msg)
        split = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e6), sources=(DataSource(rate_bps=1.0),))]
        err = refusal(() -> run_simulation(run_args(split, nothing; solver_mode=:split_imex)))
        @test err isa ArgumentError && occursin("data store needs a solver mode that integrates the full right-hand side", err.msg)
        drain = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e6), sources=(DrainSource(),))]
        err = refusal(() -> run_simulation(run_args(drain, nothing)))
        @test err isa ArgumentError && occursin("finite and nonnegative", err.msg)
    end

    @testset "one link" begin
        c, production, initial = 2e6, 1e5, 1e7
        result = run(one_link(; c, production, initial)...)
        table = result.table
        empty_at = initial / (c - production)
        # A sender holding data sends at the link's rate; once empty, it passes its production on.
        for (t, rate) in zip(table.time, table.sc1_data_rate_tx_bps)
            t < empty_at - 1.5 && @test rate == c
            t > empty_at + 1.5 && @test rate ≈ production rtol=1e-12
        end
        # The stored data follows the closed form through the event (the solution;
        # a saved row inside the step an event cuts short can hold another value).
        @test all(t -> isapprox(result.solution(t).sc[1].data_stored_bits, max(initial - (c - production) * t, 0.0);
            atol=1e-6 * initial), 0.0:0.25:60.0)
        @test run_store(result.solution, 1).empty && !run_store(result.solution, 2).empty
        # The receiver holds what was sent, and nothing is lost.
        @test table.sc2_data_received_bits[end] == table.sc1_data_sent_bits[end]
        @test table.sc2_data_stored_bits[end] ≈ initial + production * 60.0 rtol=1e-9
        @test table.sc1_data_lost_bits[end] == 0.0 && table.sc2_data_lost_bits[end] == 0.0
        @test balanced(table, 1, initial; atol=1e-6) && balanced(table, 2, 0.0; atol=1e-6)
        @test table.sc2_data_store_fill ≈ table.sc2_data_stored_bits ./ 5e7
    end

    @testset "a full receiver" begin
        c, production = 2e6, 1e5
        spacecraft = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e9, initial_bits=1e8)),
                      data_sat(2, 100e3; store=DataStore(capacity_bits=1e7), sources=(DataSource(rate_bps=production),))]
        result = run(spacecraft, data_model(spacecraft, [(1, 2, c)]))
        table = result.table
        full_at = 1e7 / (c + production)
        # The receiver fills, then takes in nothing more: the link stops, the sender
        # keeps its data, and the receiver's own production is lost.
        @test run_store(result.solution, 2).full
        @test table.sc1_data_rate_tx_bps[end] == 0.0
        @test table.sc1_data_sent_bits[end] ≈ c * full_at rtol=1e-6
        @test table.sc1_data_stored_bits[end] ≈ 1e8 - table.sc1_data_sent_bits[end] rtol=1e-12
        @test table.sc1_data_lost_bits[end] == 0.0
        @test table.sc2_data_stored_bits[end] ≈ 1e7 rtol=1e-9
        @test table.sc2_data_store_fill[end] ≈ 1.0 atol=1e-9
        @test table.sc2_data_lost_bits[end] ≈ production * (60.0 - full_at) rtol=1e-6
        @test balanced(table, 1, 1e8; atol=1e-3) && balanced(table, 2, 0.0; atol=1e-3)
    end

    @testset "a full store that drains takes data again" begin
        # Full at the start, sending on: its own production is lost only until the
        # step-end update finds room again.
        rate = 1e4
        spacecraft = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e7, initial_bits=1e7), sources=(DataSource(rate_bps=rate),)),
                      data_sat(2, 100e3; store=DataStore(capacity_bits=1e9))]
        result = run(spacecraft, data_model(spacecraft, [(1, 2, 1e5)]))
        table = result.table
        @test !run_store(result.solution, 1).full
        @test table.sc1_data_lost_bits[end] <= rate * 1.0
        @test table.sc1_data_stored_bits[end] ≈ 1e7 - (1e5 - rate) * 60.0 rtol=1e-3
        @test balanced(table, 1, 1e7; atol=1e-3)
    end

    @testset "a chain of relays" begin
        # A holds data; B, empty, relays it on to C on a faster link: B stays empty
        # and C receives what A sends.
        c1, c2 = 1e6, 3e6
        spacecraft = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e9, initial_bits=1e8)),
                      data_sat(2, 100e3; n_terminal=2, store=DataStore(capacity_bits=1e8)),
                      data_sat(3, 200e3; store=DataStore(capacity_bits=1e9))]
        result = run(spacecraft, data_model(spacecraft, [(1, (2, 1), c1), ((2, 2), 3, c2)]))
        table = result.table
        @test maximum(table.sc2_data_stored_bits) <= 1e-6 * c1
        @test table.sc2_data_rate_tx_bps[end] ≈ c1 rtol=1e-12
        @test table.sc3_data_received_bits[end] ≈ c1 * 60.0 rtol=1e-9
        @test table.sc3_data_received_bits[end] == table.sc2_data_sent_bits[end]
        @test all(i -> balanced(table, i, (1e8, 0.0, 0.0)[i]; atol=1e-3), 1:3)
    end

    @testset "production into a full store is lost" begin
        rate = 4e4
        spacecraft = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e6, initial_bits=2e5), sources=(DataSource(rate_bps=rate),))]
        table = run(spacecraft, nothing).table
        full_at = 8e5 / rate
        @test table.sc1_data_stored_bits[end] ≈ 1e6 rtol=1e-9
        @test table.sc1_data_produced_bits[end] ≈ rate * 60.0 rtol=1e-9
        @test table.sc1_data_lost_bits[end] ≈ rate * (60.0 - full_at) rtol=1e-6
        @test table.sc1_data_lost_bits ≈ table.sc1_data_produced_bits .- (table.sc1_data_stored_bits .- 2e5) atol=1e-3
        @test table.sc1_data_store_fill[end] ≈ 1.0 atol=1e-9
        # No data link: no counters.
        @test !hasproperty(table, :sc1_data_sent_bits)
    end

    @testset "a source of your own" begin
        source = WaveSource(1e5, 40.0)
        spacecraft = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e9), sources=(source, DataSource(rate_bps=1e3)))]
        table = run(spacecraft, nothing).table
        T = 60.0
        wave = source.rate_bps * T / 2 + source.rate_bps * source.period_s / (4π) * (1 - cos(2π * T / source.period_s))
        @test table.sc1_data_produced_bits[end] ≈ wave + 1e3 * T rtol=1e-6
        @test table.sc1_data_stored_bits[end] ≈ table.sc1_data_produced_bits[end] rtol=1e-12
    end

    @testset "columns" begin
        spacecraft = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e8, initial_bits=1e7)), data_sat(2, 100e3)]
        table = run(spacecraft, data_model(spacecraft, [(1, 2, 2e5)]); duration=10.0).table
        for name in (:data_stored_bits, :data_store_fill, :data_produced_bits, :data_lost_bits)
            @test hasproperty(table, Symbol("sc1_", name)) && all(isnan, getproperty(table, Symbol("sc2_", name)))
        end
        # The state has one shape: the spacecraft without a store carries the store states, at zero.
        @test run(spacecraft, data_model(spacecraft, [(1, 2, 2e5)]); duration=10.0).solution.u[end].sc[2].data_stored_bits == 0.0
        # A receiver without a store takes everything the link carries.
        @test table.sc2_data_rate_rx_bps[end] == 2e5
        @test table.sc2_data_received_bits[end] ≈ 2e6 rtol=1e-9
    end

    @testset "checkpoint state" begin
        spacecraft, model = one_link()
        args = run_args(spacecraft, model)
        state = SE._callback_checkpoint_state(args)
        @test state.configuration[1].data_store == (1e8, 1e7) && state.configuration[2].data_store == (5e7, 0.0)
        # The flags are written and restored.
        SM.data_store(spacecraft[1]).empty = true
        SM.data_store(spacecraft[2]).full = true
        saved = (runtime_state=(callbacks=SE._callback_checkpoint_state(args),),)
        @test saved.runtime_state.callbacks.spacecraft[1].data_store == (full=false, empty=true)
        fresh, fresh_model = one_link()
        parameters = (shared_buffers=(callback_resume=Ref(false),),)
        SE._restore_callback_checkpoint_state!(run_args(fresh, fresh_model), saved, parameters)
        @test SM.data_store(fresh[1]).empty && !SM.data_store(fresh[1]).full
        @test SM.data_store(fresh[2]).full && !SM.data_store(fresh[2]).empty
        @test parameters.shared_buffers.callback_resume[]
        # Another capacity or initial content refuses the checkpoint.
        for store in (DataStore(capacity_bits=2e8, initial_bits=1e7), DataStore(capacity_bits=1e8, initial_bits=2e7))
            other = [data_sat(1, 0.0; store, sources=(DataSource(rate_bps=1e5),)), data_sat(2, 100e3; store=DataStore(capacity_bits=5e7))]
            err = refusal(() -> SE._restore_callback_checkpoint_state!(run_args(other, data_model(other, [(1, 2, 2e6)])), saved,
                (shared_buffers=(callback_resume=Ref(false),),)))
            @test err isa ArgumentError && occursin("data stores differ from the checkpoint's", err.msg)
        end
        # A run without a store writes the snapshot it wrote before data stores.
        plain = [data_sat(1, 0.0), data_sat(2, 100e3)]
        add_power!(plain[1], plain[1].root, Battery(capacity_wh=10.0, initial_soc=0.5))
        plain_state = SE._callback_checkpoint_state(run_args(plain, data_model(plain, [(1, 2, 2e6)])))
        @test keys(plain_state.configuration[1]) == (:battery, :faults, :sensors, :hinges, :gimbals, :thrusters)
        @test keys(plain_state.spacecraft[1]) == (:battery, :faults, :sensors, :hinges, :gimbals, :thrusters)
    end

    @testset "a resumed run" begin
        straight = run(one_link()...).table
        mktempdir() do directory
            segments(resume) = SimulationSettings(results=resume, results_directory=directory, generate_plots=false,
                checkpoint_enabled=true, checkpoint_interval_s=10.0, checkpoint_directory=directory,
                resume_from_checkpoint=resume)
            run(one_link()...; duration=30.0, settings=segments(false))
            resumed = run(one_link()...; settings=segments(true))
            # The sender emptied in the first leg; the resumed run keeps it empty and passes its production on.
            @test run_store(resumed.solution, 1).empty
            @test resumed.table.sc1_data_rate_tx_bps[end] ≈ 1e5 rtol=1e-12
            for name in (:sc1_data_stored_bits, :sc1_data_produced_bits, :sc2_data_stored_bits, :sc2_data_received_bits)
                @test getproperty(resumed.table, name)[end] ≈ getproperty(straight, name)[end] rtol=1e-9 atol=1e-6
            end
        end
    end

    @testset "scheduling by backlog" begin
        # Three spacecraft at a point: two senders, one receiver, links from each sender.
        function backlog_fixture(; held=(8e6, 1e6), production=(0.0, 0.0), room_capacity=1e8, room_held=0.0, rate=2e6)
            spacecraft = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e9, initial_bits=held[1]),
                sources=(DataSource(rate_bps=production[1]),)),
                          data_sat(2, 50e3; store=DataStore(capacity_bits=1e9, initial_bits=held[2]),
                sources=(DataSource(rate_bps=production[2]),)),
                          data_sat(3, 100e3; store=DataStore(capacity_bits=room_capacity, initial_bits=room_held))]
            model = data_model(spacecraft, [(1, 3, rate), (2, 3, rate)])
            state = (sc=[(pos=collect(sc.initial_condition.pos), vel=collect(sc.initial_condition.vel), mass=100.0,
                          data_stored_bits=SM.data_store(sc).initial_bits) for sc in spacecraft],)
            update_availability!(model, spacecraft, state)
            return spacecraft, model, state
        end
        function scores(policy, spacecraft, model, state; context=LinkContext())
            scheduler = SchedulingPolicyModel(policy)
            score_candidates!(model, scheduler, state, PLANET.μ; spacecraft, context)
            return copy(SM.SchedulingPolicies.matching_workspace!(model, scheduler).score), select_interlinks!(model, scheduler)
        end
        # The accessors.
        spacecraft, model, state = backlog_fixture()
        ctx = SchedulingContext(PLANET.μ, LinkContext(), spacecraft, state, nothing)
        @test stored_bits(ctx, 1) === 8e6 && store_room_bits(ctx, 3) === 1e8
        plain = [data_sat(1, 0.0), data_sat(2, 50e3)]
        plain_state = (sc=[(pos=collect(sc.initial_condition.pos),) for sc in plain],)
        plain_ctx = SchedulingContext(PLANET.μ, LinkContext(), plain, plain_state, nothing)
        @test isnan(stored_bits(plain_ctx, 1)) && store_room_bits(plain_ctx, 1) == Inf
        @test isnan(stored_bits(SchedulingContext(PLANET.μ, LinkContext(), nothing, state, nothing), 1))
        # The score: the least of the rate, the held data and the room, over the horizon.
        policy = DataBacklogPolicy(horizon_s=10.0)
        values, selected = scores(policy, spacecraft, model, state)
        @test values == [min(2e6, 8e6 / 10, 1e8 / 10), min(2e6, 1e6 / 10, 1e8 / 10)]
        # The fuller store wins the receiver's one terminal.
        @test selected == [((1, 1), (3, 1))]
        # The receiver's room limits it, and a full receiver or an empty sender scores 0.
        values, _ = scores(policy, backlog_fixture(room_capacity=1e8, room_held=1e8 - 2e6)...)
        @test values == [2e6 / 10, 1e6 / 10]
        # A store within its empty or full tolerance holds nothing, or has no room.
        values, _ = scores(policy, backlog_fixture(held=(1e-3, 1e6))...)
        @test values[1] == 0.0
        values, _ = scores(policy, backlog_fixture(room_capacity=1e8, room_held=1e8 * (1 - 1e-12))...)
        @test values == [0.0, 0.0]
        values, selected = scores(policy, backlog_fixture(room_capacity=1e8, room_held=1e8)...)
        @test values == [0.0, 0.0] && isempty(selected)
        values, selected = scores(policy, backlog_fixture(held=(0.0, 3e6))...)
        @test values[1] == 0.0 && selected == [((2, 1), (3, 1))]
        # Current production adds to held-data/horizon, including at an empty store.
        values, _ = scores(policy, backlog_fixture(production=(2e5, 0.0))...)
        @test values == [1e6, 1e5]
        values, selected = scores(policy, backlog_fixture(held=(0.0, 0.0), production=(2e5, 0.0))...)
        @test values == [2e5, 0.0] && selected == [((1, 1), (3, 1))]
        values, _ = scores(policy, backlog_fixture(production=(3e6, 0.0))...)
        @test values[1] == 2e6
        values, _ = scores(policy, backlog_fixture(production=(3e6, 0.0), room_held=99e6)...)
        @test values == [1e5, 1e5]
        values, selected = scores(policy, backlog_fixture(production=(3e6, 0.0), room_held=1e8)...)
        @test values == [0.0, 0.0] && isempty(selected)

        # Directed sender 2 is sampled with its own state and both current sources.
        calls = Any[]
        reversed = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e8),
                        sources=(DataSource(rate_bps=9e6),)),
                    data_sat(2, 50e3; store=DataStore(capacity_bits=1e8),
                        sources=(BacklogProbeSource(2e5, calls), BacklogProbeSource(5e4, calls)))]
        reverse_model = data_model(reversed, [(2, 1, 2e6)])
        reverse_state = (sc=[(pos=collect(sc.initial_condition.pos), vel=collect(sc.initial_condition.vel),
                             mass=100.0, data_stored_bits=0.0) for sc in reversed],)
        update_availability!(reverse_model, reversed, reverse_state)
        context = LinkContext(time_s=2.0, planet=PLANET)
        values, selected = scores(policy, reversed, reverse_model, reverse_state; context)
        @test values == [2.5e5] && selected == [((2, 1), (1, 1))]
        @test length(calls) == 2
        @test all(call -> call[1].spacecraft_index == 2 && call[1].spacecraft === reversed[2] &&
            call[1].pos_ii == reversed[2].initial_condition.pos && call[2] === context && call[3] == 2.0, calls)
        values, selected = scores(policy, reversed, reverse_model, reverse_state;
            context=LinkContext(time_s=6.0, planet=PLANET))
        @test values == [0.0] && isempty(selected)

        # A sender without a store always holds data, and a receiver without one has room.
        open_ends = [data_sat(1, 0.0), data_sat(2, 50e3)]
        open_model = data_model(open_ends, [(1, 2, 2e6)])
        open_state = (sc=[(pos=collect(sc.initial_condition.pos), vel=collect(sc.initial_condition.vel), mass=100.0)
                          for sc in open_ends],)
        update_availability!(open_model, open_ends, open_state)
        @test first(scores(policy, open_ends, open_model, open_state)) == [2e6]
        # It shares DataThroughputPolicy's unit: a mix weighs the two directly.
        mixed = WeightedPolicy(PolicyTerm(policy; scale=1.0), PolicyTerm(DataThroughputPolicy(); scale=1.0, weight=0.5))
        values, _ = scores(mixed, backlog_fixture()...)
        @test values == [8e5 + 0.5 * 2e6, 1e5 + 0.5 * 2e6]
        @test_throws ArgumentError DataBacklogPolicy(horizon_s=0.0)
        @test_throws ArgumentError DataBacklogPolicy(horizon_s=Inf)
        # A checkpoint records the horizon.
        @test SE._scheduling_policy_entry(policy) == SE._scheduling_policy_entry(DataBacklogPolicy(horizon_s=10.0))
        @test SE._scheduling_policy_entry(policy) != SE._scheduling_policy_entry(DataBacklogPolicy(horizon_s=20.0))
        # In a run, the link closes once the sender has nothing left to send.
        spacecraft = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e8, initial_bits=1e7)),
                      data_sat(2, 100e3; store=DataStore(capacity_bits=1e9))]
        args = run_args(spacecraft, data_model(spacecraft, [(1, 2, 2e6)]))
        args = SM.SimConfig._with_configuration(args; scheduling_policy_model=DataBacklogPolicy(horizon_s=5.0))
        solution = run_simulation(args; return_solution=true)
        history = solution.prob.p.args.interlink_model.history
        @test first(history).active == [((1, 1), (2, 1))]
        @test isempty(last(history).active)
        @test solution.u[end].sc[2].data_received_bits ≈ 1e7 rtol=1e-6
    end

    @testset "backlog keeps current production connected" begin
        for initial in (0.0, 1e7)
            args = run_args(one_link(; production=1e5, initial)...)
            args = SM.SimConfig._with_configuration(args;
                scheduling_policy_model=DataBacklogPolicy(horizon_s=5.0))
            result = run_simulation(args; return_results=true, return_solution=true)
            history = result.configuration.interlink_model.history
            @test length(history) > 1
            @test all(sample -> sample.active == [((1, 1), (2, 1))], history)
            final = result.solution.u[end]
            @test final.sc[1].data_produced_bits ≈ 6e6 rtol=1e-9
            @test final.sc[1].data_sent_bits ≈ initial + 6e6 rtol=1e-9
            @test final.sc[2].data_received_bits ≈ initial + 6e6 rtol=1e-9
            @test final.sc[1].data_stored_bits ≈ 0.0 atol=1e-6
            @test final.sc[1].data_lost_bits ≈ 0.0 atol=1e-6
            @test balanced(result.table, 1, initial; atol=1e-6)
            @test balanced(result.table, 2, 0.0; atol=1e-6)
        end
    end

    @testset "backlog score units survive saved results" begin
        backlog = DataBacklogPolicy(horizon_s=5.0)
        mixed = WeightedPolicy(PolicyTerm(backlog; scale=1e6),
            PolicyTerm(DataThroughputPolicy(); scale=1e6))
        for policy in (backlog, mixed)
            mktempdir() do directory
                spacecraft, model = one_link()
                args = run_args(spacecraft, model; duration=1.0,
                    settings=SimulationSettings(results=true, generate_plots=false,
                        results_directory=directory, save_csv=false))
                args = SM.SimConfig._with_configuration(args; scheduling_policy_model=policy)
                result = run_simulation(args; return_results=true)
                expected_unit = policy === backlog ? "bit/s" : "1"
                table = link_schedule(result)
                @test !isempty(table)
                @test DataFrames.metadata(table)["score_unit"] == expected_unit
                manifest = TOML.parsefile(joinpath(directory, "simulation_results.manifest.toml"))
                @test manifest["interlinks"]["score_unit"] == expected_unit
                if policy === mixed
                    @test [term["score_unit"] for term in manifest["interlinks"]["terms"]] == ["bit/s", "bit/s"]
                end
            end
        end
    end

    @testset "without a store, as before" begin
        # Negative control: without stores a link carries its rate whatever either
        # end holds; an empty store with nothing coming in sends nothing.
        plain = [data_sat(1, 0.0), data_sat(2, 100e3)]
        table = run(plain, data_model(plain, [(1, 2, 2e6)]); duration=10.0).table
        @test table.sc1_data_sent_bits[end] ≈ 2e7 rtol=1e-9
        @test !hasproperty(table, :sc1_data_stored_bits)
        empty_sender = [data_sat(1, 0.0; store=DataStore(capacity_bits=1e8)), data_sat(2, 100e3)]
        table = run(empty_sender, data_model(empty_sender, [(1, 2, 2e6)]); duration=10.0).table
        @test table.sc1_data_sent_bits[end] == 0.0 && table.sc1_data_rate_tx_bps[end] == 0.0
    end
end

end # module InterlinkStorageTests
