const RELAY_IN = ((1, 1), (2, 1))
const DIRECT = ((1, 1), (3, 1))
const RELAY_OUT = ((2, 2), (3, 2))
const DATA_RADIUS_M = PLANET.Rp_e + 700e3

data_edge(key, rate; allowed=true, range=200e3) = (; key, rate, allowed, range)

function data_spec(kind; refined=false)
    name = string(kind) * (refined ? "_refined" : "")
    common = (; name, kind, duration=20.0, horizon=5.0, dtmax=refined ? 0.125 : 0.25)
    if kind == :producer_control
        return (; common..., capacity=[1e8, 1e9], initial=[0.0, 0.0],
            production=[2e5, 0.0], terminals=[1, 1], sink=2, relay=0,
            edges=[data_edge(((1, 1), (2, 1)), 1e6)])
    elseif kind == :production_competition
        return (; common..., capacity=[1e8, 1e8, 1e9], initial=[0.0, 8e6, 0.0],
            production=[2e5, 0.0, 0.0], terminals=[1, 1, 1], sink=3, relay=0,
            edges=[data_edge(((1, 1), (3, 1)), 1e6), data_edge(((2, 1), (3, 1)), 1e6)])
    elseif kind in (:busy_relay, :direct_route_control)
        return (; common..., capacity=[1e8, 2e7, 1e9], initial=[8e7, 2e7, 0.0],
            production=[0.0, 1e6, 0.0], terminals=[1, 2, 2], sink=3, relay=2,
            edges=[data_edge(RELAY_IN, 1e6; allowed=kind == :busy_relay),
                data_edge(DIRECT, 5e5), data_edge(RELAY_OUT, 1e6; range=50e3)])
    end
    throw(ArgumentError("Unknown data fixture: $kind"))
end

function data_configuration(spec, directory)
    spacecraft = map(eachindex(spec.capacity)) do i
        bus = Link(root=true, m=100.0)
        initial = InitialCondition(a=DATA_RADIUS_M, e=0.0, i=0.0, Ω=0.0,
            ν=rad2deg((i - 1) * 10e3 / DATA_RADIUS_M))
        sc = SpacecraftModel(root=bus, links=[bus], inertia_tensor=bus.inertia,
            initial_condition=initial, id=i, n_terminal=spec.terminals[i])
        add_data!(sc, bus, DataStore(capacity_bits=spec.capacity[i], initial_bits=spec.initial[i]))
        add_data!(sc, bus, DataSource(rate_bps=spec.production[i]))
        sc
    end
    model = InterlinkModel(record_available=true)
    for edge in spec.edges
        edge.allowed || continue
        link = DataLink(transmit_power_w=0.0, budget=FixedRateBudget(rate_bps=edge.rate),
            range_m=edge.range)
        register_link!(model, spacecraft, edge.key..., link)
    end
    return configuration(spacecraft, model, DataBacklogPolicy(horizon_s=spec.horizon);
        duration=spec.duration, dtmax=spec.dtmax, output=directory)
end

function edge_geometry(edge, state)
    a, b = state.sc[edge.key[1][1]].pos, state.sc[edge.key[2][1]].pos
    displacement = b - a
    distance = norm(displacement)
    distance > 0 || error("Fixture endpoints must be distinct.")
    fraction = clamp(-dot(a, displacement) / dot(displacement, displacement), 0.0, 1.0)
    clearance = norm(a + fraction * displacement) - PLANET.Rp_e
    return (; distance, clearance, feasible=distance <= edge.range && clearance > 0)
end

function data_candidates(spec, state, time_s, previous)
    geometries = edge_geometry.(spec.edges, Ref(state))
    outflow = zeros(length(spec.capacity))
    if spec.relay != 0
        outbound = findfirst(edge -> edge.key == RELAY_OUT, spec.edges)
        @test isapprox(state.sc[spec.relay].data_stored_bits, spec.capacity[spec.relay]; atol=BIT_ATOL, rtol=0)
        @test state.sc[spec.sink].data_stored_bits < spec.capacity[spec.sink] / 2
        if RELAY_OUT in previous && spec.edges[outbound].allowed && geometries[outbound].feasible
            # This fixture's relay stays full and its destination has ample room.
            outflow[spec.relay] = spec.edges[outbound].rate
        end
    end
    return DataFrame(map(eachindex(spec.edges)) do index
        edge, geometry = spec.edges[index], geometries[index]
        tx, rx = edge.key[1][1], edge.key[2][1]
        held = state.sc[tx].data_stored_bits
        room = spec.capacity[rx] - state.sc[rx].data_stored_bits
        held <= 1e-9 * spec.capacity[tx] && (held = 0.0)
        state.sc[rx].data_stored_bits >= (1 - 1e-9) * spec.capacity[rx] && (room = 0.0)
        hypothetical = min(edge.rate, held / spec.horizon + spec.production[tx],
            room / spec.horizon + outflow[rx])
        available = edge.allowed && geometry.feasible
        (; time_s, tx, tx_terminal=edge.key[1][2], rx, rx_terminal=edge.key[2][2],
            allowed=edge.allowed, distance_m=geometry.distance, range_m=edge.range,
            planet_clearance_m=geometry.clearance, geometrically_feasible=geometry.feasible,
            available, link_rate_bps=edge.rate, held_tx_bits=held,
            production_tx_bps=spec.production[tx], room_rx_bits=room,
            previous_feasible_outflow_rx_bps=outflow[rx],
            previous_active=edge.key in previous, hypothetical_score_bps=hypothetical,
            score_bps=available ? hypothetical : 0.0)
    end)
end

function all_data_matchings(candidates)
    allowed = candidates[candidates.allowed, :]
    return DataFrame(map(0:((1 << nrow(allowed)) - 1)) do mask
        indices = [i for i in 1:nrow(allowed) if !iszero(mask & (1 << (i - 1)))]
        rows = allowed[indices, :]
        keys = row_key.(eachrow(rows))
        terminal_feasible = allunique(collect(Iterators.flatten(keys)))
        available = all(rows.available)
        positive = all(>(0), rows.score_bps)
        (; mask, matching=matching_label(keys), terminal_feasible,
            all_available=available, positive_scores=positive,
            feasible=terminal_feasible && available && positive,
            total_score_bps=sum(rows.score_bps; init=0.0))
    end)
end

function set_replay_state!(model, spacecraft, state, previous)
    fill!(model.graph.active, false)
    for key in previous
        model.graph.active[model.graph.edge_index_lookup[key]] = true
    end
    for i in eachindex(spacecraft)
        SM.DataStorage._start_data_state!(spacecraft[i], state.sc[i].data_stored_bits)
    end
end

function replay_scores!(model, scheduler, spacecraft, state, previous, time_s)
    set_replay_state!(model, spacecraft, state, previous)
    context = LinkContext(time_s=time_s, planet=PLANET)
    update_availability!(model, spacecraft, state; context)
    score_candidates!(model, scheduler, state, PLANET.μ; spacecraft, context)
    return scheduler.matching.score
end

function outflow_probes(directory)
    spec = data_spec(:busy_relay)
    rows = NamedTuple[]
    @testset "Previous feasible outflow, not future outflow" begin
        for mode in (:no_previous, :previous_feasible, :previous_unavailable)
            args = data_configuration(spec, directory)
            state = SE.build_initial_conditions(args)
            previous = mode == :no_previous ? Set{LinkKey}() : Set([RELAY_OUT])
            if mode == :previous_unavailable
                angle = 100e3 / DATA_RADIUS_M
                state.sc[3].pos .= [DATA_RADIUS_M * cos(angle), DATA_RADIUS_M * sin(angle), 0.0]
                speed = sqrt(PLANET.μ / DATA_RADIUS_M)
                state.sc[3].vel .= [-speed * sin(angle), speed * cos(angle), 0.0]
            end
            candidates = data_candidates(spec, state, 0.0, previous)
            model, spacecraft = args.interlink_model, args.dynamics_model.spacecraft
            scheduler = SchedulingPolicyModel(DataBacklogPolicy(horizon_s=spec.horizon))
            actual_scores = replay_scores!(model, scheduler, spacecraft, state, previous, 0.0)
            selected = select_interlinks!(model, scheduler)
            expected = mode == :previous_feasible ? Set([RELAY_IN, RELAY_OUT]) :
                mode == :previous_unavailable ? Set([DIRECT]) : Set([DIRECT, RELAY_OUT])
            @test Set(selected) == expected
            for row in eachrow(candidates)
                key = row_key(row)
                actual = actual_scores[model.graph.edge_index_lookup[key]]
                @test isapprox(actual, row.score_bps; rtol=SCORE_RTOL, atol=DATA_SCORE_ATOL)
                if key == RELAY_IN
                    @test row.previous_feasible_outflow_rx_bps ==
                        (mode == :previous_feasible ? 1e6 : 0.0)
                end
                push!(rows, (; probe=string(mode), link=key_label(key),
                    row.available, row.previous_active, row.distance_m, row.range_m,
                    row.previous_feasible_outflow_rx_bps, expected_score_bps=row.score_bps,
                    actual_score_bps=actual, selected=key in selected))
            end
        end
    end
    CSV.write(joinpath(directory, "outflow_probes.csv"), DataFrame(rows))
end

function storage_accounting(spec, history, states)
    node_rows, network_rows = NamedTuple[], NamedTuple[]
    previous_counters = zeros(length(spec.capacity), 4)
    for (sample, state) in zip(history, states)
        produced = sent = received = lost = stored = 0.0
        for i in eachindex(spec.capacity)
            sc = state.sc[i]
            residual = sc.data_stored_bits - (spec.initial[i] + sc.data_produced_bits +
                sc.data_received_bits - sc.data_sent_bits - sc.data_lost_bits)
            @test abs(residual) <= BIT_ATOL
            @test -BIT_ATOL <= sc.data_stored_bits <= spec.capacity[i] + BIT_ATOL
            @test isapprox(sc.data_produced_bits, spec.production[i] * sample.time; atol=BIT_ATOL, rtol=0)
            counters = [sc.data_produced_bits, sc.data_sent_bits, sc.data_received_bits, sc.data_lost_bits]
            @test all(counters .>= -BIT_ATOL)
            @test all(counters .>= previous_counters[i, :] .- BIT_ATOL)
            previous_counters[i, :] .= counters
            allowed_overflow = spec.kind == :busy_relay && i == spec.relay
            allowed_overflow || @test abs(sc.data_lost_bits) <= BIT_ATOL
            push!(node_rows, (; time_s=sample.time, spacecraft=i, capacity_bits=spec.capacity[i],
                initial_bits=spec.initial[i], production_bps=spec.production[i],
                stored_bits=sc.data_stored_bits, produced_bits=sc.data_produced_bits,
                received_bits=sc.data_received_bits, sent_bits=sc.data_sent_bits,
                lost_bits=sc.data_lost_bits, permitted_overflow=allowed_overflow,
                balance_residual_bits=residual, x_m=sc.pos[1], y_m=sc.pos[2], z_m=sc.pos[3]))
            produced += sc.data_produced_bits
            sent += sc.data_sent_bits
            received += sc.data_received_bits
            lost += sc.data_lost_bits
            stored += sc.data_stored_bits
        end
        transfer_residual = sent - received
        network_residual = sum(spec.initial) + produced - stored - lost
        @test abs(transfer_residual) <= BIT_ATOL
        @test abs(network_residual) <= BIT_ATOL
        delivered = state.sc[spec.sink].data_received_bits
        @test isapprox(state.sc[spec.sink].data_stored_bits, delivered; atol=BIT_ATOL, rtol=0)
        push!(network_rows, (; time_s=sample.time, initial_bits=sum(spec.initial),
            produced_bits=produced, total_internal_sent_bits=sent,
            total_internal_received_bits=received, total_stored_bits=stored,
            lost_bits=lost, delivered_bits=delivered,
            transfer_residual_bits=transfer_residual, network_residual_bits=network_residual))
    end
    return DataFrame(node_rows), DataFrame(network_rows)
end

function hand_comparisons(spec, solution)
    history = solution.prob.p.args.interlink_model.history
    final = solution.u[end]
    rows = NamedTuple[]
    function compare(quantity, expected, actual)
        @test isapprox(actual, expected; atol=BIT_ATOL, rtol=0)
        push!(rows, (; case=spec.name, quantity, expected_bits=expected, actual_bits=actual,
            error_bits=actual - expected, tolerance_bits=BIT_ATOL))
    end
    if spec.kind == :producer_control
        @test all(sample -> sample.active == [((1, 1), (2, 1))], history)
        compare("delivered", spec.production[1] * spec.duration, final.sc[2].data_received_bits)
        compare("sender_stored", 0.0, final.sc[1].data_stored_bits)
    elseif spec.kind == :production_competition
        @test first(history).active == [((2, 1), (3, 1))]
        @test any(sample -> ((1, 1), (3, 1)) in sample.active, history)
        @test schedule_changes(history) > 0
        @test BIT_ATOL < final.sc[2].data_stored_bits < 1e6 + BIT_ATOL
        compare("delivered_from_initial_and_production_minus_backlog",
            sum(spec.initial) + sum(spec.production) * spec.duration -
                sum(final.sc[i].data_stored_bits for i in 1:2),
            final.sc[3].data_received_bits)
    else
        direct_duration = sum(history[i + 1].time - history[i].time
            for i in 1:(length(history) - 1) if DIRECT in history[i].active; init=0.0)
        relay_in_duration = sum(history[i + 1].time - history[i].time
            for i in 1:(length(history) - 1) if RELAY_IN in history[i].active; init=0.0)
        @test Set(first(history).active) == Set([DIRECT, RELAY_OUT])
        @test all(sample -> RELAY_OUT in sample.active, history)
        if spec.kind == :busy_relay
            @test schedule_changes(history) == 1
            @test all(sample -> Set(sample.active) == Set([RELAY_IN, RELAY_OUT]), history[2:end])
        else
            @test all(sample -> Set(sample.active) == Set([DIRECT, RELAY_OUT]), history)
            @test schedule_changes(history) == 0
        end
        compare("delivered", 1e6 * spec.duration + 5e5 * direct_duration, final.sc[3].data_received_bits)
        compare("source_sent", 5e5 * direct_duration + 1e6 * relay_in_duration, final.sc[1].data_sent_bits)
        compare("relay_produced", 1e6 * spec.duration, final.sc[2].data_produced_bits)
        compare("relay_overflow", 1e6 * relay_in_duration, final.sc[2].data_lost_bits)
        compare("relay_final_stored", spec.initial[2], final.sc[2].data_stored_bits)
    end
    return DataFrame(rows)
end

function run_data_case(spec, root)
    directory = joinpath(root, spec.name)
    mkpath(directory)
    args = data_configuration(spec, directory)
    initial = SE.build_initial_conditions(args)
    candidates = data_candidates(spec, initial, 0.0, Set{LinkKey}())
    CSV.write(joinpath(directory, "spacecraft_before_run.csv"), DataFrame(
        spacecraft=collect(eachindex(spec.capacity)), capacity_bits=spec.capacity,
        initial_bits=spec.initial, production_bps=spec.production, terminal_capacity=spec.terminals,
        arc_m=[(i - 1) * 10e3 for i in eachindex(spec.capacity)],
        destination=[i == spec.sink for i in eachindex(spec.capacity)]))
    CSV.write(joinpath(directory, "candidates_before_run.csv"), candidates)
    CSV.write(joinpath(directory, "matchings_before_run.csv"), all_data_matchings(candidates))
    write_toml(joinpath(directory, "settings.toml"), Dict(
        "duration_s" => spec.duration, "horizon_s" => spec.horizon, "dtmax_s" => spec.dtmax,
        "reltol" => 1e-10, "abstol" => 1e-8, "solver" => "Tsit5",
        "policy" => "DataBacklogPolicy", "record_available" => true,
        "occultation" => true, "attitude" => false, "virtual_terminals" => true,
        "transmit_power_w" => 0.0, "budget" => "FixedRateBudget",
        "destination_spacecraft" => spec.sink, "store_threshold_relative" => 1e-9,
        "counter_atol_bits" => BIT_ATOL, "score_rtol" => SCORE_RTOL,
        "score_atol_bps" => DATA_SCORE_ATOL))
    result = nothing
    @testset "$(spec.name)" begin
        println("Running data fixture ", spec.name)
        solution = run_checked(args)
        run_args = solution.prob.p.args
        history = run_args.interlink_model.history
        states = history_states(solution)
        allowed = Set(edge.key for edge in spec.edges if edge.allowed)
        model = deepcopy(run_args.interlink_model)
        spacecraft = deepcopy(run_args.dynamics_model.spacecraft)
        scheduler = SchedulingPolicyModel(DataBacklogPolicy(horizon_s=spec.horizon))
        audits = NamedTuple[]
        candidate_tables = DataFrame[]
        previous = Set{LinkKey}()
        maximum_score_error = 0.0
        for (sample, state) in zip(history, states)
            table = data_candidates(spec, state, sample.time, previous)
            matchings = all_data_matchings(table)
            optimum = maximum(matchings[matchings.feasible, :total_score_bps])
            check_matching(sample, allowed, run_args.dynamics_model.spacecraft)
            expected_available = Set(row_key(row) for row in eachrow(table) if row.available)
            @test Set(sample.available) == expected_available
            @test isapprox(sum(sample.scores), optimum; rtol=SCORE_RTOL, atol=DATA_SCORE_ATOL)
            replayed = replay_scores!(model, scheduler, spacecraft, state, previous, sample.time)
            replay_column = Union{Missing, Float64}[]
            selected_scores = Dict(zip(sample.active, sample.scores))
            selected = Set(sample.active)
            for row in eachrow(table)
                key = row_key(row)
                if row.allowed
                    actual = replayed[model.graph.edge_index_lookup[key]]
                    @test isapprox(actual, row.score_bps; rtol=SCORE_RTOL, atol=DATA_SCORE_ATOL)
                    maximum_score_error = max(maximum_score_error, abs(actual - row.score_bps))
                    push!(replay_column, actual)
                else
                    @test !(key in selected) && !(key in sample.available)
                    push!(replay_column, missing)
                end
                if key in selected
                    @test isapprox(selected_scores[key], row.score_bps; rtol=SCORE_RTOL, atol=DATA_SCORE_ATOL)
                end
            end
            table.replayed_score_bps = replay_column
            table.selected = [row_key(row) in selected for row in eachrow(table)]
            push!(candidate_tables, table)
            push!(audits, (; time_s=sample.time, previous_matching=matching_label(previous),
                selected_matching=matching_label(selected), allowed_count=length(allowed),
                feasible_count=length(expected_available), selected_count=length(selected),
                selected_score_bps=sum(sample.scores), exhaustive_optimum_bps=optimum,
                changed=sample.time > 0 && selected != previous,
                delivered_bits=state.sc[spec.sink].data_received_bits))
            previous = selected
        end
        candidate_audit = vcat(candidate_tables...)
        CSV.write(joinpath(directory, "candidate_scores.csv"), candidate_audit)
        CSV.write(joinpath(directory, "accepted_steps.csv"), DataFrame(audits))
        save_schedule(solution, directory)
        nodes, network = storage_accounting(spec, history, states)
        CSV.write(joinpath(directory, "node_accounting.csv"), nodes)
        CSV.write(joinpath(directory, "network_accounting.csv"), network)
        hands = hand_comparisons(spec, solution)
        CSV.write(joinpath(directory, "hand_comparisons.csv"), hands)
        if spec.kind == :production_competition
            producer = candidate_audit[candidate_audit.tx .== 1, :]
            @test all(>(0), producer.score_bps)
            @test any(.!producer.selected)
            @test any(producer.selected)
        end
        last_network = last(network)
        summary = (; case=spec.name, duration_s=spec.duration, dtmax_s=spec.dtmax,
            samples=length(history), switches=schedule_changes(history),
            delivered_bits=last_network.delivered_bits, produced_bits=last_network.produced_bits,
            lost_bits=last_network.lost_bits, stored_bits=last_network.total_stored_bits,
            maximum_node_balance_error_bits=maximum(abs, nodes.balance_residual_bits),
            maximum_network_balance_error_bits=maximum(abs, network.network_residual_bits),
            maximum_transfer_error_bits=maximum(abs, network.transfer_residual_bits),
            maximum_score_error_bps=maximum_score_error)
        result = (; summary, hands)
    end
    return result
end

function run_backlog_suite()
    directory = joinpath(OUTPUT, "backlog")
    mkpath(directory)
    outflow_probes(directory)
    results = map([
        data_spec(:producer_control), data_spec(:production_competition),
        data_spec(:busy_relay), data_spec(:direct_route_control),
        data_spec(:production_competition; refined=true), data_spec(:busy_relay; refined=true),
    ]) do spec
        run_data_case(spec, directory)
    end
    summary = DataFrame([result.summary for result in results])
    CSV.write(joinpath(directory, "summary.csv"), summary)
    CSV.write(joinpath(directory, "hand_comparisons.csv"), vcat([result.hands for result in results]...))
    row(name) = only(eachrow(summary[summary.case .== name, :]))
    busy, direct = row("busy_relay"), row("direct_route_control")
    @test direct.delivered_bits > busy.delivered_bits + 1e6
    @test busy.lost_bits > direct.lost_bits + 1e6
    refinement = DataFrame(map(("production_competition", "busy_relay")) do name
        coarse, fine = row(name), row(name * "_refined")
        (; case=name, coarse_dtmax_s=coarse.dtmax_s, refined_dtmax_s=fine.dtmax_s,
            coarse_switches=coarse.switches, refined_switches=fine.switches,
            coarse_delivered_bits=coarse.delivered_bits, refined_delivered_bits=fine.delivered_bits,
            delivered_difference_bits=fine.delivered_bits - coarse.delivered_bits,
            loss_difference_bits=fine.lost_bits - coarse.lost_bits)
    end)
    CSV.write(joinpath(directory, "step_refinement.csv"), refinement)
    return (; summary, delivery_advantage_bits=direct.delivered_bits - busy.delivered_bits,
        loss_reduction_bits=busy.lost_bits - direct.lost_bits)
end
