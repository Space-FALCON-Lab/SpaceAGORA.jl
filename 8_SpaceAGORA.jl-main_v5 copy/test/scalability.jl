using Test
using Random
using LinearAlgebra
using SpaceAGORA.SimulationModel

function scored_graph(n, edges; penalty=0.0)
    spacecraft = [(n_terminal=1, battery_energy_index=100.0, tempurature_index=100.0) for _ in 1:n]
    model = InterLinkModel(active_link_penalty=penalty)
    for (a, b, score) in edges
        key = register_candidate!(model, spacecraft, (a, 1), (b, 1))
        model.linkgraph[key].state.available = true
        model.linkgraph[key].state.score = score
    end
    return spacecraft, model
end

function exhaustive_score(model)
    candidates = sort!(collect(keys(model.linkgraph)))
    best = 0.0
    for mask in 0:((1 << length(candidates)) - 1)
        occupied = Set{Tuple{Int, Int}}()
        score = 0.0
        valid = true
        for (index, key) in enumerate(candidates)
            iszero(mask & (1 << (index - 1))) && continue
            connection = model.linkgraph[key]
            if !connection.state.available || any(endpoint in occupied for endpoint in key)
                valid = false
                break
            end
            union!(occupied, key)
            score += connection.state.score - model.active_link_penalty
        end
        valid && (best = max(best, score))
    end
    return best
end

@testset "Array-backed interlink scalability" begin
    @testset "Bulk availability and signed transitions" begin
        spacecraft, model = scored_graph(3, [(1, 2, 1.0)])
        key = only(keys(model.linkgraph))
        state = (sc=[(pos=[Float64(i), 0.0, 0.0],) for i in 1:3],)
        update_availability!(model, spacecraft, state)
        availability = model.available
        transitions = model.availability_transition
        @test availability == [true]
        @test transitions == [0] # The public connection state was already true.
        state.sc[2].pos[1] = 1e7
        update_availability!(model, spacecraft, state)
        @test availability == [false]
        @test transitions == [-1]
        update_availability!(model, spacecraft, state)
        @test transitions == [0]
        state.sc[2].pos[1] = 2.0
        update_availability!(model, spacecraft, state)
        @test transitions == [1]
        @test model.available === availability
        @test model.availability_transition === transitions
        update_availability!(model, spacecraft, state; is_active=[true, false, true])
        @test availability == [false]
        @test transitions == [-1]
        @test update_availability!(model, key, spacecraft, state)
        @test availability == [true]
        @test transitions == [1]
        spacecraft[1] = (n_terminal=1, battery_energy_index=0.0, tempurature_index=100.0)
        update_availability!(model, spacecraft, state)
        @test !only(availability)
        spacecraft[1] = (n_terminal=1, battery_energy_index=100.0, tempurature_index=100.0)
        model.eligibility = (key, spacecraft, state) -> false
        update_availability!(model, spacecraft, state)
        @test !only(availability)
        model.eligibility = nothing
        register_candidate!(model, spacecraft, (2, 1), (3, 1))
        update_availability!(model, spacecraft, state)
        @test model.available == [true, true]
        @test length(model.availability_transition) == 2
        update_availability!(model, spacecraft, state)
        @test (@allocated update_availability!(model, spacecraft, state)) <= 1024
    end

    @testset "Exact tree and cyclic matching" begin
        rng = MersenneTwister(520)
        for trial in 1:60
            edges = [(rand(rng, 1:(i - 1)), i, Float64(rand(rng, -2:10))) for i in 2:8]
            if iseven(trial)
                push!(edges, (1, 8, 4.0))
                unique!(edge -> minmax(edge[1], edge[2]), edges)
            end
            _, model = scored_graph(9, edges; penalty=1.0)
            for connection in values(model.linkgraph)
                connection.state.available = rand(rng) > 0.2
            end
            selected = select_interlinks!(model)
            @test sum((model.linkgraph[k].state.score - 1.0 for k in selected); init=0.0) == exhaustive_score(model)
            @test selected == select_interlinks!(model)
            @test length(Set(endpoint for key in selected for endpoint in key)) == 2length(selected)
        end
        edges = [(1, 2, 1.0), (2, 3, 1.0), (1, 3, 1.0),
            (4, 5, 2.0), (5, 6, 3.0), (4, 6, 3.0), (7, 8, 1.0)]
        _, model = scored_graph(8, edges)
        _, reversed_model = scored_graph(8, reverse(edges))
        @test select_interlinks!(model) == select_interlinks!(reversed_model)
        @test sum(model.linkgraph[k].state.score for k in select_interlinks!(model)) == exhaustive_score(model)
    end

    @testset "4000-terminal trees and cache invalidation" begin
        n = 4000
        for edges in (
            [(1, i, Float64(i)) for i in 2:n],
            [(i - 1, i, 1.0) for i in 2:n],
            [(i, i + 1, 1.0) for i in 1:2:(n - 1)])
            _, model = scored_graph(n, edges)
            selected = select_interlinks!(model)
            @test length(selected) == (edges[1][3] == 2.0 ? 1 : n ÷ 2)
            workspace = model.matching
            @test select_interlinks!(model) == selected
            @test model.matching === workspace
            @test length(Set(endpoint for key in selected for endpoint in key)) == 2length(selected)
            if length(selected) == 1
                @test (@allocated select_interlinks!(model)) <= 1024
            end
        end
        spacecraft, model = scored_graph(4, [(1, 2, 1.0)])
        @test length(select_interlinks!(model)) == 1
        key = register_candidate!(model, spacecraft, (3, 1), (4, 1))
        model.linkgraph[key].state.available = true
        model.linkgraph[key].state.score = 2.0
        @test length(select_interlinks!(model)) == 2
        for connection in values(model.linkgraph)
            connection.state.available = false
        end
        @test isempty(select_interlinks!(model))
        @test all(!c.state.active for c in values(model.linkgraph))
    end

    @testset "Incident-edge force equivalence" begin
        _, model = scored_graph(4, [(1, 2, 1.0), (2, 3, 3.0)])
        select_interlinks!(model)
        state = (sc=[(pos=[7e6, 1000.0 * i, 0.0],) for i in 1:4],)
        for satellite in 1:4
            reference = zeros(3)
            for (key, connection) in model.linkgraph
                connection.state.active || continue
                satellite in (key[1][1], key[2][1]) || continue
                partner = key[1][1] == satellite ? key[2][1] : key[1][1]
                reference .+= force_on_endpoint(state.sc[satellite], state.sc[partner], connection.parameters)
            end
            @test laser_force_on_spacecraft(model, state, satellite) == reference
        end
    end
end
