using Printf

function original_select_interlinks!(model, policy)
    candidates = sort!([key for (connection, key) in enumerate([edge_key(model.graph, e) for e in eachindex(model.graph.edges)])
        if model.graph.available[connection] && policy.matching.score[connection] > 0.0])
    endpoints = sort!(unique([endpoint for key in candidates for endpoint in key]))
    bits = Dict(endpoint => big(1) << (index - 1) for (index, endpoint) in enumerate(endpoints))
    neighbors = Dict(endpoint => InterLinkKey[] for endpoint in endpoints)
    for key in candidates, endpoint in key
        push!(neighbors[endpoint], key)
    end
    memo = Dict{BigInt, Tuple{Float64, Vector{InterLinkKey}}}()

    function best_matching(mask::BigInt)
        iszero(mask) && return (0.0, InterLinkKey[])
        haskey(memo, mask) && return memo[mask]
        endpoint = endpoints[findfirst(endpoint -> !iszero(mask & bits[endpoint]), endpoints)]
        remaining = mask & ~bits[endpoint]
        best_score, best_keys = best_matching(remaining)
        for key in neighbors[endpoint]
            partner = key[1] == endpoint ? key[2] : key[1]
            iszero(remaining & bits[partner]) && continue
            score, selected = best_matching(remaining & ~bits[partner])
            score += policy.matching.score[model.graph.edge_index_lookup[key]]
            if score > best_score
                best_score, best_keys = score, vcat(InterLinkKey[key], selected)
            end
        end
        memo[mask] = (best_score, best_keys)
        return memo[mask]
    end

    _, selected = best_matching((big(1) << length(endpoints)) - 1)
    selected = sort!(copy(selected))
    selected_set = Set(selected)
    for (connection, key) in enumerate([edge_key(model.graph, e) for e in eachindex(model.graph.edges)])
        model.graph.active[connection] = key in selected_set
    end
    return selected
end

function dense_fixture(count, policy)
    spacecraft = [SpacecraftModel() for _ in 1:count]
    model = ConstellationInterLinkModel()
    scores = Float64[]
    for first_satellite in 1:(count - 1), second_satellite in (first_satellite + 1):count
        key = register_allowed_link!(model, spacecraft, (first_satellite, 1), (second_satellite, 1))
        model.graph.available[model.graph.edge_index_lookup[key]] = true
        push!(scores, isodd(first_satellite) && second_satellite == first_satellite + 1 ? 2.0 : 1.0)
    end
    matching_workspace!(model, policy).score .= scores
    expected = [((satellite, 1), (satellite + 1, 1)) for satellite in 1:2:(count - 1)]
    return model, expected
end

function benchmark_worker(method, count, limit)
    policy = SchedulingPolicyModel(:gve_sma)
    selector = method == "original" ? original_select_interlinks! : select_interlinks!
    warm_model, warm_expected = dense_fixture(8, policy)
    @assert selector(warm_model, policy) == warm_expected
    model, expected = dense_fixture(count, policy)
    timings = Float64[]
    for trial in 1:5
        GC.gc()
        ccall(:alarm, Cuint, (Cuint,), limit)
        elapsed = @elapsed selected = selector(model, policy)
        ccall(:alarm, Cuint, (Cuint,), 0)
        @assert selected == expected
        @assert all(model.graph.active[connection] == (key in selected) for (connection, key) in enumerate([edge_key(model.graph, e) for e in eachindex(model.graph.edges)]))
        push!(timings, elapsed)
        @printf("%s terminals=%d trial=%d ms=%.3f\n", method, count, trial, elapsed * 1000)
        flush(stdout)
    end
    @printf("RESULT %s terminals=%d candidates=%d median_ms=%.3f\n",
        method, count, length(model.graph.edges), sort!(timings)[3] * 1000)
end

"""
Run with `julia --project=8_SpaceAGORA.jl-main_v5 --startup-file=no 8_SpaceAGORA.jl-main_v5/test/benchmark_matching.jl`.

Compare the original recursive selector with weighted Blossom on complete graphs:
adjacent odd/even terminal pairs have weight 2; every other edge has weight 1.
Five warmed calls include graph conversion and active-flag updates, but exclude
fixture construction, loading, compilation warm-up, and explicit pre-call GC.
Each child measurement has a Linux SIGALRM deadline of 10 seconds; a timeout
reports a lower bound, not a measured completion time. The production selector
is never replaced. This opt-in benchmark is not part of runtests.jl.
"""
function benchmark_main()
    Sys.islinux() || error("This benchmark uses Linux SIGALRM for process-isolated time limits.")
    project = dirname(@__DIR__)
    limit = 10
    for count in (20, 100, 200), method in ("blossom", "original")
        command = `$(Base.julia_cmd()) --project=$project --startup-file=no $(@__FILE__) --worker $method $count $limit`
        process = run(ignorestatus(command))
        if process.termsignal == 14
            println("TIMEOUT $method terminals=$count: one selector call exceeded $limit seconds.")
        elseif !success(process)
            error("Benchmark worker failed: $method, $count terminals, exit=$(process.exitcode), signal=$(process.termsignal)")
        end
    end
end

if !isempty(ARGS) && first(ARGS) == "--worker"
    include(joinpath(@__DIR__, "..", "II_examples", "gve_sma_interlinks.jl"))
    edge_key(graph::InterLinkGraph, edge::Int) =
        (graph.vertices[graph.edges[edge][1]], graph.vertices[graph.edges[edge][2]])
    using SpaceAGORA.SimulationModel.SchedulingPolicies: matching_workspace!
    benchmark_worker(ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4]))
else
    benchmark_main()
end