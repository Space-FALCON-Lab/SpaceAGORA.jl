using SpaceAGORA
using SpaceAGORA.SimulationModel
using Statistics
using CSV
using DataFrames
using Printf

function interlink_benchmark_fixture(shape, count)
    spacecraft = [SpacecraftModel(id=index, n_terminal=2) for index in 1:count]
    state = (sc=[(pos=[7e6, 10.0 * index, 0.0], vel=[0.0, -7500.0, 0.0],
        mass=200.0 + index) for index in 1:count],)
    model = ConstellationInterLinkModel()
    pairs = shape == "dense" ? [(first, second) for first in 1:(count - 1) for second in (first + 1):count] :
        shape == "star" ? [(1, partner) for partner in 2:count] : [(first, first + 1) for first in 1:(count - 1)]
    for (first, second) in pairs
        register_allowed_link!(model, spacecraft, (first, 1), (second, shape == "multi_terminal" ? 2 : 1))
    end
    return spacecraft, state, model
end

function benchmark_availability!(model, spacecraft, state)
    update_availability!(model, spacecraft, state)
    return nothing
end

function benchmark_forces(model, state)
    total = 0.0
    for satellite in eachindex(state.sc)
        total += sum(abs2, laser_force_on_spacecraft(model, state, satellite))
    end
    return total
end

function measure_interlink_stage(operation; repeats=5)
    operation()
    samples = map(1:repeats) do trial
        GC.gc()
        @timed operation()
    end
    return median(sample.time for sample in samples), median(sample.bytes for sample in samples)
end

function benchmark_interlinks(label)
    table = DataFrame(label=String[], shape=String[], satellites=Int[], candidates=Int[],
        selected=Int[], stage=String[], median_s=Float64[], median_bytes=Float64[])
    mu = make_no_gram_planet(:earth).μ
    policy = SchedulingPolicyModel(:gve_sma_first_index)
    for (shape, count) in (("star", 4000), ("chain", 4000), ("multi_terminal", 400), ("dense", 100), ("dense", 200))
        spacecraft, state, model = interlink_benchmark_fixture(shape, count)
        benchmark_availability!(model, spacecraft, state)
        score_candidates!(model, policy, state, mu; candidate_force=force_on_endpoint)
        selected = select_interlinks!(model, policy)
        endpoints = [endpoint for key in selected for endpoint in key]
        @assert !isempty(selected) && length(unique(endpoints)) == length(endpoints)
        @assert all(model.graph.available[model.graph.edge_index_lookup[key]] && policy.matching.score[model.graph.edge_index_lookup[key]] > 0 for key in selected)
        for (stage, operation) in (
                ("availability", () -> benchmark_availability!(model, spacecraft, state)),
                ("scoring", () -> score_candidates!(model, policy, state, mu; candidate_force=force_on_endpoint)),
                ("matching", () -> select_interlinks!(model, policy)),
                ("forces", () -> benchmark_forces(model, state)))
            elapsed, bytes = measure_interlink_stage(operation)
            push!(table, (label, shape, count, length(model.graph.edges), length(selected), stage, elapsed, bytes))
            @printf("%s %s N=%d E=%d %s: %.3f ms, %.0f bytes\n",
                label, shape, count, length(model.graph.edges), stage, elapsed * 1000, bytes)
            flush(stdout)
        end
    end
    output = joinpath(@__DIR__, "..", "III_output", "scalability")
    mkpath(output)
    CSV.write(joinpath(output, "interlinks_$(label).csv"), table)
    return table
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark_interlinks(isempty(ARGS) ? "optimized" : only(ARGS))
end