module BenchmarkSolveSnapshotTests
using Test
using SpaceAGORA
using TOML

const STUDIES = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies"))
for file in ("cli.jl", "modes.jl", "cases.jl", "reporting.jl", "execution.jl")
    include(joinpath(STUDIES, "parallelization_performance", file))
end

@testset "per-solve snapshot record is TOML-serializable" begin
    cache = Dict("sig" => Dict{String, Any}("mode" => "heuristic", "sweep_ns" => 2.0e8,
                                            "honoured_ns" => 1.0e9, "timings" => [Dict("ns" => 1.0)]))
    store = ppc_store_copy(cache)
    @test !haskey(store["sig"], "timings")
    @test store["sig"]["n_timings"] == 1
    policy = (rhs_plan_source=:sweep, rhs_plan_sweeps_total=1, decisions_total=7, nothing_field=nothing)
    rec = ppc_solve_record(3, false, 2.5, store, store, policy)
    back = TOML.parse(sprint(TOML.print, rec))
    @test back["seq"] == 3 && back["wall_time_s"] == 2.5
    @test back["policy"]["rhs_plan_source"] == "sweep"
    @test back["store_after"]["sig"]["sweep_ns"] == 2.0e8
    @test isempty(ppc_solve_record(1, true, 1.0, Dict(), Dict(), nothing)["policy"])
end
end
