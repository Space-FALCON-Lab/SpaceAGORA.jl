module PPBSerialBaselineTests

using Test
using DataFrames

const STUDIES = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies"))
include(joinpath(STUDIES, "parallelization_performance", "cli.jl"))
include(joinpath(STUDIES, "parallelization_performance", "modes.jl"))
include(joinpath(STUDIES, "parallelization_performance", "cases.jl"))
include(joinpath(STUDIES, "paper_parallelization_benchmarks", "cli.jl"))
include(joinpath(STUDIES, "paper_parallelization_benchmarks", "reporting.jl"))

row(phase, mode, threads, workers, wall; case = "mc_heavy", samples = 64) =
    (phase_id = phase, case = case, mode = mode, thread_count = threads,
     process_workers = workers, mc_samples = samples, wall_time_s = Float64(wall),
     throughput_samples_per_s = samples / wall, success = true)

@testset "P6ps keeps each controller's serial baseline without duplicating rows" begin
    raw = DataFrame([row("P6ps", mode, t, 32, mode == "serial" ? base : 2.0)
        for (t, base) in [(32, 20.0), (1, 10.0)]
        for mode in ["serial", "outer_threads", "outer_process", "predictive"]])
    agg = _ppb_aggregate(raw)
    @test nrow(agg) == 8
    @test all(agg.serial_median_s[agg.thread_count .== 32] .== 20.0)
    @test all(agg.serial_median_s[agg.thread_count .== 1] .== 10.0)
    @test all(agg.speedup[(agg.thread_count .== 32) .& (agg.mode .!= "serial")] .== 10.0)
    @test all(agg.speedup[agg.mode .== "serial"] .== 1.0)
end

@testset "ordinary thread ladders retain their single serial baseline" begin
    raw = DataFrame([row("P1", "serial", 1, 32, 16),
                     row("P1", "outer_threads", 1, 32, 16),
                     row("P1", "outer_threads", 32, 32, 2)])
    agg = _ppb_aggregate(raw)
    @test nrow(agg) == 3
    @test all(agg.serial_median_s .== 16.0)
    @test only(agg.speedup[agg.thread_count .== 32]) == 8.0
end

@testset "X5 splits share the full-budget serial measurement explicitly" begin
    for phase in ("X5a", "X5b")
        raw = DataFrame([row(phase, "serial", 32, 32, 24),
                         row(phase, "predictive", 32, 32, 3),
                         row(phase, "outer_threads", 4, 8, 6),
                         row(phase, "outer_process", 16, 2, 12)])
        agg = _ppb_aggregate(raw)
        @test nrow(agg) == 4
        @test all(agg.serial_median_s .== 24.0)
        @test only(agg.speedup[agg.process_workers .== 8]) == 4.0
        @test only(agg.speedup[agg.process_workers .== 2]) == 2.0
        # Unrelated phases must not borrow a serial baseline across workers.
        raw.phase_id .= "P5f"
        ordinary = _ppb_aggregate(raw)
        @test all(ismissing, ordinary.serial_median_s[ordinary.process_workers .!= 32])
    end
end

@testset "ambiguous and different-condition baselines are never combined" begin
    agg = DataFrame(phase_id = fill("X5a", 4), case = fill("mc_heavy", 4),
        mc_samples = fill(64, 4), mode = ["serial", "predictive", "serial", "predictive"],
        process_workers = [32, 8, 32, 8], thread_count = [32, 4, 32, 4],
        budget_condition = ["uncapped", "uncapped", "equal_core", "equal_core"],
        core_budget = [32, 32, 32, 4], wall_time_median_s = [24.0, 6.0, 40.0, 8.0])
    _ppb_add_serial_baseline!(agg)
    @test agg.serial_median_s == [24.0, 24.0, 40.0, 40.0]
    @test agg.core_budget == [32, 32, 32, 4]
    select!(agg, Not(:budget_condition))
    _ppb_add_serial_baseline!(agg)
    @test nrow(agg) == 4
    @test all(ismissing, agg.serial_median_s)
end

end # module
