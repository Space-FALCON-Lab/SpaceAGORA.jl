module BenchmarkSummaryColumnsTests
using Test
using CSV
using DataFrames
using Statistics

const RUNTIME = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies", "performance_runtime_analysis"))

# Select the real formatting helper without the catalog's process-policy state.
function load_definition(path::String, wanted::Symbol)
    for expr in Meta.parseall(read(path, String)).args
        expr isa Expr || continue
        node = expr.head == :macrocall ? expr.args[end] : expr
        node isa Expr && node.head == :function || continue
        signature = node.args[1]
        signature = signature.head == :(::) ? signature.args[1] : signature
        if signature.args[1] == wanted
            Core.eval(@__MODULE__, expr)
            return
        end
    end
    error("Required summary helper $wanted not found in $path")
end
load_definition(joinpath(RUNTIME, "case_catalog", "profile_types_and_helpers.jl"), :_safe_unique_join)

# The source declares an unrelated BenchmarkCase method. None of its simulation
# methods are invoked; all aggregation and statistics below are production code.
struct BenchmarkCase end
include(joinpath(RUNTIME, "reporting", "summaries.jl"))

function captured_summary(f, raw)
    before = deepcopy(raw)
    result = nothing
    caught = try
        result = f(raw)
        nothing
    catch err
        err isa InterruptException && rethrow()
        err
    end
    @test isequal(raw, before)
    return (result=result, error=caught)
end

function orbit_input(; modern=true, legacy=false, conflict=false)
    raw = DataFrame(category=fill("fixture", 2), scenario=fill("orbit_fixture", 2),
        description=fill("Synthetic summary records", 2),
        orbital_period_s=fill(100.0, 2), dt_max_orbit_s=fill(1.0, 2),
        outer_threads_safe=fill(true, 2), solve_success=fill(true, 2),
        mission_time_s=[200.0, 400.0], total_time_s=[10.0, 40.0],
        solve_time_s=[8.0, 32.0], total_bytes_mb=[2.0, 4.0],
        sim_seconds_per_wall_second=[20.0, 10.0])
    modern && (raw[!, :mission_time_multiplier] = [2, 4])
    legacy && (raw[!, :orbit_count] = conflict ? [99, 99] : [2, 4])
    return raw
end

@testset "orbit modern precedence and legacy schema" begin
    for (label, modern, legacy, conflict) in (
        ("modern only", true, false, false),
        ("legacy only", false, true, false),
        ("equal aliases", true, true, false),
        ("conflicting aliases", true, true, true),
    )
        @testset "$label" begin
            observed = captured_summary(summarize_per_orbit_results,
                orbit_input(; modern=modern, legacy=legacy, conflict=conflict))
            @test observed.error === nothing
            if observed.error === nothing
                summary = observed.result
                @test summary.mission_time_multiplier == [2, 4]
                @test summary.samples_total == [1, 1]
                @test summary.samples_success == [1, 1]
                @test summary.samples_failed == [0, 0]
                @test summary.total_time_mean_s == [10.0, 40.0]
                @test summary.time_per_baseline_period_mean_s == [5.0, 10.0]
                @test length(summary.baseline_periods_per_wall_second_mean) == 2 &&
                    summary.baseline_periods_per_wall_second_mean ≈ [0.2, 0.1]
                @test summary.time_per_orbit_mean_s == summary.time_per_baseline_period_mean_s
                @test summary.orbits_per_wall_second_mean == summary.baseline_periods_per_wall_second_mean
                # Modern grouping intentionally emits the modern key. Legacy-only
                # input retains its key and gains the human-facing modern alias.
                @test (:orbit_count in propertynames(summary)) == !modern
                if !modern
                    @test summary.orbit_count == [2, 4]
                end
                mktempdir() do dir
                    path = joinpath(dir, "orbit_summary.csv")
                    CSV.write(path, summary)
                    restored = CSV.read(path, DataFrame)
                    @test propertynames(restored) == propertynames(summary)
                    @test isequal(restored, summary)
                end
            end
        end
    end
end

@testset "orbit supplied attempts and all-failed groups" begin
    raw = orbit_input()
    raw[!, :mission_time_multiplier] = [2, 2]
    raw[!, :mission_time_s] = [200.0, 200.0]
    raw[!, :attempt] = [1, 2]
    raw[!, :solve_success] = [false, true]
    raw[!, :total_time_s] = [999.0, 8.0]
    raw[!, :scenario] = ["a_mixed", "a_mixed"]
    failed = copy(raw[1:1, :])
    failed.scenario .= "z_failed"
    failed.mission_time_multiplier .= 4
    # Reverse input group order to check production ordering as well as counts.
    observed = captured_summary(summarize_per_orbit_results, vcat(failed, raw))
    @test observed.error === nothing
    if observed.error === nothing
        summary = observed.result
        @test summary.scenario == ["a_mixed", "z_failed"]
        @test summary.mission_time_multiplier == [2, 4]
        @test summary.samples_total == [2, 1]
        @test summary.samples_success == [1, 0]
        @test summary.samples_failed == [1, 1]
        @test summary.success_rate == [0.5, 0.0]
        @test summary.total_time_mean_s[1] == 8.0
        @test summary.time_per_baseline_period_mean_s[1] == 4.0
        @test summary.baseline_periods_per_wall_second_mean[1] == 0.25
        for col in (:samples, :total_time_mean_s, :total_time_p90_s,
                    :total_time_ci95_low_s, :solve_time_mean_s,
                    :time_per_orbit_mean_s, :orbits_per_wall_second_mean,
                    :time_per_baseline_period_mean_s, :baseline_periods_per_wall_second_mean)
            @test ismissing(summary[2, col])
        end
    end
end

@testset "orbit entirely failed input retains missing statistics" begin
    raw = orbit_input()
    raw.solve_success .= false
    observed = captured_summary(summarize_per_orbit_results, raw)
    @test observed.error === nothing
    if observed.error === nothing
        summary = observed.result
        @test summary.mission_time_multiplier == [2, 4]
        @test summary.samples_total == [1, 1]
        @test summary.samples_success == [0, 0]
        @test summary.samples_failed == [1, 1]
        @test summary.success_rate == [0.0, 0.0]
        for col in (:samples, :total_time_mean_s, :total_time_p90_s,
                    :total_time_ci95_low_s, :total_time_ci95_high_s,
                    :total_time_sem_s, :total_time_cv_pct, :solve_time_mean_s,
                    :total_bytes_mean_mb, :sim_seconds_per_wall_second_mean,
                    :time_per_orbit_mean_s, :orbits_per_wall_second_mean,
                    :time_per_baseline_period_mean_s, :baseline_periods_per_wall_second_mean)
            @test all(ismissing, summary[!, col])
        end
    end
end

const BUCKETS = ["gram_point_to_point", "gram_surrogate",
    "gram_static_grid_or_cached_surrogate", "non_gram"]
const ROUTES = ["process", "threads_or_auto", "threads_or_auto", "threads_or_auto"]

@testset "density coverage and statistics use supplied records" begin
    raw = DataFrame(density_backend_bucket=[BUCKETS[1], BUCKETS[1], BUCKETS[1], BUCKETS[2], BUCKETS[4]],
        solve_success=[true, false, true, false, true],
        total_time_s=[10.0, 999.0, 20.0, 500.0, 8.0],
        sim_seconds_per_wall_second=[4.0, 999.0, 6.0, 500.0, 2.0],
        density_family=["z_family", "a_family", "z_family", "surrogate", "constant"],
        outer_route=["process", "threads", "process", "none", "none"],
        attempt=[1, 1, 2, 1, 1])
    observed = captured_summary(summarize_density_backend_breakdown, raw)
    @test observed.error === nothing
    if observed.error === nothing
        summary = observed.result
        @test summary.density_backend_bucket == BUCKETS
        @test summary.recommended_route == ROUTES
        @test summary.covered == [true, true, false, true]
        @test summary.samples_total == [3, 1, 0, 1]
        @test summary.samples_success == [2, 0, 0, 1]
        @test summary.samples_failed == [1, 1, 0, 0]
        @test summary.success_rate_pct[1] ≈ 200 / 3
        @test summary.success_rate_pct[2] == 0.0
        @test ismissing(summary.success_rate_pct[3])
        @test summary.success_rate_pct[4] == 100.0
        @test summary.total_time_mean_s[1] == 15.0
        @test summary.total_time_p90_s[1] ≈ 19.0
        @test summary.sim_seconds_per_wall_second_mean[1] == 5.0
        @test summary.density_families[1] == "a_family,z_family"
        @test summary.outer_routes[1] == "process,threads"
        for row in (2, 3), col in (:total_time_mean_s, :total_time_p90_s,
                                   :sim_seconds_per_wall_second_mean)
            @test ismissing(summary[row, col])
        end
        @test ismissing(summary.density_families[3])
        @test ismissing(summary.outer_routes[3])
        @test summary.total_time_mean_s[4] == 8.0
        mktempdir() do dir
            path = joinpath(dir, "density_summary.csv")
            CSV.write(path, summary)
            @test isequal(CSV.read(path, DataFrame), summary)
        end
    end
end

@testset "density absent-schema and empty fallbacks" begin
    for raw in (DataFrame(), DataFrame(density_backend_bucket=String[]),
                DataFrame(unrelated=["not a density schema"]))
        observed = captured_summary(summarize_density_backend_breakdown, raw)
        @test observed.error === nothing
        if observed.error === nothing
            summary = observed.result
            @test summary.density_backend_bucket == BUCKETS
            @test summary.recommended_route == ROUTES
            @test all(.!summary.covered)
            @test all(==(0), summary.samples_total)
            @test all(==(0), summary.samples_success)
            @test all(==(0), summary.samples_failed)
            @test summary.density_families == fill("", 4)
            @test summary.outer_routes == fill("", 4)
            for col in (:success_rate_pct, :total_time_mean_s, :total_time_p90_s,
                        :sim_seconds_per_wall_second_mean)
                @test all(ismissing, summary[!, col])
            end
        end
    end
end

@testset "summary fixture has no native or simulation imports" begin
    for name in (:SpaceAGORA, :SimulationModel, :GRAMSuite, :SPICE, :Plots)
        @test !isdefined(@__MODULE__, name)
    end
end
end # module
