module RuntimeAnalysisCopyOverheadGate
using Test
using DataFrames
using Dates
using Statistics

const RUNTIME = normpath(joinpath(@__DIR__, "..", "..", "benchmarks", "studies", "performance_runtime_analysis"))

include(joinpath(RUNTIME, "reporting", "runtime_reports.jl"))
using .RuntimeReports: ProfileSpec, summarize_per_orbit_results, summarize_entry_duration_results,
    summarize_density_backend_breakdown
const PERF_BASELINE_SCENARIO = "baseline_fixture"
const FALLBACK = (hardware_class="fallback-class", machine_label="fallback-machine",
    host_name="fallback-host", cpu_name="fallback-cpu", cpu_threads=3,
    julia_threads=2, os="fallback-os", arch="fallback-arch")
const SPEC = ProfileSpec(name="quick", repeats=1, warmup=0, max_attempts=2,
    mission_short_s=120.0, mission_long_s=600.0, montecarlo_samples=0,
    montecarlo_mission_s=120.0)

include(joinpath(@__DIR__, "../support", "runtime_report_fixtures.jl"))
isolation_section(text) = first(split(last(split(text, "## State Isolation Overhead")), "\n## "))

const COPY_METRICS = (:copy_time_mean_s, :solve_time_mean_s, :total_time_mean_s,
    :copy_bytes_mean_mb, :solve_bytes_mean_mb, :total_bytes_mean_mb,
    :copy_alloc_mean, :solve_alloc_mean, :copy_time_share, :copy_bytes_share,
    :copy_alloc_share)

@testset "copy-overhead successful-run means and rendered shares" begin
    # Unequal denominators distinguish ratio-of-means from mean-of-ratios.
    # The failed attempt deliberately dwarfs all accepted measurements.
    raw = DataFrame([raw_row(),
        raw_row(total_time_s=6.0, copy_time_s=3.5, solve_time_s=2.5,
            total_bytes_mb=28.0, copy_bytes_mb=21.0, solve_bytes_mb=7.0,
            copy_alloc_calls=700, solve_alloc_calls=100),
        raw_row(solve_success=false, total_time_s=1000.0, copy_time_s=999.0,
            solve_time_s=1.0, total_bytes_mb=5000.0, copy_bytes_mb=4999.0,
            solve_bytes_mb=1.0, copy_alloc_calls=9999, solve_alloc_calls=1),
        raw_row(scenario="high_copy_fixture", total_time_s=4.0,
            copy_time_s=3.0, solve_time_s=1.0, total_bytes_mb=40.0,
            copy_bytes_mb=10.0, solve_bytes_mb=30.0,
            copy_alloc_calls=1000, solve_alloc_calls=1000)])
    before = deepcopy(raw)
    summary = summarize_results(raw)
    @test isequal(raw, before)
    baseline = only(eachrow(summary[summary.scenario .== PERF_BASELINE_SCENARIO, :]))
    @test (baseline.samples_total, baseline.samples_success, baseline.samples_failed) == (3, 2, 1)
    for (metric, expected) in ((:copy_time_mean_s, 2.0), (:solve_time_mean_s, 2.0),
            (:total_time_mean_s, 4.0), (:copy_bytes_mean_mb, 12.0),
            (:solve_bytes_mean_mb, 8.0), (:total_bytes_mean_mb, 20.0),
            (:copy_alloc_mean, 400.0), (:solve_alloc_mean, 200.0),
            (:copy_time_share, 0.5), (:copy_bytes_share, 0.6),
            (:copy_alloc_share, 2 / 3))
        @test baseline[metric] == expected
    end
    @test baseline.setup_share == baseline.copy_time_share
    section = isolation_section(render_report(raw, summary))
    @test occursin("isolate_state=false", section)
    @test occursin("| baseline_fixture | baseline | 2.0 | 2.0 | 4.0 | 0.5 | 12.0 | 0.6 | 400.0 | 200.0 | 0.667 |", section)
    @test occursin("| high_copy_fixture | baseline | 3.0 | 1.0 | 4.0 | 0.75 | 10.0 | 0.25 | 1000.0 | 1000.0 | 0.5 |", section)
    @test findfirst("| high_copy_fixture |", section).start < findfirst("| baseline_fixture |", section).start

    # Input order must not affect successful-run means or rendered ranking.
    for order in ([4, 3, 2, 1], [2, 4, 1, 3])
        reordered = raw[order, :]
        reordered_before = deepcopy(reordered)
        reordered_summary = summarize_results(reordered)
        @test isequal(reordered, reordered_before)
        for scenario in (PERF_BASELINE_SCENARIO, "high_copy_fixture")
            expected = only(eachrow(summary[summary.scenario .== scenario, :]))
            actual = only(eachrow(reordered_summary[reordered_summary.scenario .== scenario, :]))
            for metric in COPY_METRICS
                @test isequal(actual[metric], expected[metric])
            end
        end
        @test isolation_section(render_report(reordered, reordered_summary)) == section
    end
end

@testset "failure-only copy metrics stay unavailable" begin
    raw = DataFrame([raw_row(solve_success=false),
        raw_row(solve_success=false, attempt=2, copy_time_s=100.0)])
    summary = summarize_results(raw)
    @test only(summary.samples_failed) == 2
    for metric in COPY_METRICS
        @test ismissing(only(summary[!, metric]))
    end
    section = isolation_section(render_report(raw, summary))
    @test occursin("No successful scenarios produced copy-overhead metrics.", section)
    @test !occursin("| baseline_fixture |", section)
end

@testset "zero, missing and empty copy metrics" begin
    # A real measured zero has a zero share when its denominator is positive.
    zero_copy = DataFrame([raw_row(copy_time_s=0.0, solve_time_s=2.0,
        copy_bytes_mb=0.0, solve_bytes_mb=12.0,
        copy_alloc_calls=0, solve_alloc_calls=400)])
    zero_copy_summary = summarize_results(zero_copy)
    for metric in (:copy_time_share, :copy_bytes_share, :copy_alloc_share)
        @test only(zero_copy_summary[!, metric]) == 0.0
    end
    @test occursin("| baseline_fixture | baseline | 0.0 | 2.0 | 2.0 | 0.0 | 0.0 | 0.0 | 0.0 | 400.0 | 0.0 |",
        isolation_section(render_report(zero_copy, zero_copy_summary)))

    # Undefined denominators and absent copy measurements must not become zero.
    for raw in (DataFrame([raw_row(total_time_s=0.0, copy_time_s=0.0,
                    solve_time_s=0.0, total_bytes_mb=0.0, copy_bytes_mb=0.0,
                    solve_bytes_mb=0.0, copy_alloc_calls=0, solve_alloc_calls=0)]),
                DataFrame([raw_row(copy_time_s=missing, copy_bytes_mb=missing,
                    copy_alloc_calls=missing)]))
        summary = summarize_results(raw)
        for metric in (:copy_time_share, :copy_bytes_share, :copy_alloc_share)
            @test ismissing(only(summary[!, metric]))
        end
        section = isolation_section(render_report(raw, summary))
        @test occursin("No successful scenarios produced copy-overhead metrics.", section)
        @test !occursin("| baseline_fixture |", section)
    end

    # Preserve the schema while supplying no attempts at all.
    empty_raw = zero_copy[1:0, :]
    empty_summary = summarize_results(empty_raw)
    @test nrow(empty_summary) == 0
    for metric in COPY_METRICS
        @test metric in propertynames(empty_summary)
    end
    section = isolation_section(render_report(empty_raw, empty_summary))
    @test occursin("No successful scenarios produced copy-overhead metrics.", section)
end

@testset "copy-overhead reporting remains independent of simulation and native loading" begin
    for name in (:SpaceAGORA, :SimulationModel, :GRAMSuite, :SPICE, :Plots)
        @test !isdefined(@__MODULE__, name)
    end
end

end # module
