module RuntimeAnalysisCopyOverheadGate
using Test
using DataFrames
using Dates
using Statistics

const RUNTIME = normpath(joinpath(@__DIR__, "..", "..", "benchmarks", "studies", "performance_runtime_analysis"))

# Load actual small definitions without the catalog's process-policy state or
# main.jl's unconditional GRAMSuite import. No simulation or plotting is run.
function load_definition(file::String, wanted::Symbol)
    for expr in Meta.parseall(read(file, String)).args
        expr isa Expr || continue
        node = expr.head == :macrocall ? expr.args[end] : expr
        node isa Expr || continue
        name = if node.head == :struct
            node.args[2]
        elseif node.head == :function
            signature = node.args[1]
            signature = signature.head == :(::) ? signature.args[1] : signature
            signature.args[1]
        else
            nothing
        end
        if name == wanted
            Core.eval(@__MODULE__, expr)
            return
        end
    end
    error("Required reporting definition $wanted not found in $file")
end
const HELPERS = joinpath(RUNTIME, "case_catalog", "profile_types_and_helpers.jl")
load_definition(HELPERS, :ProfileSpec)
load_definition(HELPERS, :_safe_unique_join)
for name in (:_fmt, :_scenario_metric)
    load_definition(joinpath(RUNTIME, "reporting", "plots_and_reporting.jl"), name)
end

# The summaries file declares an unused BenchmarkCase method. Its simulation
# type is deliberately absent here. Host/settings providers are fixed context
# so recorded provenance cannot accidentally equal the fallback machine.
struct BenchmarkCase end
const PERF_BASELINE_SCENARIO = "baseline_fixture"
const FALLBACK = (hardware_class="fallback-class", machine_label="fallback-machine",
    host_name="fallback-host", cpu_name="fallback-cpu", cpu_threads=3,
    julia_threads=2, os="fallback-os", arch="fallback-arch")
_runtime_hardware_snapshot() = FALLBACK
_perf_default_solver_mode(::String) = "fixture-default"
_perf_solver_mode_env(::String) = "fixture-effective"
_split_rollout_enabled() = false
_split_rollout_enforce() = false
_split_rollout_case_names() = String[]
_split_rollout_solver_variants() = String[]
_multirate_rollout_enabled() = false
_multirate_rollout_enforce() = false
_multirate_rollout_case_names() = String[]
_multirate_rollout_max_slowdown_ratio() = 1.0
include(joinpath(RUNTIME, "reporting", "summaries.jl"))
include(joinpath(RUNTIME, "reporting", "report_writing.jl"))
const SPEC = ProfileSpec(name="quick", repeats=1, warmup=0, max_attempts=2,
    mission_short_s=120.0, mission_long_s=600.0, montecarlo_samples=0,
    montecarlo_mission_s=120.0)

# Supplied attempt rows check the real summary/report boundary. They do not run
# the collector or establish simulation, copy performance or retry-retention claims.
function raw_row(; kwargs...)
    return merge((category="baseline", scenario=PERF_BASELINE_SCENARIO,
        description="Synthetic reporting fixture", satellites=1, orientation=false,
        mission_time_s=120.0, outer_route="serial", outer_threads_safe=true,
        density_parallel_mode="serial", control_parallel_mode="serial",
        multibody_parallel_mode="serial", dt_max_orbit_s=10.0,
        dynamic_effectors="fixture", control_effectors="none", solve_success=true,
        total_time_s=2.0, solve_time_s=1.5, copy_time_s=0.5,
        copy_compile_time_s=0.0, solve_compile_time_s=0.0,
        copy_gctime_s=0.0, solve_gctime_s=0.0,
        total_bytes_mb=12.0, copy_bytes_mb=3.0, solve_bytes_mb=9.0,
        copy_alloc_calls=100, solve_alloc_calls=300, saved_points=10,
        accepted_steps=20, rejected_steps=0, solver_mode="auto",
        solver_sequence="Tsit5", solver_fallback_used=false,
        solver_fallback_count=0, solver_fallback_trigger=missing,
        policy_threads_enabled_total=0, policy_decisions_total=0,
        policy_density_threads_enabled=0, policy_control_threads_enabled=0,
        policy_multibody_threads_enabled=0, policy_other_threads_enabled=0,
        sim_seconds_per_wall_second=60.0, satellite_sim_seconds_per_wall_second=60.0,
        nbody_spkpos_runtime_calls=0, nbody_spkpos_cache_build_calls=0,
        nbody_spkpos_total_calls=0, srp_spkpos_runtime_calls=0,
        srp_spkpos_cache_build_calls=0, srp_spkpos_total_calls=0,
        planet_pxform_runtime_calls=0, planet_pxform_cache_build_calls=0,
        planet_pxform_total_calls=0, attempt=1), (; kwargs...))
end
function render_report(raw, summary=summarize_results(raw), entry=nothing; kwargs...)
    # Writing must not reorder or alter any supplied table.
    tables = Any[raw, summary, entry, values(kwargs)...]
    before = deepcopy(tables)
    text = mktempdir() do dir
        path = joinpath(dir, "report.md")
        write_report(path, SPEC, raw, summary, DataFrame(), entry; kwargs...)
        read(path, String)
    end
    @test isequal(tables, before)
    return text
end
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
