module BenchmarkReportTests
using Test
using CSV
using DataFrames
using Dates
using Statistics

const RUNTIME = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies", "performance_runtime_analysis"))

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

# Synthetic supplied attempt rows, following the existing copy-overhead gate's
# schema. This checks reporting, not whether an upstream collector retains retries.
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
entry_section(text) = last(split(text, "## Entry-Duration Sweep Results"))

@testset "recorded report provenance and optional stage columns" begin
    raw = DataFrame([raw_row(), raw_row(), raw_row()])
    summary = summarize_results(raw)
    recorded = (hardware_class="recorded-class", machine_label="recorded-machine",
        host_name="recorded-host", cpu_name="recorded-cpu", cpu_threads=64,
        julia_threads=8, os="recorded-os", arch="recorded-arch")
    for (key, value) in pairs(recorded)
        later = value isa Int ? value + 1 : "later-$value"
        raw[!, key] = [missing, value, later]
    end
    stages = DataFrame(stage=["total", "run_benchmarks", "run_split_rollout_gate",
        "run_multirate_rollout_gate", "run_per_orbit", "run_entry_duration_sweep"],
        elapsed_s=[21.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    text = render_report(raw, summary; stage_timing_df=stages)
    for expected in ("Machine label: `recorded-machine`", "Hardware class: `recorded-class`",
        "Hostname: `recorded-host`", "CPU: `recorded-cpu` (`64` system threads)",
        "Julia threads in process: `8`", "OS/Arch: `recorded-os` / `recorded-arch`")
        @test occursin(expected, text)
    end
    @test !occursin("fallback-host", text)
    @test !occursin("later-recorded", text)
    @test occursin("Stage elapsed [s]: run_benchmarks=`1.0`, split_gate=`2.0`, multirate_gate=`3.0`, mission_time_sweep=`4.0`, entry_duration_sweep=`5.0`, total=`21.0`", text)
    raw[!, :machine_label] = fill(missing, nrow(raw))
    select!(raw, Not(:host_name))
    fallback_text = render_report(raw, summary)
    @test occursin("Machine label: `fallback-machine`", fallback_text)
    @test occursin("Hostname: `fallback-host`", fallback_text)
    @test occursin("CPU: `recorded-cpu`", fallback_text)
    partial = render_report(raw, summary; stage_timing_df=DataFrame(stage=["total"], elapsed_s=[7.0]))
    @test occursin("run_benchmarks=`n/a`", partial)
    @test occursin("total=`7.0`", partial)
    for table in (DataFrame(stage=["total"]), DataFrame(elapsed_s=[7.0]),
                  DataFrame(stage=String[], elapsed_s=Float64[]))
        @test !occursin("Stage elapsed [s]", render_report(raw, summary; stage_timing_df=table))
    end
end

@testset "gate counts describe supplied rows without enforcing a gate" begin
    raw = DataFrame([raw_row()])
    text = render_report(raw; split_gate_df=DataFrame(pass_all=[true, false]),
                         multirate_gate_df=DataFrame(pass_all=[false]))
    @test occursin("Split rollout guardrail: `1/2` pass.", text)
    @test occursin("Multirate rollout guardrail: `0/1` pass.", text)
    @test occursin("Any gate failure: `true`", text)
    @test occursin("Multirate any gate failure: `true`", text)
    @test occursin("Split rollout verification rows: `2` (`1` pass).", text)
    @test occursin("Multirate rollout verification rows: `1` (`0` pass).", text)
    passing = render_report(raw; split_gate_df=DataFrame(pass_all=[true]),
                            multirate_gate_df=DataFrame(pass_all=[true, true]))
    @test occursin("Any gate failure: `false`", passing)
    @test occursin("Multirate any gate failure: `false`", passing)
    for table in (nothing, DataFrame(pass_all=Bool[]), DataFrame(other=[true]))
        empty_text = render_report(raw; split_gate_df=table, multirate_gate_df=table)
        @test occursin("Split rollout guardrail: disabled or no gate rows.", empty_text)
        @test occursin("Multirate rollout guardrail: disabled or no gate rows.", empty_text)
    end
end

@testset "failure and retry accounting for supplied rows" begin
    raw = DataFrame([raw_row(solve_success=false, total_time_s=100.0, solver_fallback_count=1),
        raw_row(attempt=2, solver_fallback_count=2), raw_row(total_time_s=4.0),
        raw_row(scenario="failed_fixture", solve_success=false, total_time_s=9.0,
                solver_fallback_count=missing)])
    before = deepcopy(raw)
    summary = summarize_results(raw)
    @test isequal(raw, before)
    @test summary.scenario == [PERF_BASELINE_SCENARIO, "failed_fixture"]
    @test summary.samples_total == [3, 1]
    @test summary.samples_success == [2, 0]
    @test summary.samples_failed == [1, 1]
    @test summary.requested_runs == [2, 1]
    @test summary.retries_total == [1, 0]
    @test isequal(summary.total_time_mean_s, [3.0, missing])
    @test summary.penalized_expected_wall_time_s == [53.0, 9.0]
    text = render_report(raw, summary)
    @test occursin("Successful samples: `2/4`", text)
    @test occursin("Failed attempts: `2/4`", text)
    @test occursin("Retry overhead: `1` retries across `3` requested runs", text)
    @test occursin("Robustness-adjusted expected wall time (all attempts): `38.333 s/requested run`", text)
    @test occursin("Mean fallback count across all attempts: `0.75`", text)
    @test occursin("Solver failures detected in `2` scenario groups", text)
    @test occursin("| failed_fixture | baseline | 0/1 |", text)
    failed = DataFrame([raw_row(solve_success=false, total_time_s=8.0),
        raw_row(solve_success=false, attempt=2, total_time_s=missing)])
    failed_summary = summarize_results(failed)
    @test only(failed_summary.samples_total) == 2
    @test only(failed_summary.samples_failed) == 2
    @test ismissing(only(failed_summary.total_time_mean_s))
    @test only(failed_summary.penalized_expected_wall_time_s) == 8.0
    failed_text = render_report(failed, failed_summary)
    @test occursin("No successful runs were recorded.", failed_text)
    @test !occursin("Fastest successful scenario", failed_text)
    @test occursin("Failed attempts: `2/2`", failed_text)
    @test occursin("Retry overhead: `1` retries across `1` requested runs", failed_text)
    @test occursin("[n/a, n/a]", failed_text)
    empty_text = render_report(raw[1:0, :], summary[1:0, :])
    @test occursin("Successful samples: `0/0` (`n/a%`)", empty_text)
    @test occursin("No successful runs were recorded.", empty_text)
    legacy = select(raw, Not([:attempt, :solver_fallback_count]))
    legacy_text = render_report(legacy, summary)
    @test occursin("Retry overhead: `0` retries across `4` requested runs", legacy_text)
    @test occursin("Mean fallback count across all attempts: `n/a`", legacy_text)
end

function entry_row(; kwargs...)
    merge((category="entry", scenario="entry_fixture", description="Entry fixture",
        outer_threads_safe=true, entry_run_role="measured",
        entry_atmospheric_interface_count=2, entry_passage_duration_s=10.0,
        entry_wall_time_per_passage_s=3.0, entry_event_time_abs_error_s=0.25,
        entry_reference_terminal_time_s=200.0, terminal_time_s=200.25,
        solve_success=true, total_time_s=6.0, sim_seconds_per_wall_second=2.0), (; kwargs...))
end
@testset "entry producer preserves roles and metrics through CSV and Markdown" begin
    entry_raw = DataFrame([entry_row(entry_run_role="reference", total_time_s=100.0,
                                    entry_passage_duration_s=80.0),
        entry_row(solve_success=false, total_time_s=999.0, entry_passage_duration_s=999.0),
        entry_row(), entry_row(scenario="z_failed_entry", solve_success=false),
        entry_row(entry_atmospheric_interface_count=3, entry_passage_duration_s=15.0)])
    before = deepcopy(entry_raw)
    supplied_columns = collect(eachcol(entry_raw))
    entry = summarize_entry_duration_results(entry_raw)
    @test isequal(entry_raw, before)
    @test all(a === b for (a, b) in zip(eachcol(entry_raw), supplied_columns))
    @test entry.scenario == ["entry_fixture", "entry_fixture", "entry_fixture", "z_failed_entry"]
    @test entry.entry_run_role == ["measured", "reference", "measured", "measured"]
    @test entry.entry_atmospheric_interface_count == [2, 2, 3, 2]
    @test entry.samples_total == [2, 1, 1, 1]
    @test entry.samples_success == [1, 1, 1, 0]
    @test entry.samples_failed == [1, 0, 0, 1]
    @test entry.passage_duration_mean_s[1] == 10.0
    @test entry.wall_time_per_passage_mean_s[1] == 3.0
    @test entry.event_time_abs_error_mean_s[1] == 0.25
    @test entry.reference_event_time_mean_s[1] == 200.0
    @test entry.event_time_mean_s[1] == 200.25
    @test entry.total_time_mean_s[1] == 6.0
    @test entry.passage_duration_mean_s[3] == 15.0
    @test ismissing(entry.total_time_mean_s[4])
    raw = DataFrame([raw_row()])
    summary = summarize_results(raw)
    section = entry_section(render_report(raw, summary, entry))
    @test occursin("| entry_fixture | 2 | 1/2 | 10.0 | 0.25 | 3.0 | 6.0 |", section)
    @test occursin("| entry_fixture | 3 | 1/1 | 15.0 |", section)
    @test !occursin("80.0", section)
    @test occursin("| z_failed_entry | 2 | 0/1 | n/a | n/a | n/a | n/a |", section)
    @test findfirst("| entry_fixture |", section).start < findfirst("| z_failed_entry |", section).start
    reference = entry[entry.entry_run_role .== "reference", :]
    @test occursin("no measured rows were available", entry_section(render_report(raw, summary, reference)))
    legacy_summary = select(reference, Not(:entry_run_role))
    @test occursin("80.0", entry_section(render_report(raw, summary, legacy_summary)))
    for absent in (nothing, DataFrame())
        @test occursin("No entry-duration sweep rows were produced.", entry_section(render_report(raw, summary, absent)))
    end
    mktempdir() do dir
        csv = joinpath(dir, "entry_summary.csv")
        CSV.write(csv, entry)
        restored = CSV.read(csv, DataFrame)
        @test isequal(restored, entry)
        @test entry_section(render_report(raw, summary, restored)) == section
    end
    legacy_raw = select(DataFrame([entry_row()]), Not([:entry_run_role,
        :entry_atmospheric_interface_count, :entry_passage_duration_s,
        :entry_wall_time_per_passage_s, :entry_event_time_abs_error_s,
        :entry_reference_terminal_time_s, :terminal_time_s]))
    legacy_before = deepcopy(legacy_raw)
    legacy_entry = summarize_entry_duration_results(legacy_raw)
    @test isequal(select(legacy_raw, names(legacy_before)), legacy_before)
    @test legacy_raw.entry_run_role == ["measured"]
    @test ismissing(only(legacy_entry.entry_atmospheric_interface_count))
    for col in (:passage_duration_mean_s, :wall_time_per_passage_mean_s,
                :event_time_abs_error_mean_s, :reference_event_time_mean_s, :event_time_mean_s)
        @test ismissing(only(legacy_entry[!, col]))
    end
    @test only(legacy_entry.total_time_mean_s) == 6.0
    @test occursin("| entry_fixture | missing | 1/1 | n/a | n/a | n/a | 6.0 |",
                   entry_section(render_report(raw, summary, legacy_entry)))
    failed_entry = summarize_entry_duration_results(DataFrame([entry_row(solve_success=false)]))
    @test only(failed_entry.samples_failed) == 1
    @test ismissing(only(failed_entry.event_time_mean_s))
    @test occursin("| entry_fixture | 2 | 0/1 | n/a | n/a | n/a | n/a |",
                   entry_section(render_report(raw, summary, failed_entry)))
    @test isempty(summarize_entry_duration_results(DataFrame()))
end

@testset "report fixture does not import native or simulation modules" begin
    for name in (:SpaceAGORA, :SimulationModel, :GRAMSuite, :SPICE, :Plots)
        @test !isdefined(@__MODULE__, name)
    end
end
end # module
