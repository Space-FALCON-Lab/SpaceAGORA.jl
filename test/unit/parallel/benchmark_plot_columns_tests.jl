module BenchmarkPlotColumnsTests
using Test
using DataFrames
using Statistics
using Logging

const RUNTIME = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies", "performance_runtime_analysis"))

# Only rendering and renderer availability are replaced. Selection, sorting,
# metric aliases, artifact names and the full plot function are production code.
module Plots
const RECORDING_BACKEND = true
const mm = 1
mutable struct RecordedPlot
    kind::Symbol
    args::Tuple
    attrs::NamedTuple
    series::Vector{Any}
end
const SAVED = Dict{String,RecordedPlot}()
theme(args...) = nothing
default(; kwargs...) = nothing
font(args...) = args
cgrad(args...) = args
for name in (:plot, :bar, :heatmap, :histogram)
    @eval $name(args...; kwargs...) = RecordedPlot($(QuoteNode(name)), deepcopy(args), (; kwargs...), Any[])
end
for name in (:plot!, :scatter!, :hline!, :vline!)
    @eval function $name(plt::RecordedPlot, args...; kwargs...)
        push!(plt.series, (kind=$(QuoteNode(name)), args=deepcopy(args), attrs=(; kwargs...)))
        return plt
    end
end
function savefig(plt::RecordedPlot, path::String)
    SAVED[path] = deepcopy(plt)
    return nothing
end
end # recording renderer, no native plotting import

function load_definition(path::String, wanted::Symbol)
    nodes = collect(Meta.parseall(read(path, String)).args)
    for expr in nodes
        expr isa Expr || continue
        if expr.head == :module
            append!(nodes, expr.args[end].args)
            continue
        end
        node = expr.head == :macrocall ? expr.args[end] : expr
        node isa Expr || continue
        name = if node.head == :struct
            node.args[2]
        elseif node.head in (:function, :(=))
            signature = node.args[1]
            signature isa Expr || continue
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
    error("Required plotting definition $wanted not found in $path")
end

const HELPERS = joinpath(RUNTIME, "case_catalog", "profile_types_and_helpers.jl")
for name in (:ProfileSpec, :_parse_bool_env, :_safe_unique_join)
    load_definition(HELPERS, name)
end
for name in (:_plot_label, :_plot_axis_label, :_plot_number)
    load_definition(joinpath(RUNTIME, "reporting", "rollout_gates.jl"), name)
end
struct BenchmarkCase end
const PERF_BASELINE_SCENARIO = "baseline_fixture"
include(joinpath(RUNTIME, "reporting", "summaries.jl"))
load_definition(joinpath(@__DIR__, "benchmark_report_tests.jl"), :raw_row)
_plot_ready() = true
const _runtime_plot_theme_applied = Ref(false)
for name in (:_plot_wrapped_label, :_ensure_runtime_plot_theme!, :_plot_margins,
             :_plot_metric_pairs, :_has_row_fields, :_save_runtime_plot!,
             :_paper_figure_pack_enabled, :_sorted_orbit_groups, :generate_runtime_plots)
    load_definition(joinpath(RUNTIME, "reporting", "plots_and_reporting.jl"), name)
end
const SPEC = ProfileSpec(name="quick", repeats=1, warmup=0, max_attempts=1,
    mission_short_s=120.0, mission_long_s=600.0, montecarlo_samples=0,
    montecarlo_mission_s=120.0)
const ORBIT_NAMES = ["runtime_plot_per_orbit_scaling", "runtime_plot_per_orbit_efficiency",
    "runtime_plot_per_orbit_heatmap"]
const ENTRY_NAMES = ["runtime_plot_entry_duration_passage", "runtime_plot_entry_duration_wall_per_passage",
    "runtime_plot_entry_duration_event_time_error"]

function record_core(orbit=DataFrame(), entry=DataFrame())
    raw = DataFrame([raw_row(seed=missing)])
    summary = summarize_results(raw)
    before = deepcopy((raw, summary, orbit, entry))
    empty!(Plots.SAVED)
    paths = String[]
    err = nothing
    mktempdir() do dir
        try
            paths = withenv("SPACEAGORA_PERF_PAPER_FIGURE_PACK" => "false") do
                @test_logs min_level=Logging.Warn generate_runtime_plots(dir, SPEC, "fixture", raw, summary, orbit, entry)
            end
        catch caught
            caught isa InterruptException && rethrow()
            err = caught
        end
        @test all(path -> dirname(path) == dir && endswith(path, "_quick_fixture.png"), paths)
    end
    @test err === nothing
    @test isequal((raw, summary, orbit, entry), before)
    @test Set(paths) == Set(keys(Plots.SAVED))
    plots = Dict(replace(basename(path), "_quick_fixture.png" => "") => plt for (path, plt) in Plots.SAVED)
    return (plots=plots, paths=paths, error=err)
end

function orbit_summary(; schema=:modern)
    df = DataFrame(scenario=["slow", "fast", "slow", "fast", "failed"],
        samples_success=[1, 1, 1, 1, 0], total_time_mean_s=[80.0, 16.0, 30.0, 8.0, 100.0],
        time_per_orbit_mean_s=[20.0, 4.0, 15.0, 4.0, 50.0],
        orbits_per_wall_second_mean=[0.05, 0.25, 2/30, 0.25, 0.02])
    if schema != :legacy
        df[!, :mission_time_multiplier] = [4, 4, 2, 2, 2]
        df[!, :time_per_baseline_period_mean_s] = copy(df.time_per_orbit_mean_s)
        df[!, :baseline_periods_per_wall_second_mean] = copy(df.orbits_per_wall_second_mean)
    end
    if schema != :modern
        df[!, :orbit_count] = schema == :conflict ? [99, 98, 97, 96, 95] : [4, 4, 2, 2, 2]
    end
    if schema == :conflict
        # Keep legacy columns nonmissing so only alias precedence determines
        # displayed values. The modern aliases are deliberately different.
        df.time_per_orbit_mean_s .= 999.0
        df.orbits_per_wall_second_mean .= 999.0
    end
    return df
end

@testset "core orbit plots select modern and legacy columns" begin
    for schema in (:modern, :legacy, :dual, :conflict)
        @testset "$schema" begin
            observed = record_core(orbit_summary(; schema=schema))
            @test all(name -> haskey(observed.plots, name), ORBIT_NAMES)
            @test all(name -> !haskey(observed.plots, name), ENTRY_NAMES)
            if all(name -> haskey(observed.plots, name), ORBIT_NAMES)
                scaling = observed.plots[ORBIT_NAMES[1]]
                @test scaling.attrs.xlabel == "Mission-Time Multiplier [x baseline period]"
                @test scaling.attrs.ylabel == "Mean Time per Baseline-Period Unit [s]"
                @test [s.attrs.label for s in scaling.series] == ["slow", "fast"]
                @test [s.args[1] for s in scaling.series] == [[2, 4], [2, 4]]
                @test [s.args[2] for s in scaling.series] == [[15.0, 20.0], [4.0, 4.0]]
                efficiency = observed.plots[ORBIT_NAMES[2]]
                @test efficiency.attrs.ylabel == "Baseline-Period Units / Wall-sec"
                @test [s.attrs.label for s in efficiency.series] == ["fast", "slow"]
                @test [s.args[1] for s in efficiency.series] == [[2, 4], [2, 4]]
                @test efficiency.series[1].args[2] == [0.25, 0.25]
                @test efficiency.series[2].args[2] ≈ [2/30, 0.05]
                heat = observed.plots[ORBIT_NAMES[3]]
                @test heat.kind == :heatmap
                @test heat.args[1] == [2, 4]
                @test collect(heat.args[2]) == [1, 2]
                @test heat.args[3] == [15.0 20.0; 4.0 4.0]
                @test heat.attrs.yticks[2] == ["slow", "fast"]
                @test heat.attrs.colorbar_title == "s / baseline-period unit"
            end
        end
    end
end

function entry_summary()
    return DataFrame(scenario=["entry", "reference", "entry", "entry", "failed"],
        entry_run_role=["measured", "serial_reference", "measured", "measured", "measured"],
        entry_atmospheric_interface_count=[3, 1, 1, missing, 2],
        passage_duration_mean_s=[30.0, 999.0, 10.0, 777.0, missing],
        wall_time_per_passage_mean_s=[6.0, 999.0, 2.0, 777.0, missing],
        event_time_abs_error_mean_s=[0.3, 999.0, 0.1, 777.0, missing])
end

@testset "entry plots select measured complete records" begin
    observed = record_core(DataFrame(), entry_summary())
    @test all(name -> haskey(observed.plots, name), ENTRY_NAMES)
    @test all(name -> !haskey(observed.plots, name), ORBIT_NAMES)
    for (name, expected, ylabel) in zip(ENTRY_NAMES,
            ([10.0, 30.0], [2.0, 6.0], [0.1, 0.3]),
            ("Passage Duration [s]", "Wall Time per Passage [s]", "|Event-Time Error| [s]"))
        if haskey(observed.plots, name)
            plt = observed.plots[name]
            @test plt.attrs.xlabel == "Atmospheric-Interface Count"
            @test plt.attrs.ylabel == ylabel
            @test length(plt.series) == 1
            @test plt.series[1].attrs.label == "entry"
            @test plt.series[1].args[1] == [1, 3]
            @test plt.series[1].args[2] == expected
        end
    end
    for (col, absent_name) in zip((:passage_duration_mean_s, :wall_time_per_passage_mean_s,
                                   :event_time_abs_error_mean_s), ENTRY_NAMES)
        observed = record_core(DataFrame(), select(entry_summary(), Not(col)))
        @test !haskey(observed.plots, absent_name)
        @test all(name -> haskey(observed.plots, name), filter(!=(absent_name), ENTRY_NAMES))
    end
end

@testset "core plots skip empty incomplete and unsuccessful inputs" begin
    for col in (:samples_success, :total_time_mean_s, :orbits_per_wall_second_mean, :time_per_orbit_mean_s)
        observed = record_core(select(orbit_summary(), Not(col)))
        @test all(name -> !haskey(observed.plots, name), ORBIT_NAMES)
    end
    failed = orbit_summary()
    failed.samples_success .= 0
    no_multiplier = orbit_summary()
    no_multiplier[!, :mission_time_multiplier] = fill(missing, nrow(no_multiplier))
    for orbit in (DataFrame(), orbit_summary()[1:0, :], failed, no_multiplier)
        observed = record_core(orbit)
        @test all(name -> !haskey(observed.plots, name), ORBIT_NAMES)
    end
    entry = entry_summary()
    for absent in (:entry_run_role, :entry_atmospheric_interface_count)
        observed = record_core(DataFrame(), select(entry, Not(absent)))
        @test all(name -> !haskey(observed.plots, name), ENTRY_NAMES)
    end
    for entry in (DataFrame(), entry_summary()[1:0, :], entry_summary()[2:2, :], entry_summary()[5:5, :])
        observed = record_core(DataFrame(), entry)
        @test all(name -> !haskey(observed.plots, name), ENTRY_NAMES)
    end
end

@testset "recording fixture avoids native workloads" begin
    @test Plots.RECORDING_BACKEND
    for name in (:SpaceAGORA, :SimulationModel, :GRAMSuite, :SPICE)
        @test !isdefined(@__MODULE__, name)
    end
end
end # module
