module BenchmarkPaperPlotColumnsTests
using Test
using CSV
using DataFrames
using Statistics
using Logging

const RUNTIME = normpath(joinpath(@__DIR__, "..", "..", "..", "benchmarks", "studies", "performance_runtime_analysis"))
const CORE_FIXTURE = joinpath(@__DIR__, "benchmark_plot_columns_tests.jl")

# Reuse only the recording renderer and definition loader. No core testset,
# setup invocation or other fixture statement is evaluated here.
function load_core_definition(wanted::Symbol)
    parsed = Meta.parseall(read(CORE_FIXTURE, String))
    outer = only(filter(x -> x isa Expr && x.head == :module, parsed.args))
    for expr in outer.args[end].args
        expr isa Expr || continue
        name = expr.head == :module ? expr.args[2] :
            expr.head == :function ? expr.args[1].args[1] : nothing
        if name == wanted
            Core.eval(@__MODULE__, expr)
            return
        end
    end
    error("Required recording definition $wanted not found")
end
load_core_definition(:Plots)
load_core_definition(:load_definition)

const HELPERS = joinpath(RUNTIME, "case_catalog", "profile_types_and_helpers.jl")
for name in (:ProfileSpec, :_parse_bool_env, :_safe_unique_join, :_mode_token, :_is_adaptive_mode_token)
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
const DEFAULT_OUTPUT_DIR = "unused_fixture_default"
const REPO_ROOT = dirname(RUNTIME)
for name in (:_plot_wrapped_label, :_ensure_runtime_plot_theme!, :_plot_margins,
             :_plot_metric_pairs, :_has_row_fields, :_save_runtime_plot!,
             :_paper_figure_pack_enabled, :_paper_ladder_outdir, :_paper_cross_machine_outdir,
             :_latest_profile_artifact_optional, :_read_optional_dataframe,
             :_paper_figure_external_data, :_sorted_orbit_groups, :generate_runtime_plots)
    load_definition(joinpath(RUNTIME, "reporting", "plots_and_reporting.jl"), name)
end
const SPEC = ProfileSpec(name="quick", repeats=1, warmup=0, max_attempts=1,
    mission_short_s=120.0, mission_long_s=600.0, montecarlo_samples=0,
    montecarlo_mission_s=120.0)
const EXTERNAL_FILES = (
    (:layer_attribution_speedup_df, "smart_parallel_ladder_layer_attribution_speedup", :ladder),
    (:deep_accuracy_df, "smart_parallel_ladder_deep_accuracy_parity", :ladder),
    (:montecarlo_parity_df, "smart_parallel_ladder_montecarlo_distribution_parity", :ladder),
    (:route_mix_df, "smart_parallel_ladder_route_mix", :ladder),
    (:cross_speedup_summary_df, "smart_parallel_ladder_cross_machine_speedup_summary", :cross),
    (:cross_adaptive_regret_summary_df, "smart_parallel_ladder_cross_machine_adaptive_regret_summary", :cross),
    (:cross_route_mix_summary_df, "smart_parallel_ladder_cross_machine_route_mix_summary", :cross),
)

empty_raw() = DataFrame([raw_row(seed=missing)])[1:0, :]
function route_raw(routes; markers=nothing)
    raw = DataFrame([raw_row(seed=missing, outer_route=route) for route in routes])
    markers === nothing || (raw[!, :is_terminal_attempt] = markers)
    return raw
end

function record_paper(; external=(;), inner=DataFrame(), density=DataFrame(),
                      raw=empty_raw(), split=nothing, multirate=nothing, enabled=true)
    # Baseline plots use real summaries of the existing synthetic reporting row.
    # Optional inputs are independently supplied artifacts, as in production.
    summary = summarize_results(DataFrame([raw_row(seed=missing)]))
    before = deepcopy((external, inner, density, raw, split, multirate, summary))
    empty!(Plots.SAVED)
    paths = String[]
    err = nothing
    mktempdir() do dir
        ladder = joinpath(dir, "ladder"); cross = joinpath(dir, "cross")
        mkpath(ladder); mkpath(cross)
        for (field, prefix, location) in EXTERNAL_FILES
            df = get(external, field, nothing)
            df === nothing && continue
            CSV.write(joinpath(location == :ladder ? ladder : cross, "$(prefix)_quick_fixture.csv"), df)
        end
        try
            paths = withenv("SPACEAGORA_PERF_PAPER_FIGURE_PACK" => string(enabled),
                           "SPACEAGORA_PERF_PAPER_LADDER_OUTDIR" => ladder,
                           "SPACEAGORA_PERF_PAPER_CROSS_MACHINE_OUTDIR" => cross) do
                @test_logs min_level=Logging.Warn generate_runtime_plots(dir, SPEC, "fixture", raw,
                    summary, DataFrame(), DataFrame(); inner_hint_layer_df=inner,
                    density_backend_breakdown_df=density, split_gate_df=split, multirate_gate_df=multirate)
            end
        catch caught
            caught isa InterruptException && rethrow()
            err = caught
        end
        @test all(path -> dirname(path) == dir && endswith(path, "_quick_fixture.png"), paths)
    end
    @test err === nothing
    @test isequal((external, inner, density, raw, split, multirate, summary), before)
    @test Set(paths) == Set(keys(Plots.SAVED))
    plots = Dict(replace(basename(path), "_quick_fixture.png" => "") => plt for (path, plt) in Plots.SAVED)
    return plots
end

function titled_panels(plots, artifact)
    found = Dict{String,Plots.RecordedPlot}()
    function visit(plt)
        plt isa Plots.RecordedPlot || return
        haskey(plt.attrs, :title) && (found[plt.attrs.title] = plt)
        foreach(visit, plt.args)
    end
    haskey(plots, artifact) && visit(plots[artifact])
    return found
end
const ADAPT = "runtime_plot_paper_adaptivity_behavior"
const LAYER = "runtime_plot_paper_layer_attribution"
const DENSITY = "runtime_plot_paper_density_backend"
const ACCURACY = "runtime_plot_paper_accuracy_suite"
const CROSS = "runtime_plot_paper_cross_machine"
const ARTIFACTS = (ADAPT, LAYER, DENSITY, ACCURACY, CROSS)

route_table() = DataFrame(mode=["fixed", "outer_inner_adaptive", "outer_inner_full_smart"],
    rung=["r0", "r2", "r1"], none_pct=[90.0, 10.0, missing],
    threads_pct=[5.0, 20.0, 50.0], process_pct=[5.0, 70.0, 50.0])
regret_table() = DataFrame(adaptive_mode=["mode_b", "mode_a"],
    mean_time_regret_pct=[5.0, missing], win_rate_pct=[80.0, 90.0])
hint_table() = DataFrame(layer=["density", "control", "skip"],
    confidence_mean=[0.2, 0.4, missing], regret_mean_ns=[2e6, 4e6, 8e6])

@testset "adaptive panels preserve selection and unit conversions" begin
    plots = record_paper(external=(route_mix_df=route_table(), cross_adaptive_regret_summary_df=regret_table()),
        inner=hint_table())
    panels = titled_panels(plots, ADAPT)
    route_title = "Adaptive Route-Choice Distribution"
    hint_title = "Inner Adaptive Hint Confidence/Regret by Layer"
    regret_title = "Adaptive Regret vs Best Fixed (Cross-Machine)"
    @test Set(keys(panels)) == Set((route_title, hint_title, regret_title))
    if haskey(panels, route_title)
        @test panels[route_title].args[1] == ["r2", "r1"]
        @test panels[route_title].args[2] == [10.0 20.0 70.0; 0.0 50.0 50.0]
        @test panels[route_title].attrs.ylabel == "Route Share [%]"
    end
    if haskey(panels, hint_title)
        @test panels[hint_title].args == (["density", "control"], [0.2, 0.4])
        @test panels[hint_title].series[1].args == (["density", "control"], [2.0, 4.0])
        @test panels[hint_title].series[1].attrs.label == "Mean regret [ms]"
    end
    if haskey(panels, regret_title)
        @test panels[regret_title].args[1] == ["mode\nb", "mode\na"]
        @test isequal(panels[regret_title].args[2], [5.0, NaN])
        @test panels[regret_title].series[1].args[2] == [80.0, 90.0]
    end
    for (df, expected) in ((select(route_table(), Not(:rung)), ["outer\ninner\nadaptive", "outer\ninner\nfull\nsmart"]),
                          (select(route_table(), Not([:rung, :mode])), ["adaptive\n1", "adaptive\n2", "adaptive\n3"]))
        panels = titled_panels(record_paper(external=(route_mix_df=df,)), ADAPT)
        @test haskey(panels, route_title)
        haskey(panels, route_title) && (@test panels[route_title].args[1] == expected)
    end
end

@testset "raw route fallback preserves terminal outcomes and legacy rows" begin
    title = "Observed Outer-Route Distribution"
    for (markers, expected) in ((nothing, [0.0, 25.0, 25.0, 50.0]),
            ([false, true, missing, true], [0.0, 0.0, 100/3, 200/3]),
            ([false, false, false, false], nothing))
        raw = route_raw(["threads", " PROCESS ", "legacy_unknown", "other"]; markers=markers)
        panels = titled_panels(record_paper(raw=raw), ADAPT)
        @test haskey(panels, title) == (expected !== nothing)
        if expected !== nothing && haskey(panels, title)
            @test panels[title].args[1] == ["none", "threads", "process", "other"]
            @test panels[title].args[2] ≈ expected
            @test panels[title].attrs.ylabel == "Share of Runs [%]"
        end
    end
    # Only the Boolean false marker is excluded. Integer zero and missing retain
    # the compatibility behavior of a row that is not explicitly intermediate.
    panels = titled_panels(record_paper(raw=route_raw(["none", "process", "threads"];
        markers=Any[0, missing, false])), ADAPT)
    @test haskey(panels, title)
    haskey(panels, title) && (@test panels[title].args[2] == [50.0, 0.0, 50.0, 0.0])
    # A supported external route table keeps precedence over the raw fallback.
    panels = titled_panels(record_paper(external=(route_mix_df=route_table(),),
        raw=route_raw(["none"])), ADAPT)
    @test !haskey(panels, title)
    @test haskey(panels, "Adaptive Route-Choice Distribution")
end

layer_table() = DataFrame(layer_set=["density", "outer_only", " Density ", "thermal", "control", "effector", "unknown"],
    total_speedup_vs_outer_only=[1.3, 1.0, 1.7, Inf, missing, 1.1, 999.0])
density_table() = DataFrame(density_backend_bucket=["non_gram", "gram_point_to_point", "gram_surrogate"],
    total_time_mean_s=[8.0, 20.0, missing], sim_seconds_per_wall_second_mean=[2.0, 5.0, missing],
    success_rate_pct=[100.0, 50.0, 0.0])

@testset "layer and density panels retain existing metrics and ordering" begin
    plots = record_paper(external=(layer_attribution_speedup_df=layer_table(),), density=density_table())
    @test haskey(plots, LAYER)
    if haskey(plots, LAYER)
        @test plots[LAYER].args == (["outer\nonly", "density", "effector"], [1.0, 1.7, 1.1])
        @test plots[LAYER].series[1].args == ([1.0],)
    end
    panels = titled_panels(plots, DENSITY)
    titles = ("Density Backend Benchmark: Runtime", "Density Backend Benchmark: Throughput", "Density Backend Benchmark: Solve Success")
    @test Set(keys(panels)) == Set(titles)
    for (title, expected) in zip(titles, ([20.0, NaN, 8.0], [5.0, NaN, 2.0], [50.0, 0.0, 100.0]))
        if haskey(panels, title)
            @test panels[title].args[1] == ["gram\npoint\nto\npoint", "gram\nsurrogate", "non\ngram"]
            @test isequal(panels[title].args[2], expected)
        end
    end
    failed = density_table()
    failed.total_time_mean_s .= missing
    failed.sim_seconds_per_wall_second_mean .= missing
    failed.success_rate_pct .= 0.0
    panels = titled_panels(record_paper(density=failed), DENSITY)
    @test Set(keys(panels)) == Set(titles)
    if haskey(panels, titles[1])
        @test all(isnan, panels[titles[1]].args[2])
        @test panels[titles[3]].args[2] == [0.0, 0.0, 0.0]
    end
end

deep_table() = DataFrame(mode=["z", "a"], rung=["r_z", "r_a"],
    traj_pos_rel_rms_median_pct=[2.0, missing], traj_vel_rel_rms_median_pct=[4.0, 3.0],
    periapsis_time_abs_err_p90_s=[10.0, missing], interface_time_abs_err_p90_s=[20.0, 15.0],
    propellant_rel_err_p90_pct=[0.8, 0.2], control_impulse_rel_err_p90_pct=[0.9, 0.3],
    callback_exact_match_pct=[99.0, 100.0])
mc_table() = DataFrame(mode=["z", "a", "a"], rung=["r_z", "r_a", "r_a"],
    event_time_ks_distance=[0.8, 0.1, 0.5], pos_ks_distance=[0.9, 0.2, 0.6],
    vel_ks_distance=[0.7, 0.3, 0.5], mass_ks_distance=[0.5, missing, 0.7])
const DEEP_TITLES = ("Accuracy Parity: Trajectory RMS Relative Error", "Accuracy Parity: Event-Time Error",
    "Accuracy Parity: Control/Propellant", "Accuracy Parity: Callback-State Exact Match")
const MC_TITLE = "Accuracy Parity: Monte Carlo Distribution KS Distance"
const GATE_TITLE = "Accuracy Gate Fallback: Trajectory Relative Error Maxima"

@testset "accuracy panels preserve mode ordering medians and fallbacks" begin
    split = DataFrame(scenario=["split_case"], pos_rel_max=[0.1], vel_rel_max=[0.2])
    multirate = DataFrame(scenario=["failed_case"], pos_rel_max=[missing], vel_rel_max=[missing])
    panels = titled_panels(record_paper(external=(deep_accuracy_df=deep_table(), montecarlo_parity_df=mc_table()),
        split=split), ACCURACY)
    @test Set(keys(panels)) == Set((DEEP_TITLES..., MC_TITLE))
    expected = ([NaN 3.0; 2.0 4.0], [NaN 15.0; 10.0 20.0], [0.2 0.3; 0.8 0.9], [100.0, 99.0])
    for (title, values) in zip(DEEP_TITLES, expected)
        if haskey(panels, title)
            @test panels[title].args[1] == ["r\na", "r\nz"]
            @test isequal(panels[title].args[2], values)
        end
    end
    if haskey(panels, MC_TITLE)
        @test panels[MC_TITLE].args[1] == ["r\na", "r\nz"]
        @test panels[MC_TITLE].args[2] ≈ [0.3 0.4 0.4 0.7; 0.8 0.9 0.7 0.5]
        @test panels[MC_TITLE].attrs.ylabel == "KS Distance"
    end
    @test !haskey(panels, GATE_TITLE)
    panels = titled_panels(record_paper(split=split, multirate=multirate), ACCURACY)
    @test Set(keys(panels)) == Set((GATE_TITLE,))
    if haskey(panels, GATE_TITLE)
        @test panels[GATE_TITLE].args[1] == ["split\ncase", "failed\ncase"]
        @test isequal(panels[GATE_TITLE].args[2], [0.1 0.2; NaN NaN])
        @test panels[GATE_TITLE].attrs.ylabel == "Relative Error"
    end
    # Mode is optional; without it the external row ordering remains unchanged.
    panels = titled_panels(record_paper(external=(deep_accuracy_df=select(deep_table(), Not(:mode)),)), ACCURACY)
    @test haskey(panels, DEEP_TITLES[1])
    haskey(panels, DEEP_TITLES[1]) && (@test panels[DEEP_TITLES[1]].args[1] == ["r\nz", "r\na"])
    for (col, absent) in ((:periapsis_time_abs_err_p90_s, DEEP_TITLES[2]),
                          (:interface_time_abs_err_p90_s, DEEP_TITLES[2]),
                          (:propellant_rel_err_p90_pct, DEEP_TITLES[3]),
                          (:control_impulse_rel_err_p90_pct, DEEP_TITLES[3]),
                          (:callback_exact_match_pct, DEEP_TITLES[4]))
        panels = titled_panels(record_paper(external=(deep_accuracy_df=select(deep_table(), Not(col)),)), ACCURACY)
        @test !haskey(panels, absent)
        @test haskey(panels, DEEP_TITLES[1])
    end
end

speed_table() = DataFrame(rung=["r_slow", "r_fast"], median_speedup_vs_r0=[1.2, 2.4])
cross_route_table() = rename(route_table(), :none_pct=>:none_pct_mean,
    :threads_pct=>:threads_pct_mean, :process_pct=>:process_pct_mean)
const CROSS_TITLES = ("Cross-Machine Median Speedup vs R0", "Cross-Machine Adaptive Regret vs Best Fixed", "Cross-Machine Adaptive Route Mix")

@testset "cross-machine panels preserve independent artifact statistics" begin
    plots = record_paper(external=(cross_speedup_summary_df=speed_table(),
        cross_adaptive_regret_summary_df=regret_table(), cross_route_mix_summary_df=cross_route_table()))
    panels = titled_panels(plots, CROSS)
    @test Set(keys(panels)) == Set(CROSS_TITLES)
    if haskey(panels, CROSS_TITLES[1])
        @test panels[CROSS_TITLES[1]].args == (["r\nfast", "r\nslow"], [2.4, 1.2])
        @test panels[CROSS_TITLES[1]].series[1].args == ([1.0],)
    end
    if haskey(panels, CROSS_TITLES[2])
        @test isequal(panels[CROSS_TITLES[2]].args, (["mode\nb", "mode\na"], [5.0, NaN]))
        @test panels[CROSS_TITLES[2]].series[1].args[2] == [80.0, 90.0]
    end
    if haskey(panels, CROSS_TITLES[3])
        @test panels[CROSS_TITLES[3]].args == (["r2", "r1"], [10.0 20.0 70.0; 0.0 50.0 50.0])
    end
    panels = titled_panels(record_paper(external=(cross_adaptive_regret_summary_df=select(regret_table(), Not(:win_rate_pct)),)), CROSS)
    @test haskey(panels, CROSS_TITLES[2])
    haskey(panels, CROSS_TITLES[2]) && (@test all(isnan, panels[CROSS_TITLES[2]].series[1].args[2]))
    for (df, expected) in ((select(cross_route_table(), Not(:rung)), ["outer\ninner\nadaptive", "outer\ninner\nfull\nsmart"]),
                          (select(cross_route_table(), Not([:rung, :mode])), ["adaptive\n1", "adaptive\n2", "adaptive\n3"]))
        panels = titled_panels(record_paper(external=(cross_route_mix_summary_df=df,)), CROSS)
        @test haskey(panels, CROSS_TITLES[3])
        haskey(panels, CROSS_TITLES[3]) && (@test panels[CROSS_TITLES[3]].args[1] == expected)
    end
end

@testset "optional panels skip absent empty and incomplete supported schemas" begin
    @test all(name -> !haskey(record_paper(), name), ARTIFACTS)
    families = ((:route_mix_df, route_table(), (:none_pct, :threads_pct, :process_pct), ADAPT),
        (:layer_attribution_speedup_df, layer_table(), (:layer_set, :total_speedup_vs_outer_only), LAYER),
        (:deep_accuracy_df, deep_table(), (:rung, :traj_pos_rel_rms_median_pct, :traj_vel_rel_rms_median_pct), ACCURACY),
        (:montecarlo_parity_df, mc_table(), (:mode, :rung, :event_time_ks_distance), ACCURACY),
        (:cross_speedup_summary_df, speed_table(), (:rung, :median_speedup_vs_r0), CROSS),
        (:cross_adaptive_regret_summary_df, regret_table(), (:adaptive_mode, :mean_time_regret_pct), CROSS),
        (:cross_route_mix_summary_df, cross_route_table(), (:none_pct_mean, :threads_pct_mean, :process_pct_mean), CROSS))
    for (field, full, required, artifact) in families
        for df in (full[1:0, :], (select(full, Not(col)) for col in required)...)
            plots = record_paper(external=NamedTuple{(field,)}((df,)))
            @test !haskey(plots, artifact)
        end
    end
    for df in (hint_table()[1:0, :], (select(hint_table(), Not(col)) for col in (:layer, :confidence_mean, :regret_mean_ns))...)
        @test !haskey(record_paper(inner=df), ADAPT)
    end
    for df in (density_table()[1:0, :], (select(density_table(), Not(col)) for col in
            (:density_backend_bucket, :total_time_mean_s, :sim_seconds_per_wall_second_mean))...)
        @test !haskey(record_paper(density=df), DENSITY)
    end
    # Disabled pack keeps every optional artifact absent despite complete input.
    plots = record_paper(external=(route_mix_df=route_table(), layer_attribution_speedup_df=layer_table(),
        deep_accuracy_df=deep_table(), montecarlo_parity_df=mc_table(), cross_speedup_summary_df=speed_table(),
        cross_adaptive_regret_summary_df=regret_table(), cross_route_mix_summary_df=cross_route_table()),
        inner=hint_table(), density=density_table(), enabled=false)
    @test all(name -> !haskey(plots, name), ARTIFACTS)
end

@testset "paper plotting fixture stays native-free" begin
    @test Plots.RECORDING_BACKEND
    for name in (:SpaceAGORA, :SimulationModel, :GRAMSuite, :SPICE)
        @test !isdefined(@__MODULE__, name)
    end
end
end # module
