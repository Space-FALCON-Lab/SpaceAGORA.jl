# Child process of rpo_comparison_report_data_tests.jl: the RPO Plotly figures in
# one load order, `core-first` (SpaceAGORA, then PlotlyJS) or `plotlyjs-first`.
# Core-first also checks every refusal before PlotlyJS loads.
using Test
const ORDER = only(ARGS)
ORDER in ("core-first", "plotlyjs-first") || error("Unknown load order $(ORDER).")
ORDER == "plotlyjs-first" && @eval import PlotlyJS
using SpaceAGORA
const G = SpaceAGORA.SimulationModel.GuidanceHooks
include(joinpath(@__DIR__, "rpo_comparison_report_fixture.jl"))

extension() = Base.get_extension(SpaceAGORA, :SpaceAGORAPlotlyJSExt)
const FIGURES = (
    () -> SpaceAGORA.rpo_path_plot(zeros(3, 2)),
    () -> SpaceAGORA.rpo_tracking_plot([0.0, 1.0], [1.0, 0.5]),
    () -> G.rpo_comparison_metric_summary_plot(batch),
    () -> G.rpo_comparison_station_trace(batch),
    () -> G.rpo_comparison_station_mesh_trace(Float64[0 1 0; 0 0 1; 0 0 0]),
    () -> G.rpo_comparison_path_family_plot(batch),
    () -> G.rpo_comparison_single_path_plot(batch, batch.plans_by_planner[:chomp][1]; planner=:chomp),
    () -> G.rpo_comparison_failed_paths_plot(batch; planner=:rrt_connect, failed_indices=[2]),
    () -> G.rpo_comparison_cost_iteration_plot(batch; planner=:chomp),
)

if ORDER == "core-first"
    # The planner comparison's entry needs the HYPR companion before it checks the plot flags.
    @eval using HYPR
    @testset "Without PlotlyJS the figures refuse and the CSV is still written" begin
        @test extension() === nothing
        for figure in FIGURES
            @test_throws PlotlyJSUnavailableError figure()
        end
        message = sprint(showerror, PlotlyJSUnavailableError())
        for remedy in ("`Pkg.add(\"PlotlyJS\")`", "`using PlotlyJS`",
                       "`write_plotly_outputs=false`", "`write_failed_path_outputs=false`")
            @test occursin(remedy, message)
        end
        mktempdir() do dir
            # Asked for figures, the writer refuses before it writes anything.
            for (plots, failed) in ((true, false), (false, true))
                out = joinpath(dir, "figures_$(plots)_$(failed)")
                @test_throws PlotlyJSUnavailableError G.rpo_write_planner_comparison_outputs(batch;
                    output_dir=out, write_plotly_outputs=plots, write_failed_path_outputs=failed)
                @test !ispath(out)
            end
            out = joinpath(dir, "csv_only")
            outputs = G.rpo_write_planner_comparison_outputs(batch; output_dir=out,
                write_plotly_outputs=false, write_failed_path_outputs=false)
            @test readdir(out) == ["rpo_planner_comparison_results.csv"]
            @test outputs.metrics_plot === nothing
            @test isempty(outputs.failed_path_plots)
        end
        # A comparison asked for figures refuses before it plans anything: no
        # case or geometry is read.
        for (plots, failed) in ((true, false), (false, true))
            config = G.RPOPlannerComparisonConfig(write_plotly_outputs=plots,
                write_failed_path_outputs=failed, show_progress=false)
            @test_throws PlotlyJSUnavailableError G.rpo_run_planner_comparison_batch([], nothing, config)
        end
    end
    @eval import PlotlyJS
end

@testset "With PlotlyJS loaded ($(ORDER)) the figures are made" begin
    @test extension() !== nothing
    for figure in FIGURES
        @test figure() isa Union{PlotlyJS.Plot, PlotlyJS.GenericTrace}
    end
    mktempdir() do dir
        outputs = G.rpo_write_planner_comparison_outputs(batch; output_dir=dir,
            write_plotly_outputs=true, write_failed_path_outputs=true)
        @test isfile(outputs.metrics_plot)
        @test all(isfile, values(outputs.method_path_plots))
        @test collect(keys(outputs.failed_path_plots)) == [:rrt_connect]
    end
end
println("PLOTLYJS_EXTENSION_CHILD_OK ", ORDER)
