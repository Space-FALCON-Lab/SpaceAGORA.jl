module RPOComparisonReportDataTests
using Test, TOML, SpaceAGORA
# The figures come from the PlotlyJS package extension; the tests that make them
# load it here. The child processes below test the refusals without it.
import PlotlyJS
const G = SpaceAGORA.SimulationModel.GuidanceHooks
const Reports = G.RPOComparisonReports
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))

# Every symbol in a file's code. Docstrings are strings and comments are not
# parsed, so neither counts.
function _code_symbols(path)
    symbols = Set{Symbol}()
    walk(x::Symbol) = push!(symbols, x)
    walk(x::Expr) = foreach(walk, x.args)
    walk(x) = nothing
    walk(Meta.parseall(read(path, String); filename=path))
    return symbols
end

@testset "PlotlyJS is a weak dependency and only its extension uses it" begin
    project = TOML.parsefile(joinpath(ROOT, "Project.toml"))
    @test !haskey(project["deps"], "PlotlyJS")
    @test haskey(project["weakdeps"], "PlotlyJS")
    @test project["extensions"]["SpaceAGORAPlotlyJSExt"] == "PlotlyJS"
    users = String[]
    for (dir, _, files) in walkdir(joinpath(ROOT, "src")), file in files
        path = joinpath(dir, file)
        endswith(file, ".jl") && :PlotlyJS in _code_symbols(path) && push!(users, relpath(path, ROOT))
    end
    @test isempty(users)
    @test :PlotlyJS in _code_symbols(joinpath(ROOT, "ext", "SpaceAGORAPlotlyJSExt.jl"))
    # The 0.5 qualified bindings are gone.
    @test !isdefined(G, :PlotlyJS)
    @test !isdefined(Reports, :PlotlyJS)
    @test !isdefined(SpaceAGORA.RPOVisualization, :PlotlyJS)
    @test Base.get_extension(SpaceAGORA, :SpaceAGORAPlotlyJSExt) !== nothing
end

include(joinpath(@__DIR__, "rpo_comparison_report_fixture.jl"))

@testset "Each figure shows the data its preparation returns" begin
    series = Reports._rpo_comparison_metric_series(batch)
    @test [s.metric for s in series] == first.(G.rpo_comparison_metric_specs())
    fuel = only(filter(s -> s.metric === :fuel_used_pct, series))
    @test fuel.planners == ["RRT-Connect", "CHOMP"] # the failed run is left out of the fuel aggregate
    @test fuel.values == [1.0, 3.0]
    @test fuel.case_labels == ["case a", "case a"]
    success = only(filter(s -> s.metric === :success_pct, series))
    @test success.values == [100.0, 0.0, 100.0]
    boxes = filter(t -> get(t.fields, :type, "") == "box", G.rpo_comparison_metric_summary_plot(batch).data)
    @test [b.fields[:y] for b in boxes] == [s.values for s in series]
    @test [b.fields[:x] for b in boxes] == [s.planners for s in series]

    @test Reports._rpo_comparison_failed_indices(batch.results_by_planner[:rrt_connect]) == [2]
    @test Reports._rpo_comparison_failed_indices(batch.results_by_planner[:chomp]) == Int[]

    costs = Reports._rpo_comparison_cost_series(batch, :rrt_connect)
    @test [(c.case_label, c.iterations, c.costs) for c in costs] == [("a", [1, 2], [5.0, 3.0])]
    single = Reports._rpo_comparison_cost_series(batch, :chomp)
    @test [(c.case_label, c.iterations, c.costs) for c in single] == [("a", [6], [4.0])]
    plot = G.rpo_comparison_cost_iteration_plot(batch; planner=:rrt_connect)
    @test [(t.fields[:x], t.fields[:y]) for t in plot.data] == [([1, 2], [5.0, 3.0])]

    with_raw = Reports._rpo_comparison_case_path(batch.plans_by_planner[:rrt_connect][2])
    @test with_raw.has_raw && with_raw.raw == [0 2; 0 1; 0 0]
    @test with_raw.tracked == batch.plans_by_planner[:rrt_connect][2].tracking.x_hist[1:3, :]
    without_raw = Reports._rpo_comparison_case_path(batch.plans_by_planner[:chomp][1])
    @test !without_raw.has_raw && without_raw.raw === nothing
    failed_plot = G.rpo_comparison_failed_paths_plot(batch; planner=:rrt_connect, failed_indices=[2])
    names = [t.fields[:name] for t in failed_plot.data]
    @test names == ["station point cloud", "case b desired", "case b tracked", "case b raw planner path",
                    "case b start", "case b goal"]

    @test Reports._rpo_comparison_station_point_indices(batch.geometry.station.points_body, 2) == 1:2:4
    @test !Reports._rpo_comparison_has_station_mesh(batch)
    @test !Reports._rpo_comparison_has_station_mesh(merge(batch, (station_triangles=nothing,)))
    mesh = Reports._rpo_comparison_station_mesh(Float64[0 1 0; 0 0 1; 0 0 0])
    @test (mesh.i, mesh.j, mesh.k) == ([0], [1], [2])
    @test Reports._rpo_comparison_has_station_mesh(merge(batch, (station_triangles=mesh.vertices,)))
end

@testset "The figures need PlotlyJS, in either load order" begin
    # Fresh processes in this environment. Their load path is this environment,
    # the one that provides the HYPR companion (the comparison entry needs it)
    # and the standard library, not this process's whole load path: an earlier
    # test may have put another project first, as the native GRAM tests do with
    # GRAMSuite's, and the children would resolve packages there. GRAMSuite's
    # manifest has StaticArrays without ChainRulesCore, so StaticArrays'
    # ChainRulesCore extension then failed to precompile in a coverage run.
    child = joinpath(@__DIR__, "rpo_plotlyjs_extension_child.jl")
    project = dirname(Base.active_project())
    hypr = Base.identify_package_env("HYPR")
    hypr === nothing && error("HYPR is not on the load path; load test/helpers/load_hypr.jl first.")
    load_path = join(unique(["@", last(hypr), "@stdlib"]), Sys.iswindows() ? ';' : ':')
    for order in ("core-first", "plotlyjs-first")
        output = IOBuffer()
        command = addenv(`$(Base.julia_cmd()) --startup-file=no --project=$(project) $(child) $(order)`,
            "JULIA_LOAD_PATH" => load_path)
        process = run(pipeline(ignorestatus(command); stdout=output, stderr=output))
        text = String(take!(output))
        @test success(process)
        @test occursin("PLOTLYJS_EXTENSION_CHILD_OK $(order)", text)
        success(process) || print(text)
    end
end

end # module RPOComparisonReportDataTests
