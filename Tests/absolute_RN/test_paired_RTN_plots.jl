import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..", "9_SpaceAGORA.jl-main_full integration")))
using Test, CSV, Plots
include(joinpath(@__DIR__, "..", "..", "Simulation_absolute_RN", "Code", "plot_paired_dv_RTN.jl"))
const PP = PairedRTNPlots
const OE = PP.OE
const ROOT = joinpath(OE.StudyConfig.DEFAULT_OUTPUT, "images", "delta_v_RTN")
rows = vcat([collect(CSV.File(joinpath(ROOT, "3_positive_negative_RTN",
    "cumulative_dv_RTN_branches_$(window)orbits.csv"))) for window in (80, 100)]...)
thresholds = collect(CSV.File(joinpath(ROOT, "continuous_access_thresholds.csv")))

@testset "Paired branches, styles, signs and separate scheduler outputs" begin
    @test length(PP.LINE_STYLES) == 36
    @test PP.PLOT_ALTITUDES == (1000, 1050, 1100)
    mktempdir() do output
        seen = Dict{Symbol,Tuple{Float64,Float64}}()
        for scheduler in OE.POLICIES, window in (80, 100)
            figure = PP.make_plot(rows, thresholds, scheduler, window, output)
            @test length(figure.subplots) == 10
            for (r, component) in enumerate(PP.COMPONENTS), (c, altitude) in enumerate(PP.PLOT_ALTITUDES)
                panel = figure.subplots[(r-1)*3+c]
                @test panel[:title] == (r == 1 ? "h_t = $altitude km" : "")
                limits = ylims(panel)
                @test limits == get!(seen, component, limits)
                @test limits == PP.component_limits(rows, component)
                @test limits[1] < 0 < limits[2]
                curves = filter(s -> s[:seriestype] == :path, panel.series_list)
                @test length(curves) == 6
                for (k, inclination) in enumerate(OE.INCLINATIONS), (b, branch) in enumerate(PP.BRANCHES)
                    curve = curves[(k-1)*2+b]
                    field = Symbol("dv_$(component)_$(branch)_mps")
                    expected = sort(filter(row -> row.orbits == window && row.scheduler == scheduler &&
                        row.target_altitude_km == altitude && row.delta_i_deg == inclination, rows);
                        by=row -> row.helpers)
                    @test collect(curve[:x]) == collect(OE.COUNTS)
                    @test collect(curve[:y]) == [getproperty(row, field) for row in expected]
                    @test all(y -> limits[1] <= y <= limits[2], curve[:y])
                    style = PP.LINE_STYLES[(scheduler, inclination, component, branch)]
                    @test curve[:linecolor] == Plots.plot_color(style.color)
                    @test curve[:linestyle] == style.linestyle
                    @test curve[:markercolor] == Plots.plot_color(style.markercolor)
                    @test curve[:markerstrokecolor] == Plots.plot_color(style.markerstrokecolor)
                end
                @test length(filter(s -> s[:seriestype] == :straightline, panel.series_list)) == 4
            end
            stem = "dv_RTN_paired_$(scheduler)_$(window)orbits_h_cols"
            for extension in ("png", "pdf")
                @test filesize(joinpath(output, "$stem.$extension")) > 10_000
            end
            records = CSV.File(joinpath(output, "$stem.csv"))
            @test length(records) == 81
            @test all(row -> row.scheduler == scheduler && row.orbits == window, records)
            @test Set(row.target_altitude_km for row in records) == Set(PP.PLOT_ALTITUDES)
        end
        styles = copy(PP.LINE_STYLES)
        key = ("gve_sma", 0.5, :R, :negative)
        styles[key] = merge(styles[key], (color="#ff6600", linestyle=:dashdot, markershape=:diamond))
        figure = PP.make_plot(rows, thresholds, "gve_sma", 100, output; styles)
        curves = filter(s -> s[:seriestype] == :path, figure.subplots[1].series_list)
        @test curves[4][:linecolor] == Plots.plot_color("#ff6600")
        @test curves[4][:linestyle] == :dashdot
        @test curves[4][:markershape] == :diamond
        legend_lines = filter(s -> s[:seriestype] == :path, figure.subplots[end].series_list)
        @test legend_lines[4][:linecolor] == Plots.plot_color("#ff6600")
        @test legend_lines[4][:linestyle] == :dashdot
        @test_throws ArgumentError PP.make_plot(rows, thresholds, "unknown", 80, output)
        retained = filter(row -> row.target_altitude_km in PP.PLOT_ALTITUDES, rows)
        @test_throws ErrorException PP.make_plot(retained[2:end], thresholds,
            first(retained).scheduler, first(retained).orbits, output)
    end
end
