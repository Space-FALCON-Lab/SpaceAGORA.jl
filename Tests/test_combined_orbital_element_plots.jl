import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "9_SpaceAGORA.jl-main_full integration")))
using Test, CSV, Plots
include(joinpath(@__DIR__, "..", "Simulation", "Code", "plot_orbital_elements.jl"))
const OE = OrbitalElementPlots

@testset "Combined scheduler figures from saved endpoints" begin
    gr()
    thresholds = NamedTuple.(CSV.File(joinpath(OE.OUTPUT, "continuous_access_thresholds.csv")))
    colors = ("#652884", "#4E5689", "#018B38", "#B46DA9", "#41B9C1", "#57AF37")
    fields = (:da_m, :de, :di_deg, :draan_deg)
    mktempdir() do output
        for window in (80, 100)
            rows = NamedTuple.(CSV.File(joinpath(OE.OUTPUT, "orbital_element_deltas_$(window)orbits.csv")))
            @test length(rows) == 462
            figure = OE.make_plot(rows, thresholds, "combined", window, output)
            @test length(figure.subplots) == 29
            for (index, panel) in enumerate(figure.subplots[1:28])
                @test panel[:xaxis][:grid]
                @test panel[:yaxis][:grid]
                row, column = divrem(index - 1, 7) .+ 1
                field = fields[row]
                curves = filter(s -> s[:seriestype] == :path, panel.series_list)
                @test length(curves) == (row == 4 ? 4 : 6)
                curve_index = 0
                for (j, policy) in enumerate(("gve_sma", "gve_inc")),
                        (k, inclination) in enumerate(OE.INCLINATIONS)
                    row == 4 && inclination == 0 && continue
                    curve_index += 1
                    expected = sort(filter(r -> r.scheduler == policy &&
                        r.delta_i_deg == inclination &&
                        r.target_altitude_km == OE.ALTITUDES[column], rows); by=r -> r.helpers)
                    series = curves[curve_index]
                    @test series[:linecolor] == Plots.plot_color(colors[(j - 1) * 3 + k])
                    @test series[:markercolor] == Plots.plot_color(colors[(j - 1) * 3 + k])
                    @test series[:linestyle] == (policy == "gve_inc" ? :dot : :solid)
                    @test collect(series[:x]) == collect(OE.COUNTS)
                    @test collect(series[:y]) == [getproperty(r, field) for r in expected]
                end
                lines = filter(s -> s[:seriestype] == :straightline, panel.series_list)
                expected_lines = sort([t.minimum_helpers for t in thresholds
                    if t.target_altitude_km == OE.ALTITUDES[column] && t.minimum_helpers <= 800])
                @test sort([first(s[:x]) for s in lines]) == expected_lines
            end
            legend_curves = filter(s -> s[:seriestype] == :path, figure.subplots[end].series_list)
            @test [s[:linestyle] for s in legend_curves[1:6]] ==
                [:solid, :solid, :solid, :dot, :dot, :dot]
            for extension in ("png", "pdf")
                path = joinpath(output, "OE_from_feather_combined_e0.0000_nu1.0000_$(window)orbits_h_cols.$extension")
                @test filesize(path) > 10_000
            end
            @test_throws ErrorException OE.make_plot(rows[2:end], thresholds, "combined", window, output)
        end
    end
end
