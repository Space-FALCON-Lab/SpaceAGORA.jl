import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "9_SpaceAGORA.jl-main_full integration")))
using Test, Plots
include(joinpath(@__DIR__, "..", "Simulation", "Code", "plot_pooled_orbital_elements.jl"))
const P = PooledOrbitalElementPlots
const OE = P.OE

rows = [(scheduler=policy, target_altitude_km=altitude, delta_i_deg=inclination,
    helpers=helpers, orbits=window, da_m=Float64(helpers + altitude + window),
    de=helpers * 1e-6, di_deg=inclination * 1e-3,
    draan_deg=inclination == 0 ? NaN : -inclination * window,
    raan_defined=inclination != 0)
    for window in (80, 100) for policy in OE.POLICIES for altitude in OE.ALTITUDES
    for inclination in OE.INCLINATIONS for helpers in OE.COUNTS]

@testset "Pooled scatter grouping and colors" begin
    gr()
    for window in (80, 100), xfield in (:helpers, :target_altitude_km)
        figure = P.make_plot(rows, window, xfield)
        @test length(figure.subplots) == 4
        groupfield = xfield == :helpers ? :scheduler : :delta_i_deg
        groups = xfield == :helpers ? ("gve_sma", "gve_inc") : OE.INCLINATIONS
        colors = xfield == :helpers ? (:red, :blue) : (:orange, :green, :purple)
        for (index, field) in enumerate(P.FIELDS)
            panel = figure.subplots[index]
            @test all(s -> s[:seriestype] == :scatter, panel.series_list)
            @test sum(length(s[:x]) for s in panel.series_list) == (index == 4 ? 308 : 462)
            expected_groups = index == 4 && xfield == :target_altitude_km ? 2 : length(groups)
            @test length(panel.series_list) == expected_groups
            series_index = 0
            for (group, color) in zip(groups, colors)
                points = filter(row -> row.orbits == window &&
                    getproperty(row, groupfield) == group &&
                    isfinite(getproperty(row, field)), rows)
                isempty(points) && continue
                series_index += 1
                series = panel.series_list[series_index]
                @test collect(series[:x]) == [getproperty(row, xfield) for row in points]
                @test collect(series[:y]) == [getproperty(row, field) for row in points]
                @test series[:markercolor] == Plots.plot_color(color)
            end
        end
    end
    @test_throws ArgumentError P.make_plot(rows, 80, :unknown)
    @test_throws ErrorException P.make_plot(rows[2:end], 80, :helpers)
    @test_throws ErrorException P.make_plot(vcat(rows, [first(rows)]), 80, :helpers)
    @test_throws ErrorException P.make_plot(rows, 90, :helpers)
end
