#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

module PooledOrbitalElementPlots

include("plot_orbital_elements.jl")
using .OrbitalElementPlots
using Plots

const OE = OrbitalElementPlots
const OUTPUT = joinpath(OE.OUTPUT, "2_separate_plots_by_N-sceduler_h-i")
const FIELDS = (:da_m, :de, :di_deg, :draan_deg)
const LABELS = ("\u0394a (m)", "\u0394e", "\u0394i (deg)", "\u0394\u03a9 (deg)")

function make_plot(rows, window, xfield)
    if xfield == :helpers
        groupfield = :scheduler
        groups = ("gve_sma", "gve_inc")
        colors = (:red, :blue)
        legends = ("SMA scheduler", "Inclination scheduler")
        xlabel = "Number of helpers"
        ticks = collect(OE.COUNTS)
        limits = (-15, maximum(OE.COUNTS) + 25)
        title = "OE vs number of helpers"
        subtitle = "All target altitudes and initial inclinations pooled"
    elseif xfield == :target_altitude_km
        groupfield = :delta_i_deg
        groups = OE.INCLINATIONS
        colors = (:orange, :green, :purple)
        legends = ("Initial inclination = 0 deg", "Initial inclination = 0.5 deg",
            "Initial inclination = 1 deg")
        xlabel = "Target altitude (km)"
        ticks = sort(collect(OE.ALTITUDES))
        limits = (minimum(ticks) - 10, maximum(ticks) + 10)
        title = "OE vs target altitude"
        subtitle = "All helper counts and both schedulers pooled"
    else
        throw(ArgumentError("Unsupported x-axis: $xfield"))
    end
    selected = filter(row -> row.orbits == window, rows)
    expected = Set((policy, altitude, inclination, helpers)
        for policy in OE.POLICIES for altitude in OE.ALTITUDES
        for inclination in OE.INCLINATIONS for helpers in OE.COUNTS)
    actual = Set((row.scheduler, row.target_altitude_km, row.delta_i_deg, row.helpers)
        for row in selected)
    length(selected) == OE.CASE_COUNT && actual == expected ||
        error("Expected the complete $(OE.CASE_COUNT)-case grid for $window orbits")
    all(row -> all(isfinite, (row.da_m, row.de, row.di_deg)) &&
        (row.delta_i_deg == 0 ? !row.raan_defined && isnan(row.draan_deg) :
            row.raan_defined && isfinite(row.draan_deg)), selected) ||
        error("Invalid orbital-element values for $window orbits")

    panels = Plots.Plot[]
    for (index, field) in enumerate(FIELDS)
        panel = plot(; title=index == 1 ? "$title - $window orbits\n$subtitle" :
                index == 4 ? "Unwrapped RAAN; equatorial cases omitted (undefined)" : "",
            titlefontsize=index == 4 ? 9 : 13, ylabel=LABELS[index],
            xlabel=index == 4 ? xlabel : "", xlims=limits, xticks=ticks,
            xrotation=xfield == :helpers ? 45 : 0,
            ylims=OE.ylimits([getproperty(row, field) for row in rows], field),
            legend=index == 1 ? :outertop : false, legend_columns=length(groups),
            legendfontsize=9, legend_foreground_color=nothing,
            framestyle=:axes, grid=true, gridalpha=0.18, tick_direction=:out,
            tickfontsize=10, guidefontsize=12)
        for (group, color, label) in zip(groups, colors, legends)
            points = filter(row -> getproperty(row, groupfield) == group &&
                isfinite(getproperty(row, field)), selected)
            isempty(points) && continue
            scatter!(panel, [getproperty(row, xfield) for row in points],
                [getproperty(row, field) for row in points];
                color, label, markershape=:circle, markersize=3,
                markerstrokewidth=0, markeralpha=0.65)
        end
        push!(panels, panel)
    end
    plot(panels...; layout=(4, 1), size=(1000, 1300), dpi=160,
        fontfamily="Computer Modern", left_margin=8Plots.mm, right_margin=5Plots.mm,
        top_margin=5Plots.mm, bottom_margin=5Plots.mm)
end

function main(args=ARGS)
    options = OE.parse_options(vcat(["--output", OUTPUT], args))
    rows, _ = OE.load_grid(options.input, options.windows)
    mkpath(options.output)
    gr()
    for window in options.windows, xfield in (:helpers, :target_altitude_km)
        figure = make_plot(rows, window, xfield)
        stem = xfield == :helpers ? "OE_vs_helpers_by_scheduler" : "OE_vs_target_altitude_by_inclination"
        for extension in ("png", "pdf")
            path = joinpath(options.output, "$(stem)_$(window)orbits.$extension")
            savefig(figure, path)
            println("Saved $path")
        end
    end
    println("Plotted $(length(rows)) endpoints from $(OE.CASE_COUNT) verified trajectories; no averaging or jitter.")
    nothing
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    PooledOrbitalElementPlots.main()
end
