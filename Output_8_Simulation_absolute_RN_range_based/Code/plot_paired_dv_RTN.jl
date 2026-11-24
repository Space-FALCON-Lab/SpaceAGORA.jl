#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..", "9_SpaceAGORA.jl-main_full integration")))
end

module PairedRTNPlots

# MANUAL APPEARANCE SETTINGS
# Each key is (scheduler, initial inclination in degrees, component, branch).
# Styles apply to that line in every altitude column and both orbit windows.
# Supported line styles include :solid, :dash, :dot, :dashdot.
const CASE_COLORS = (
    gve_sma=("#652884", "#4E5689", "#018B38"),
    gve_inc=("#B46DA9", "#41B9C1", "#57AF37"),
)
const LINE_STYLES = Dict(
    (scheduler, inclination, component, branch) => (
        color=getproperty(CASE_COLORS, Symbol(scheduler))[k],
        linestyle=scheduler == "gve_sma" ? :solid : :dot,
        linewidth=1.6, markershape=:circle, markersize=3.5,
        markercolor=branch == :positive ? getproperty(CASE_COLORS, Symbol(scheduler))[k] : "white",
        markerstrokecolor=getproperty(CASE_COLORS, Symbol(scheduler))[k],
        markerstrokewidth=0.9,
    )
    for scheduler in ("gve_sma", "gve_inc")
    for (k, inclination) in enumerate((0.0, 0.5, 1.0))
    for component in (:R, :T, :N) for branch in (:positive, :negative)
)
# To customize one line, uncomment and edit this example:
# LINE_STYLES[("gve_sma", 0.5, :R, :positive)] = merge(
#     LINE_STYLES[("gve_sma", 0.5, :R, :positive)],
#     (color="#ff6600", markercolor="#ff6600", markerstrokecolor="#ff6600", linestyle=:dash))
const PLOT_ALTITUDES = (1000, 1050, 1100)
const THRESHOLD_COLORS = ("#652884", "#4E5689", "#018B38")
const THRESHOLD_STYLE = :dash
const ZERO_LINE_COLOR = "#555555"

include("plot_orbital_elements.jl")
using .OrbitalElementPlots
using CSV, Plots
const OE = OrbitalElementPlots
const OUTPUT = joinpath(OE.StudyConfig.DEFAULT_OUTPUT, "images", "delta_v_RTN",
    "8_paired_RTN_by_scheduler")
const COMPONENTS = (:R, :T, :N)
const BRANCHES = (:positive, :negative)
const FIELDS = Tuple(Symbol("dv_$(component)_$(branch)_mps")
    for component in COMPONENTS for branch in BRANCHES)
const CSV_FIELDS = (:scheduler, :target_altitude_km, :delta_i_deg, :helpers,
    :orbits, :time_s, :row, :source_trajectory, FIELDS...)

function component_limits(rows, component)
    values = [getproperty(row, Symbol("dv_$(component)_$(branch)_mps"))
        for row in rows for branch in BRANCHES]
    OE.dv_ylimits(values, Symbol("dv_$(component)_mps"))
end

function make_plot(rows, thresholds, scheduler, window, output; styles=LINE_STYLES)
    scheduler in OE.POLICIES || throw(ArgumentError("Unsupported scheduler: $scheduler"))
    !isempty(PLOT_ALTITUDES) && length(unique(PLOT_ALTITUDES)) == length(PLOT_ALTITUDES) &&
        all(in(OE.ALTITUDES), PLOT_ALTITUDES) || error("Invalid PLOT_ALTITUDES")
    selected = filter(r -> r.scheduler == scheduler && r.orbits == window &&
        r.target_altitude_km in PLOT_ALTITUDES, rows)
    expected = Set((h, i, n) for h in PLOT_ALTITUDES for i in OE.INCLINATIONS for n in OE.COUNTS)
    length(selected) == length(expected) &&
        Set((r.target_altitude_km, r.delta_i_deg, r.helpers) for r in selected) == expected ||
        error("Incomplete or duplicate grid for $scheduler, $window orbits")
    panels = Plots.Plot[]
    for (r, component) in enumerate(COMPONENTS), (c, altitude) in enumerate(PLOT_ALTITUDES)
        limits = component_limits(rows, component)
        panel = plot(; legend=false, grid=true, gridalpha=0.18, framestyle=:axes,
            title=r == 1 ? "h_t = $altitude km" : "", titlefontsize=12,
            ylabel=c == 1 ? "\u0394v_$(component)+ / \u0394v_$(component)- (m/s)" : "",
            xlabel=r == 3 ? "N_helpers" : "", ylims=limits,
            xlims=(-8, maximum(OE.COUNTS) + 15), xticks=collect(OE.COUNTS), xrotation=60,
            guidefontsize=12, tickfontsize=10, tick_direction=:out,
            left_margin=3Plots.mm, right_margin=2Plots.mm,
            top_margin=2Plots.mm, bottom_margin=3Plots.mm)
        hline!(panel, [0.0]; color=ZERO_LINE_COLOR, linewidth=0.8, label="")
        for (k, inclination) in enumerate(OE.INCLINATIONS)
            threshold = only(filter(t -> t.target_altitude_km == altitude &&
                t.delta_i_deg == inclination, thresholds))
            if !ismissing(threshold.minimum_helpers) && threshold.minimum_helpers <= maximum(OE.COUNTS)
                vline!(panel, [threshold.minimum_helpers]; color=THRESHOLD_COLORS[k],
                    linestyle=THRESHOLD_STYLE, linewidth=1.2, alpha=0.75, label="")
            elseif r == 1
                label = ismissing(threshold.minimum_helpers) ? "No finite N_cont" :
                    "N_cont = $(threshold.minimum_helpers) > $(maximum(OE.COUNTS))"
                annotate!(panel, maximum(OE.COUNTS), limits[2] - k * 0.08 * (limits[2] - limits[1]),
                    text("$inclination deg: $label", 9, THRESHOLD_COLORS[k], :right))
            end
            curve = sort(filter(row -> row.target_altitude_km == altitude &&
                row.delta_i_deg == inclination, selected); by=row -> row.helpers)
            for branch in BRANCHES
                field = Symbol("dv_$(component)_$(branch)_mps")
                values = [getproperty(row, field) for row in curve]
                all(isfinite, values) || error("Nonfinite $field for $scheduler, $altitude km")
                plot!(panel, [row.helpers for row in curve], values;
                    styles[(scheduler, inclination, component, branch)]..., label="")
            end
        end
        push!(panels, panel)
    end
    footer = plot(; framestyle=:none, ticks=nothing, legend=false, grid=false,
        xlims=(0, 1), ylims=(0, 1), margin=0Plots.mm)
    for (column, component) in enumerate(COMPONENTS)
        x = 0.03 + (column - 1) * 0.33
        annotate!(footer, x, 0.96, text("$component branches", 11, :left))
        for (k, inclination) in enumerate(OE.INCLINATIONS), (b, branch) in enumerate(BRANCHES)
            y = 0.84 - ((k - 1) * 2 + b - 1) * 0.095
            style = styles[(scheduler, inclination, component, branch)]
            plot!(footer, [x, x + 0.035], [y, y]; style..., label="")
            sign = branch == :positive ? "+" : "-"
            annotate!(footer, x + 0.045, y, text("i = $inclination deg, $component$sign", 10, :left))
        end
    end
    for (k, inclination) in enumerate(OE.INCLINATIONS)
        x = 0.03 + (k - 1) * 0.33
        plot!(footer, [x, x + 0.035], [0.22, 0.22]; color=THRESHOLD_COLORS[k],
            linestyle=THRESHOLD_STYLE, linewidth=1.2, label="")
        annotate!(footer, x + 0.045, 0.22, text("N_cont, i = $inclination deg", 10, :left))
    end
    annotate!(footer, 0.5, 0.105,
        text("$scheduler | $window initial target periods | N_cont = ceil(\u03c0 / \u03b1_min), fixed circular range-only estimate", 10))
    annotate!(footer, 0.5, 0.01,
        text("Separate laser-only positive and negative integrals in instantaneous RTN; NOT net delta-v. Negative signs retained.", 10))
    ncols = length(PLOT_ALTITUDES)
    layout = @layout [grid(3, ncols); footer{0.25h}]
    figure = plot(panels..., footer; layout, size=(1600, 1400), dpi=160,
        fontfamily="Computer Modern", left_margin=10Plots.mm, right_margin=3Plots.mm)
    plot!(figure; subplot=length(panels) + 1, top_margin=10Plots.mm)
    mkpath(output)
    stem = "dv_RTN_paired_$(scheduler)_$(window)orbits_h_cols"
    for extension in ("png", "pdf")
        path = joinpath(output, "$stem.$extension")
        savefig(figure, path)
        println("Saved $path")
    end
    CSV.write(joinpath(output, "$stem.csv"),
        [NamedTuple{CSV_FIELDS}(Tuple(getproperty(row, field) for field in CSV_FIELDS)) for row in selected])
    figure
end

function render_all(rows, thresholds, windows, output=OUTPUT)
    mkpath(output)
    selected_thresholds = filter(t -> t.target_altitude_km in PLOT_ALTITUDES, thresholds)
    length(selected_thresholds) == length(PLOT_ALTITUDES) * length(OE.INCLINATIONS) ||
        error("Incomplete continuous-access thresholds")
    gr()
    for scheduler in OE.POLICIES, window in windows
        make_plot(rows, selected_thresholds, scheduler, window, output)
    end
    CSV.write(joinpath(output, "continuous_access_thresholds.csv"), selected_thresholds)
    nothing
end

function main(args=ARGS)
    options = OE.parse_options(vcat(["--output", OUTPUT], args))
    rows, thresholds = OE.load_grid(options.input, options.windows; include_rtn=true)
    render_all(rows, thresholds, options.windows, options.output)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    PairedRTNPlots.main()
end
