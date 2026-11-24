#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

module PaperPlots

# EDIT COLORS AND LINE STYLES HERE: entries correspond to i = 0, 0.5, 1 degree.
const LINE_COLORS = (
    # gve_sma=("#652884", "#D9A421", "#018B38"),
    # gve_inc=("#B46DA9", "#F5A216", "#57AF37"),
    gve_sma=("#D9A421", "#652884", "#018B38"),
    gve_inc=("#F5A216", "#B46DA9", "#57AF37"),
    # gve_sma=("#0076B9", "#D9A421", "#018B38"),
    # gve_inc=("#699ECA", "#F5A216", "#57AF37"),
)
const LINE_STYLES = (
    gve_sma=(:solid, :solid, :solid),
    gve_inc=(:dot, :dot, :dot),
)
const LINE_WIDTH = 1.0
const MARKER_SIZE = 1.8
const ACCESS_LINE_WIDTH = 0.8
const FONT_FAMILY = "Computer Modern"
const FONT_SIZE = 8
const TICK_FONT_SIZE = FONT_SIZE - 1
const GUIDE_FONT_SIZE = FONT_SIZE + 1
const X_TICK_ROTATION = 0
const LABELED_HELPERS = (1, 200, 400, 600, 800)

# EDIT SIZES HERE (Plots canvas units, before DPI scaling).
# With figure_size=nothing, width = 3*subplot_width and
# height = 4*subplot_height + footer_height. Subplot sizes include axes/margins.
# Set figure_size=(width, height) to override the overall size; this rescales
# the subplot regions, so their absolute dimensions then follow the override.
const OE_SIZE = (subplot_width=170, subplot_height=135, footer_height=46, figure_size=nothing)
const DV_SIZE = (subplot_width=170, subplot_height=135, footer_height=32, figure_size=nothing)
const DPI = 300
# These are the final margins; no later figure-level margins override them.
# Increase y_label_x_mm to move that figure's y-axis titles right.
# The -1 margins partly offset GR's built-in 2 mm padding.
const OE_SPACING = (
    panel_margins_mm=(left=1, right=1, top=1, bottom=1),
    y_label_x_mm=5,
    y_label_margin_mm=7,
    bottom_row_margin_mm=-1,
    legend_top_margin_mm=-1,
)
const DV_SPACING = (
    panel_margins_mm=(left=1, right=1, top=1, bottom=1),
    y_label_x_mm=5,
    y_label_margin_mm=7,
    bottom_row_margin_mm=-1,
    legend_top_margin_mm=-1,
)
# Offsets are fractions of the legend footer: negative x moves left, positive y up.
# row_spacing multiplies the gaps between rows, keeping the first row anchored.
const OE_LEGEND = (x_offset=-0.05, y_offset=0.0, row_spacing=1.2)
const DV_LEGEND = (x_offset=-0.05, y_offset=0.0, row_spacing=1.2)
const SHOW_OE_ACCESS = true
const SHOW_DV_ACCESS = false

# Source units (m, dimensionless, deg, m/s); RAAN is fitted to the 100-orbit CSV.
const Y_LIMITS = (
    da_m=(-5000.0, 20000.0), de=(0.0, 0.001), di_deg=(-0.005, 0.025),
    dv_R_mps=(-2, 8), dv_T_mps=(-2.0, 10.0),
    dv_N_positive_mps=(-1.0, 3.0), dv_N_negative_mps=(-3.0, 1.0),
)

using CSV
get!(ENV, "GKSwstype", "100")
using Plots

const ALTITUDES = (1000, 1050, 1100)
const INCLINATIONS = (0.0, 0.5, 1.0)
const HELPERS = (1, 100, 200, 300, 400, 500, 600, 700, 800)
const POLICIES = ("gve_sma", "gve_inc")
const SPECS = (
    oe=(csv="orbital_element_deltas_100orbits.csv",
        stem="OE_from_feather_combined_e0.0000_nu1.0000_100orbits_h_cols",
        fields=(:da_m, :de, :di_deg, :draan_deg),
        scales=(1e-3, 1e4, 1e3, 1.0),
        labels=("\$\\Delta a\\, (\\mathrm{km})\$", "\$\\Delta e\\quad (\\,\\!\\times10^{-4})\$",
            "\$\\Delta i\\quad (\\,\\!\\times10^{-3}\\,\\mathrm{deg})\$",
            "\$\\Delta\\Omega\\; (\\mathrm{deg})\$")),
    dv=(csv="cumulative_dv_RT_N_branches_100orbits.csv",
        stem="dv_RTN_split_N_from_feather_combined_e0.0000_nu1.0000_100orbits_h_cols",
        fields=(:dv_R_mps, :dv_T_mps, :dv_N_positive_mps, :dv_N_negative_mps),
        scales=(1.0, 1.0, 1.0, 1.0),
        labels=("\$\\Delta v_R\\; (\\mathrm{m/s})\$", "\$\\Delta v_T\\; (\\mathrm{m/s})\$",
            "\$\\Delta v_{N+}\\; (\\mathrm{m/s})\$", "\$\\Delta v_{N-}\\; (\\mathrm{m/s})\$")),
)

function canvas_geometry(settings)
    all(x -> isfinite(x) && x > 0,
        (settings.subplot_width, settings.subplot_height, settings.footer_height)) ||
        error("Subplot dimensions and footer height must be positive and finite")
    grid_height = 4 * settings.subplot_height
    automatic = (3 * settings.subplot_width, grid_height + settings.footer_height)
    size = isnothing(settings.figure_size) ? automatic : settings.figure_size
    length(size) == 2 && all(x -> isfinite(x) && x > 0, size) ||
        error("figure_size must contain two positive, finite dimensions")
    (; size, footer_fraction=settings.footer_height / automatic[2])
end

function raan_limits(rows)
    values = [row.draan_deg for row in rows if row.delta_i_deg != 0]
    !isempty(values) && all(isfinite, values) ||
        error("Expected finite non-equatorial RAAN changes for the 100-orbit range")
    lo, hi = extrema(values)
    if lo == hi
        padding = max(abs(lo)*0.05, 0.1)
        lo, hi = lo-padding, hi+padding
    end
    _, lower, upper = Plots.PlotUtils.optimize_ticks(lo, hi; strict_span=false, k_ideal=5)
    (lower, upper)
end

# Adapted from Simulation_absolute_RN/Code/plot_orbital_elements.jl.
# Rendering is local; this script does not import or modify the original code.
function make_plot(kind; directory=@__DIR__, output=directory,
        settings=kind == :oe ? OE_SIZE : DV_SIZE,
        spacing=kind == :oe ? OE_SPACING : DV_SPACING,
        legend_settings=kind == :oe ? OE_LEGEND : DV_LEGEND,
        show_access=kind == :oe ? SHOW_OE_ACCESS : SHOW_DV_ACCESS,
        colors=LINE_COLORS, styles=LINE_STYLES, policies=POLICIES,
        policy_labels=(gve_sma="\$\\Delta a\$-max", gve_inc="\$\\Delta i\$-max"),
        specification=nothing, orbit_window=100, endpoint_row=1001, source_limits=nothing,
        averaging_orbits::Int=0, rows=nothing, prediction_curve=nothing)
    kind in (:oe, :dv) || throw(ArgumentError("Expected :oe or :dv"))
    all(isfinite, (legend_settings.x_offset, legend_settings.y_offset,
        legend_settings.row_spacing)) && legend_settings.row_spacing > 0 ||
        error("Legend offsets must be finite and row_spacing must be positive and finite")
    spec = specification === nothing ? getproperty(SPECS, kind) : specification
    length(spec.fields) == length(spec.scales) == length(spec.labels) == 4 ||
        error("Paper figures require four quantity rows")
    !isempty(policies) && length(unique(policies)) == length(policies) && length(policies) <= 2 ||
        error("Paper figures require one or two distinct schedulers")
    rows = isnothing(rows) ? collect(CSV.File(joinpath(directory, spec.csv))) : rows
    expected = Set((p, h, i, n) for p in policies for h in ALTITUDES
        for i in INCLINATIONS for n in HELPERS)
    length(rows) == length(expected) &&
        Set((r.scheduler, r.target_altitude_km, r.delta_i_deg, r.helpers) for r in rows) == expected ||
        error("Expected the complete $(length(expected))-case grid in $(spec.csv)")
    all(r -> r.orbits == orbit_window && r.row == endpoint_row, rows) ||
        error("Expected $orbit_window-orbit endpoints at row $endpoint_row")
    averaging_orbits in (0, 2) || error("Only endpoint or last-two-orbit figures are supported")
    if averaging_orbits > 0
        all(r -> r.aggregation == "time_weighted_mean" &&
            r.averaging_orbits == averaging_orbits &&
            r.window_start_orbits == orbit_window - averaging_orbits &&
            r.window_end_orbits == orbit_window && r.window_end_row == endpoint_row &&
            r.window_start_row == round(Int, 10(orbit_window - averaging_orbits)) + 1 &&
            r.sample_count == r.window_end_row - r.window_start_row + 1 &&
            r.window_end_time_s > r.window_start_time_s, rows) ||
            error("Invalid last-two-orbit averaging metadata")
    end
    ranges = source_limits === nothing ?
        (kind == :oe ? merge(Y_LIMITS, (draan_deg=raan_limits(rows),)) : Y_LIMITS) : source_limits
    thresholds = show_access ? collect(CSV.File(joinpath(directory, "continuous_access_thresholds.csv"))) : []
    simulation_marker(policy) = !isnothing(prediction_curve) && policy == "gve_inc" ? :utriangle : :circle
    panels = Plots.Plot[]
    for (r, field) in enumerate(spec.fields), (c, altitude) in enumerate(ALTITUDES)
        scale = spec.scales[r]
        limits = getproperty(ranges, field) .* scale
        panel = plot(; title=r == 1 ? "Target Alt. = $altitude km" : "",
            titlefontsize=FONT_SIZE, titlelocation=:center,
            legend=false, grid=true, gridalpha=0.18, framestyle=:axes,
            ylabel="",
            xlabel=r == 4 ? "\$N_{\\mathrm{helpers}}\$" : "",
            xlims=(-8, maximum(HELPERS) + 15),
            xticks=(collect(HELPERS), [h in LABELED_HELPERS ? string(h) : "" for h in HELPERS]),
            xrotation=X_TICK_ROTATION, ylims=limits, yformatter=:plain,
            tickfontsize=TICK_FONT_SIZE,
            guidefontsize=GUIDE_FONT_SIZE, tick_direction=:out,
            left_margin=(spacing.panel_margins_mm.left + (c == 1 ? spacing.y_label_margin_mm : 0))*Plots.mm,
            right_margin=spacing.panel_margins_mm.right*Plots.mm,
            top_margin=spacing.panel_margins_mm.top*Plots.mm,
            bottom_margin=(r == 4 ? spacing.bottom_row_margin_mm : spacing.panel_margins_mm.bottom)*Plots.mm)
        for policy in policies, (k, inclination) in enumerate(INCLINATIONS)
            curve = sort(filter(row -> row.scheduler == policy &&
                row.target_altitude_km == altitude && row.delta_i_deg == inclination, rows);
                by=row -> row.helpers)
            values = [getproperty(row, field)*scale for row in curve]
            if field == :draan_deg && inclination == 0
                all(isnan, values) || error("Equatorial RAAN must be undefined")
                continue
            end
            all(v -> isfinite(v) && limits[1] <= v <= limits[2], values) ||
                error("Invalid or clipped data for $field, $policy, $altitude km, $inclination deg")
            color = getproperty(colors, Symbol(policy))[k]
            linestyle = getproperty(styles, Symbol(policy))[k]
            plot!(panel, [row.helpers for row in curve], values; color, linestyle,
                seriestype=isnothing(prediction_curve) ? :path : :scatter,
                lw=LINE_WIDTH, marker=simulation_marker(policy), ms=MARKER_SIZE, markerstrokewidth=0)
            if !isnothing(prediction_curve)
                prediction = prediction_curve(field, policy, altitude, inclination)
                if !isnothing(prediction)
                    helpers, source_values = prediction
                    predicted = source_values .* scale
                    all(v -> isnan(v) || limits[1] <= v <= limits[2], predicted) ||
                        error("Prediction outside $field axis: $policy, $altitude km, $inclination deg")
                    plot!(panel, helpers, predicted; color, linestyle, lw=LINE_WIDTH, marker=:none)
                end
            end
        end
        if show_access
            for (k, inclination) in enumerate(INCLINATIONS)
                threshold = only(filter(t -> t.target_altitude_km == altitude &&
                    t.delta_i_deg == inclination, thresholds))
                !ismissing(threshold.minimum_helpers) && threshold.minimum_helpers <= maximum(HELPERS) ||
                    error("Continuous-access threshold is outside this paper figure's helper range")
                vline!(panel, [threshold.minimum_helpers];
                    color=getproperty(colors, Symbol(first(policies)))[k],
                    linestyle=:dash, lw=ACCESS_LINE_WIDTH, alpha=0.75)
            end
        end
        push!(panels, panel)
    end
    footer = plot(; title="", framestyle=:none, ticks=nothing, legend=false, grid=false,
        xlims=(0, 1), ylims=(0, 1), margin=1Plots.mm,
        top_margin=spacing.legend_top_margin_mm*Plots.mm, bottom_margin=0Plots.mm)
    base_y = if isnothing(prediction_curve)
        show_access ? (0.85, 0.85-0.33, 0.17) : (0.78, 0.78-0.50)
    else
        row_count = 2*length(policies) + show_access
        extra_rows = 2*(length(policies)-1)
        footer_height = settings.footer_height + 16*extra_rows + 5*row_count
        settings = merge(settings, (; footer_height))
        [1 - (r-0.5)/(row_count+1) for r in 1:row_count]
    end
    row_y = map(y -> first(base_y) + legend_settings.y_offset +
        (y-first(base_y))*legend_settings.row_spacing, base_y)
    geometry = canvas_geometry(settings)
    layout = grid(2, 1; heights=[1-geometry.footer_fraction, geometry.footer_fraction])
    layout[1, 1] = grid(4, 3)
    average_title = averaging_orbits == 0 ? (;) :
        (; plot_title="Time-weighted mean: reference orbits $(orbit_window - averaging_orbits)-$orbit_window",
            plot_titlefontsize=FONT_SIZE + 1)
    canvas_size = averaging_orbits == 0 ? geometry.size :
        (geometry.size[1] + 90, geometry.size[2])
    figure = plot(panels..., footer; layout, size=canvas_size, dpi=DPI,
        fontfamily=FONT_FAMILY, average_title...)
    footer = last(figure.subplots)
    left = isnothing(prediction_curve) ? 0.06 + legend_settings.x_offset : 0.06
    if isnothing(prediction_curve)
        for (j, policy) in enumerate(policies), (k, inclination) in enumerate(INCLINATIONS)
            x = 0.06 + legend_settings.x_offset + (k-1)/3
            y = row_y[j]
            plot!(footer, [x, x+0.065], [y, y];
                color=getproperty(colors, Symbol(policy))[k],
                linestyle=getproperty(styles, Symbol(policy))[k], lw=LINE_WIDTH,
                marker=simulation_marker(policy), ms=MARKER_SIZE, markerstrokewidth=0)
            label = "$(getproperty(policy_labels, Symbol(policy))), $inclination\u00b0"
            annotate!(footer, x+0.08, y, text(label, FONT_SIZE, :left, FONT_FAMILY))
        end
    else
        # The grid footer inherits the y-axis gutter. Keep its allocated size and
        # height, but anchor the actual legend directly to the full canvas.
        Plots.prepare_output(figure)
        area = Plots.plotarea(footer)
        inset = Plots.bbox(legend_settings.x_offset*Plots.width(area), Plots.top(area),
            Plots.width(area), Plots.height(area); h_anchor=:center)
        plot!(figure; inset_subplots=inset, subplot=length(figure.subplots)+1,
            framestyle=:none, ticks=nothing, legend=false, grid=false, margin=0Plots.mm)
        footer = last(figure.subplots)
        pitch = row_y[1] - row_y[2]
        top = first(row_y) + 0.6*pitch
        bottom = last(row_y) - (show_access ? 0.95 : 0.6)*pitch
        x0, x1 = left - 0.012, left + 0.98
        source_edge, policy_edge = left + 0.146, left + 0.292
        boundary(i) = (row_y[i] + row_y[i+1])/2
        sim_bottom = boundary(length(policies))
        policy_bottom = show_access ? boundary(2*length(policies)) : bottom
        section_edges = (top, sim_bottom, policy_bottom)
        cell_width = (x1-policy_edge)/length(INCLINATIONS)
        function draw_key!(k, y, color; marker=:none, linestyle=:solid, lw=LINE_WIDTH, alpha=1)
            # Reserve the same sample slot in every row so the angle labels align.
            inset = 0.13
            sample_width = 0.048
            sample_left = policy_edge + (k-1+inset)*cell_width
            # GR's right anchor includes the degree glyph's trailing side bearing.
            text_right = policy_edge + (k-inset)*cell_width + 0.004
            if marker == :none
                plot!(footer, [sample_left, sample_left+sample_width], [y, y]; color, linestyle, lw, alpha, marker=:none)
            else
                scatter!(footer, [sample_left+sample_width/2], [y]; color, marker, ms=MARKER_SIZE, markerstrokewidth=0)
            end
            label = "\$\\Delta i=$(INCLINATIONS[k])^\\circ\$"
            annotate!(footer, text_right, y, text(label, FONT_SIZE, :right, FONT_FAMILY))
        end
        for (i, source) in enumerate(("Simulation", "Prediction"))
            upper, lower = section_edges[i], section_edges[i+1]
            annotate!(footer, (x0+source_edge)/2, (upper+lower)/2,
                text(source, FONT_SIZE, :center, FONT_FAMILY))
            for (j, policy) in enumerate(policies)
                y = upper - (j-0.5)*(upper-lower)/length(policies)
                key = Symbol(policy)
                annotate!(footer, (source_edge+policy_edge)/2, y,
                    text(policy_labels[key], FONT_SIZE, :center, FONT_FAMILY))
                for k in eachindex(INCLINATIONS)
                    draw_key!(k, y, colors[key][k]; linestyle=styles[key][k],
                        marker=i == 1 ? simulation_marker(policy) : :none)
                end
            end
        end
        if show_access
            y = (policy_bottom+bottom)/2
            # The mixed text/math label has the same right bearing, at half weight.
            annotate!(footer, (x0+policy_edge)/2 + 0.002, y,
                text("Threshold, \$N_{\\mathrm{cont}}\$", FONT_SIZE, :center, FONT_FAMILY))
            for k in eachindex(INCLINATIONS)
                draw_key!(k, y, colors[Symbol(first(policies))][k]; linestyle=:dash, lw=ACCESS_LINE_WIDTH, alpha=0.75)
            end
        end
    end
    if show_access && isnothing(prediction_curve)
        y = last(row_y)
        for (k, inclination) in enumerate(INCLINATIONS)
            x = left + (k-1)/3
            plot!(footer, [x, x+0.065], [y, y];
                color=getproperty(colors, Symbol(first(policies)))[k],
                linestyle=:dash, lw=ACCESS_LINE_WIDTH, alpha=0.75)
            label = "\$N_{\\mathrm{cont}},\\; $inclination^\\circ\$"
            annotate!(footer, x+0.08, y, text(label, FONT_SIZE, :left, FONT_FAMILY))
        end
    end
    if !isnothing(prediction_curve)
        xlims!(footer, (x0-0.007, x1+0.007))
        ylims!(footer, (min(0, bottom-0.01), max(1, top+0.01)))
        rule!(xs, ys) = plot!(footer, xs, ys; color=:black, lw=0.6, marker=:none)
        rule!([x0, x1, x1, x0, x0], [top, top, bottom, bottom, top])
        rule!([x0, x1], fill(sim_bottom, 2))
        # The threshold label spans the source and policy columns.
        show_access && rule!([x0, x1], fill(policy_bottom, 2))
        rule!(fill(source_edge, 2), [top, policy_bottom])
        rule!(fill(policy_edge, 2), [top, bottom])
    end
    # Anchor titles to the same canvas coordinate, not each row's tick-label width.
    Plots.prepare_output(figure)
    for (r, label) in enumerate(spec.labels)
        panel = figure.subplots[(r-1)*3+1]
        area = Plots.plotarea(panel)
        x = (spacing.y_label_x_mm*Plots.mm - Plots.left(area)) / Plots.width(area)
        annotate!(panel, ((x, 0.5),
            text(label, GUIDE_FONT_SIZE, :center, FONT_FAMILY; rotation=90)))
    end
    mkpath(output)
    for extension in ("png", "pdf")
        path = joinpath(output, "$(spec.stem).$extension")
        savefig(figure, path)
        println("Saved $path")
    end
    figure
end

function main()
    gr()
    make_plot(:oe)
    make_plot(:dv)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    PaperPlots.main()
end
