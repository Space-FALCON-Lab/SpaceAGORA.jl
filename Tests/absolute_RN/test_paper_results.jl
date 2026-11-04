import Pkg
const PROJECT = normpath(joinpath(@__DIR__, "..", ".."))
Pkg.activate(joinpath(PROJECT, "9_SpaceAGORA.jl-main_full integration"))
using Test, CSV, Plots
const PAPER = joinpath(PROJECT, "Simulation_absolute_RN", "Output", "images", "for_paper")
include(joinpath(PAPER, "plot_paper_results.jl"))
const PP = PaperPlots

function check_ytitle_alignment(figure, spec, spacing)
    for r in 1:4
        panel = figure.subplots[(r-1)*3+1]
        position, label = only(panel[:annotations])
        area = Plots.plotarea(panel)
        @test isapprox((Plots.left(area) + position[1]*Plots.width(area)) / Plots.mm,
            spacing.y_label_x_mm; atol=1e-8)
        @test position[2] == 0.5
        @test label.str == spec.labels[r]
        @test label.font.rotation == 90
        @test label.font.pointsize == PP.GUIDE_FONT_SIZE
    end
end

function check_spacing(figure, spacing)
    for r in 1:4, c in 1:3
        panel = figure.subplots[(r-1)*3+c]
        @test panel[:left_margin] ==
            (spacing.panel_margins_mm.left + (c == 1 ? spacing.y_label_margin_mm : 0))*Plots.mm
        @test panel[:right_margin] == spacing.panel_margins_mm.right*Plots.mm
        @test panel[:top_margin] == spacing.panel_margins_mm.top*Plots.mm
        @test panel[:bottom_margin] ==
            (r == 4 ? spacing.bottom_row_margin_mm : spacing.panel_margins_mm.bottom)*Plots.mm
    end
    @test figure.subplots[end][:top_margin] == spacing.legend_top_margin_mm*Plots.mm
end

function check_legend_geometry(figure, show_access, settings)
    footer = figure.subplots[end]
    baseline = show_access ? (0.85, 0.52, 0.17) : (0.78, 0.28)
    for row in eachindex(baseline), col in 1:3
        index = (row-1)*3+col
        sample = footer.series_list[index]
        annotation = footer[:annotations][index]
        x = 0.06 + (col-1)/3 + settings.x_offset
        y = baseline[1] + settings.y_offset -
            (baseline[1]-baseline[row])*settings.row_spacing
        @test collect(sample[:x]) ≈ [x, x+0.065]
        @test collect(sample[:y]) ≈ [y, y]
        @test annotation[1] ≈ x+0.08
        @test annotation[2] ≈ y
    end
end

@testset "Paper source data are exact copies" begin
    sources = (
        "orbital_element_deltas_100orbits.csv" => "orbital_elements/4_combined_1000_1050_1100km",
        "cumulative_dv_RT_N_branches_100orbits.csv" => "delta_v_RTN/6_signed_RT_positive_negative_N_1000_1050_1100km",
        "continuous_access_thresholds.csv" => "orbital_elements/4_combined_1000_1050_1100km",
    )
    for (file, directory) in sources
        @test read(joinpath(PAPER, file)) == read(joinpath(dirname(PAPER), directory, file))
    end
end

@testset "Paper rendering, altitude headers and configurable dimensions" begin
    mktempdir() do output
        for kind in (:oe, :dv)
            spec = getproperty(PP.SPECS, kind)
            rows = collect(CSV.File(joinpath(PAPER, spec.csv)))
            figure = PP.make_plot(kind; output)
            @test length(figure.subplots) == 13
            @test all(panel -> panel[:title] == "", figure.subplots[4:end])
            settings = kind == :oe ? PP.OE_SIZE : PP.DV_SIZE
            spacing = kind == :oe ? PP.OE_SPACING : PP.DV_SPACING
            @test figure[:size] == PP.canvas_geometry(settings).size
            check_ytitle_alignment(figure, spec, spacing)
            check_spacing(figure, spacing)
            check_legend_geometry(figure, kind == :oe, kind == :oe ? PP.OE_LEGEND : PP.DV_LEGEND)
            for (r, field) in enumerate(spec.fields), (c, altitude) in enumerate(PP.ALTITUDES)
                panel = figure.subplots[(r-1)*3+c]
                @test panel[:title] == (r == 1 ? "Target Alt.: $altitude km" : "")
                @test panel[:titlefontsize] == PP.FONT_SIZE
                scale = field == :da_m ? 1e-3 : field == :de ? 1e4 : field == :di_deg ? 1e3 : 1.0
                if field == :draan_deg
                    actual = [row.draan_deg for row in rows if row.delta_i_deg != 0]
                    lower, upper = ylims(panel)
                    @test lower <= minimum(actual) <= maximum(actual) <= upper
                    @test upper - lower <= 4
                    @test upper < -39
                    @test ylims(panel) == ylims(figure.subplots[10])
                else
                    @test ylims(panel) == getproperty(PP.Y_LIMITS, field) .* scale
                end
                @test panel[:xaxis][:guide] == (r == 4 ? "\$N_{\\mathrm{helpers}}\$" : "")
                @test panel[:yaxis][:guide] == ""
                @test length(panel[:annotations]) == (c == 1 ? 1 : 0)
                @test panel[:yaxis][:formatter] == :plain
                @test panel[:xaxis][:rotation] == 0
                @test panel[:xaxis][:ticks] ==
                    (collect(PP.HELPERS), ["1", "", "200", "", "400", "", "600", "", "800"])
                @test panel[:xaxis][:tickfontsize] == PP.TICK_FONT_SIZE
                @test panel[:yaxis][:guidefontsize] == PP.GUIDE_FONT_SIZE
                curves = filter(s -> s[:seriestype] == :path, panel.series_list)
                @test length(curves) == (field == :draan_deg ? 4 : 6)
                index = 0
                for policy in PP.POLICIES, (k, inclination) in enumerate(PP.INCLINATIONS)
                    field == :draan_deg && inclination == 0 && continue
                    index += 1
                    curve = curves[index]
                    expected = sort(filter(row -> row.scheduler == policy &&
                        row.target_altitude_km == altitude && row.delta_i_deg == inclination, rows);
                        by=row -> row.helpers)
                    @test collect(curve[:x]) == collect(PP.HELPERS)
                    @test collect(curve[:y]) == [getproperty(row, field)*scale for row in expected]
                    @test curve[:linecolor] == Plots.plot_color(getproperty(PP.LINE_COLORS, Symbol(policy))[k])
                    @test curve[:linestyle] == (policy == "gve_sma" ? :solid : :dot)
                    @test curve[:linewidth] == PP.LINE_WIDTH
                    @test curve[:markersize] == PP.MARKER_SIZE
                end
                @test length(filter(s -> s[:seriestype] == :straightline, panel.series_list)) ==
                    (kind == :oe ? 3 : 0)
            end
            footer = figure.subplots[end]
            @test length(footer[:annotations]) == (kind == :oe ? 9 : 6)
            @test length(footer.series_list) == (kind == :oe ? 9 : 6)
            labels = [a[3].str for a in footer[:annotations]]
            @test labels[1:6] == ["$policy, $inclination\u00b0"
                for policy in ("\$\\Delta a\$-max", "\$\\Delta i\$-max") for inclination in PP.INCLINATIONS]
            @test all(a -> a[3].font.pointsize == PP.FONT_SIZE, footer[:annotations])
            if kind == :oe
                @test labels[7:9] == ["\$N_{\\mathrm{cont}},\\; $i^\\circ\$" for i in PP.INCLINATIONS]
                for k in 1:3
                    sample = footer.series_list[6+k]
                    @test sample[:linecolor] == Plots.plot_color(PP.LINE_COLORS.gve_sma[k])
                    @test sample[:linestyle] == :dash
                    @test sample[:linewidth] == PP.ACCESS_LINE_WIDTH
                    @test sample[:linealpha] == 0.75
                    for panel in figure.subplots[1:12]
                        threshold = filter(s -> s[:seriestype] == :straightline, panel.series_list)[k]
                        @test sample[:linecolor] == threshold[:linecolor]
                    end
                end
                @test occursin(raw"\mathrm{km}", PP.SPECS.oe.labels[1])
                @test occursin("10^{-4}", PP.SPECS.oe.labels[2])
                @test occursin("10^{-3}", PP.SPECS.oe.labels[3])
                @test occursin(raw"\mathrm{deg}", PP.SPECS.oe.labels[3])
            end
            for ext in ("png", "pdf")
                @test filesize(joinpath(output, "$(spec.stem).$ext")) > 10_000
            end
        end
        changed = (subplot_width=600, subplot_height=300, footer_height=200, figure_size=nothing)
        @test PP.canvas_geometry(changed).size == (1800, 1400)
        dv_spacing = merge(PP.DV_SPACING, (
            panel_margins_mm=(left=2, right=3, top=4, bottom=2),
            y_label_x_mm=6, y_label_margin_mm=7,
            bottom_row_margin_mm=0, legend_top_margin_mm=1))
        dv_legend = (x_offset=-0.03, y_offset=0.04, row_spacing=0.8)
        customized = PP.make_plot(:dv; output, settings=changed, spacing=dv_spacing,
            legend_settings=dv_legend,
            colors=merge(PP.LINE_COLORS, (gve_sma=("#ff6600", "#4E5689", "#018B38"),)),
            styles=merge(PP.LINE_STYLES, (gve_sma=(:dash, :solid, :solid),)))
        @test customized[:size] == (1800, 1400)
        check_ytitle_alignment(customized, PP.SPECS.dv, dv_spacing)
        check_spacing(customized, dv_spacing)
        check_legend_geometry(customized, false, dv_legend)
        @test customized.subplots[1].series_list[1][:linecolor] == Plots.plot_color("#ff6600")
        @test customized.subplots[1].series_list[1][:linestyle] == :dash
        @test customized.subplots[end].series_list[1][:linestyle] == :dash
        override = merge(changed, (figure_size=(1500, 1300),))
        oe_spacing = merge(PP.OE_SPACING, (
            panel_margins_mm=(left=3, right=2, top=1, bottom=3),
            y_label_x_mm=4, y_label_margin_mm=6,
            bottom_row_margin_mm=1, legend_top_margin_mm=2))
        oe_legend = (x_offset=0.02, y_offset=-0.02, row_spacing=0.9)
        overridden = PP.make_plot(:oe; output, settings=override, spacing=oe_spacing,
            legend_settings=oe_legend,
            colors=merge(PP.LINE_COLORS, (gve_sma=("#ff6600", "#4E5689", "#018B38"),)))
        @test overridden[:size] == (1500, 1300)
        check_ytitle_alignment(overridden, PP.SPECS.oe, oe_spacing)
        check_spacing(overridden, oe_spacing)
        check_legend_geometry(overridden, true, oe_legend)
        check_ytitle_alignment(customized, PP.SPECS.dv, dv_spacing)
        check_spacing(customized, dv_spacing)
        check_legend_geometry(customized, false, dv_legend)
        @test overridden.subplots[end].series_list[7][:linecolor] == Plots.plot_color("#ff6600")
        @test_throws ErrorException PP.raan_limits([(delta_i_deg=0.5, draan_deg=NaN)])
        @test_throws ErrorException PP.raan_limits([(delta_i_deg=0.0, draan_deg=NaN)])
        lo, hi = PP.raan_limits([(delta_i_deg=0.5, draan_deg=-42.0)])
        @test lo < -42 < hi
        @test_throws ErrorException PP.canvas_geometry(merge(changed, (subplot_height=-1,)))
        @test_throws ErrorException PP.canvas_geometry(merge(changed, (figure_size=(0, 500),)))
        for invalid in (
                (x_offset=NaN, y_offset=0.0, row_spacing=1.0),
                (x_offset=0.0, y_offset=Inf, row_spacing=1.0),
                (x_offset=0.0, y_offset=0.0, row_spacing=0.0),
                (x_offset=0.0, y_offset=0.0, row_spacing=-1.0),
                (x_offset=0.0, y_offset=0.0, row_spacing=NaN))
            @test_throws ErrorException PP.make_plot(:oe; output, legend_settings=invalid)
        end
    end
end
