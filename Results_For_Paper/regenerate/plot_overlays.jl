import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..", "..", "9_SpaceAGORA.jl-main_full integration")))
include("../plot_paper_results.jl")
using CSV

function paper_prediction_curve(predictions)
    rules = (gve_sma="amax", gve_inc="imax")
    columns = (da_m=:da_km, di_deg=:di_deg, draan_deg=:dOmega_deg, dv_R_mps=:dvR,
        dv_T_mps=:dvT, dv_N_positive_mps=:dvNp, dv_N_negative_mps=:dvNm)
    function prediction_curve(field, policy, altitude, inclination)
        field == :de && return nothing
        rule = getproperty(rules, Symbol(policy))
        curve = sort(filter(r -> r.rule == rule && r.h_t_km == altitude &&
            r.inc_deg == inclination, predictions); by=r -> r.N_h)
        helpers = [r.N_h for r in curve]
        helpers == [1; 25:25:800] || error("Incomplete prediction: $rule, $altitude km, $inclination deg")
        # The renderer takes source units (metres for da_m); NaNs break unavailable lines.
        factor = field == :da_m ? 1000.0 : 1.0
        helpers, [factor * coalesce(getproperty(r, columns[field]), NaN) for r in curve]
    end
    return prediction_curve
end

function plot_oe_overlay(prediction_curve; output=@__DIR__, specification=PaperPlots.SPECS.oe,
        legend_x_offset=0.0)
    paper = dirname(@__DIR__)
    key(r) = (r.target_altitude_km, r.delta_i_deg, r.helpers, r.orbits, r.row)
    file = specification.csv
    off = Dict(key(r) => r.draan_deg for r in CSV.File(joinpath(dirname(paper),
        "Output_4_no_laser", "Output", "images", "for_paper", file)))
    oe_rows = [merge(NamedTuple(r), (; draan_deg=r.draan_deg - off[key(r)]))
        for r in CSV.File(joinpath(paper, file))]
    # The links-off subtraction shortens the node tick labels. Retain their original
    # 1.277 mm of extra gutter so the aligned columns and legend do not shift left.
    margins = PaperPlots.OE_SPACING.panel_margins_mm
    spacing = merge(PaperPlots.OE_SPACING,
        (; panel_margins_mm=merge(margins, (; left=margins.left + 1.276991666))))
    PaperPlots.make_plot(:oe; output, specification, rows=oe_rows, prediction_curve, spacing,
        legend_settings=merge(PaperPlots.OE_LEGEND, (; x_offset=legend_x_offset)))
end

function plot_overlays(; legend_x_offset=-0.0)
    prediction_curve = paper_prediction_curve(collect(CSV.File(joinpath(@__DIR__, "plotted_predictions.csv"))))
    plot_oe_overlay(prediction_curve; legend_x_offset)
    PaperPlots.make_plot(:dv; output=@__DIR__, prediction_curve,
        legend_settings=merge(PaperPlots.DV_LEGEND, (; x_offset=legend_x_offset)))
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_overlays()
end
