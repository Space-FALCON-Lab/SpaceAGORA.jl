include("oracle_closed_forms.jl")
include("../Output_2_Simulation_absolute_RN/Code/StudyConfig.jl")
include("../Results_For_Paper/regenerate/plot_overlays.jl")

function plot_original_leapfrog(kind=:dv)
    kind in (:oe, :dv) || throw(ArgumentError("Expected :oe or :dv"))
    p = Params(R_E=StudyConfig.EARTH_RADIUS_M, mu=StudyConfig.MU)
    grid = joinpath(@__DIR__, "..", "Results_For_Paper", "regenerate", "oracle_predictions.csv")
    predictions = map(CSV.File(grid)) do row
        r = NamedTuple(row)
        if r.rule == "amax" && r.h_t_km == 1000 && r.N_h < 232
            b0 = r_h(p) * deg2rad(mod(1.0, 360 / r.N_h))
            # Without inc_deg, leapfrog returns only tangential impulse and shell gain.
            scalar = leapfrog(r.N_h, b0; n_orbits=100, p)
            return merge(r, (; dvR=missing, dvT=scalar.dvT, dvNp=missing, dvNm=missing,
                da_km=scalar.da_km, di_deg=missing, dOmega_deg=missing))
        end
        return r
    end
    prediction_curve = paper_prediction_curve(predictions)
    if kind == :oe
        return plot_oe_overlay(prediction_curve; output=@__DIR__,
            specification=merge(PaperPlots.SPECS.oe, (; stem="OE_original_leapfrog_100orbits")))
    end
    PaperPlots.make_plot(:dv; output=@__DIR__, prediction_curve,
        specification=merge(PaperPlots.SPECS.dv, (; stem="dv_RTN_original_leapfrog_100orbits")),
        legend_settings=merge(PaperPlots.DV_LEGEND, (; x_offset=0.0)))
end

if abspath(PROGRAM_FILE) == @__FILE__
    plot_original_leapfrog(:dv)
    plot_original_leapfrog(:oe)
end
