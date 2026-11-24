include("oracle_closed_forms.jl")
include("../../Output_2_Simulation_absolute_RN/Code/StudyConfig.jl")
include("plot_overlays.jl")

cd(@__DIR__) do
    p = Params(R_E = StudyConfig.EARTH_RADIUS_M, mu = StudyConfig.MU)
    prediction_table("oracle_predictions.csv"; p = p)

    # Keep the unmodified grid as well as the exact table passed to the renderer.
    open("plotted_predictions.csv", "w") do io
        lines = readlines("oracle_predictions.csv")
        println(io, lines[1], ",method,b0_m")
        for line in lines[2:end]
            fields = String.(split(line, ','))
            rule, altitude, inclination, helpers = fields
            nh = parse(Int, helpers)
            if rule == "amax" && parse(Float64, altitude) == 1000 && nh < 232
                # Target at +1 degree; helpers start at multiples of 360/N_h.
                u0_deg = 1.0
                b0 = r_h(p) * deg2rad(mod(u0_deg, 360 / nh))
                result = leapfrog(nh, b0; n_orbits = 100, p, inc_deg=parse(Float64, inclination), u0_deg)
                fields[5:9] = string.([result.dvR, result.dvT, result.dvNp, result.dvNm, result.da_km])
                fields[10:11] .= ""  # Fixed-plane model: no inclination or node prediction.
                println(io, join(fields, ','), ",leapfrog_rtn,", b0)
            else
                println(io, line, ",drift,")
            end
        end
    end
    plot_overlays()
end
