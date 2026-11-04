#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

include("plot_orbital_elements.jl")
using CSV

function main(args=ARGS)
    output = joinpath(OrbitalElementPlots.StudyConfig.DEFAULT_OUTPUT, "images", "delta_v_RTN")
    options = OrbitalElementPlots.parse_options(vcat(["--output", output], args))
    rows, thresholds = OrbitalElementPlots.load_grid(options.input, options.windows; include_rtn=true)
    mkpath(options.output)
    OrbitalElementPlots.gr()
    for window in options.windows
        endpoints = [(scheduler=row.scheduler, target_altitude_km=row.target_altitude_km,
            delta_i_deg=row.delta_i_deg, helpers=row.helpers, orbits=row.orbits,
            time_s=row.time_s, row=row.row, dv_R_mps=row.dv_R_mps,
            dv_T_mps=row.dv_T_mps, dv_N_mps=row.dv_N_mps,
            source_trajectory=row.source_trajectory) for row in rows if row.orbits == window]
        CSV.write(joinpath(options.output, "cumulative_dv_RTN_$(window)orbits.csv"), endpoints)
        OrbitalElementPlots.make_plot(rows, thresholds, "combined", window, options.output; quantity=:rtn)
    end
    CSV.write(joinpath(options.output, "continuous_access_thresholds.csv"), thresholds)
    println("Plotted cumulative RTN delta-v at $(length(rows)) endpoints from verified trajectories.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
