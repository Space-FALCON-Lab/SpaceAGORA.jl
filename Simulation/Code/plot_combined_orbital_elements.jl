#!/usr/bin/env julia

if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(normpath(joinpath(@__DIR__, "..", "..",
        "9_SpaceAGORA.jl-main_full integration")))
end

include("plot_orbital_elements.jl")

function main(args=ARGS)
    output = joinpath(OrbitalElementPlots.OUTPUT, "3_combined_plots")
    options = OrbitalElementPlots.parse_options(vcat(["--output", output], args))
    rows, thresholds = OrbitalElementPlots.load_grid(options.input, options.windows)
    mkpath(options.output)
    OrbitalElementPlots.gr()
    for window in options.windows
        OrbitalElementPlots.make_plot(rows, thresholds, "combined", window, options.output)
    end
    println("Combined both schedulers for $(length(rows)) verified endpoints; existing separate figures retained.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
