include(joinpath(@__DIR__, "..", "..", "..", "scripts", "setup_visualizer.jl"))
length(ARGS) == 1 || error("Usage: setup.jl <scratch-directory>")
for kind in ("core", "visualizer")
    VisualizerInstallation.setup(joinpath(abspath(only(ARGS)), kind); with_visualizer=kind == "visualizer")
end
