# Plotting examples require CSV 0.10 / JSON 0.21 through PlotlyJS.
# Keep their generated manifest separate from the newer core dependency set.
include(joinpath(@__DIR__, "..", "scripts", "prepare_plotting_environment.jl"))
PlottingEnvironmentSetup.prepare(@__DIR__)
