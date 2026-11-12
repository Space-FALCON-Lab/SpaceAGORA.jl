# The full suite includes PlotlyJS and retains CSV 0.10 / JSON 0.21.
# The root environment independently exercises the newer CSV / JSON versions.
include(joinpath(@__DIR__, "..", "scripts", "prepare_plotting_environment.jl"))
PlottingEnvironmentSetup.prepare(@__DIR__)
