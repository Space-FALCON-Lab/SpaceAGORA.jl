include(joinpath(@__DIR__, "..", "..", "..", "scripts", "setup_core.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "optional_packages.jl"))
length(ARGS) == 1 || error("Usage: setup.jl <scratch-directory>")
directory = abspath(only(ARGS))
CoreEnvironment.setup(joinpath(directory, "core"))
OptionalPackages.setup("hypr", joinpath(directory, "hypr"))
# The integration campaign (optional.jl) also runs the RPO comparison figure
# tests, whose figures come from SpaceAGORA's PlotlyJS extension: the
# environment that runs them adds PlotlyJS, as a user's project would.
using Pkg
Pkg.activate(joinpath(directory, "hypr"))
Pkg.add("PlotlyJS"; preserve=Pkg.PRESERVE_ALL)
