include(joinpath(@__DIR__, "..", "..", "..", "scripts", "setup_core.jl"))
include(joinpath(@__DIR__, "..", "..", "..", "scripts", "optional_packages.jl"))
length(ARGS) == 1 || error("Usage: setup.jl <scratch-directory>")
directory = abspath(only(ARGS))
CoreEnvironment.setup(joinpath(directory, "core"))
OptionalPackages.setup("hypr", joinpath(directory, "hypr"))
# The integration campaign directly imports these test tools. Declare them in
# its environment rather than relying on the old copied core project to expose
# transitive dependencies. PlotlyJS activates the comparison-figure extension.
using Pkg
Pkg.activate(joinpath(directory, "hypr"))
# Keep direct package versions, but let transitive dependencies resolve with
# the newly requested plotting extension in this fresh test environment.
Pkg.add(["ComponentArrays", "DiffEqCallbacks", "JSON", "PlotlyJS"]; preserve=Pkg.PRESERVE_DIRECT)
