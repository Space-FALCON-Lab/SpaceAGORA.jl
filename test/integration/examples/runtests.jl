using Pkg
# Choose the plotting dependency family before companion helpers import core.
Pkg.activate(joinpath(@__DIR__, "..", "..", "..", "examples"); io=devnull)
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "..", "helpers", "load_edg.jl"))
Base.include(@__MODULE__, joinpath(@__DIR__, "..", "..", "helpers", "load_hypr.jl"))
include(joinpath(@__DIR__, "..", "..", "smoke", "ci_examples_regression.jl"))
