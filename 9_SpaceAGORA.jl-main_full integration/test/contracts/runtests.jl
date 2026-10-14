# Local alias for the PR contract gates. Every contract gate that applies to a
# pull request is wired into pr_runtests.jl; this file keeps the historical
# local entry point and runs exactly that set.
include(joinpath(@__DIR__, "pr_runtests.jl"))
