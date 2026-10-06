include(joinpath(@__DIR__, "..", "..", "scripts", "setup_hypr.jl"))
length(ARGS) == 1 || error("Usage: setup.jl <scratch-directory>")
for kind in ("core", "hypr")
    HYPRInstallation.setup(joinpath(abspath(only(ARGS)), kind); with_hypr=kind == "hypr")
end
