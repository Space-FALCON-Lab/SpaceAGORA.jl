include(joinpath(@__DIR__, "..", "..", "scripts", "setup_hypr.jl"))
# Generated installation metadata stays local; the tracked template is immutable.
# The examples need no optional package, so HYPR is not installed.
project = joinpath(@__DIR__, "Project.toml")
isfile(project) || cp(joinpath(@__DIR__, "Project.template.toml"), project)
HYPRInstallation.setup(@__DIR__; with_hypr=false)
