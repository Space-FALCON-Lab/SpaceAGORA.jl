include(joinpath(@__DIR__, "..", "..", "..", "scripts", "setup_edg.jl"))
length(ARGS) == 1 || error("Usage: setup.jl <scratch-directory>")
directory = abspath(only(ARGS))
EDGTestInstallation.CoreEnvironment.setup(joinpath(directory, "core"))
EDGTestInstallation.setup(joinpath(directory, "edg"))

