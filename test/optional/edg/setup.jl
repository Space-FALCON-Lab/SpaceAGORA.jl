include(joinpath(@__DIR__, "..", "..", "..", "scripts", "setup_edg.jl"))
length(ARGS) == 1 || error("Usage: setup.jl <scratch-directory>")
directory = abspath(only(ARGS))
EDGTestInstallation.setup(joinpath(directory, "core"); with_edg=false)
EDGTestInstallation.setup(joinpath(directory, "edg"))

