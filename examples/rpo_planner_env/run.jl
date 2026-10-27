include(joinpath(@__DIR__, "environment.jl"))
environment = RPOExampleEnvironment.directory()
isfile(joinpath(environment, "Project.toml")) || error(
    "Prepare the RPO environment first: julia --startup-file=no examples/rpo_planner_env/setup.jl")
file = isempty(ARGS) ? joinpath(@__DIR__, "smoke.jl") : abspath(first(ARGS))
arguments = isempty(ARGS) ? String[] : ARGS[2:end]
isfile(file) || error("Example file does not exist: $file")
run(`$(Base.julia_cmd()) --startup-file=no --project=$environment $file $arguments`)
