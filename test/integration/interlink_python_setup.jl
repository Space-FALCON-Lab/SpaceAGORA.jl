# Prepare the explicit cyclic-matching test environment without changing core inputs.
using Pkg, TOML
length(ARGS) == 1 || error("Usage: interlink_python_setup.jl <new-environment>")
root = normpath(joinpath(@__DIR__, "..", ".."))
environment = abspath(only(ARGS))
ispath(environment) && !isempty(readdir(environment)) && error("Use an empty test environment.")
mkpath(environment)
project = Dict("deps" => copy(TOML.parsefile(joinpath(root, "Project.toml"))["deps"]))
project["deps"]["SpaceAGORA"] = "afbfb69f-5c0b-4832-b760-43725dff8540"
open(joinpath(environment, "Project.toml"), "w") do io
    TOML.print(io, project; sorted=true)
end
cp(joinpath(root, "Manifest.toml"), joinpath(environment, "Manifest.toml"))
Pkg.activate(environment)
Pkg.develop(PackageSpec(path=root); preserve=Pkg.PRESERVE_ALL)
Pkg.add(PackageSpec(name="PythonCall", version="0.9.36"); preserve=Pkg.PRESERVE_ALL)
Pkg.instantiate()
