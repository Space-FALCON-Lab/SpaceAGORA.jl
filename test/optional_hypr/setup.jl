# Prepare independent installation projects outside the source tree.
using Pkg, TOML
root = normpath(joinpath(@__DIR__, "..", ".."))
length(ARGS) == 1 || error("Usage: setup.jl <scratch-directory>")
for kind in ("core", "hypr")
    env = joinpath(abspath(only(ARGS)), kind)
    mkpath(env)
    open(joinpath(env,"Project.toml"), "w") do io
        TOML.print(io, Dict("deps" => Dict("StaticArrays" => "90137ffa-7385-5640-81b9-e52037218182",
            "JSON" => "682c06a0-de6a-54ab-a142-c8b1cf79cde6",
            "ComponentArrays" => "b0b7db55-cfe3-40fc-9ded-d10e2dbeff66",
            "DiffEqCallbacks" => "459566f4-90b8-5000-8ac3-15dfb0a30def")))
    end
    cp(joinpath(root,"Manifest.toml"), joinpath(env,"Manifest.toml");force=true)
    Pkg.activate(env)
    specs = [PackageSpec(path=root)]
    kind == "hypr" && push!(specs, PackageSpec(path=joinpath(root,"packages","SpaceAGORAHYPR")))
    Pkg.develop(specs; preserve=Pkg.PRESERVE_ALL)
end
