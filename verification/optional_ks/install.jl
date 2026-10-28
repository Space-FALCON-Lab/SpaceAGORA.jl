# From this candidate checkout: julia --startup-file=no verification/optional_ks/install.jl NEW_ENVIRONMENT [--without-ks]
using Pkg
1 <= length(ARGS) <= 2 || error("usage: install.jl NEW_ENVIRONMENT [--without-ks]")
without_ks=length(ARGS)==2
without_ks && ARGS[2]!="--without-ks" && error("unknown option: $(ARGS[2])")
directory=abspath(ARGS[1])
isdir(directory) && !isempty(readdir(directory)) && error("Use a new environment directory: $directory")
Pkg.activate(directory)
Pkg.develop(PackageSpec(path=normpath(joinpath(@__DIR__,"../.."))))
if !without_ks
    Pkg.add(PackageSpec(url="https://github.com/gfalcon2/KSPropagator.jl",rev="8cfc936fe2ba403996ed7682aa700b56d6e04c5d"))
end
Pkg.add(["Test","LinearAlgebra","TOML","SHA","CSV","DataFrames","Serialization","DiffEqCallbacks","StaticArrays"])
Pkg.instantiate()
using SpaceAGORA
if without_ks
    isnothing(Base.find_package("KSPropagator")) || error("KS must be absent from this environment; use JULIA_LOAD_PATH=@:@stdlib")
else
    @eval using KSPropagator
end
println("OPTIONAL_KS_INSTALL_SUCCESS without_ks=",without_ks," julia=",VERSION," spaceagora=",pathof(SpaceAGORA))
