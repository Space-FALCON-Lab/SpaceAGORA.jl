using SatelliteToolboxSgp4, SatelliteToolboxTle, TOML
const root = normpath(joinpath(@__DIR__, "..", ".."))
@assert Base.active_project() == joinpath(root,"Project.toml")
@assert !haskey(TOML.parsefile(Base.active_project())["deps"],"SatelliteToolboxSgp4")
@assert pkgversion(SatelliteToolboxSgp4) == v"2.3.0"
@assert pkgversion(SatelliteToolboxTle) == v"1.1.1"
println("Development dependencies available; core project unchanged")

include(joinpath(@__DIR__, "test_dev_launcher.jl"))
