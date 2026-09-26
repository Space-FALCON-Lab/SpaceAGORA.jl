using Pkg
Pkg.activate(@__DIR__)
# The manifest is not committed. When an earlier setup left one, update GRAMSuite so that
# a changed pin in Project.toml (a new rev after a pull) replaces the recorded revision;
# Pkg.instantiate alone keeps the old one.
isfile(joinpath(@__DIR__, "Manifest.toml")) && Pkg.update("GRAMSuite")
Pkg.instantiate()
println("Odyssey surrogate environment is ready. No native GRAM setup is required.")
