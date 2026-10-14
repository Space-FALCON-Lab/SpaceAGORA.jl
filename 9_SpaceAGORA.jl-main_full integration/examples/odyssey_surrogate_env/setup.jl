using Pkg
Pkg.activate(@__DIR__)
# The manifest is not committed. When an earlier setup left one, update GRAMSuite so that
# a changed pin in Project.toml (a new rev after a pull) replaces the recorded revision;
# Pkg.instantiate alone keeps the old one.
isfile(joinpath(@__DIR__, "Manifest.toml")) && Pkg.update("GRAMSuite")
Pkg.instantiate()
# EDG is a separately installed provider even for the frozen-atmosphere example.
# Keep its environment outside both source checkouts; persist only the selected
# dependency in this example's already explicit environment.
include(joinpath(@__DIR__, "..", "..", "scripts", "optional_packages.jl"))
provider = get(ENV, "SPACEAGORA_EDG_PATH", "")
if isempty(provider)
    provider = OptionalPackages.checkout("edg")
else
    include(joinpath(@__DIR__, "..", "..", "scripts", "setup_edg.jl"))
    provider = EDGTestInstallation.provider_path(EDGTestInstallation.ROOT; local_path=provider)
end
Pkg.develop(PackageSpec(path=provider); preserve=Pkg.PRESERVE_ALL)
println("Odyssey surrogate environment is ready with the EDG provider. No native GRAM setup is required.")
