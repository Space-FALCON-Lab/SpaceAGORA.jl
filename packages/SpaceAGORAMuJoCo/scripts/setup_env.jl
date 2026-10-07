# Build an isolated environment with SpaceAGORA (from this checkout, root manifest versions) and
# SpaceAGORAMuJoCo, then precompile it:
#   julia --startup-file=no packages/SpaceAGORAMuJoCo/scripts/setup_env.jl <project-directory>
# Test with: julia --project=<project-directory> packages/SpaceAGORAMuJoCo/test/runtests.jl
using Pkg
const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
length(ARGS) == 1 || error("Usage: julia packages/SpaceAGORAMuJoCo/scripts/setup_env.jl <project-directory>")
include(joinpath(ROOT, "scripts", "setup_hypr.jl"))
env = abspath(only(ARGS))
HYPRInstallation.setup(env; with_hypr=false)
previous = Base.active_project()
try
    Pkg.activate(env)
    Pkg.develop(PackageSpec(path=joinpath(ROOT, "packages", "SpaceAGORAMuJoCo")); preserve=Pkg.PRESERVE_ALL)
    Pkg.add(["Test", "LazyArtifacts", "Libdl", "Artifacts", "LinearAlgebra", "StaticArrays"]; preserve=Pkg.PRESERVE_ALL)
    Pkg.instantiate()
    Pkg.precompile(["SpaceAGORA", "SpaceAGORAMuJoCo"]; strict=true)
finally
    previous === nothing ? Pkg.activate() : Pkg.activate(dirname(previous))
end
