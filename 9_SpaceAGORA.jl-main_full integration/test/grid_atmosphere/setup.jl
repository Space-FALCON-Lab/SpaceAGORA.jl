using Pkg

# Resolve the two local packages together. Stacking independently resolved
# projects can shadow StaticArrays' optional dependencies during precompilation.
repo = normpath(joinpath(@__DIR__, "..", ".."))
wrapper = joinpath(repo, "data", "GRAMSuite.jl")
isfile(joinpath(wrapper, "Project.toml")) || error(
    "The pinned GRAMSuite checkout is missing; retrieve it before preparing grid tests.")
Pkg.activate(@__DIR__)
Pkg.develop([Pkg.PackageSpec(path=repo), Pkg.PackageSpec(path=wrapper)])
Pkg.instantiate()
