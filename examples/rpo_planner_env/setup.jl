using Pkg
root=normpath(joinpath(@__DIR__, "..", ".."))
# Start from the reviewed dependency versions on a clean installation.
manifest=joinpath(@__DIR__, "Manifest.toml")
isfile(manifest) || cp(joinpath(root,"Manifest.toml"),manifest)
Pkg.develop([PackageSpec(path=root), PackageSpec(path=joinpath(root,"packages","SpaceAGORAHYPR"))];preserve=Pkg.PRESERVE_ALL)
Pkg.instantiate()
