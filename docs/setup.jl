using Pkg

# Prepare the documentation environment: Documenter and this checkout's
# SpaceAGORA (docs/Project.toml). Its manifest is generated, not committed, and
# is seeded with the root manifest on every run. Resolution starts from that
# seed: the packages the root manifest records keep its versions and package
# images, the documentation tools are added, and the records of the libraries
# bundled with Julia, such as OpenSSL_jll, follow the installed Julia.
docs = @__DIR__
cp(joinpath(dirname(docs), "Manifest.toml"), joinpath(docs, "Manifest.toml"); force=true)
Pkg.activate(docs)
Pkg.resolve()
Pkg.instantiate()
