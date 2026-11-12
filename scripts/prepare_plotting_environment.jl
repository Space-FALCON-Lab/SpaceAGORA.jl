module PlottingEnvironmentSetup

using Pkg, TOML

"""Prepare a test/example project with PlotlyJS without changing the core lockfile."""
function prepare(directory)
    root = normpath(joinpath(@__DIR__, ".."))
    project = TOML.parsefile(joinpath(directory, "Project.toml"))
    core = TOML.parsefile(joinpath(root, "Project.toml"))
    undeclared = setdiff(keys(core["deps"]), keys(project["deps"]))
    isempty(undeclared) || error(
        "$(joinpath(directory, "Project.toml")) must also declare these dependencies of the root project: $(join(sort(collect(undeclared)), ", ")).")

    # Seed from the core versions, relocating its path for the child project.
    # Only this generated manifest is written; the core lockfile stays intact.
    manifest = TOML.parsefile(joinpath(root, "Manifest.toml"))
    manifest["deps"]["SpaceAGORA"][1]["path"] = relpath(root, directory)
    open(joinpath(directory, "Manifest.toml"), "w") do io
        TOML.print(io, manifest; sorted=true)
    end
    Pkg.activate(directory)
    if VersionNumber(manifest["deps"]["CSV"][1]["version"]) >= v"1" ||
       VersionNumber(manifest["deps"]["JSON"][1]["version"]) >= v"1"
        # PlotlyJS requires JSON 0.21 and cannot share CSV 1's Parsers 3.
        # The child projects intentionally retain CSV 0.10 / JSON 0.21.
        # Relax only that family (SpaceIndices also depends on JSON), preserving
        # all unrelated core versions and their reusable package images.
        Pkg.update(["CSV", "JSON", "Parsers", "SpaceIndices"];
            mode=Pkg.PKGMODE_MANIFEST, preserve=Pkg.PRESERVE_ALL, update_registry=false)
    else
        Pkg.resolve()
    end
    Pkg.instantiate()
    return nothing
end

end
