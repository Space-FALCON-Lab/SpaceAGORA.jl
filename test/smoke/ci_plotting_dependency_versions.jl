using Test, TOML, CSV, JSON, PlotlyJS, SpaceAGORA

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const CORE = TOML.parsefile(joinpath(ROOT, "Manifest.toml"))["deps"]
const CHILD = TOML.parsefile(joinpath(dirname(Base.active_project()), "Manifest.toml"))["deps"]
const ADJUSTED = Set(["CSV", "JSON", "Parsers", "SpaceIndices"])

@testset "plotting environment remains compatible and isolated" begin
    @test v"0.10" <= pkgversion(CSV) < v"0.11"
    @test v"0.21" <= pkgversion(JSON) < v"0.22"
    @test Base.get_extension(SpaceAGORA, :SpaceAGORAPlotlyJSExt) !== nothing
    for name in intersect(keys(CORE), keys(CHILD))
        name in ADJUSTED && continue
        a, b = only(CORE[name]), only(CHILD[name])
        # Registry package trees are portable; Julia-bundled stdlibs can follow
        # the installed patch release and the local core path is relocated.
        haskey(a, "git-tree-sha1") || continue
        @test get(b, "git-tree-sha1", nothing) == a["git-tree-sha1"]
    end
end
