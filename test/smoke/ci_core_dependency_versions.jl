using Test, TOML, CSV, JSON

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const MANIFEST = TOML.parsefile(joinpath(ROOT, "Manifest.toml"))["deps"]

@testset "core CSV and JSON use the committed versions" begin
    @test realpath(Base.active_project()) == realpath(joinpath(ROOT, "Project.toml"))
    @test pkgversion(CSV) == VersionNumber(MANIFEST["CSV"][1]["version"])
    @test pkgversion(JSON) == VersionNumber(MANIFEST["JSON"][1]["version"])
end

# The full suite uses PlotlyJS's older dependency set. Exercise the core's
# committed CSV/JSON versions independently, including saved and streamed data.
include(joinpath(ROOT, "test", "unit", "io", "long_layout_tests.jl"))
include(joinpath(ROOT, "test", "unit", "io", "long_streaming_tests.jl"))
