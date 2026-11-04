using SatelliteToolboxSgp4, SatelliteToolboxTle, TOML
const root = normpath(joinpath(@__DIR__, "..", ".."))
@assert Base.active_project() == joinpath(root,"Project.toml")
@assert !haskey(TOML.parsefile(Base.active_project())["deps"],"SatelliteToolboxSgp4")
@assert pkgversion(SatelliteToolboxSgp4) == v"3.0.0"
@assert pkgversion(SatelliteToolboxTle) == v"2.0.0"
println("Development catalogue dependencies available outside the core project")

# Public catalogue fixture, frozen with SGP4 2.3.0 before the coordinated
# solver/SatelliteToolbox upgrade. Keep the same WGS84 model and km units.
using Test
@testset "Catalogue propagation survives the dependency upgrade" begin
    lines = readlines(joinpath(@__DIR__, "viewer_demos", "cygnss_historical.tle"))
    idx = findfirst(line -> startswith(line, "1 "), lines)
    tle = only(read_tles("baseline\n" * lines[idx] * "\n" * lines[idx + 1]))
    for (minutes, expected) in (
        (0.0, [163.57392958010647, -6822.351745492144, 0.014209575539262971,
               6.2624251378079325, 0.1617064427133884, 4.383221189871168]),
        (5760.0, [-4395.353711354704, 3558.9576542199457, -3825.1306305112666,
                  -4.282870831425125, -6.258311984271359, -0.9140800846104301]),
    )
        prop = sgp4_init(tle; sgp4c=SGP4C_WGS84)
        r, v = sgp4!(prop, minutes)
        @test r ≈ expected[1:3] atol=1e-8 rtol=0
        @test v ≈ expected[4:6] atol=1e-11 rtol=0
    end
end

include(joinpath(@__DIR__, "test_dev_launcher.jl"))
