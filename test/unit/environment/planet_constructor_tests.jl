@testset "Typed Planet Constructors + Topography Workspace" begin
    mars = Mars("", SPICE_PATH)
    venus = Venus("", SPICE_PATH)
    titan = Titan("", SPICE_PATH)

    @test mars.name == "Mars"
    @test venus.name == "Venus"
    @test titan.name == "Titan"
    @test mars.μ > 0.0
    @test venus.μ > 0.0
    @test titan.μ > 0.0

    earth = Earth("", SPICE_PATH)
    moon = Moon("", SPICE_PATH)
    earth_radii_km = lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
        bodvrd("EARTH", "RADII")
    end
    moon_radii_km = lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
        bodvrd("MOON", "RADII")
    end
    @test isapprox(earth.Rp_e, earth_radii_km[1] * 1e3; atol=1e-6, rtol=0.0)
    @test isapprox(earth.Rp_p, earth_radii_km[3] * 1e3; atol=1e-6, rtol=0.0)
    @test isapprox(moon.Rp_e, moon_radii_km[1] * 1e3; atol=1e-6, rtol=0.0)
    @test isapprox(moon.Rp_p, moon_radii_km[3] * 1e3; atol=1e-6, rtol=0.0)

    mktempdir() do tmp
        topo_file = joinpath(tmp, "topo_harmonics.csv")
        write(topo_file, "degree,order,C,S\n0,0,1.0,0.0\n1,0,0.1,0.0\n1,1,0.05,0.02\n")
        SimulationModel.Planets.TopographyHarmonicsWorkspace!(topo_file, earth)
        @test size(earth.topography_workspace.Clm) == (2, 2)
        @test size(earth.topography_workspace.Slm) == (2, 2)
        @test isapprox(earth.topography_workspace.Clm[2, 2], 0.05; atol=0.0, rtol=0.0)
        @test isapprox(earth.topography_workspace.Slm[2, 2], 0.02; atol=0.0, rtol=0.0)
    end
end

