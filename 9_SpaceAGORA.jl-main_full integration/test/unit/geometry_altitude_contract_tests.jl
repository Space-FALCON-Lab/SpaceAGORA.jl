module GeometryAltitudeContractTests
using Test
using LinearAlgebra
using StaticArrays
using SPICE: recgeo
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const G = SM.Geodesy
const FT = SM.FrameTransforms
const REPO = normpath(joinpath(@__DIR__, "..", ".."))

module RawGeometry
const args = :geometry_contract_includer
include(joinpath(@__DIR__, "..", "helpers", "reference_system_sandbox.jl"))
end

@testset "Shared Bowring expressions retain caller contracts" begin
    @test FT.Geodesy === G
    @test parentmodule(G._bowring_latitude) === G
    @test parentmodule(G._bowring_altitude) === G
    for planet in (SM.Earth(), SM.Mars(), SM.Venus()), latitude in (-pi/2, -pi/3, 0.0, pi/4, pi/2)
        for height in (0.0, 150e3, 500e3)
            longitude = -0.7
            r = SVector{3, Float64}(FT.latlongtor([latitude, longitude, height], planet, 0.0, 0.0, 0.0))
            result = FT.rtolatlong(r, planet)
            spice_lon, spice_lat, spice_height = recgeo(collect(r), planet.Rp_e,
                (planet.Rp_e - planet.Rp_p) / planet.Rp_e)
            @test result isa SVector{3, Float64}
            @test result[1] ≈ spice_height atol=1e-3 rtol=0
            # One-pass Bowring latitude is approximate at finite height. On
            # this 0--500 km grid the pinned original differs from CSPICE by
            # at most 4.452e-9 rad (Mars, 500 km); retain the independent
            # 1 mm altitude check and bound latitude here to 1e-8 rad.
            @test result[2] ≈ spice_lat atol=1e-8 rtol=0
            @test result[3] ≈ spice_lon atol=1e-14 rtol=0
            @test result[1] ≈ height atol=1e-3 rtol=0
            @test RawGeometry.rtolatlong(r, planet) == result
            @test FT.rtolatlong(r, planet, SM.SimpleEphemeridesModel()) == result
            direction = r / norm(r)
            @test G.geodetic_altitude(norm(r), direction, planet) ≈ spice_height atol=1e-3 rtol=0
        end
    end

    # Spherical geometry supplies an analytic large-height and below-surface
    # reference without assuming Bowring is exact at arbitrary ellipsoid heights.
    sphere = (Rp_e=6.0e6, Rp_p=6.0e6)
    for height in (-500.0, 0.0, 1.0e9)
        r = SVector(0.0, 0.0, sphere.Rp_p + height)
        @test FT.rtolatlong(r, sphere)[1] ≈ height atol=2e-7 rtol=0
    end
    # Public wrappers retain their existing validation and conversion policies.
    direction = SVector(1.0, 0.0, 0.0)
    @test isnan(G.radius_for_geodetic_altitude(-1.0, direction, SM.Earth()))
    @test_throws ArgumentError SM.InitialCondition(SM.Earth(); ra=500e3, hp=-1.0)
    @test_throws MethodError FT.rtolatlong(Float32[7e6, 0, 0], sphere)
    @test_throws MethodError G.geodetic_altitude(7.0f6, direction, sphere)
    for T in (Float32, BigFloat)
        p = (Rp_e=T(6378137), Rp_p=T(6356752))
        @test FT._planet_flattening(p) isa Float64
        @test FT.rtolatlong(SVector(7e6, 2e6, -1e6), p) isa SVector{3, Float64}
        @test G.geodetic_altitude(7e6, direction, p) isa Float64
    end
end

@testset "Topography stays in the reference-system caller" begin
    seen = Ref{Any}(nothing)
    p = (Rp_e=6378137.0, Rp_p=6356752.314245,
        Clm_topo=:clm, Slm_topo=:slm, A_topo=6e6,
        topography_function=(args, clm, slm, lat, lon, a) -> begin
            seen[] = (args, clm, slm, lat, lon, a)
            a
        end)
    r = SVector(6.8e6, -2e6, 3e6)
    ellipsoid = FT.rtolatlong(r, p)
    @test seen[] === nothing
    actual = FT.rtolatlong(r, p, true)
    @test actual[1] == norm(r) - p.A_topo
    @test actual[2:3] == ellipsoid[2:3]
    @test seen[] == (nothing, :clm, :slm, ellipsoid[2], ellipsoid[3], p.A_topo)
    raw = RawGeometry.rtolatlong(r, p, true)
    @test raw == actual
    @test seen[][1] === RawGeometry.args
    bad = merge(p, (topography_function=(args...) -> error("topography witness"),))
    @test FT.rtolatlong(r, bad) == ellipsoid
    @test_throws ErrorException FT.rtolatlong(r, bad, true)
    @test_throws ErrorException RawGeometry.rtolatlong(r, bad, true)
end
end
