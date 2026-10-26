module LegacyVelocityFrameContractTests
using Test
using LinearAlgebra
using StaticArrays
using SPICE: georec
using SpaceAGORA

const FT = SpaceAGORA.SimulationModel.FrameTransforms
const IDENTITY = SMatrix{3, 3, Float64}(I)
const TILTED = @SMatrix [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
const BASE = (Rp_e=6378137.0, Rp_p=6356752.314245, μ=3.986004418e14,
    ω=SVector(0.0, 0.0, 7.292115e-5), L_PI=IDENTITY)

function _legacy_state(location, planet, γ, α, speed)
    elements = redirect_stdout(devnull) do
        FT.latlongtoOE(location, planet, γ, α, speed)
    end
    r, v = FT.orbitalelemtorv(SVector{7, Float64}(vcat(elements, 0.0)), planet)
    return elements, SVector{3, Float64}(r), SVector{3, Float64}(v)
end

@testset "Legacy speed is inertial, not rotating-planet relative" begin
    radius = BASE.Rp_e + 400e3
    speed = 7700.0
    for (azimuth, sign) in ((pi/2, 1.0), (-pi/2, -1.0))
        oe, r, v = _legacy_state([0.0, 0.0, 400e3], BASE, 0.0, azimuth, speed)
        expected_v = SVector(0.0, sign*speed, 0.0)
        @test r ≈ SVector(radius, 0.0, 0.0) atol=1e-6 rtol=0
        @test v ≈ expected_v atol=1e-8 rtol=0
        @test norm(v) ≈ speed atol=1e-8 rtol=0
        @test oe[1] ≈ inv(2/radius - speed^2/BASE.μ) rtol=1e-13
        @test oe[2] ≈ abs(speed^2*radius/BASE.μ - 1) atol=1e-13
        _, relative_interpretation = FT.r_pintor_i(SVector(radius, 0.0, 0.0), expected_v, BASE)
        @test relative_interpretation ≈ expected_v + SVector(0.0, BASE.ω[3]*radius, 0.0) atol=1e-10
        @test norm(relative_interpretation - v) ≈ BASE.ω[3]*radius atol=1e-8
        @test norm(relative_interpretation) ≈ speed + sign*BASE.ω[3]*radius atol=1e-8
        zero_spin = merge(BASE, (ω=zeros(SVector{3, Float64}),))
        zero_oe, zero_r, zero_v = _legacy_state([0.0, 0.0, 400e3], zero_spin, 0.0, azimuth, speed)
        @test isequal(zero_oe, oe)
        @test zero_r == r
        @test zero_v == v
        _, zero_relative = FT.r_pintor_i(SVector(radius, 0.0, 0.0), expected_v, zero_spin)
        @test zero_relative ≈ v atol=1e-8
    end
end

@testset "Legacy local angles and caller-selected epoch orientation" begin
    for latitude in (-0.6, 0.7), azimuth in (0.0, pi/2, -pi/2), γ in (-0.15, 0.2)
        longitude, height, speed = 0.4, 250e3, 7900.0
        location = [latitude, longitude, height]
        planet = merge(BASE, (L_PI=TILTED,))
        # CSPICE constructs the ellipsoidal location independently, without kernels.
        r_p = SVector{3, Float64}(georec(longitude, latitude, height,
            planet.Rp_e, (planet.Rp_e-planet.Rp_p)/planet.Rp_e))
        up = SVector(cos(latitude)*cos(longitude), cos(latitude)*sin(longitude), sin(latitude))
        east = normalize(cross(SVector(0.0, 0.0, 1.0), up))
        north = cross(up, east)
        resolved = speed * (cos(γ)*cos(azimuth)*north + cos(γ)*sin(azimuth)*east + sin(γ)*up)
        _, r, v = _legacy_state(location, planet, γ, azimuth, speed)
        @test r ≈ TILTED' * r_p atol=1e-5 rtol=0
        @test v ≈ TILTED' * resolved atol=1e-8 rtol=0
        @test norm(v) ≈ speed atol=1e-8 rtol=0
        @test dot(TILTED*v, up)/speed ≈ sin(γ) atol=2e-14
        @test dot(TILTED*v, east)/speed ≈ cos(γ)*sin(azimuth) atol=2e-14
        @test dot(TILTED*v, north)/speed ≈ cos(γ)*cos(azimuth) atol=2e-14
        _, with_transport = FT.r_pintor_i(r_p, resolved, planet)
        @test with_transport - v ≈ TILTED' * cross(planet.ω, r_p) atol=1e-8
        @test norm(with_transport - v) > 100.0
        identity_oe, identity_r, identity_v = _legacy_state(location, BASE, γ, azimuth, speed)
        @test r ≈ TILTED' * identity_r atol=1e-5 rtol=0
        @test v ≈ TILTED' * identity_v atol=1e-8 rtol=0
    end
end
end
