module ApsisRotationContractTests
using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const OG = SM.OrbitalGeometry
const FT = SM.FrameTransforms
const SC = SM.SpacecraftModels

# Independent active elementary rotations: apply argument of periapsis, then
# inclination, then ascending node to a perifocal vector.
_zturn(v, a) = SVector(cos(a)*v[1] - sin(a)*v[2], sin(a)*v[1] + cos(a)*v[2], v[3])
_xturn(v, a) = SVector(v[1], cos(a)*v[2] - sin(a)*v[3], sin(a)*v[2] + cos(a)*v[3])
_inertial(v, i, ω, Ω) = _zturn(_xturn(_zturn(v, ω), i), Ω)

@testset "Shared orbital matrix retains direction and caller scaling" begin
    @test FT.OrbitalGeometry === OG
    @test SC.OrbitalGeometry === OG
    @test parentmodule(OG._perifocal_from_inertial) === OG
    x = SVector(1.0, 0.0, 0.0)
    y = SVector(0.0, 1.0, 0.0)
    @test OG._perifocal_from_inertial(0.0, 0.0, pi/2)' * x ≈ y atol=1e-15
    @test OG._perifocal_from_inertial(pi/2, 0.0, 0.0)' * y ≈ SVector(0.0, 0.0, 1.0) atol=1e-15
    @test_throws MethodError OG._perifocal_from_inertial(0.0f0, 0.0f0, 0.0f0)
    planet = SM.Earth()
    for (i, ω, Ω) in ((0.0, 0.0, 0.0), (pi/2, pi/2, 0.7), (2.4, -1.2, 3.1), (0.8, 7.1, -4.3))
        q = OG._perifocal_from_inertial(i, ω, Ω)
        @test q isa SMatrix{3, 3, Float64, 9}
        @test q * q' ≈ SMatrix{3, 3, Float64}(I) atol=1e-14
        @test det(q) ≈ 1.0 atol=1e-14
        for ν in (0.0, pi/2, Float64(pi), -0.6)
            perifocal = SVector(cos(ν), sin(ν), 0.0)
            independent = _inertial(perifocal, i, ω, Ω)
            direction = SC._initial_condition_apsis_direction_ii(i, ω, Ω, ν)
            @test direction isa SVector{3, Float64}
            @test direction ≈ independent atol=1e-14
            @test norm(direction) ≈ 1.0 atol=5e-16
            a, e = 8e6, 0.2
            p = a * (1 - e^2)
            oe = SVector(a, e, i, Ω, ω, ν, 0.0)
            position, velocity = FT.orbitalelemtorv(oe, planet)
            @test position isa Vector{Float64}
            @test velocity isa Vector{Float64}
            @test position ≈ p/(1+e*cos(ν)) * independent rtol=2e-15
            independent_velocity = sqrt(planet.μ/p) * _inertial(SVector(-sin(ν), e+cos(ν), 0.0), i, ω, Ω)
            @test velocity ≈ independent_velocity rtol=2e-15
            @test norm(cross(position, velocity)) ≈ sqrt(planet.μ*p) rtol=2e-15
            @test dot(velocity, velocity)/2 - planet.μ/norm(position) ≈ -planet.μ/(2a) rtol=2e-15
        end
    end
end
end
