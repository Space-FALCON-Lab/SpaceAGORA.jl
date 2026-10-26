module J2SecularRateBoundaryTests

using Test
using SpaceAGORA
const rate = SpaceAGORA.SimulationModel.DynamicEffectors.GravityEffectors.j2_secular_rates
const UNIT_PLANET = (μ=1.0, Rp_e=1.0, J2=0.001)

@testset "J2 secular rates: inclination symmetries and scaling" begin
    # Dimensionless μ=a=R=1 makes n=1. Special-angle predictions are
    # independent closed-form values, including zero nodal/periapsis drift.
    @test collect(rate(1.0, 0.0, 0.0, UNIT_PLANET)) ≈ [-0.0015, 0.003]
    @test rate(1.0, 0.0, π/2, UNIT_PLANET)[1] ≈ 0.0 atol=1e-18
    @test rate(1.0, 0.0, π/2, UNIT_PLANET)[2] ≈ -0.00075
    @test collect(rate(1.0, 0.0, Float64(π), UNIT_PLANET)) ≈ [0.0015, 0.003]
    critical = acos(inv(sqrt(5.0)))
    @test rate(1.0, 0.0, critical, UNIT_PLANET)[1] ≈ -0.0015/sqrt(5.0)
    @test rate(1.0, 0.0, critical, UNIT_PLANET)[2] ≈ 0.0 atol=1e-18
    @test collect(rate(1.0, 0.0, π/3, UNIT_PLANET)) ≈ [-0.00075, 0.0001875]

    expected = [-0.00075, 0.0001875]
    @test collect(rate(4.0, 0.0, π/3, UNIT_PLANET)) ≈ expected/128 # a^(-7/2)
    @test collect(rate(1.0, 0.6, π/3, UNIT_PLANET)) ≈ expected*(625/256) # (1-e²)^(-2)
    @test collect(rate(1.0, 0.0, π/3, merge(UNIT_PLANET, (μ=4.0,)))) ≈ 2expected
    @test collect(rate(1.0, 0.0, π/3, merge(UNIT_PLANET, (Rp_e=2.0,)))) ≈ 4expected
    @test collect(rate(1.0, 0.0, π/3, merge(UNIT_PLANET, (J2=-0.001,)))) ≈ -expected
    for inclination in (0.1, 0.7, 1.2)
        prograde = rate(1.0, 0.2, inclination, UNIT_PLANET)
        retrograde = rate(1.0, 0.2, π-inclination, UNIT_PLANET)
        @test prograde[1] ≈ -retrograde[1]
        @test prograde[2] ≈ retrograde[2]
    end
end

@testset "J2 secular rates: invalid and degenerate elliptic inputs" begin
    for j2 in (0.0, NaN, Inf, -Inf)
        @test rate(1.0, 0.1, 0.3, merge(UNIT_PLANET, (J2=j2,))) === (0.0,0.0)
    end
    for a in (NaN, Inf, -Inf, 0.0, -1.0)
        @test rate(a, 0.1, 0.3, UNIT_PLANET) === (0.0,0.0)
    end
    for e in (NaN, Inf, -Inf, -0.1, 1.0, 1.1)
        @test rate(1.0, e, 0.3, UNIT_PLANET) === (0.0,0.0)
    end
    # Positive a and elliptic e can still underflow the semilatus rectum.
    @test rate(nextfloat(0.0), 0.9, 0.3, UNIT_PLANET) === (0.0,0.0)
    # A zero/undefined mean motion has no finite secular prediction.
    @test rate(1.0e-110, 0.1, 0.3, UNIT_PLANET) === (0.0,0.0)
    @test rate(1.0e110, 0.1, 0.3, UNIT_PLANET) === (0.0,0.0)
    for μ in (0.0, NaN, Inf)
        @test rate(1.0, 0.1, 0.3, merge(UNIT_PLANET, (μ=μ,))) === (0.0,0.0)
    end
end

end # module J2SecularRateBoundaryTests
