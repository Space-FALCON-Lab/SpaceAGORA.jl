module FacetSRPAndThermalModelTests
# The per-facet SRP force model and the Sutton-Graves and tabular thermal
# models, checked against their closed forms. No SPICE or native GRAM: the
# SRP wrench gets a hand-built solar sample.
using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine

const AU = 149_597_870_700.0
const P1AU = 4.56e-6

function facet_spacecraft(facets::Vector{SM.Facet})
    root = SM.Link(root=true, m=500.0, ref_area=4.0)
    SM.add_facet!(root, facets)
    return SM.SpacecraftModel(links=SM.Link[root], root=root)
end

function srp_wrench(sc, pos; q_ib=SVector(0.0, 0.0, 0.0, 1.0), sun=SVector(AU, 0.0, 0.0))
    x = SM.StateSample(pos, SVector(0.0, 7.5e3, 0.0), 500.0; q_ib=q_ib, spacecraft=sc)
    env = SM.EnvironmentSample(SM.Earth(); solar=SM.SolarEphemerisSample(sun))
    return SM.wrench(SM.FacetSolarRadiationPressureModel(), x, env, 0.0)
end

# Sunlit, off the Sun line: the spacecraft sits on +y, the Sun on +x.
const POS = SVector(0.0, 7.0e6, 0.0)
s_hat(pos=POS, sun=SVector(AU, 0.0, 0.0)) = normalize(sun - pos)
pressure(pos=POS, sun=SVector(AU, 0.0, 0.0)) = P1AU * (AU / norm(sun - pos))^2

@testset "FacetSolarRadiationPressureModel" begin
    s = s_hat(); P = pressure(); n = SVector(1.0, 0.0, 0.0); c = dot(n, s)

    @testset "absorbing facet facing the Sun: F = -P A cos(theta) s" begin
        sc = facet_spacecraft([SM.Facet(area=2.0, normal_vector=[1.0, 0.0, 0.0], δ=0.0, ρ=0.0)])
        F, τ = srp_wrench(sc, POS)
        @test F ≈ -P * 2.0 * c * s rtol=1e-12
        @test τ == SVector(0.0, 0.0, 0.0)              # centre of pressure at the CoM
    end

    @testset "fully specular facet: F = -2 P A cos^2(theta) n" begin
        sc = facet_spacecraft([SM.Facet(area=2.0, normal_vector=[1.0, 0.0, 0.0], δ=1.0, ρ=0.0)])
        F, _ = srp_wrench(sc, POS)
        @test F ≈ -2.0 * P * 2.0 * c^2 * n rtol=1e-12
    end

    @testset "diffuse term: 2 (rho/3) n" begin
        sc = facet_spacecraft([SM.Facet(area=1.0, normal_vector=[1.0, 0.0, 0.0], δ=0.0, ρ=0.6)])
        F, _ = srp_wrench(sc, POS)
        @test F ≈ -P * c * (s + 2.0 * (0.6 / 3.0) * n) rtol=1e-12
    end

    @testset "facets facing away, and spacecraft in shadow, feel nothing" begin
        away = facet_spacecraft([SM.Facet(area=2.0, normal_vector=[-1.0, 0.0, 0.0])])
        @test srp_wrench(away, POS) == (SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0))
        lit = facet_spacecraft([SM.Facet(area=2.0, normal_vector=[1.0, 0.0, 0.0])])
        behind_earth = SVector(-7.0e6, 0.0, 0.0)
        @test srp_wrench(lit, behind_earth) == (SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0))
        @test srp_wrench(facet_spacecraft(SM.Facet[]), POS) == (SVector(0.0, 0.0, 0.0), SVector(0.0, 0.0, 0.0))
    end

    @testset "offset centre of pressure gives tau = r x F (body = J2000 here)" begin
        sc = facet_spacecraft([SM.Facet(area=2.0, normal_vector=[1.0, 0.0, 0.0], cp=[0.0, 0.5, 0.2])])
        F, τ = srp_wrench(sc, POS)
        @test τ ≈ cross(SVector(0.0, 0.5, 0.2), F) rtol=1e-12
    end

    @testset "the attitude orients the facet" begin
        # Rotate the body 90 deg about z: body x then points along inertial -y
        # under the scalar-last, inertial-to-body `rot` convention, so the
        # facet turns away from the Sun (which is along +x).
        q = SVector(0.0, 0.0, sin(pi / 4), cos(pi / 4))
        n_ii = SM.rot(q)' * SVector(1.0, 0.0, 0.0)
        sc = facet_spacecraft([SM.Facet(area=2.0, normal_vector=[1.0, 0.0, 0.0])])
        F, _ = srp_wrench(sc, POS; q_ib=q)
        cq = dot(n_ii, s)
        @test F ≈ (cq > 0 ? -P * 2.0 * cq * s : SVector(0.0, 0.0, 0.0)) atol=1e-20
    end

    @testset "without a propagated attitude: force, no torque" begin
        sc = facet_spacecraft([SM.Facet(area=2.0, normal_vector=[1.0, 0.0, 0.0], cp=[0.0, 0.5, 0.0])])
        F, τ = srp_wrench(sc, POS; q_ib=nothing)
        @test F ≈ -P * 2.0 * c * s rtol=1e-12
        @test τ == SVector(0.0, 0.0, 0.0)
    end

    @testset "construction, environment request, engine wiring" begin
        @test_throws ArgumentError SM.FacetSolarRadiationPressureModel(AU_m=0.0)
        @test_throws ArgumentError SM.FacetSolarRadiationPressureModel(p_srp_1au=-1.0)
        @test SM.environment_requirements(SM.FacetSolarRadiationPressureModel()).solar
        @test SE._wrench_method_available(SM.FacetSolarRadiationPressureModel())
    end
end

@testset "SuttonGravesHeat" begin
    earth = SM.Earth()
    m = SM.SuttonGravesHeat(planet=earth, nose_radius_m=0.5)
    @test m.k == earth.k
    ρ, v = 1.0e-6, 7.5e3
    @test SM.getHeatRate(m, 20.0, 200.0, ρ, v, 0.3) ≈ earth.k * sqrt(ρ / 0.5) * v^3 * 1e-4 rtol=1e-14
    @test SM.getHeatRate(m, 20.0, 200.0, 0.0, v, 0.3) == 0.0
    @test SM.SuttonGravesHeat(planet=earth, nose_radius_m=0.5, k=2.0e-4).k == 2.0e-4
    @test_throws ArgumentError SM.SuttonGravesHeat(planet=earth, nose_radius_m=0.0)
end

@testset "TabularHeat" begin
    vs = [5000.0, 7000.0]
    ds = [1.0e-8, 1.0e-6]
    q = [1.0 3.0; 2.0 6.0]                       # q[i, j] at (vs[i], ds[j])
    m = SM.TabularHeat(vs, ds, q)
    h(ρ, v) = SM.getHeatRate(m, 20.0, 200.0, ρ, v, 0.3)

    @test h(1.0e-8, 5000.0) == 1.0               # grid points are exact
    @test h(1.0e-6, 7000.0) == 6.0
    @test h(1.0e-8, 6000.0) ≈ 1.5                # linear in velocity
    @test h(1.0e-7, 5000.0) ≈ 2.0                # linear in log density (1e-7 is mid-way)
    @test h(1.0e-6, 9000.0) == 6.0               # held at the velocity edge
    @test h(1.0e-4, 7000.0) == 6.0               # held above the highest density
    @test h(1.0e-9, 5000.0) ≈ 0.1                # proportional to density below the grid
    @test h(0.0, 5000.0) == 0.0

    @testset "CSV database" begin
        mktempdir() do dir
            path = joinpath(dir, "aerothermal.csv")
            open(path, "w") do io
                println(io, "velocity_m_s,density_kg_m3,heat_rate_W_cm2")
                for (j, d) in enumerate(ds), (i, v) in reverse(collect(enumerate(vs)))
                    println(io, v, ",", d, ",", q[i, j])
                end
            end
            mf = SM.TabularHeat(path)
            @test mf.heat_rates == q && mf.velocities == vs && mf.densities == ds
            open(io -> println(io, "velocity_m_s,density_kg_m3,heat_rate_W_cm2\n5000,1e-8,1.0\n7000,1e-6,6.0"), path, "w")
            @test_throws ArgumentError SM.TabularHeat(path)        # incomplete grid
            open(io -> println(io, "velocity_m_s,density_kg_m3\n5000,1e-8"), path, "w")
            @test_throws ArgumentError SM.TabularHeat(path)        # missing column
        end
    end

    @test_throws ArgumentError SM.TabularHeat(vs, ds, q[:, 1:1])
    @test_throws ArgumentError SM.TabularHeat([7000.0, 5000.0], ds, q)
    @test_throws ArgumentError SM.TabularHeat(vs, [0.0, 1.0e-6], q)
end

@testset "thermal models pass the engine's support check" begin
    earth = SM.Earth()
    for model in (SM.SuttonGravesHeat(planet=earth, nose_radius_m=0.5),
                  SM.TabularHeat([5000.0, 7000.0], [1.0e-8, 1.0e-6], [1.0 3.0; 2.0 6.0]))
        @test hasmethod(SM.getHeatRate, Tuple{typeof(model), Float64, Float64, Float64, Float64, Float64})
    end
end

end # module
