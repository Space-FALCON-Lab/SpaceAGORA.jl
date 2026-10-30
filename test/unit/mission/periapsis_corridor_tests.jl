using Test
using LinearAlgebra
using StaticArrays

include(joinpath(@__DIR__, "../../../examples/periapsis_corridor_support.jl"))
using .PeriapsisCorridorSupport

const PC = PeriapsisCorridorSupport

struct SyntheticDenseSolution
    t::Vector{Float64}
    u::Vector{Any}
    radius_reference_m::Float64
end

function (solution::SyntheticDenseSolution)(t::Real)
    radius = solution.radius_reference_m + 100_000.0 + (Float64(t) - 5.0)^2
    velocity = 2.0 * (Float64(t) - 5.0)
    return (sc=[(pos=SVector(radius, 0.0, 0.0), vel=SVector(velocity, 0.0, 0.0))],)
end

@testset "Earth periapsis corridor configuration and diagnostics" begin
    config = load_corridor_config()
    mission = config.values["mission"]
    @test mission["horizon_s"] == 5_184_000.0
    @test mission["corridor_min_altitude_m"] == 95_000.0
    @test mission["corridor_max_altitude_m"] == 100_000.0
    @test mission["altitude_reference"] == "geodetic height above the planet ellipsoid"
    @test mission["safety_margin_m"] == 250.0
    vehicle = config.values["vehicle"]
    reference_wet_mass = vehicle["dry_mass_kg"] + vehicle["propellant_mass_kg"]
    @test reference_wet_mass / vehicle["ballistic_coefficient_kg_m2"] == 0.82
    @test (reference_wet_mass / vehicle["ballistic_coefficient_kg_m2"]) /
        vehicle["drag_coefficient"] == 0.41
    @test config.values["comparison"]["max_burn_duration_s"] == 1800.0
    @test config.values["integration"]["solver_mode"] == "tsit5"

    @test config.values["calibration"]["tolerance_m"] == 1.0
    @test config.values["calibration"]["max_iterations"] == 3

    @test corridor_coordinate(95_000.0, 95_000.0, 100_000.0) == 0.0
    @test corridor_coordinate(97_500.0, 95_000.0, 100_000.0) == 0.5
    @test corridor_coordinate(100_000.0, 95_000.0, 100_000.0) == 1.0
    @test_throws ArgumentError corridor_coordinate(1.0, 3.0, 3.0)
end

@testset "Corridor piecewise exponential atmosphere fit" begin
    fit = fitted_atmosphere_parameters(load_corridor_config())
    nodes = fit.node_altitudes_m
    # corridor-band node densities are the PDAS 1976 table values
    @test corridor_density(fit, 95_000.0) == 1.3873e-6
    @test corridor_density(fit, 100_000.0) == 5.6044e-7
    # every node reproduces its table density
    for (h, rho) in zip(nodes, fit.node_density_kg_m3)
        @test corridor_density(fit, h) ≈ rho rtol = 1e-12
    end
    @test nodes[1] <= 80_000.0 && nodes[end] >= 150_000.0
    # continuous across each interior boundary
    for h in nodes[2:end-1]
        @test corridor_density(fit, prevfloat(h)) ≈ corridor_density(fit, nextfloat(h)) rtol = 1e-9
    end
    # strictly decreasing with altitude over the covered range
    grid = range(nodes[1], nodes[end]; length=2001)
    densities = corridor_density.(Ref(fit), grid)
    @test all(diff(densities) .< 0.0)
    # held-out table rows: the worst-case relative error is bounded and reported
    errors = atmosphere_fit_relative_errors(fit)
    held_out = [!(h in nodes) for h in fit.table_altitudes_m]
    @test count(held_out) == 7
    worst = maximum(abs, errors[held_out])
    println("atmosphere fit worst-case relative error at held-out table rows: ", worst)
    @test worst < 0.09
    @test all(==(0.0), errors[.!held_out]) || all(abs.(errors[.!held_out]) .< 1e-12)
    # explicit out-of-range behavior
    @test corridor_density(fit, nodes[end] + 10_000.0) < corridor_density(fit, nodes[end])
    @test corridor_density(fit, nodes[1] - 10_000.0) > corridor_density(fit, nodes[1])
    zero_fit = CorridorAtmosphereFit(fit.node_altitudes_m, fit.node_density_kg_m3,
        fit.scale_heights_m, :zero, :zero, fit.temperature_k, fit.source,
        fit.table_altitudes_m, fit.table_density_kg_m3)
    @test corridor_density(zero_fit, nodes[end] + 1.0) == 0.0
    @test corridor_density(zero_fit, nodes[1] - 1.0) == 0.0
    @test corridor_density(zero_fit, nodes[1]) == fit.node_density_kg_m3[1]
end

@testset "Earth corridor context preserves signed J3 and fixed bounds" begin
    bounds = corridor_eccentricity_bounds(36.0e6, 6_378_136.3 + 95_000.0, 6_378_136.3 + 100_000.0)
    @test bounds.minimum_eccentricity < bounds.maximum_eccentricity

    ctx = dimensionless_context(36.0e6, 0.82, deg2rad(50), 1.0826e-3,
                                -2.532e-6, 6_378_136.3, 6_378_136.3,
                                6_378_136.3 + 95_000.0, 6_378_136.3 + 100_000.0)
    @test ctx.Pi_J2 > 0.0
    @test ctx.Pi_J3 > 0.0
    @test ctx.signed_Pi_J3 < 0.0
    @test ctx.J3_sign == -1.0
    @test ctx.third_body_context == 0.0
    @test ctx.Gamma_3B_J2 == 0.0
    @test ctx.Pi_3B == 0.0
    @test ctx.reference_radius_over_a ≈ 6_378_136.3 / 36.0e6
    drag = PC.drag_context_ratio(5.6e-7, 7_800.0, 1_000.0,
                                 3.986004415e14, 6_475_000.0)
    @test drag.Pi_drag > 0.0
    @test drag.dynamic_pressure_pa > 0.0
    @test_throws ArgumentError corridor_eccentricity_bounds(6.4e6, 6.378e6 + 95_000.0, 6.378e6 + 100_000.0)
end

@testset "corridor eccentricity bounds allow osculating a below the corridor top" begin
    R, lo, hi = 6_378_136.6, 687_922.0, 697_922.0
    a = R + 697_000.0
    b = corridor_eccentricity_bounds(a, R + lo, R + hi)
    @test b.minimum_eccentricity == 0.0
    @test b.maximum_eccentricity ≈ 1 - (R + lo) / a
end

@testset "physical periapsis is located at a radial-velocity root" begin
    radius_reference = 6_378_136.3
    solution = SyntheticDenseSolution(
        [0.0, 4.0, 6.0, 10.0],
        Any[(sc=[(pos=SVector(radius_reference + 125_000.0, 0.0, 0.0),
                         vel=SVector(-10.0, 0.0, 0.0))],),
            (sc=[(pos=SVector(radius_reference + 101_000.0, 0.0, 0.0),
                         vel=SVector(-2.0, 0.0, 0.0))],),
            (sc=[(pos=SVector(radius_reference + 101_000.0, 0.0, 0.0),
                         vel=SVector(2.0, 0.0, 0.0))],),
            (sc=[(pos=SVector(radius_reference + 125_000.0, 0.0, 0.0),
                         vel=SVector(10.0, 0.0, 0.0))],)],
        radius_reference)
    # Spherical stand-in for the geodetic conversion: the root finding is what is tested here.
    geodetic_of(t, r) = (norm(r) - radius_reference, 0.0)
    events = physical_apsis_events(solution, geodetic_of)
    @test length(events) == 1
    @test events[1].kind == "periapsis"
    @test events[1].time_s ≈ 5.0 atol=1e-7
    @test events[1].physical_altitude_m ≈ 100_000.0 atol=1e-7
    @test abs(events[1].radial_velocity_residual_m2_s) < 1e-3

    # Solver callbacks may save the apsis itself as a segment endpoint. Dense
    # interpolation then has a zero endpoint rather than a strict sign bracket.
    endpoint_solution = SyntheticDenseSolution(
        [0.0, 5.0, 10.0],
        Any[(sc=[(pos=SVector(radius_reference + 125_000.0, 0.0, 0.0),
                         vel=SVector(-10.0, 0.0, 0.0))],),
            (sc=[(pos=SVector(radius_reference + 100_000.0, 0.0, 0.0),
                         vel=SVector(0.0, 0.0, 0.0))],),
            (sc=[(pos=SVector(radius_reference + 125_000.0, 0.0, 0.0),
                         vel=SVector(10.0, 0.0, 0.0))],)],
        radius_reference)
    endpoint_events = physical_apsis_events(endpoint_solution, geodetic_of)
    @test length(endpoint_events) == 1
    @test endpoint_events[1].kind == "periapsis"
    @test endpoint_events[1].time_s == 5.0
end

# Regression: the solver steps exactly onto the event, so the saved knot holds the
# root while the dense interpolant at that knot is off (same sign as its neighbor).
struct OffsetDenseSolution
    t::Vector{Float64}
    u::Vector{Any}
end
function (solution::OffsetDenseSolution)(t::Real)
    return (sc=[(pos=SVector(6.5e6, 0.0, 0.0), vel=SVector(0.5 + Float64(t), 0.0, 0.0))],)
end

@testset "periapsis at a saved knot survives a mismatched dense interpolant" begin
    saved(v) = (sc=[(pos=SVector(6.5e6, 0.0, 0.0), vel=SVector(v, 0.0, 0.0))],)
    solution = OffsetDenseSolution([0.0, 5.0, 10.0], Any[saved(-5.0), saved(-1e-12), saved(5.0)])
    events = physical_apsis_events(solution, (t, r) -> (norm(r) - 6.378e6, 0.0))
    @test length(events) == 1
    @test events[1].kind == "periapsis"
    @test events[1].time_s == 5.0
end
