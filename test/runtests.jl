# test/runtests.jl
#
# Minimal Julia test suite for SpaceAGORA.jl.
#
# Run with either:
#   julia --project=.SpaceAGORA -e 'using Pkg; Pkg.test()'
# or directly:
#   julia --project=.SpaceAGORA test/runtests.jl
# or, if a sysimage has been built:
#   julia --project=.SpaceAGORA --sysimage SpaceAGORA.so test/runtests.jl
#
# These tests use only built-in density and gravity models (no GRAM, no SPICE)
# so they can run in CI without external data files.

using Test
using StaticArrays

# ---------------------------------------------------------------------------
# Load source files once, at the top level, to avoid redefinition warnings
# when running multiple testsets.
# ---------------------------------------------------------------------------
const SRC = joinpath(@__DIR__, "..", "src")

include(joinpath(SRC, "physical_models", "Planet_data.jl"))
include(joinpath(SRC, "config.jl"))
import .config

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
@testset "SpaceAGORA.jl" begin

    @testset "Planet data" begin
        earth = planet_data(0)
        @test earth.μ > 0
        @test earth.Rp_e > earth.Rp_p          # equatorial > polar radius
        mars  = planet_data(1)
        @test mars.μ > 0
        @test mars.Rp_e > 0
        venus = planet_data(2)
        @test venus.μ > 0
    end

    @testset "Spacecraft model" begin
        sc  = config.SpacecraftModel()
        bus = config.Link(
            root     = true,
            r        = SVector{3, Float64}(0.0, 0.0, 0.0),
            ṙ        = SVector{3, Float64}([0, 0, 0]),
            dims     = SVector{3, Float64}([3.7, 2.05, 2.8]),
            ref_area = 10.0,
            m        = 690.0,
        )
        config.add_body!(sc, bus, prop_mass = 50.0)
        @test length(sc.links) == 1
        com = config.get_COM(sc, bus)
        @test length(com) == 3
        I   = config.get_inertia_tensor(sc, bus)
        @test size(I) == (3, 3)
        # Inertia tensor diagonal elements must be non-negative
        @test all(diag(I) .>= 0)
    end

end
