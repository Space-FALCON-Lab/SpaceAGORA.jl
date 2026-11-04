using Test
using SpaceAGORA
using StaticArrays
using LinearAlgebra
using Dates
using SPICE

const SM_SUN = SpaceAGORA.SimulationModel
const SUN_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SUN_SPICE_PATH = joinpath(SUN_REPO, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")

# Seconds past the simple model's J2000 epoch (2000-01-01T12:00 UTC).
_sun_et(dt::DateTime) = Dates.value(dt - DateTime(2000, 1, 1, 12, 0, 0)) / 1000.0

# Declination of an equatorial unit vector, degrees.
_sun_declination_deg(d) = rad2deg(asin(clamp(d[3], -1.0, 1.0)))

@testset "SunDirection" begin
    @testset "simple ephemerides at Earth" begin
        earth = make_no_gram_planet(:earth)
        simple = SM_SUN.SimpleEphemeridesModel()

        # June solstice 2020: the Sun is over the northern tropic and its right
        # ascension is 6 h, so the direction is +y-ish with the obliquity in z.
        june = SM_SUN.ephemerides_sun_direction_ii(earth, _sun_et(DateTime(2020, 6, 20, 21, 44, 0)), simple)
        @test june !== nothing
        @test isapprox(norm(june), 1.0; atol=1e-12)
        @test _sun_declination_deg(june) ≈ 23.44 atol = 0.2
        @test june[2] > 0.9
        @test abs(june[1]) < 0.05

        # December solstice: the mirror image, the Sun south of the equator.
        december = SM_SUN.ephemerides_sun_direction_ii(earth, _sun_et(DateTime(2020, 12, 21, 10, 2, 0)), simple)
        @test _sun_declination_deg(december) ≈ -23.44 atol = 0.2
        @test december[2] < -0.9

        # March equinox: the Sun crosses the equator at the vernal equinox
        # itself, the +x axis of the frame.
        march = SM_SUN.ephemerides_sun_direction_ii(earth, _sun_et(DateTime(2020, 3, 20, 3, 50, 0)), simple)
        @test abs(_sun_declination_deg(march)) < 0.3
        @test march[1] > 0.99

        # September equinox: the opposite crossing, the Sun near -x.
        september = SM_SUN.ephemerides_sun_direction_ii(earth, _sun_et(DateTime(2020, 9, 22, 13, 31, 0)), simple)
        @test abs(_sun_declination_deg(september)) < 0.3
        @test september[1] < -0.99

        # The simple model carries no planetary ephemeris, so it declines to
        # invent a Sun for any other body.
        @test SM_SUN.ephemerides_sun_direction_ii(make_no_gram_planet(:mars), 0.0, simple) === nothing
        @test SM_SUN.ephemerides_sun_direction_ii(make_no_gram_planet(:venus), 0.0, simple) === nothing
    end

    # Apollo 11 landed with the Sun 10.8 degrees above the horizon at Tranquility
    # Base, low enough for the crew to read the surface relief by its shadows.
    @testset "SPICE ephemerides at the Moon" begin
        if !isdir(SUN_SPICE_PATH)
            @info "SPICE kernels absent; skipping the lunar sun elevation test" SUN_SPICE_PATH
            @test true
        else
            moon = SM_SUN.Moon("", SUN_SPICE_PATH)
            spice = SM_SUN.SpiceEphemeridesModel()
            et = utc2et("1969-07-20T20:05:05")
            direction = SM_SUN.ephemerides_sun_direction_ii(moon, et, spice)
            @test direction !== nothing
            @test isapprox(norm(direction), 1.0; atol=1e-12)

            # Rotate into the body-fixed frame and take the elevation over the
            # site's local vertical.
            lpi = SM_SUN.planet_frame_lpi(moon, et, spice)
            sun_pcpf = lpi * direction
            lat, lon = deg2rad(0.67416), deg2rad(23.47314)
            up = SVector{3, Float64}(cos(lat) * cos(lon), cos(lat) * sin(lon), sin(lat))
            elevation_deg = rad2deg(asin(clamp(dot(up, sun_pcpf), -1.0, 1.0)))
            @test elevation_deg ≈ 10.8 atol = 0.6
        end
    end
end
