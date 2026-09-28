module EpochConversionTests
# The engine's InitialTime -> ephemeris-time conversion must land on the
# requested UTC instant. Before this test, the SPICE route went through an
# AstroTime TAI epoch and back to a string; before 1972 the fractional TAI-UTC
# offset made that round trip inexact and the string came out as
# "1969-07-20T20:05:04.1000", which SPICE read 0.9 s early (the Apollo 11 PDI
# epoch). The cases below are the ones that reproduced it plus their neighbours.
using Test
using SpaceAGORA
using Dates
using SPICE: str2et

const SM_EPOCH = SpaceAGORA.SimulationModel
const EPH = SM_EPOCH.EphemeridesModels
const EPOCH_REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const EPOCH_SPICE_PATH = get(ENV, "SPACEAGORA_SPICE_PATH",
    joinpath(EPOCH_REPO, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE"))
const EPOCH_LSK = joinpath(EPOCH_SPICE_PATH, "lsk", "naif0012.tls")
const EPOCH_SPICE_READY = isfile(EPOCH_LSK)

# A model that carries the SPICE model's name but is not the package type: the
# engine's flexible resolver has no `ephemerides_time_seconds` method for it and
# takes its own SPICE branch, which must agree with the package method.
module DuckTypes
struct SpiceEphemeridesModel end
end

_it(y, mo, d, h, mi, s) = SM_EPOCH.InitialTime(year=y, month=mo, day=d, hour=h, minute=mi, second=Float32(s))

# Requested UTC strings and the InitialTime fields that name the same instant.
const EPOCH_CASES = (
    (_it(1969, 7, 20, 20, 5, 5.0), "1969-07-20T20:05:05"),      # Apollo 11 PDI: read as 20:05:04.100 before the fix
    (_it(1969, 7, 20, 20, 5, 59.0), "1969-07-20T20:05:59"),     # same failure one second before a minute rollover
    (_it(1969, 7, 20, 17, 44, 0.0), "1969-07-20T17:44:00"),     # an exact minute was unaffected
    (_it(1969, 7, 20, 20, 5, 5.5), "1969-07-20T20:05:05.5"),    # a half second was unaffected
    (_it(1971, 12, 31, 23, 59, 30.0), "1971-12-31T23:59:30"),   # last day of the fractional-offset era
    (_it(1971, 12, 31, 23, 59, 59.0), "1971-12-31T23:59:59"),   # fractional-offset boundary also changed the old round trip
    (_it(1972, 1, 1, 0, 0, 5.0), "1972-01-01T00:00:05"),        # first day of integer leap seconds
    (_it(2025, 6, 6, 0, 0, 0.0), "2025-06-06T00:00:00"),
)

@testset "EpochConversion" begin
    @testset "millisecond UTC string" begin
        @test EPH._initial_time_utc_string(_it(1969, 7, 20, 20, 5, 5.0)) == "1969-07-20T20:05:05.000"
        @test EPH._initial_time_utc_string(_it(1969, 7, 20, 17, 44, 0.0)) == "1969-07-20T17:44:00.000"
        @test EPH._initial_time_utc_string(_it(2025, 6, 6, 0, 0, 0.0)) == "2025-06-06T00:00:00.000"
        @test EPH._initial_time_utc_string(_it(1969, 7, 20, 20, 5, 5.5)) == "1969-07-20T20:05:05.500"
        # Sub-millisecond seconds round to the engine's millisecond clock, and a
        # rounded minute rolls over instead of printing second 60.
        @test EPH._initial_time_utc_string(_it(1969, 7, 20, 20, 5, 59.9996)) == "1969-07-20T20:06:00.000"
        @test EPH._initial_time_utc_string(_it(1999, 12, 31, 23, 59, 59.9999)) == "2000-01-01T00:00:00.000"
        # The simple model's clock reads the same DateTime.
        @test SM_EPOCH.ephemerides_time_seconds(_it(2000, 1, 1, 12, 0, 5.0), SM_EPOCH.SimpleEphemeridesModel()) == 5.0
    end

    if EPOCH_SPICE_READY
        @testset "SPICE ET equals the requested UTC" begin
            SM_EPOCH.Planets._furnsh_once(EPOCH_LSK)   # UTC conversion needs only the leap-second kernel
            spice = SM_EPOCH.SpiceEphemeridesModel()
            for (it, utc) in EPOCH_CASES
                requested = lock(SpaceAGORA.RuntimeServices.SPICE_LOCK) do
                    str2et(utc)
                end
                @test isapprox(SM_EPOCH.ephemerides_time_seconds(it, spice), requested; atol=1e-6, rtol=0.0)
                @test isapprox(SpaceAGORA.SimulationEngine._ephemerides_time_seconds_flexible(it, DuckTypes.SpiceEphemeridesModel()),
                               requested; atol=1e-6, rtol=0.0)
            end
        end
    else
        @info "Set SPACEAGORA_SPICE_PATH to a directory containing lsk/naif0012.tls to run the SPICE epoch regression" EPOCH_SPICE_PATH
        @test_skip EPOCH_SPICE_READY
    end
end
end # module
