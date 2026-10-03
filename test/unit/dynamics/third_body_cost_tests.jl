module ThirdBodyCostTests

# Guards for the per-evaluation (time-only) environment lookups that the
# third-body and solar-radiation-pressure effectors depend on.
#
# The Sun's and the Moon's positions, and the planet's inertial-to-fixed
# rotation, depend only on the evaluation time: every spacecraft in a
# constellation shares them within one derivative evaluation. The cost of
# getting that wrong is measured in docs/architecture/third_body_cost.md --
# at 256 spacecraft, one degree-50 harmonics effector plus a Sun/Moon third
# body, name resolution alone was ~10% of the constellation RHS because it ran
# once per spacecraft per evaluation instead of once per run.
#
# The tests below pin the two invariants that keep that work off the
# per-spacecraft path: the canonical SPICE query name is resolved once per
# distinct name and never allocates again, and the planet-frame lookup asks the
# ephemerides backend once per evaluation no matter how many spacecraft consume
# it. Both check the shared value against the uncached one, so a cache that
# returns something different fails rather than being quietly faster.

using Test
using LinearAlgebra
using StaticArrays
using SpaceAGORA
using SpaceAGORA.SimulationModel

const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PE = SM.DynamicEffectors.PerturbationEffectors
const EM = SM.EphemeridesModels

# Stand-in ephemerides backend that counts how many times the planet-frame
# lookup is actually performed. It answers from the analytic simple model, so
# the test needs no SPICE kernels.
struct CountingEphemeridesModel <: SM.AbstractEphemeridesModel
    inner::SimpleEphemeridesModel
    calls::Base.Threads.Atomic{Int}
end

CountingEphemeridesModel() = CountingEphemeridesModel(SimpleEphemeridesModel(), Base.Threads.Atomic{Int}(0))

EM.ephemerides_requires_spice(::CountingEphemeridesModel) = false
EM.ephemerides_cache_key(::CountingEphemeridesModel) = (:counting,)
EM.ephemerides_time_seconds(initial_time, model::CountingEphemeridesModel) =
    EM.ephemerides_time_seconds(initial_time, model.inner)

function EM.planet_frame_lpi(planet, et::Float64, model::CountingEphemeridesModel)
    Base.Threads.atomic_add!(model.calls, 1)
    return EM.planet_frame_lpi(planet, et, model.inner)
end

EM.ephemerides_sun_direction_ii(planet, et::Float64, model::CountingEphemeridesModel) =
    EM.ephemerides_sun_direction_ii(planet, et, model.inner)

# A four-coefficient degree-2 field, written to a temporary file because the
# harmonics model only has a file-backed constructor. The values are a J2-only
# Earth field; nothing here depends on their magnitude.
function write_degree2_coefficients(path::String)
    open(path, "w") do io
        println(io, "degree,order,C,S")
        println(io, "0,0,1.0,0.0")
        println(io, "1,0,0.0,0.0")
        println(io, "1,1,0.0,0.0")
        println(io, "2,0,-1.08262668e-3,0.0")
        println(io, "2,1,0.0,0.0")
        println(io, "2,2,1.5745e-6,-9.0387e-7")
    end
    return path
end

function counting_configuration(planet, harmonics, n_sats::Int, ephemerides)
    spacecraft = SpacecraftModel[]
    for i in 1:n_sats
        root = Link(root=true, m=500.0, ref_area=12.0)
        ic = InitialCondition(
            ra=planet.Rp_e + 550e3 + 1e3 * i,
            rp=planet.Rp_e + 500e3 + 1e3 * i,
            i=53.0, ω=24.0, Ω=10.0, ν=40.0 * i
        )
        push!(spacecraft, SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, i))
    end
    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime, keplerian=true, number_of_orbits=1,
            mission_time=10.0, orientation_sim=false, num_steps_to_save=10, data_rate=1.0
        ),
        environment_model=EnvironmentModel(
            planet=planet, EI=120.0, density_model=NoAtmosphereModel(),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false, wind=false, ephemerides_model=ephemerides
        ),
        dynamics_model=DynamicsModel(spacecraft, (harmonics,)),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(
            reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=20.0
        )
    )
end

@testset "third body and SRP per-evaluation lookups" begin

    @testset "canonical SPICE query names are resolved once, not per spacecraft" begin
        names = ["Sun", "Moon", "Earth", "sun", " Mars ", "jupiter", "PLUTO",
                 "saturn_barycenter", "Solar System Barycenter"]

        # The interned value is exactly what the uncached expression produces.
        for name in names
            @test PE._spice_query_name(name) == PE._spice_query_name_uncached(name)
        end

        # Bodies whose own SPK is usually absent still resolve to the
        # barycenter form, and an already-resolved name is a fixed point.
        @test PE._spice_query_name("Mars") == "mars_barycenter"
        @test PE._spice_query_name("mars_barycenter") == "mars_barycenter"
        @test PE._spice_query_name("Solar System Barycenter") == "solar_system_barycenter"

        # Repeat resolution returns the same object, so the third-body sample
        # of spacecraft 4096 costs no more than that of spacecraft 1.
        first_call = PE._spice_query_name("Moon")
        @test PE._spice_query_name("Moon") === first_call

        # ...and costs nothing at all: the per-spacecraft path must not
        # allocate. 4096 resolutions stand in for a full constellation's worth
        # of per-spacecraft third-body samples in one derivative evaluation.
        resolve_many(name, n) = begin
            last = PE._spice_query_name(name)
            for _ in 2:n
                last = PE._spice_query_name(name)
            end
            last
        end
        resolve_many("Sun", 2)  # warm up compilation and the intern table
        @test @allocated(resolve_many("Sun", 4096)) == 0

        # Resolution runs on every RHS worker, so concurrent callers must all
        # get the correct answer (single-threaded sessions still exercise the
        # loop, just without the concurrency).
        concurrent = Vector{String}(undef, length(names) * 8)
        Threads.@threads for idx in eachindex(concurrent)
            concurrent[idx] = PE._spice_query_name(names[1 + (idx - 1) % length(names)])
        end
        for idx in eachindex(concurrent)
            @test concurrent[idx] == PE._spice_query_name_uncached(names[1 + (idx - 1) % length(names)])
        end
    end

    @testset "planet-frame lookup is consulted once per evaluation, not once per spacecraft" begin
        planet = Earth()
        mktempdir() do dir
            coefficients = write_degree2_coefficients(joinpath(dir, "degree2.csv"))
            harmonics = GravitationalHarmonicsModel(2, 2, coefficients, planet)
            n_sats = 64
            ephemerides = CountingEphemeridesModel()
            args = counting_configuration(planet, harmonics, n_sats, ephemerides)
            p = SM.ODEParams(n_sats=n_sats, args=args)

            et = 1.2345e6
            reference = EM.planet_frame_lpi(planet, et, SimpleEphemeridesModel())
            baseline_calls = ephemerides.calls[]

            # One derivative evaluation: every spacecraft asks for the frame at
            # the same et.
            lpi_values = [PE._harmonics_lpi_at!(harmonics, p, et) for _ in 1:n_sats]

            @test ephemerides.calls[] - baseline_calls == 1
            # The shared value is the uncached one, bit for bit.
            for lpi in lpi_values
                @test lpi == reference
            end

            # The next evaluation is a different time, so it costs exactly one
            # more lookup however many spacecraft consume it.
            et2 = et + 20.0
            reference2 = EM.planet_frame_lpi(planet, et2, SimpleEphemeridesModel())
            for _ in 1:n_sats
                @test PE._harmonics_lpi_at!(harmonics, p, et2) == reference2
            end
            @test ephemerides.calls[] - baseline_calls == 2
            @test reference2 != reference
        end
    end

    @testset "third-body positions read from the ephemeris cache are shareable across the constellation" begin
        # The interpolated cache is what a GRAM/SPICE-backed run reads inside
        # the RHS. Its result depends on (et, body) alone, so a constellation
        # that looks it up once per spacecraft is doing the same arithmetic N
        # times for the same answer -- the property this test pins is that the
        # answer really is identical, which is what makes hoisting it legal.
        ets = collect(0.0:600.0:6000.0)
        body_names = ["sun", "moon"]
        positions = Matrix{SVector{3, Float64}}(undef, length(ets), length(body_names))
        for (i, et) in pairs(ets), (k, _) in pairs(body_names)
            angle = 1e-6 * et + 0.7 * k
            radius = k == 1 ? 1.49e11 : 3.84e8
            positions[i, k] = SVector{3, Float64}(
                radius * cos(angle), radius * sin(angle), 0.1 * radius * sin(0.5 * angle)
            )
        end
        cache = SM.NBodyEphemerisCache(
            "earth",
            body_names,
            Dict(name => k for (k, name) in pairs(body_names)),
            ets,
            positions,
        )

        et = 2750.0  # between samples, so the interpolation actually runs
        reference = PE._nbody_body_position_from_cache_j2000_m(cache, et, "sun", "earth")
        @test reference !== nothing
        for _ in 1:512
            @test PE._nbody_body_position_from_cache_j2000_m(cache, et, "sun", "earth") == reference
        end

        # A cache built for another primary body must not answer at all, rather
        # than answering with the wrong frame.
        @test PE._nbody_body_position_from_cache_j2000_m(cache, et, "sun", "mars_barycenter") === nothing
        # Outside the cached span the caller has to fall back to the backend.
        @test PE._nbody_body_position_from_cache_j2000_m(cache, last(ets) + 1.0, "sun", "earth") === nothing
    end
end

end # module
