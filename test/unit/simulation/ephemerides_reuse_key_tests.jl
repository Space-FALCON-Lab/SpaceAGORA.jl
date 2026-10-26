module EphemeridesReuseKeyTests

# The planet-frame cache built at the start of a run is kept for later runs,
# keyed in part by the ephemerides model's `ephemerides_cache_key`. The engine
# used to turn every key other than `(:spice,)` into "simple:<2nd>:<3rd>", so a
# model of your own with a key shorter than three elements stopped the run at
# setup with a BoundsError, and two models whose keys differ only in the tag
# shared one cache. These tests pin the built-in strings and the whole-key
# route for keys of your own.

using Test, SpaceAGORA, StaticArrays
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const EM = SM.EphemeridesModels

# Answers from the analytic model under a key of its own. `turn_rad` turns the
# planet frame about the pole on top of the analytic one, so models with
# different turns give different frames.
struct KeyedEphemeridesModel <: SM.AbstractEphemeridesModel
    key::Tuple
    turn_rad::Float64
end
KeyedEphemeridesModel(key::Tuple) = KeyedEphemeridesModel(key, 0.0)

const ANALYTIC = SM.SimpleEphemeridesModel()

EM.ephemerides_requires_spice(::KeyedEphemeridesModel) = false
EM.ephemerides_cache_key(model::KeyedEphemeridesModel) = model.key
EM.ephemerides_epoch_seconds(::KeyedEphemeridesModel, initial_time) =
    EM.ephemerides_epoch_seconds(ANALYTIC, initial_time)
function EM.ephemerides_frame_lpi(model::KeyedEphemeridesModel, planet, et::Float64)
    c, s = cos(model.turn_rad), sin(model.turn_rad)
    turn = SMatrix{3, 3, Float64}(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0)
    return turn * EM.ephemerides_frame_lpi(ANALYTIC, planet, et)
end
EM.ephemerides_sun_ii(::KeyedEphemeridesModel, planet, et::Float64) =
    EM.ephemerides_sun_ii(ANALYTIC, planet, et)

@testset "Ephemerides reuse key: built-in keys keep their strings" begin
    @test SE._ephemerides_model_reuse_key(SM.SpiceEphemeridesModel()) == "spice"
    @test SE._ephemerides_model_reuse_key(SM.SimpleEphemeridesModel()) ==
        "simple:0:-9223372036854775808"
    @test SE._ephemerides_model_reuse_key(SM.SimpleEphemeridesModel(
        reference_epoch_seconds=1.5, prime_meridian_at_reference_rad=0.25)) ==
        "simple:1500000:250000000000"
    # A key of your own equal to a built-in one names the same cache: models
    # with equal keys give the same frames by contract.
    @test SE._ephemerides_model_reuse_key(KeyedEphemeridesModel((:spice,))) == "spice"
    @test SE._ephemerides_model_reuse_key(KeyedEphemeridesModel((:simple, 0, typemin(Int64)))) ==
        SE._ephemerides_model_reuse_key(SM.SimpleEphemeridesModel())
end

@testset "Ephemerides reuse key: keys of your own" begin
    reuse_key(key) = SE._ephemerides_model_reuse_key(KeyedEphemeridesModel(key))
    # Any length, the empty tuple included.
    @test reuse_key((:mine,)) == "(:mine,)"
    @test reuse_key((:mine, 1)) == "(:mine, 1)"
    @test reuse_key(()) == "()"
    @test reuse_key((:mine, 1, 2, 3.5, "x")) == "(:mine, 1, 2, 3.5, \"x\")"
    @test reuse_key((missing,)) == "(missing,)"
    # The tag counts: keys that differ only in the tag do not share, and a key
    # that ends like the analytic model's does not take the analytic model's.
    @test reuse_key((:tag_a, 1, 2)) != reuse_key((:tag_b, 1, 2))
    @test reuse_key((:tag_a, 1, 2)) != "simple:1:2"
    simple_key = EM.ephemerides_cache_key(ANALYTIC)
    @test reuse_key((:other, simple_key[2], simple_key[3])) != SE._ephemerides_model_reuse_key(ANALYTIC)
    # Distinct keys give distinct reuse keys.
    distinct_keys = (
        (:spice,), simple_key, ("spice",), ("simple", simple_key[2], simple_key[3]),
        (:other, simple_key[2], simple_key[3]), (:simple, simple_key[2], simple_key[3], 0),
        (:mine,), (:mine, 1), (:other, 1), (:mine, "a:b"), (:mine, :a, :b), (),
    )
    @test allunique(distinct_keys)
    @test allunique(map(reuse_key, distinct_keys))
    # A three-element key tagged :simple with other element types is printed
    # whole, so its elements cannot run into one another.
    @test reuse_key((:simple, "1:2", 3)) != reuse_key((:simple, 1, "2:3"))
    # No key of your own other than the built-in ones prints as a built-in string.
    for key in ((:mine,), (:simple,), (:simple, 1), (:simple, 1, 2, 3), (:spice, 1), ("spice",))
        @test reuse_key(key) != "spice"
        @test !startswith(reuse_key(key), "simple:")
    end
    planet = SM.Earth()
    @test SE._planet_frame_ephemeris_reuse_key(planet, KeyedEphemeridesModel((:mine,)), 0.0, 10.0, 1.0) isa
        SE.PlanetFrameEphemerisReuseKey
    models = (KeyedEphemeridesModel((:mine,)), KeyedEphemeridesModel((:tag_a, 1, 2)),
        KeyedEphemeridesModel((:tag_b, 1, 2)), SM.SimpleEphemeridesModel(), SM.SpiceEphemeridesModel())
    @test allunique([SE._planet_frame_ephemeris_reuse_key(planet, model, 0.0, 10.0, 1.0) for model in models])
end

# A stand-in for the run parameters with the fields the planet-frame cache
# reads and writes.
frame_cache_params(planet, model) = (
    args=(environment_model=(planet=planet, ephemerides_model=model),
        simulation_settings=(verbose=false,)),
    shared_buffers=(planet_frame_ephemeris_cache=Ref{Any}(nothing),
        spice_runtime_counters=(planet_pxform_cache_build_calls=Threads.Atomic{Int64}(0),)),
)

const FRAME_PLANET = SM.Earth()
# A start time no other test uses, so the shared reuse store holds no entry for
# this span before the tests; each test removes the entries for this span.
const FRAME_SPAN = (et_start=4_321_000.5, mission_end_s=90.0, dt_s=30.0)

# The run-start planet-frame cache for `model`, built or reused as at the start
# of a run.
function frame_cache_for(model)
    p = frame_cache_params(FRAME_PLANET, model)
    SE._initialize_planet_frame_ephemeris_cache!(p, FRAME_SPAN.et_start, FRAME_SPAN.mission_end_s)
    return p.shared_buffers.planet_frame_ephemeris_cache[]
end

frames_match(cache, model) = all(
    cache.quaternions[i] == SM.dcm_to_quaternion(EM.ephemerides_frame_lpi(model, FRAME_PLANET, cache.ets[i]))
    for i in eachindex(cache.ets))

# Runs `body` with the planet-frame cache and its reuse on, then removes this
# span's entries from the shared reuse store.
function with_frame_cache(body)
    span_key = (SE._cache_time_key(FRAME_SPAN.et_start), SE._cache_time_key(FRAME_SPAN.mission_end_s),
        SE._cache_time_key(FRAME_SPAN.dt_s))
    try
        withenv(body, "SPACEAGORA_PLANET_FRAME_CACHE" => "1",
            "SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "1",
            "SPACEAGORA_PLANET_FRAME_CACHE_DT_S" => string(FRAME_SPAN.dt_s),
            "SPACEAGORA_PLANET_FRAME_CACHE_MAX_SAMPLES" => "100",
            "SPACEAGORA_EPHEMERIS_CACHE_REUSE_MAX_ENTRIES" => "1000")
    finally
        lock(SE._EPHEMERIS_REUSE_LOCK) do
            filter!(entry -> entry.first[4:6] != span_key, SE._PLANET_FRAME_EPHEMERIS_REUSE_CACHE)
        end
    end
end

@testset "Planet-frame cache with a one-element key" begin
    one_element = KeyedEphemeridesModel((:one_element,))
    with_frame_cache() do
        cache = frame_cache_for(one_element)
        @test cache isa SM.PlanetFrameEphemerisCache
        @test length(cache.ets) == 4
        @test frames_match(cache, one_element)
        # A later run with the same key reuses it.
        @test frame_cache_for(one_element) === cache
    end
end

@testset "Planet-frame cache with keys that differ only in the tag" begin
    tag_a = KeyedEphemeridesModel((:tag_a, 1, 2), 0.0)
    tag_b = KeyedEphemeridesModel((:tag_b, 1, 2), 0.3)
    with_frame_cache() do
        cache_a = frame_cache_for(tag_a)
        cache_b = frame_cache_for(tag_b)
        @test cache_b !== cache_a
        @test frames_match(cache_a, tag_a)
        @test frames_match(cache_b, tag_b)
        @test cache_b.quaternions != cache_a.quaternions
        @test frame_cache_for(tag_a) === cache_a
        @test frame_cache_for(tag_b) === cache_b
        # The built-in analytic model keeps reusing its own cache beside them.
        analytic = frame_cache_for(ANALYTIC)
        @test analytic !== cache_a
        @test frames_match(analytic, ANALYTIC)
        @test frame_cache_for(ANALYTIC) === analytic
        # A model whose key ends like the analytic model's builds its own frames
        # instead of reading the analytic model's cache.
        simple_key = EM.ephemerides_cache_key(ANALYTIC)
        other = KeyedEphemeridesModel((:other, simple_key[2], simple_key[3]), 0.3)
        cache_other = frame_cache_for(other)
        @test cache_other !== analytic
        @test frames_match(cache_other, other)
        @test frame_cache_for(ANALYTIC) === analytic
    end
end

function simulation_args(ephemerides)
    planet = SM.Earth()
    root = SM.Link(root=true, m=500.0, ref_area=2.0)
    ic = SM.InitialCondition(ra=planet.Rp_e + 500e3, rp=planet.Rp_e + 500e3,
        i=35.0, ω=40.0, Ω=10.0, ν=120.0)
    sc = SM.SpacecraftModel(SM.Joint[], [root], root, true, 500.0, 0.0,
        root.inertia, 0, 0, ic, 1)
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime,
            keplerian=true, number_of_orbits=1, mission_time=30.0,
            orientation_sim=false, num_steps_to_save=10),
        environment_model=SM.EnvironmentModel(planet=planet, EI=120.0,
            density_model=SM.NoAtmosphereModel(),
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            ephemerides_model=ephemerides, topography=false, wind=false),
        dynamics_model=SM.DynamicsModel([sc], (SM.InverseSquaredGravityModel(),)),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances(),
        solver_config=SE.SolverConfig(solver_mode=:tsit5),
    )
end

@testset "Public simulation with a one-element ephemerides key" begin
    # The planet-frame cache and its reuse are on by default, so setup builds
    # the reuse key for every run.
    final(ephemerides) = withenv("SPACEAGORA_PLANET_FRAME_CACHE" => "1",
            "SPACEAGORA_EPHEMERIS_CACHE_REUSE" => "1") do
        collect(SE.run_simulation(simulation_args(ephemerides); return_solution=true).u[end])
    end
    # Frames and times equal to the analytic model's, so the run matches it.
    @test final(KeyedEphemeridesModel((:one_element,))) == final(ANALYTIC)
end

end
