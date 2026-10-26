module ThirdBodyEphemerisRequirementsTests

# A force or torque model asks for third-body positions through the public
# declaration `environment_requirements(model) =
# EffectorEnvironmentRequirements(third_body_names=...)` and reads them from
# `env.third_bodies` in `wrench`. The engine used to sample third bodies only
# for `NBodyGravityModel`, from that model's own `body_names`, so any other
# model following the declaration failed at its first sample, and the setup
# checks that pick the ephemerides backend and build the N-body ephemeris cache
# did not see it either. These tests pin the declaration route for a model the
# engine knows nothing about.

using Test, SpaceAGORA, StaticArrays, LinearAlgebra
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine
const PE = SM.DynamicEffectors.PerturbationEffectors

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const SPICE_PATH = joinpath(REPO, "data/GRAMSuite.jl/GRAM Suite 2.0", "SPICE")

const DECLARED_BODIES = ("Moon", "Sun")
# Gravitational parameters passed explicitly, so building an
# `NBodyGravityModel` for comparison needs no kernels.
const DECLARED_MUS = (4.9028e12, 1.32712440018e20)

# Declares two third bodies and reads them in `wrench`. The force depends on
# both positions and on their order, so a missing, reordered or wrongly
# referenced sample changes the result.
struct DeclaredThirdBodyEffector <: SM.AbstractForceTorqueModel end
SpaceAGORA.environment_requirements(::DeclaredThirdBodyEffector) =
    SpaceAGORA.EffectorEnvironmentRequirements(third_body_names=DECLARED_BODIES)
function SpaceAGORA.wrench(::DeclaredThirdBodyEffector, x::SpaceAGORA.StateSample,
        env::SpaceAGORA.EnvironmentSample, t::Float64)
    third_bodies = env.third_bodies
    third_bodies === nothing && throw(ArgumentError("env.third_bodies was not sampled"))
    moon_ii, sun_ii = third_bodies.positions_ii
    return 1e-12 * (moon_ii - x.pos_ii) + 1e-15 * sun_ii, SVector(0.0, 0.0, 0.0)
end

# The Moon's tidal acceleration on the spacecraft, computed by the built-in
# n-body kernel from the declared sample: a custom model that reproduces
# `NBodyGravityModel(body_names=("Moon",))` through the public route.
struct DeclaredMoonGravity <: SM.AbstractForceTorqueModel
    inner::SM.NBodyGravityModel
end
DeclaredMoonGravity() = DeclaredMoonGravity(
    SM.NBodyGravityModel(body_names=("Moon",), primary_body_name="Earth",
        planet=SM.Earth(), body_mus=(DECLARED_MUS[1],)))
SpaceAGORA.environment_requirements(::DeclaredMoonGravity) =
    SpaceAGORA.EffectorEnvironmentRequirements(third_body_names=("Moon",))
SpaceAGORA.wrench(model::DeclaredMoonGravity, x::SpaceAGORA.StateSample,
        env::SpaceAGORA.EnvironmentSample, t::Float64) = SpaceAGORA.wrench(model.inner, x, env, t)

function spacecraft(planet)
    root = SM.Link(root=true, m=500.0, ref_area=2.0)
    ic = SM.InitialCondition(ra=planet.Rp_e + 700e3, rp=planet.Rp_e + 500e3,
        i=35.0, ω=40.0, Ω=10.0, ν=120.0)
    return SM.SpacecraftModel(SM.Joint[], [root], root, true, 500.0, 0.0,
        root.inertia, 0, 0, ic, 1)
end

function simulation_args(effectors::Tuple; ephemerides=SM.SimpleEphemeridesModel(), mission_time=30.0,
        planet=SM.Earth())
    return SM.SimulationConfiguration(
        simulation_settings=SM.SimulationSettings(results=false, verbose=false,
            generate_plots=false, normalize=false, save_csv=false),
        mission_configuration=SM.MissionConfiguration(mission_type=SM.MissionTime,
            keplerian=true, number_of_orbits=1, mission_time=mission_time,
            orientation_sim=false, num_steps_to_save=10),
        environment_model=SM.EnvironmentModel(planet=planet, EI=120.0,
            density_model=SM.NoAtmosphereModel(),
            thermal_model=SM.MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            ephemerides_model=ephemerides, topography=false, wind=false),
        dynamics_model=SM.DynamicsModel([spacecraft(planet)], effectors),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances(reltol_orbit=1e-12, abstol_orbit=1e-12),
        solver_config=SE.SolverConfig(solver_mode=:tsit5),
    )
end

# An interpolated third-body cache measured from the Earth, built from an
# analytic curve so the samplers below read positions without kernels.
function seeded_cache(body_names)
    query_names = [PE._spice_query_name(name) for name in body_names]
    ets = collect(0.0:600.0:6000.0)
    positions = Matrix{SVector{3, Float64}}(undef, length(ets), length(query_names))
    for (i, et) in pairs(ets), k in eachindex(query_names)
        angle = 1e-6 * et + 0.7 * k
        radius = query_names[k] == "sun" ? 1.49e11 : 3.84e8
        positions[i, k] = SVector(radius * cos(angle), radius * sin(angle), 0.1 * radius * sin(0.5 * angle))
    end
    return SM.NBodyEphemerisCache(PE._spice_query_name("Earth"), query_names,
        Dict(name => k for (k, name) in pairs(query_names)), ets, positions)
end

function seeded_params(effectors::Tuple)
    args = simulation_args(effectors)
    p = SM.ODEParams(n_sats=1, args=args)
    p.shared_buffers.nbody_ephemeris_cache[] = seeded_cache(DECLARED_BODIES)
    p.shared_buffers.et_start[] = 0.0
    return p, SE.build_initial_conditions(args)
end

bits(v::SVector{3, Float64}) = ntuple(i -> reinterpret(UInt64, v[i]), 3)

@testset "Declared third bodies are sampled for any model" begin
    model = DeclaredThirdBodyEffector()
    p, u = seeded_params((SM.InverseSquaredGravityModel(), model))
    req = SpaceAGORA.environment_requirements(model)
    cache = p.shared_buffers.nbody_ephemeris_cache[]
    calls = p.shared_buffers.spice_runtime_counters.nbody_spkpos_runtime_calls
    state = SE.build_state_sample(u.sc[1], p.args.dynamics_model.spacecraft[1], false, 1)
    for t in (0.0, 1234.5, 2750.0)
        expected = map(DECLARED_BODIES) do name
            PE._nbody_body_position_from_cache_j2000_m(cache, t,
                PE._spice_query_name(name), PE._spice_query_name("Earth"))
        end
        for env in (SE.sample_environment(req, model, u.sc[1], p, 1, t),
                    SE.sample_environment_with_reusable_buffers(req, model, u.sc[1], p, 1, t))
            @test env.third_bodies isa SpaceAGORA.ThirdBodyEphemerisSample{2}
            @test env.third_bodies.names === DECLARED_BODIES
            for k in 1:2
                @test bits(env.third_bodies.positions_ii[k]) == bits(expected[k])
            end
        end
        # Through the engine's effector evaluation, as the right-hand side
        # calls it: the wrench reads the sampled positions.
        force, torque = SE._evaluate_dynamic_effector(model, u.sc[1], state, p, 1, t)
        @test force == 1e-12 * (expected[1] - state.pos_ii) + 1e-15 * expected[2]
        @test torque == SVector(0.0, 0.0, 0.0)
    end
    # Served from the cache measured from the run's planet, not from SPICE.
    @test calls[] == 0
end

@testset "NBodyGravityModel and a declaring model read the same sample" begin
    nbody = SM.NBodyGravityModel(body_names=DECLARED_BODIES, primary_body_name="Earth",
        planet=SM.Earth(), body_mus=DECLARED_MUS)
    model = DeclaredThirdBodyEffector()
    p, u = seeded_params((SM.InverseSquaredGravityModel(), nbody, model))
    for t in (0.0, 2750.0)
        from_nbody = SE.sample_third_body_ephemerides(nbody, u.sc[1], p, 1, t)
        from_declared = SE.sample_environment(SpaceAGORA.environment_requirements(model),
            model, u.sc[1], p, 1, t).third_bodies
        nbody_env = SE.sample_environment_with_reusable_buffers(
            SpaceAGORA.environment_requirements(nbody), nbody, u.sc[1], p, 1, t).third_bodies
        @test from_nbody.names === nbody.body_names
        @test from_declared.names == from_nbody.names
        for k in 1:2
            @test bits(from_declared.positions_ii[k]) == bits(from_nbody.positions_ii[k])
            @test bits(nbody_env.positions_ii[k]) == bits(from_nbody.positions_ii[k])
        end
    end
end

validation_args(effectors, backend=SM.SimpleEphemeridesModel()) = (
    dynamics_model=(dynamic_effectors=effectors,),
    environment_model=(planet=SM.Earth(), ephemerides_model=backend),
    simulation_settings=(verbose=false,),
)

@testset "Setup sees declared third bodies" begin
    declared = DeclaredThirdBodyEffector()
    nbody = SM.NBodyGravityModel(body_names=("Sun", "Jupiter"), primary_body_name="Earth",
        planet=SM.Earth(), body_mus=(DECLARED_MUS[2], 1.26686534e17))
    @test SE._has_active_nbody_effector((SM.InverseSquaredGravityModel(), declared))
    @test !SE._has_active_nbody_effector((SM.InverseSquaredGravityModel(),))
    @test SE._collect_nbody_query_names((declared,)) == ["moon", "sun"]
    # NBodyGravityModel's bodies keep their order; declared ones join once.
    @test SE._collect_nbody_query_names((nbody, declared)) ==
        [PE._spice_query_name("Sun"), PE._spice_query_name("Jupiter"), "moon"]

    # The analytic backend has no planetary ephemeris: refuse the run before
    # integration, as for NBodyGravityModel.
    err = try
        SE._validate_ephemerides_support!(validation_args((declared,)))
        nothing
    catch ex
        ex
    end
    @test err isa ArgumentError
    @test occursin("third_body_names", sprint(showerror, err))
    @test SE._validate_ephemerides_support!(validation_args((declared,), SM.SpiceEphemeridesModel())) === nothing

    # The run-start cache is built for the declared bodies as well. Seed the
    # prewarmed store under the exact key, so nothing is read from SPICE, and
    # restore it afterwards.
    et_start, duration, dt = 1_234_567.125, 60.0, 30.0
    primary = PE._spice_query_name("Earth")
    names = SE._collect_nbody_query_names((declared,))
    seed = SE._nbody_ephemeris_cache_from_samples(primary, names,
        [et_start, et_start + 30.0, et_start + 60.0], fill(SVector(1.0, 2.0, 3.0), 3, 2))
    key = SE._nbody_ephemeris_reuse_key(primary, names, et_start, duration, dt)
    previous = lock(SE._EPHEMERIS_REUSE_LOCK) do
        old = get(SE._NBODY_EPHEMERIS_PREWARMED_CACHE, key, nothing)
        SE._NBODY_EPHEMERIS_PREWARMED_CACHE[key] = seed
        old
    end
    try
        withenv("SPACEAGORA_NBODY_EPHEMERIS_CACHE" => "1",
                "SPACEAGORA_NBODY_EPHEMERIS_CACHE_DT_S" => "30.0",
                "SPACEAGORA_NBODY_EPHEMERIS_CACHE_MAX_SAMPLES" => "100") do
            p = SM.ODEParams(n_sats=1, args=simulation_args((declared,); ephemerides=SM.SpiceEphemeridesModel()))
            SE._initialize_nbody_ephemeris_cache!(p, et_start, duration)
            @test p.shared_buffers.nbody_ephemeris_cache[] === seed
        end
    finally
        lock(SE._EPHEMERIS_REUSE_LOCK) do
            if previous === nothing
                delete!(SE._NBODY_EPHEMERIS_PREWARMED_CACHE, key)
            else
                SE._NBODY_EPHEMERIS_PREWARMED_CACHE[key] = previous
            end
        end
    end
end

@testset "Public simulation with a declaring model" begin
    # The analytic backend refuses before integration.
    err = try
        SE.run_simulation(simulation_args((SM.InverseSquaredGravityModel(), DeclaredThirdBodyEffector()));
            return_solution=true)
        nothing
    catch ex
        ex
    end
    @test err isa ArgumentError
    @test occursin("SimpleEphemeridesModel does not support third-body ephemerides", sprint(showerror, err))

    # With SPICE, a model that computes the Moon's tidal term from its declared
    # sample follows NBodyGravityModel's run.
    if !isdir(SPICE_PATH)
        @info "SPICE kernels absent; skipping the declared third-body run" SPICE_PATH
        @test true
    else
        planet = SM.Earth("", SPICE_PATH)
        final(effectors) = collect(SE.run_simulation(
            simulation_args(effectors; ephemerides=SM.SpiceEphemeridesModel(), mission_time=600.0,
                planet=planet);
            return_solution=true).u[end])
        declared = DeclaredMoonGravity()
        keplerian = final((SM.InverseSquaredGravityModel(),))
        via_declaration = final((SM.InverseSquaredGravityModel(), declared))
        via_nbody = final((SM.InverseSquaredGravityModel(), declared.inner))
        moon_effect = via_nbody[1:3] - keplerian[1:3]
        @test norm(moon_effect) > 1e-3
        @test norm(via_declaration[1:3] - via_nbody[1:3]) <= 1e-6 * norm(moon_effect)
    end
end

end
