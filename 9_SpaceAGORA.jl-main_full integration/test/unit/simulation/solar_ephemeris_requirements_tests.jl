module SolarEphemerisRequirementsTests

using Test, SpaceAGORA, StaticArrays
const SM = SpaceAGORA.SimulationModel
const SE = SpaceAGORA.SimulationEngine

struct DeclaredSolarEffector <: SM.AbstractForceTorqueModel end
SM.environment_requirements(::DeclaredSolarEffector) = SM.EffectorEnvironmentRequirements(solar=true)
struct LegacyEffector <: SM.AbstractForceTorqueModel end

validation_args(effectors, backend=SM.SimpleEphemeridesModel()) = (
    dynamics_model=(dynamic_effectors=effectors,),
    environment_model=(planet=SM.Earth(), ephemerides_model=backend),
    simulation_settings=(verbose=false,),
)

const ACTIVE = (
    SM.FacetSolarRadiationPressureModel(),
    DeclaredSolarEffector(),
    SM.SolarRadiationPressureModel(1.2, 12.0),
    SM.SolarRadiationPressureModel(1.2, 12.0; direct=false, albedo=true),
)
const INACTIVE = (
    SM.SolarRadiationPressureModel(1.2, 0.0),
    SM.SolarRadiationPressureModel(1.2, 12.0; direct=false, albedo=false, ir=true),
    SM.SolarRadiationPressureModel(1.2, 12.0; direct=false, albedo=false, ir=false),
    LegacyEffector(),
    SM.InverseSquaredGravityModel(),
)

@testset "Solar requirements determine ephemeris compatibility" begin
    for effector in ACTIVE
        args = validation_args((effector,))
        @test_throws ArgumentError SE._validate_ephemerides_support!(args)
        @test SE._validate_ephemerides_support!(
            validation_args((effector,), SM.SpiceEphemeridesModel())) === nothing
    end
    for effector in INACTIVE
        @test SE._validate_ephemerides_support!(validation_args((effector,))) === nothing
    end
    @test SE._validate_ephemerides_support!(validation_args(())) === nothing
    # An inactive first entry must not hide a later solar consumer.
    @test_throws ArgumentError SE._validate_ephemerides_support!(
        validation_args((INACTIVE[1], ACTIVE[1])))
end

function simulation_args(effector)
    planet = SM.Earth()
    root = SM.Link(root=true, m=500.0, ref_area=2.0)
    SM.add_facet!(root, SM.Facet(area=2.0, normal_vector=[1.0, 0.0, 0.0]))
    ic = SM.InitialCondition(ra=planet.Rp_e+500e3, rp=planet.Rp_e+500e3,
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
            ephemerides_model=SM.SimpleEphemeridesModel(), topography=false, wind=false),
        dynamics_model=SM.DynamicsModel([sc], (SM.InverseSquaredGravityModel(), effector)),
        guidance_model=SM.GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=SM.NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=SM.ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=SM.InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=SM.IntegrationTolerances(),
        solver_config=SE.SolverConfig(solver_mode=:tsit5),
    )
end

@testset "Public simulation rejects facet SRP before integration" begin
    for effector in (ACTIVE[1], ACTIVE[3])
        err = try
            SE.run_simulation(simulation_args(effector); return_solution=true)
            nothing
        catch ex
            ex
        end
        @test err isa ArgumentError
        @test occursin("SimpleEphemeridesModel does not support solar-radiation-pressure ephemerides",
            sprint(showerror, err))
    end
end

@testset "Declared solar consumers initialize and sample the Sun cache" begin
    et_start, duration, dt = 1_234_567.125, 30.0, 30.0
    sun0 = SVector(149_597_870_700.0, 1000.0, 0.0)
    sun1 = sun0 + SVector(30.0, 60.0, 0.0)
    seed = SM.SRPSunEphemerisCache([et_start, et_start+duration], [sun0, sun1])
    key = SE._srp_ephemeris_reuse_key("earth", et_start, duration, dt)
    # Seed the real reuse store to exercise initialization without loading
    # kernels or replacing the production sampler. Restore only this key.
    previous = lock(SE._EPHEMERIS_REUSE_LOCK) do
        old = get(SE._SRP_EPHEMERIS_REUSE_CACHE, key, nothing)
        SE._SRP_EPHEMERIS_REUSE_CACHE[key] = seed
        old
    end
    try
        withenv("SPACEAGORA_SRP_EPHEMERIS_CACHE"=>"1",
                "SPACEAGORA_EPHEMERIS_CACHE_REUSE"=>"1",
                "SPACEAGORA_SRP_EPHEMERIS_CACHE_DT_S"=>"30.0",
                "SPACEAGORA_SRP_EPHEMERIS_CACHE_MAX_SAMPLES"=>"100") do
            for effector in (ACTIVE..., INACTIVE...)
                counters = (srp_spkpos_runtime_calls=Threads.Atomic{Int64}(0),
                            srp_spkpos_cache_build_calls=Threads.Atomic{Int64}(0))
                p = (args=validation_args((effector,), SM.SpiceEphemeridesModel()),
                     shared_buffers=(srp_sun_ephemeris_cache=Ref{Any}(nothing),
                        et_start=Ref(et_start), spice_rhs_memo_enabled=Ref(false),
                        spice_rhs_memo=SM.SpiceRhsMemo(), spice_runtime_counters=counters))
                SE._initialize_srp_sun_ephemeris_cache!(p, et_start, duration)
                if effector in ACTIVE
                    @test p.shared_buffers.srp_sun_ephemeris_cache[] === seed
                    # Avoid native fallback on a regression; the missing cache
                    # assertion above remains the useful failure.
                    if p.shared_buffers.srp_sun_ephemeris_cache[] === seed
                        sample = SE.sample_solar_ephemeris(nothing, p, 1, 15.0)
                        @test sample.sun_pos_ii ≈ (sun0+sun1)/2
                    end
                else
                    @test p.shared_buffers.srp_sun_ephemeris_cache[] === nothing
                end
                @test counters.srp_spkpos_runtime_calls[] == 0
                @test counters.srp_spkpos_cache_build_calls[] == 0
            end
        end
    finally
        lock(SE._EPHEMERIS_REUSE_LOCK) do
            if previous === nothing
                delete!(SE._SRP_EPHEMERIS_REUSE_CACHE, key)
            else
                SE._SRP_EPHEMERIS_REUSE_CACHE[key] = previous
            end
        end
    end
end

end
