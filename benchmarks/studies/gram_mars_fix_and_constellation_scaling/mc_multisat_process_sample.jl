# Shared per-sample logic for mc_multisat_process_backend.jl. Included on the
# driver process (serial route) and, via `remotecall_wait`, on every Distributed
# worker (outer_process route). Each process independently constructs its own
# config/model -- safe even for the GRAM surrogate because process isolation
# means no shared native GRAM/CSPICE global state to race on, unlike
# Threads.@spawn-based outer parallelism (see
# project_gram_concurrent_construction_crash memory / Finding 4 in
# THREAD_ALLOCATION_AND_GRAM_CONCURRENCY_HANDOFF.md).

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO_ROOT, "examples", "common.jl"))
ensure_gramsuite_loaded!()
const GRAMAtmosphereModelSurrogate = SimulationModel.GRAMAtmosphereModelSurrogate

const ALT_M = 150e3
const MISSION_TIME_S = 600.0

function mc_env_pairs(mode::String)
    common = [
        "SPACEAGORA_RHS_EXECUTION_MODE" => "flat",
        "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "1",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_RHS_CALIBRATE" => "off",
        "SPACEAGORA_INNER_THREAD_BUDGET" => "1",
    ]
    mode == "surrogate" || return common
    return vcat(common, ["SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on"])
end

function mc_build_sample_config(mode::String, n_sats::Int, seed::Int)
    planet = Earth("", SPICE_PATH)
    harmonics = GravitationalHarmonicsModel(
        20, 20, joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv"), planet
    )
    density_model = mode == "surrogate" ? GRAMAtmosphereModelSurrogate(planet_name="earth") : NoAtmosphereModel()
    effectors = mode == "no_gram" ? (harmonics,) : (harmonics, AerodynamicCoefficientfM())

    # Deterministic per-seed jitter, same scheme as
    # mc_multisat_thread_allocation_worker.jl, so each sample is a distinct
    # constellation. alt_jitter is cycled through the same 1-8 magnitude range
    # that study empirically validated as non-pathological, rather than
    # growing unbounded with `seed`, since N_SAMPLES here is much larger and an
    # unbounded altitude offset risks pushing some sample into a slow/stiff
    # integration regime (see mc_multisat_process_backend.jl header).
    alt_jitter_m = 200.0 * (mod(seed - 1, 8) + 1)
    raan_jitter_deg = 0.75 * seed

    spacecraft = SpacecraftModel[]
    for i in 1:n_sats
        root = Link(root=true, m=500.0, ref_area=12.0)
        phase = 50.0 * (i - 1) / max(n_sats, 1)
        ic = InitialCondition(
            ra=planet.Rp_e + ALT_M + phase + alt_jitter_m,
            rp=planet.Rp_e + ALT_M + phase + alt_jitter_m,
            i=53.0,
            ω=0.0,
            Ω=10.0 + raan_jitter_deg,
            ν=360.0 * (i - 1) / max(n_sats, 1)
        )
        push!(spacecraft, SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, i))
    end

    return SimulationConfiguration(
        simulation_settings=SimulationSettings(
            results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
        ),
        mission_configuration=MissionConfiguration(
            mission_type=MissionTime,
            keplerian=true,
            number_of_orbits=1,
            mission_time=MISSION_TIME_S,
            orientation_sim=false,
            num_steps_to_save=20
        ),
        environment_model=EnvironmentModel(
            planet=planet,
            EI=300.0,
            density_model=density_model,
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel(spacecraft, effectors),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=2.0)
    )
end

function mc_run_sample(mode::String, n_sats::Int, index::Int, seed::Int)
    start_ns = time_ns()
    success = true
    err_str = ""
    try
        withenv(mc_env_pairs(mode)...) do
            args = mc_build_sample_config(mode, n_sats, seed)
            result = SpaceAGORA.run_simulation(
                args; isolate_state=false, return_solution=true, return_solver_metadata=true
            )
            success = string(result.solution.retcode) == "Success"
        end
    catch e
        success = false
        err_str = sprint(showerror, e)
    end
    elapsed_s = (time_ns() - start_ns) / 1.0e9
    return (index=index, seed=seed, success=success, elapsed_s=elapsed_s, error=err_str)
end
