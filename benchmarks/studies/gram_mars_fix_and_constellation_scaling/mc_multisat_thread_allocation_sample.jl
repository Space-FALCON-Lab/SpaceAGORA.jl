# Shared per-sample logic for mc_multisat_thread_allocation_worker.jl. Included on
# the entry-point process (both outer_backend="threads" and "process" routes) and,
# for outer_backend="process", via `remotecall_wait` on every SpaceAGORA.ParallelProcess
# pool worker too -- mirroring mc_multisat_process_backend.jl / mc_multisat_process_sample.jl's
# existing split.
#
# Deliberately ARGS-independent (no top-level ARGS parsing): a pool worker's own ARGS
# is empty (Distributed workers aren't launched with this study's command-line
# arguments), so `build_sample_config`/`env_pairs` take `mode`/`outer_backend`/etc. as
# explicit parameters rather than reading process-global constants derived from ARGS.
# The entry-point script's run_monte_carlo closure then captures those values BY VALUE
# (ordinary Julia closure semantics), which Distributed ships alongside the closure --
# no separate step is needed to make them "visible" on the worker.

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO_ROOT, "examples", "common.jl"))
ensure_gramsuite_loaded!()
const GRAMAtmosphereModel = SimulationModel.GRAMAtmosphereModel
const GRAMAtmosphereModelSurrogate = SimulationModel.GRAMAtmosphereModelSurrogate

const N_SATS_PER_SAMPLE = 4
const N_SAMPLES = 8
const ALT_M = 150e3
const MISSION_TIME_S = 600.0

function build_sample_config(mode::String, seed::Int)
    planet = Earth("", SPICE_PATH)
    harmonics = GravitationalHarmonicsModel(
        20, 20, joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv"), planet
    )
    density_model = if mode == "surrogate"
        GRAMAtmosphereModelSurrogate(planet_name="earth")
    elseif mode == "standard"
        GRAMAtmosphereModel(planet_name="earth")
    else
        NoAtmosphereModel()
    end
    effectors = mode == "no_gram" ? (harmonics,) : (harmonics, AerodynamicCoefficientfM())

    # Deterministic per-seed jitter (distinct from the fixed inter-satellite phasing
    # below) so each of the 8 Monte Carlo samples is a different constellation, not
    # 8 repeats of the same one.
    alt_jitter_m = 200.0 * seed
    raan_jitter_deg = 0.75 * seed

    spacecraft = SpacecraftModel[]
    for i in 1:N_SATS_PER_SAMPLE
        root = Link(root=true, m=500.0, ref_area=12.0)
        phase = 50.0 * (i - 1) / max(N_SATS_PER_SAMPLE, 1)
        ic = InitialCondition(
            ra=planet.Rp_e + ALT_M + phase + alt_jitter_m,
            rp=planet.Rp_e + ALT_M + phase + alt_jitter_m,
            i=53.0,
            ω=0.0,
            Ω=10.0 + raan_jitter_deg,
            ν=360.0 * (i - 1) / max(N_SATS_PER_SAMPLE, 1)
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

function env_pairs(mode::String, outer_backend::String, outer_workers::Int, inner_threads::Int)
    common = [
        "SPACEAGORA_RHS_EXECUTION_MODE" => "flat",
        "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "1",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
        # Only meaningful for outer_backend="threads": tells inner thread policies
        # they're already running under an outer split so they yield/share budget.
        # The "process" backend's adaptive threads=:auto dispatch checks this SAME
        # flag itself (ParallelPolicy.outer_parallel_active()) to decide whether it's
        # nested inside an *enclosing* outer split and should route serially -- pre-
        # setting it to "1" here would make every "process" point silently collapse
        # to serial before ever reaching the process pool.
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => (outer_backend == "threads" && outer_workers > 1 ? "1" : "0"),
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_RHS_CALIBRATE" => "off",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(inner_threads),
    ]
    if mode == "standard"
        # Same rationale as leo_constellation_size_scaling_point.jl: without the
        # vacuum-predicted cache, concurrent native GRAM calls are unsafe under
        # outer_backend="threads" (now locked, but still serialized -- see the
        # worker script's header), and the cache itself hangs on rebuild with 2+
        # satellites -- freeze per-step and disable the cache so only the
        # proven-safe initial build runs.
        common = vcat(common, [
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "1",
            "SPACEAGORA_VACUUM_GRAM_CACHE" => "0",
            "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => string(MISSION_TIME_S + 500.0),
            "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1e8",
        ])
    end
    mode == "surrogate" || return common
    # See leo_constellation_size_scaling_worker.jl: the generic density-callback
    # auto-gate defaults to a 16-thread minimum regardless of density-model type, so
    # "auto" would silently leave density sampling serial for every inner_threads
    # value tested here (max 8). Force "on" so the lock-free surrogate's inner
    # parallelism is actually exercised.
    return vcat(common, ["SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "on"])
end
