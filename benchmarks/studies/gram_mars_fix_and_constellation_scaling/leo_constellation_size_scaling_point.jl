# Shared per-(n_sats, mode) simulation logic for leo_constellation_size_scaling.jl
# (in-process sweep) and leo_constellation_size_scaling_worker.jl (standalone
# single-point CLI). Factored out of the worker's former top-level-ARGS-only
# form so multiple points can run back-to-back in ONE process.
#
# This used to require one subprocess per point specifically because ODEParams
# was parameterized on N_sats as well as the density-model type, so every
# distinct satellite count forced its own fresh JIT specialization of the
# whole RHS/solver pipeline -- accumulating compiled code across a 33-point
# sweep in one process. That's no longer true: N_sats is now a runtime field
# (src/core/types/runtime_types.jl's ODEParams/SharedBuffers), so every n_sats
# value for a given mode shares one compiled specialization, and only the
# density-model type (3 modes) still varies the type. Subprocess-per-point is
# no longer necessary to bound compiled-code growth to a sane amount.

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
include(joinpath(REPO_ROOT, "examples", "common.jl"))
include(joinpath(@__DIR__, "resource_monitor.jl"))
ensure_gramsuite_loaded!()
const GRAMAtmosphereModel = SimulationModel.GRAMAtmosphereModel
const GRAMAtmosphereModelSurrogate = SimulationModel.GRAMAtmosphereModelSurrogate

const ALT_M = 150e3
const MISSION_TIME_S = parse(Float64, get(ENV, "LEO_SCALING_MISSION_TIME_S", "600.0"))
const N_REPEATS = 1

function build_constellation_config(n_sats::Int, mode::String)
    planet = Earth("", SPICE_PATH)
    harmonics = GravitationalHarmonicsModel(
        20, 20, joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv"), planet
    )

    density_model = if mode == "standard"
        GRAMAtmosphereModel(planet_name="earth")
    elseif mode == "surrogate"
        GRAMAtmosphereModelSurrogate(planet_name="earth")
    else
        NoAtmosphereModel()
    end

    # no_gram drops the aerodynamic effector entirely (matching this repo's
    # own no-GRAM baseline convention, e.g. examples/AGORA_Earth_NoGRAM.jl) so
    # the comparison isolates atmosphere-model overhead rather than measuring
    # a drag effector that always computes zero force against a zero-density
    # model.
    effectors = mode == "no_gram" ? (harmonics,) : (harmonics, AerodynamicCoefficientfM())

    spacecraft = SpacecraftModel[]
    for i in 1:n_sats
        root = Link{0}(root=true, m=500.0, ref_area=12.0)
        jitter = 50.0 * (i - 1) / max(n_sats, 1)
        ic = InitialCondition(
            ra=planet.Rp_e + ALT_M + jitter,
            rp=planet.Rp_e + ALT_M + jitter,
            i=53.0,
            ω=0.0,
            Ω=10.0,
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

function env_pairs_for(mode::String, route::String)
    if route == "monolithic"
        common = [
            "SPACEAGORA_RHS_EXECUTION_MODE" => "flat",
            "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "1",
            "SPACEAGORA_EFFECTOR_PARALLEL" => "auto",
            "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
            "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
            "SPACEAGORA_RHS_CALIBRATE" => "off",
            "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
        ]
        if mode == "standard"
            return vcat(common, [
                "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
                "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "1",
                "SPACEAGORA_VACUUM_GRAM_CACHE" => "0",
            ])
        elseif mode == "surrogate"
            # Surrogate has no native lock, but "auto" still defers to the generic
            # density-callback gate (auto_thread_min_budget(:density_callback), default
            # 16 threads) regardless of density-model type -- override via this env var
            # to test whether forcing it "on" below that threshold helps a lock-free model.
            density_callback_mode = get(ENV, "SPACEAGORA_TS_SURROGATE_DENSITY_CALLBACK_OVERRIDE", "auto")
            return vcat(common, ["SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => density_callback_mode])
        else
            return common
        end
    end
    # route == "ensemble": each member is a single satellite, so there is no
    # multi-satellite batching to gain from rhs_mode="flat" -- "serial" is the
    # natural per-member routing (matches leo_2048_constellation_gram_scaling.jl's
    # own ensemble mode). Parallelism comes entirely from run_constellation_ensemble's
    # outer worker-task split (SPACEAGORA_OUTER_PARALLEL_ACTIVE is set to "1" internally
    # by that runner while it dispatches members, overriding the "0" below for the
    # duration of each member's solve).
    common = [
        "SPACEAGORA_RHS_EXECUTION_MODE" => "serial",
        "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "0",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_RHS_CALIBRATE" => "off",
    ]
    if mode == "standard"
        # Same rationale as leo_2048_constellation_gram_scaling.jl: without the
        # vacuum-predicted cache, concurrent native GRAM calls across members are
        # unsafe, and the cache itself hangs on rebuild with 2+ satellites -- freeze
        # per-step and disable the cache so only the proven-safe initial build runs.
        return vcat(common, [
            "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "1",
            "SPACEAGORA_VACUUM_GRAM_CACHE" => "0",
            "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => string(MISSION_TIME_S + 500.0),
            "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1e8",
        ])
    end
    return common
end

function run_once_monolithic(args)
    result = SpaceAGORA.run_simulation(args; isolate_state=false, return_solution=true, return_solver_metadata=true)
    return result.solution
end

function run_once_ensemble(args)
    return SpaceAGORA.run_constellation_ensemble(
        args; threads=Threads.nthreads(), return_solution=true, return_solver_metadata=true
    )
end

function report_outcome(label::String, result, route::String, n_sats::Int)
    if route == "monolithic"
        println("$(label) retcode=$(result.retcode)")
        flush(stdout)
        return string(result.retcode) == "Success"
    end
    n_ok = length(result.successful)
    println("$(label): $(n_ok)/$(n_sats) members succeeded")
    flush(stdout)
    return n_ok == n_sats
end

# Runs one (n_sats, mode) point (warmup + N_REPEATS timed repeats) and returns
# its median wall time plus a ResourceUsage (peak RSS, mean/peak CPU%) sampled
# from this process (resource_monitor.jl) across the timed repeats only --
# warmup is excluded since it also pays one-time JIT/specialization cost that
# isn't representative of steady-state execution at this point. A
# physics/retcode failure only warns (matching the original worker's
# behavior: the timing is still meaningful and reported); only a thrown
# exception propagates to the caller.
function run_scaling_point(n_sats::Int, mode::String, route::String="monolithic")
    args = build_constellation_config(n_sats, mode)
    pairs = env_pairs_for(mode, route)
    run_once = route == "monolithic" ? run_once_monolithic : run_once_ensemble

    warmup_result = withenv(pairs...) do
        run_once(args)
    end
    report_outcome("warmup", warmup_result, route, n_sats) ||
        @warn "warmup run did not report full success -- results below may not reflect a full propagation." n_sats mode route

    GC.gc()
    times = Float64[]
    local last_result
    _, usage = measure_resource_usage() do
        for r in 1:N_REPEATS
            GC.gc()
            t = @elapsed begin
                last_result = withenv(pairs...) do
                    run_once(args)
                end
            end
            push!(times, t)
            report_outcome("repeat $(r) ($(round(t; digits=4)) s)", last_result, route, n_sats)
        end
    end
    median_t = sort(times)[cld(length(times), 2)]
    println("median wall time (mode=$(mode), route=$(route), n_sat=$(n_sats), $(Threads.nthreads()) threads): $(round(median_t; digits=4)) s")
    println("resource usage (mode=$(mode), n_sat=$(n_sats), $(Threads.nthreads()) threads): peak_rss=$(round(usage.peak_rss_mb; digits=1)) MB, mean_cpu=$(round(usage.mean_cpu_pct; digits=1))%, peak_cpu=$(round(usage.peak_cpu_pct; digits=1))% ($(usage.n_samples) samples)")
    flush(stdout)
    return median_t, usage
end
