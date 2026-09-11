# Prototype: process-pool (Distributed) backend for propagating an uncoupled LEO
# constellation with real, uncached ("standard") GRAM, as an alternative to
# leo_2048_constellation_gram_scaling.jl's "ensemble" mode -- which does the same
# uncoupled per-satellite propagation via `run_constellation_ensemble`'s Julia-thread
# worker pool.
#
# Motivation: real (native) GRAM is not thread-safe -- concurrent calls from multiple
# Julia threads in one process corrupt or hang on the vendored library's shared
# Fortran/SPICE global state (see [[project_gram_isolated_pool_unsafe]] in project
# memory: SPACEAGORA_GRAM_ISOLATED_POOL=on, which deep-copies model *objects* within one
# process, still hangs or crashes with a SPICE BADSUBSCRIPT error, because deepcopy
# cannot isolate the native library's own global state). Separate OS *processes* each
# independently load their own copy of the native library with no shared global state --
# this codebase already relies on exactly that fact for Monte Carlo campaigns
# (aerobraking_perturbation_mc/main.jl, performance_mc_thread_scaling_worker.jl, whose
# header states plainly: "no shared Fortran global state between concurrent cases"). This
# script asks whether the same fix applies to a LEO constellation's satellites.
#
# Each satellite here is fully independent (no coupling, no GNC effectors), matching
# run_constellation_ensemble's own precondition -- so it maps cleanly onto one
# `run_simulation` call per worker *process* instead of per worker *thread*: same orbit
# parameters, mission length, and env config as leo_2048_constellation_gram_scaling.jl's
# "ensemble" mode, so the two are a direct, apples-to-apples comparison.
#
# Usage (compare directly against the threaded baseline):
#   julia --project=. --threads=8 benchmarks/studies/leo_2048_constellation_gram_scaling.jl ensemble
#   julia --project=. benchmarks/studies/leo_ensemble_gram_process_pool.jl
#
# Env knobs (defaults chosen to match a specific already-measured comparison point):
#   LEO_SCALING_N_SATS          number of satellites (default 1024)
#   LEO_SCALING_MISSION_TIME_S  mission length in seconds (default 10.0)
#   LEO_POOL_N_WORKERS          number of Distributed worker processes (default 8)

const _LEPP_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using Distributed
using Statistics
using Printf

include(joinpath(_LEPP_REPO_ROOT, "examples", "common.jl"))

const N_SATS = parse(Int, get(ENV, "LEO_SCALING_N_SATS", "1024"))
const ALT_M = 150e3
const MISSION_TIME_S = parse(Float64, get(ENV, "LEO_SCALING_MISSION_TIME_S", "10.0"))
const N_WORKERS = parse(Int, get(ENV, "LEO_POOL_N_WORKERS", "8"))

# Returns an Expr that, when eval'd on a worker that has already loaded common.jl and
# GRAMSuite, defines _lepp_worker_run(sat_idx) -> Nothing (timing is done by the caller,
# wrapping the remotecall_fetch itself in @elapsed). Each worker builds its own
# single-satellite config locally -- same orbit geometry formula as
# leo_2048_constellation_gram_scaling.jl's build_constellation_config -- and propagates
# it with real, uncached GRAM (freeze-per-step, no vacuum look-ahead cache).
function _lepp_worker_fn_expr(n_sats::Int, alt_m::Float64, mission_time_s::Float64)::Expr
    return quote
        function _lepp_worker_run(sat_idx::Int)::Nothing
            planet = Earth("", SPICE_PATH)
            harmonics = GravitationalHarmonicsModel(
                20, 20, joinpath($_LEPP_REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv"), planet
            )
            root = Link{0}(root=true, m=500.0, ref_area=12.0)
            jitter = 50.0 * (sat_idx - 1) / $n_sats
            ic = InitialCondition(
                ra=planet.Rp_e + $alt_m + jitter,
                rp=planet.Rp_e + $alt_m + jitter,
                i=53.0,
                ω=0.0,
                Ω=10.0,
                ν=360.0 * (sat_idx - 1) / $n_sats
            )
            spacecraft = SpacecraftModel(Joint[], [root], root, true, 500.0, 0.0, root.inertia, 0, 0, ic, sat_idx)

            cfg = SimulationConfiguration(
                simulation_settings=SimulationSettings(
                    results=false, verbose=false, generate_plots=false, normalize=false, save_csv=false
                ),
                mission_configuration=MissionConfiguration(
                    mission_type=MissionTime,
                    keplerian=true,
                    number_of_orbits=1,
                    mission_time=$mission_time_s,
                    orientation_sim=false,
                    num_steps_to_save=20
                ),
                environment_model=EnvironmentModel(
                    planet=planet,
                    EI=300.0,
                    density_model=SimulationModel.GRAMAtmosphereModel(planet_name="earth"),
                    thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
                    topography=false,
                    wind=false
                ),
                dynamics_model=DynamicsModel([spacecraft], (harmonics, AerodynamicCoefficientfM())),
                guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
                navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
                control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
                initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
                integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=2.0)
            )

            result = SpaceAGORA.run_simulation(cfg; isolate_state=false, return_solution=true, return_solver_metadata=true)
            string(result.solution.retcode) == "Success" || error("sat_idx=$(sat_idx) did not report Success (retcode=$(result.solution.retcode))")
            return nothing
        end
        nothing
    end
end

function main()
    @printf("[leo-ensemble-pool] n_sat=%d mission_s=%.1f n_workers=%d gram_method=standard (no vacuum cache)\n",
            N_SATS, MISSION_TIME_S, N_WORKERS)
    flush(stdout)

    worker_exeflags = String[
        "--project=$(_LEPP_REPO_ROOT)",
        "--compiled-modules=existing",
        "--threads=1",
        "--gcthreads=1,1",
    ]
    # Matches env_pairs_for(rhs_mode="serial", gram_method="standard") in
    # leo_2048_constellation_gram_scaling.jl, so this is a direct apples-to-apples
    # comparison against that script's threaded "ensemble" mode -- just swapping the
    # worker-pool backend from Julia threads to OS processes.
    worker_env = [
        "SPACEAGORA_RHS_EXECUTION_MODE" => "serial",
        "SPACEAGORA_HARMONICS_BATCH_ENABLED" => "0",
        "SPACEAGORA_EFFECTOR_PARALLEL" => "off",
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => "off",
        "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "1",
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "off",
        "SPACEAGORA_RHS_CALIBRATE" => "off",
        "SPACEAGORA_INNER_THREAD_BUDGET" => "1",
        "SPACEAGORA_VACUUM_GRAM_CACHE" => "0",
    ]

    addprocs(N_WORKERS; exeflags=worker_exeflags, env=worker_env)

    println("[leo-ensemble-pool] initializing $(length(workers())) worker(s) (loading packages + GRAM + defining run fn)...")
    flush(stdout)

    pkg_init_ex = quote
        import Pkg
        Pkg.activate($_LEPP_REPO_ROOT; io=devnull)
        include($(joinpath(_LEPP_REPO_ROOT, "examples", "common.jl")))
        ensure_gramsuite_loaded!()
        nothing
    end
    fn_init_ex = _lepp_worker_fn_expr(N_SATS, ALT_M, MISSION_TIME_S)

    init_tasks = map(workers()) do w
        @async begin
            remotecall_wait(Core.eval, w, Main, pkg_init_ex)
            remotecall_wait(Core.eval, w, Main, fn_init_ex)
        end
    end
    foreach(wait, init_tasks)

    println("[leo-ensemble-pool] warming up $(length(workers())) worker(s) in parallel (satellite 1 on each)...")
    flush(stdout)
    warmup_tasks = [
        @async remotecall_fetch(Core.eval, w, Main, :(_lepp_worker_run(1)))
        for w in workers()
    ]
    foreach(wait, warmup_tasks)

    println("[leo-ensemble-pool] launching $(N_SATS) satellite(s) on $(length(workers())) distributed worker(s)...")
    flush(stdout)

    sat_times_s = Vector{Float64}(undef, N_SATS)
    next_sat = Threads.Atomic{Int}(1)

    wall_s = @elapsed begin
        worker_tasks = Task[]
        for w in workers()
            task = @async begin
                while true
                    idx = Threads.atomic_add!(next_sat, 1)
                    idx > N_SATS && break
                    t_i = @elapsed remotecall_fetch(Core.eval, w, Main, :(_lepp_worker_run($idx)))
                    sat_times_s[idx] = t_i
                end
            end
            push!(worker_tasks, task)
        end
        foreach(wait, worker_tasks)
    end

    @printf(
        "[leo-ensemble-pool] done  wall=%.3f s  mean_sat_s=%.4f  min_sat_s=%.4f  max_sat_s=%.4f\n",
        wall_s, mean(sat_times_s), minimum(sat_times_s), maximum(sat_times_s),
    )
    println("median wall time (ensemble-process-pool, standard GRAM, $(length(workers())) worker processes): $(round(wall_s; digits=3)) s")
end

main()
