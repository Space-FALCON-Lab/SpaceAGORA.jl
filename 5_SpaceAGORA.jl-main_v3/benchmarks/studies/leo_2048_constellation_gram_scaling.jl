# Quick serial-vs-threaded scaling check for a large LEO constellation with a
# real (native) GRAM atmosphere and harmonics gravity, both applied per satellite.
#
# Julia's thread count is fixed at process startup, so "serial" and "8 threads"
# must be two separate invocations -- run both and compare the printed median:
#
#   julia --project=. --threads=1 benchmarks/studies/leo_2048_constellation_gram_scaling.jl serial
#   julia --project=. --threads=8 benchmarks/studies/leo_2048_constellation_gram_scaling.jl parallel
#
# "serial"   forces SPACEAGORA_RHS_EXECUTION_MODE=serial (no Polyester batching,
#            no harmonics SIMD batch kernel) AND standard/direct native GRAM
#            calls (vacuum-predicted trajectory cache disabled): every
#            satellite's atmosphere sample is a fresh per-step native GRAM
#            query, the pre-look-ahead-cache baseline.
# "parallel" forces SPACEAGORA_RHS_EXECUTION_MODE=flat (batched harmonics SIMD
#            kernel, flat constellation effector queue, all using the
#            available threads) AND the GRAM trajectory look-ahead cache
#            (drag-free propagation + cubic-spline interpolation).
#
# A third mode compares the two GRAM-calling methods for trajectory accuracy,
# independent of threading/runtime:
#
#   julia --project=. --threads=8 benchmarks/studies/leo_2048_constellation_gram_scaling.jl accuracy
#
# "accuracy" runs the same constellation twice in one process -- once with
# standard/direct GRAM (the reference), once with the trajectory look-ahead
# cache -- holding execution routing fixed between the two runs so any
# resulting state difference is attributable only to the GRAM-calling method,
# then reports position/velocity/attitude parity via
# parallelization_performance/trajectory_parity.jl's ppc_compare_trajectories.
#
# A fourth mode propagates the (uncoupled) constellation via
# run_constellation_ensemble instead of one coupled state vector -- each
# satellite becomes an independent single-satellite run_simulation call,
# dispatched to a worker once for its entire propagation instead of paying
# per-timestep thread dispatch across satellites, and each keeps its own
# adaptive step size instead of the whole constellation sharing the global
# minimum:
#
#   julia --project=. --threads=8 benchmarks/studies/leo_2048_constellation_gram_scaling.jl ensemble
#
# KNOWN BUG (worked around below, not fixed) in the vendored GRAMSuite native
# atmosphere model -- see [[project_gram_multisat_rebuild_bug]] in memory: it
# hangs when 2+ satellites need a vacuum-predicted-cache *rebuild* (mission
# longer than the cache horizon). A single satellite handles rebuilds fine over
# a full hour; any 2+ satellites needing one (tested at 2/4/8) hang
# indefinitely. Reproduces under fully serial, single-threaded execution, so
# it's sequential-call state corruption in the native library, not a Julia
# `Threads` race. Workaround: cache horizon + deviation threshold both set well
# beyond what a full mission can reach, so only the initial build (proven safe
# for any satellite count, and safe under real concurrent threaded access too)
# ever runs.
#
# A second issue -- a genuine data race in this codebase's own aero
# force-caching (drag_cache/lift_cache/cross_cache in SaveCache, resized
# lazily and unsafely from concurrent worker threads) -- was found and fixed
# in setup.jl/execution.jl (_initialize_save_cache_buffers!). It was initially
# mistaken for a GRAM thread-safety problem because it only surfaced once
# gravity was threaded across 2+ satellites (ConcurrencyViolationError at
# small N, a GC segfault from corrupted memory at larger N); it had nothing to
# do with GRAM or density-callback threading, both of which are safe once this
# fix is in place.

mode = length(ARGS) >= 1 ? ARGS[1] : ""
mode in ("serial", "parallel", "accuracy", "ensemble") || error(
    "Usage: julia --project=. --threads=<N> $(basename(@__FILE__)) <serial|parallel|accuracy|ensemble>\n" *
    "  serial   : standard/direct GRAM, run with --threads=1\n" *
    "  parallel : GRAM trajectory look-ahead cache, run with --threads=8\n" *
    "  accuracy : compare standard vs. look-ahead trajectories in one process\n" *
    "  ensemble : run_constellation_ensemble (one worker per satellite), run with --threads=8"
)

include(joinpath(@__DIR__, "..", "..", "examples", "common.jl"))
ensure_gramsuite_loaded!()
const GRAMAtmosphereModel = SimulationModel.GRAMAtmosphereModel

using LinearAlgebra
using StaticArrays
using Statistics

include(joinpath(@__DIR__, "parallelization_performance", "trajectory_parity.jl"))

# Overridable via LEO_SCALING_N_SATS so leo_constellation_size_scaling.jl can sweep
# constellation size across subprocess invocations of this same script (guaranteeing
# identical orbit parameters/config at every size).
const N_SATS = parse(Int, get(ENV, "LEO_SCALING_N_SATS", "1024"))
const ALT_M = 150e3
# ~1/9 of one orbital period at 150 km altitude (~10 min) by default; also overridable
# so ad-hoc process-pool-vs-thread comparisons can hold this fixed at a short duration.
const MISSION_TIME_S = parse(Float64, get(ENV, "LEO_SCALING_MISSION_TIME_S", "600.0"))
const N_REPEATS = 1

function build_constellation_config()
    planet = Earth("", SPICE_PATH)
    harmonics = GravitationalHarmonicsModel(
        20, 20, joinpath(REPO_ROOT, "data", "Gravity_harmonics_data", "EarthGGM05C.csv"), planet
    )

    # EI is set well above the orbital altitude (see EnvironmentModel below), and
    # in_atmosphere[] is now initialized from each satellite's actual starting
    # altitude (see _initialize_in_atmosphere_flags! in setup.jl), so a plain
    # circular orbit starts -- and stays -- correctly flagged as in-atmosphere
    # for the whole mission without needing any eccentricity to force a crossing.
    spacecraft = SpacecraftModel[]
    for i in 1:N_SATS
        root = Link{0}(root=true, m=500.0, ref_area=12.0)
        jitter = 50.0 * (i - 1) / N_SATS  # sub-km jitter: avoids exactly duplicate states
        ic = InitialCondition(
            ra=planet.Rp_e + ALT_M + jitter,
            rp=planet.Rp_e + ALT_M + jitter,
            i=53.0,
            ω=0.0,
            Ω=10.0,
            ν=360.0 * (i - 1) / N_SATS
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
            density_model=GRAMAtmosphereModel(planet_name="earth"),
            thermal_model=MaxwellianHeat(thermal_accomodation_factor=1.0, planet=planet),
            topography=false,
            wind=false
        ),
        dynamics_model=DynamicsModel(spacecraft, (harmonics, AerodynamicCoefficientfM())),
        guidance_model=GuidanceModel(guidance_effectors=(), guidance_rates=Float64[]),
        navigation_model=NavigationModel(navigation_effectors=(), navigation_rates=Float64[]),
        control_model=ControlModel(control_effectors=(), control_rates=Float64[]),
        initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
        integration_tolerances=IntegrationTolerances(reltol_orbit=1e-9, abstol_orbit=1e-9, dt_max_orbit=2.0)
    )
end

# gram_method="standard"  : density is sampled once per accepted solver step
#                           (SPACEAGORA_DENSITY_FREEZE_PER_STEP=1) via real,
#                           uncached native GRAM queries (SPACEAGORA_VACUUM_GRAM_CACHE=0)
#                           -- the pre-look-ahead-cache baseline, and the
#                           reference trajectory for the "accuracy" comparison.
#                           Freeze-per-step is required, not optional: without
#                           it, real GRAM's per-call perturbation-model noise is
#                           resampled at every RK stage, which an adaptive
#                           solver's step-size controller reacts to as if it
#                           were stiffness, collapsing dt (measured: a 2-sat,
#                           1s-mission run without freeze took 604-613s and
#                           2.42M solver steps; with freeze, 35.7s and 37
#                           steps -- matching the look-ahead cache's own
#                           timing). Freezing is accurate here because a LEO
#                           satellite's altitude -- the dominant driver of the
#                           smooth mean density -- changes negligibly over one
#                           integration step; it just also removes GRAM's
#                           small-scale perturbation noise from within a step,
#                           which is what actually fixes the solver.
# gram_method="lookahead" : the vacuum-predicted trajectory cache (drag-free
#                           propagation + cubic-spline interpolation).
#
# Horizon is set above the mission length, AND the deviation threshold is set
# very large, as a workaround for the known rebuild bug documented at the top
# of this file: a rebuild triggers either when t exceeds the cache horizon OR
# when the real (drag-perturbed) trajectory deviates from the drag-free vacuum
# prediction by more than the deviation threshold -- over a full hour with real
# drag, that deviation is reached well before 3600s regardless of horizon, so
# both knobs need loosening to guarantee only the (proven-safe) initial build
# ever runs and no rebuild is ever triggered.
function env_pairs_for(; rhs_mode::String, gram_method::String)
    gram_method in ("standard", "lookahead") || error("gram_method must be \"standard\" or \"lookahead\"")
    vacuum_enabled = gram_method == "lookahead" ? "1" : "0"
    freeze_per_step = gram_method == "standard" ? "1" : "0"
    # Native GRAM is documented (see header) as unsafe under concurrent/interleaved
    # per-satellite queries -- that hazard is exactly what the vacuum-predicted
    # cache exists to avoid. Without the cache (gram_method="standard"), density
    # sampling must be forced serial regardless of rhs_mode/threading, or
    # concurrent direct native calls across satellites hang (reproduced: a
    # flat/threaded standard-GRAM run hung indefinitely at 2048 sats, killed by
    # a 1800s external timeout with zero further output after the initial GRAM
    # data load). Gravity/harmonics can still batch/thread fine either way.
    density_parallel = gram_method == "standard" ? "off" : (rhs_mode == "serial" ? "off" : "auto")
    return [
        "SPACEAGORA_RHS_EXECUTION_MODE" => rhs_mode,
        "SPACEAGORA_HARMONICS_BATCH_ENABLED" => (rhs_mode == "serial" ? "0" : "1"),
        "SPACEAGORA_EFFECTOR_PARALLEL" => (rhs_mode == "serial" ? "off" : "auto"),
        "SPACEAGORA_DENSITY_CALLBACK_PARALLEL" => density_parallel,
        "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => freeze_per_step,
        "SPACEAGORA_OUTER_PARALLEL_ACTIVE" => "0",
        "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE" => "0",
        "SPACEAGORA_RHS_CALIBRATE" => "off",
        "SPACEAGORA_INNER_THREAD_BUDGET" => string(Threads.nthreads()),
        "SPACEAGORA_VACUUM_GRAM_CACHE" => vacuum_enabled,
        "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "20",
        "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => string(MISSION_TIME_S + 500.0),
        "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1e8",
    ]
end

function run_once(args)
    result = SpaceAGORA.run_simulation(args; isolate_state=false, return_solution=true, return_solver_metadata=true)
    return result.solution
end

args = build_constellation_config()

function timed_run(label::String, pairs)
    sol_warmup = withenv(pairs...) do
        run_once(args)
    end
    println("$(label) warmup retcode=$(sol_warmup.retcode)")
    flush(stdout)
    string(sol_warmup.retcode) == "Success" ||
        @warn "$(label) warmup run did not report Success -- results below may not reflect a full propagation."

    GC.gc()
    times = Float64[]
    local last_sol
    for r in 1:N_REPEATS
        GC.gc()
        t = @elapsed begin
            last_sol = withenv(pairs...) do
                run_once(args)
            end
        end
        push!(times, t)
        println("$(label) repeat $(r): $(round(t; digits=3)) s  retcode=$(last_sol.retcode)")
        flush(stdout)
    end
    median_t = sort(times)[cld(length(times), 2)]
    return (solution=last_sol, median_s=median_t)
end

println("mode=$(mode)  julia_threads=$(Threads.nthreads())  n_sat=$(N_SATS)  alt_km=$(ALT_M/1e3)  mission_s=$(MISSION_TIME_S)")
flush(stdout)

function timed_ensemble_run(label::String, pairs)
    ensemble_once = () -> SpaceAGORA.run_constellation_ensemble(
        args; threads=Threads.nthreads(), return_solution=true, return_solver_metadata=true
    )

    warmup = withenv(pairs...) do
        ensemble_once()
    end
    n_ok = length(warmup.successful)
    println("$(label) warmup: $(n_ok)/$(N_SATS) members succeeded")
    flush(stdout)
    n_ok == N_SATS || @warn "$(label) warmup: $(length(warmup.failed)) member(s) failed -- results below may not reflect a full propagation."

    GC.gc()
    times = Float64[]
    local last_result
    for r in 1:N_REPEATS
        GC.gc()
        t = @elapsed begin
            last_result = withenv(pairs...) do
                ensemble_once()
            end
        end
        push!(times, t)
        println("$(label) repeat $(r): $(round(t; digits=3)) s  $(length(last_result.successful))/$(N_SATS) succeeded")
        flush(stdout)
    end
    median_t = sort(times)[cld(length(times), 2)]
    return (result=last_result, median_s=median_t)
end

if mode in ("serial", "parallel")
    rhs_mode = mode == "serial" ? "serial" : "flat"
    gram_method = mode == "serial" ? "standard" : "lookahead"
    pairs = env_pairs_for(rhs_mode=rhs_mode, gram_method=gram_method)
    println("gram_method=$(gram_method)")
    flush(stdout)
    result = timed_run(mode, pairs)
    println("median wall time ($(mode), $(gram_method) GRAM, $(Threads.nthreads()) threads): $(round(result.median_s; digits=3)) s")
elseif mode == "ensemble"
    # Each ensemble member is a single satellite, so there is no multi-satellite
    # batching to gain from rhs_mode="flat" -- "serial" is the natural per-member
    # routing. Freeze-per-step still applies (real, uncached GRAM per member).
    pairs = env_pairs_for(rhs_mode="serial", gram_method="standard")
    result = timed_ensemble_run("ensemble", pairs)
    println("median wall time (ensemble, standard GRAM, $(Threads.nthreads()) outer workers): $(round(result.median_s; digits=3)) s")
else # mode == "accuracy"
    # Routing is held fixed (flat/threaded) between the two runs so any resulting
    # state difference is attributable only to the GRAM-calling method, not to a
    # different execution path.
    pairs_standard = env_pairs_for(rhs_mode="flat", gram_method="standard")
    pairs_lookahead = env_pairs_for(rhs_mode="flat", gram_method="lookahead")

    result_ref = timed_run("standard (reference)", pairs_standard)
    result_cmp = timed_run("lookahead", pairs_lookahead)

    println()
    println("median wall time (standard):  $(round(result_ref.median_s; digits=3)) s")
    println("median wall time (lookahead): $(round(result_cmp.median_s; digits=3)) s")

    comparison = ppc_compare_trajectories(result_ref.solution, result_cmp.solution, args; sample_count=31)
    println()
    println("Trajectory parity (standard GRAM = reference, look-ahead = candidate), $(N_SATS) satellites, $(comparison.samples) time samples:")
    println("  position error:  rms=$(comparison.pos_rel_rms)  p90=$(comparison.pos_rel_p90)  max=$(comparison.pos_rel_max)  (relative)")
    println("  velocity error:  rms=$(comparison.vel_rel_rms)  p90=$(comparison.vel_rel_p90)  max=$(comparison.vel_rel_max)  (relative)")
    println("  quaternion angle max: $(comparison.q_angle_max_rad) rad")
    println("  omega error max: $(comparison.omega_rel_max) (relative)")
    println("  mass error max:  $(comparison.mass_rel_max) (relative)")
    println("  periapsis crossings: ref=$(comparison.ref_periapsis_count) cmp=$(comparison.cmp_periapsis_count)")
    println("  EI crossings:        ref=$(comparison.ref_interface_count) cmp=$(comparison.cmp_interface_count)")
    println("  event time delta max: $(comparison.event_time_abs_max_s) s")
    println("  PASS (within parallelization_performance's parity thresholds): $(comparison.pass)")
end
