# Worker subprocess for performance_mc_thread_scaling.jl.
#
# Invoked by the orchestrator. Starts SPACEAGORA_MC_DIST_NWORKERS distributed worker
# processes via addprocs — mirroring aerobraking_perturbation_mc/main.jl — then dispatches
# one simulation per worker via remotecall_fetch using the same work-stealing @async loop
# used by _run_samples_distributed in study.jl. Wall-clock time and per-case elapsed times
# are recorded and written as a one-row CSV to SPACEAGORA_MC_THREAD_OUTDIR.
#
# Each distributed worker builds its own simulation config (SPICE calls, struct construction)
# locally in its own process, so there is no shared Fortran global state between concurrent
# cases. This differs from the old Threads.@spawn approach where SPICE calls had to be
# serialised to the main thread.
#
# Environment variables (set by orchestrator):
#   SPACEAGORA_MC_THREAD_OUTDIR      output directory for the result CSV
#   SPACEAGORA_MC_THREAD_NORBITS     orbital periods per case (default: 30)
#   SPACEAGORA_MC_DIST_NWORKERS      number of distributed worker processes to spawn
#   SPACEAGORA_MC_THREAD_SYSIMAGE    path to sysimage .so (forwarded to workers)

const _MCTSW_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using CSV
using DataFrames
using Distributed
using Printf
using Statistics

include(joinpath(_MCTSW_REPO_ROOT, "examples", "common.jl"))

const _MCTSW_RP_ALT_M = 125e3
const _MCTSW_RA_ALT_M = 2_000e3
const _MCTSW_CASE_ID  = "mars_twobody_rp125km_ra2000km"

# Returns an Expr that, when eval'd on a worker that has already loaded common.jl,
# defines _mctsw_worker_run(n_orbits) -> elapsed_seconds. Each worker builds its own
# config locally (no serialisation of SPICE-dependent structs across processes).
function _mctsw_worker_fn_expr(rp_alt_m::Float64, ra_alt_m::Float64)::Expr
    return quote
        function _mctsw_worker_run(n_orbits::Int)::Float64
            planet = Mars("", SPICE_PATH)
            spacecraft = make_three_body_spacecraft(
                bus_dims=(2.2, 2.6, 1.7),
                panel_dims=(0.01, 3.0, 1.5),
                bus_mass=350.0,
                panel_mass_each=10.0,
                panel_offset_y=2.6,
                ic=InitialCondition(
                    ra=planet.Rp_e + $ra_alt_m,
                    rp=planet.Rp_e + $rp_alt_m,
                    i=93.0,
                    ω=80.0,
                    Ω=30.0,
                    ν=180.0,
                ),
                reflection_coefficient=0.9,
                prop_mass=50.0,
                id=1,
            )
            a = planet.Rp_e + 0.5 * ($rp_alt_m + $ra_alt_m)
            mission_time = n_orbits * 2π * sqrt(a^3 / planet.μ)
            cfg = make_example_config(
                planet=planet,
                spacecraft=spacecraft,
                mission_time=mission_time,
                initial_time=InitialTime(year=2020, month=1, day=1, hour=0, minute=0, second=0.0),
                dynamic_effectors=(InverseSquaredGravityModel(),),
                density_model=NoAtmosphereModel(),
                orientation_sim=false,
                keplerian=false,
                EI_km=220.0,
                verbose=false,
                results=false,
                results_directory="",
            )
            return @elapsed run_simulation(cfg; return_solution=false, isolate_state=false)
        end
        nothing
    end
end

function main()
    n_workers = parse(Int, get(ENV, "SPACEAGORA_MC_DIST_NWORKERS",
                               string(Threads.nthreads())))
    n_orbits  = parse(Int, get(ENV, "SPACEAGORA_MC_THREAD_NORBITS", "30"))
    outdir    = strip(get(ENV, "SPACEAGORA_MC_THREAD_OUTDIR", "."))
    sysimage  = strip(get(ENV, "SPACEAGORA_MC_THREAD_SYSIMAGE", ""))
    project   = _MCTSW_REPO_ROOT
    common_jl = joinpath(project, "examples", "common.jl")

    @printf("[mc-ts-worker n=%d] norbits=%d  starting %d distributed worker(s)...\n",
            n_workers, n_orbits, n_workers)
    flush(stdout)

    # Mirror the exeflags pattern from aerobraking_perturbation_mc/main.jl.
    # Each worker runs single-threaded; inner parallelism is disabled via env.
    worker_exeflags = String[
        "--project=$(project)",
        "--compiled-modules=existing",
        "--threads=1",
        "--gcthreads=1,1",
    ]
    if !isempty(sysimage) && isfile(sysimage)
        push!(worker_exeflags, "--sysimage=$(sysimage)")
    end

    addprocs(n_workers;
        exeflags=worker_exeflags,
        env=["SPACEAGORA_INNER_THREAD_BUDGET"      => "1",
             "SPACEAGORA_PARALLEL_POLICY_ADAPTIVE"  => "off"],
    )

    @printf("[mc-ts-worker n=%d] initializing %d worker(s) (loading packages + defining run fn)...\n",
            n_workers, length(workers()))
    flush(stdout)

    # Step 1: activate project and load simulation code on each worker.
    pkg_init_ex = quote
        import Pkg
        Pkg.activate($project; io=devnull)
        include($common_jl)
        nothing
    end
    # Step 2: define _mctsw_worker_run on each worker.
    fn_init_ex = _mctsw_worker_fn_expr(_MCTSW_RP_ALT_M, _MCTSW_RA_ALT_M)

    init_tasks = map(workers()) do w
        @async begin
            remotecall_wait(Core.eval, w, Main, pkg_init_ex)
            remotecall_wait(Core.eval, w, Main, fn_init_ex)
        end
    end
    foreach(wait, init_tasks)

    @printf("[mc-ts-worker n=%d] warming up %d worker(s) in parallel...\n", n_workers, n_workers)
    flush(stdout)

    warmup_tasks = [
        @async remotecall_fetch(Core.eval, w, Main, :(_mctsw_worker_run($n_orbits)))
        for w in workers()
    ]
    foreach(wait, warmup_tasks)

    @printf("[mc-ts-worker n=%d] warmup done\n", n_workers)
    flush(stdout)

    @printf("[mc-ts-worker n=%d] launching %d concurrent case(s) on distributed workers...\n",
            n_workers, n_workers)
    flush(stdout)

    case_times_s = Vector{Float64}(undef, n_workers)
    next_case    = Threads.Atomic{Int}(1)

    # Work-stealing dispatch loop: mirrors _run_samples_distributed in study.jl.
    wall_s = @elapsed begin
        worker_tasks = Task[]
        for w in workers()
            task = @async begin
                while true
                    idx = Threads.atomic_add!(next_case, 1)
                    idx > n_workers && break
                    case_times_s[idx] = remotecall_fetch(
                        Core.eval, w, Main, :(_mctsw_worker_run($n_orbits)),
                    )
                end
            end
            push!(worker_tasks, task)
        end
        foreach(wait, worker_tasks)
    end

    time_per_case_s = wall_s / n_workers
    throughput      = n_workers / wall_s

    @printf(
        "[mc-ts-worker n=%d] done  wall=%.3f s  time_per_case=%.3f s  throughput=%.3f cases/s\n",
        n_workers, wall_s, time_per_case_s, throughput,
    )
    flush(stdout)

    row = (
        julia_threads       = n_workers,
        n_concurrent        = n_workers,
        wall_time_s         = wall_s,
        time_per_case_s     = time_per_case_s,
        throughput_cases_s  = throughput,
        mean_case_s         = mean(case_times_s),
        min_case_s          = minimum(case_times_s),
        max_case_s          = maximum(case_times_s),
        imbalance_pct       = 100.0 * (maximum(case_times_s) - minimum(case_times_s)) / mean(case_times_s),
        n_orbits            = n_orbits,
        case_id             = _MCTSW_CASE_ID,
    )

    outfile = joinpath(outdir, "mc_ts_worker_$(n_workers)t.csv")
    CSV.write(outfile, DataFrame([row]))
    @printf("[mc-ts-worker n=%d] wrote %s\n", n_workers, outfile)
    flush(stdout)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
