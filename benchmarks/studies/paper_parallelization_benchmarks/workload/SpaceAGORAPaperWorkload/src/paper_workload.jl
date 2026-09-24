"""
    SpaceAGORAPaperWorkload

Precompile workload for the paper benchmark harness (the P1-P6p phases of
`benchmarks/studies/paper_parallelization_benchmarks`). Loading it changes no
behavior: at precompile time it runs every (case, mode) point those phases name,
at the case's real constellation size and under the mode's route environment,
through the harness's own worker entry points with the `test` profile's 10 s
mission, so Julia caches the native code for the solver and right-hand-side
specializations in this package's image. A harness worker that loads it skips
most of its first-solve compilation.

The harness code runs in the worker's `Main`, not from here: this module
includes its own copy of the harness files only to drive the workload. What a
worker reuses from the image are the specializations of SpaceAGORA and its
dependencies (`run_simulation` on a given `SimulationConfiguration` type, the
ODE solver on its state type, CSV writing of the row type); the harness's own
functions and closures, being `Main`'s, still compile at run time.

ComponentArrays puts the satellite count into the state vector's type, so each
constellation size is its own specialization and every size the phases name is
run. GRAM-backed cases are left out; see `ppb_workload_excluded`.

See `../README.md` for building, how the harness picks the image up, and the
measured effect.
"""
module SpaceAGORAPaperWorkload

using PrecompileTools
using Distributed

const PPBW_STUDIES_DIR = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const PPBW_PPC_DIR = joinpath(PPBW_STUDIES_DIR, "parallelization_performance")
const PPBW_PPB_DIR = joinpath(PPBW_STUDIES_DIR, "paper_parallelization_benchmarks")

# The same files, in the same order, as benchmarks/studies/paper_parallelization_benchmarks.jl
# (minus the paper reporting/main, which run only in the controller). `include`
# records each file as a source dependency, so editing any of them invalidates
# the image and the harness falls back to running without it.
include(joinpath(PPBW_PPC_DIR, "cli.jl"))
include(joinpath(PPBW_PPC_DIR, "modes.jl"))
include(joinpath(PPBW_PPC_DIR, "cases.jl"))
include(joinpath(PPBW_PPC_DIR, "trajectory_parity.jl"))
include(joinpath(PPBW_PPC_DIR, "reporting.jl"))
include(joinpath(PPBW_PPC_DIR, "execution.jl"))
include(joinpath(PPBW_PPB_DIR, "cli.jl"))
include("points.jl")

# [assumed] The workload runs each point on the harness's `test` profile, whose
# mission is 10 s for every P-phase case (see the ppc_mission_time calls in
# cases.jl). Long enough to take several steps at every case's step cap (1-20 s),
# short enough that the 4096-satellite cases stay cheap. Only the types matter.
const PPBW_PROFILE = "test"

# Seed of the workload's samples; any value works, it does not enter a type.
const PPBW_SEED = 20260615

# Every state file the runtime would otherwise read or write under
# pwd()/output/parallel_policy_state is pointed into the workload's scratch
# directory (and the workload also runs with that directory as its working
# directory), so precompiling never reads this machine's learned routing state
# and never writes into it.
function _ppbw_state_env(scratch::String)
    return [
        "SPACEAGORA_OUTER_ROUTE_STATE_PATH" => joinpath(scratch, "outer_route_state.toml"),
        "SPACEAGORA_PARALLEL_POLICY_STATE_PATH" => joinpath(scratch, "policy_state.toml"),
        "SPACEAGORA_COST_CONSTANTS_PATH" => joinpath(scratch, "cost_constants.toml"),
        "SPACEAGORA_RHS_CALIBRATION_PATH" => joinpath(scratch, "rhs_calibration.toml"),
        "SPACEAGORA_CAMPAIGN_CORRECTIONS_PATH" => joinpath(scratch, "campaign_corrections.toml"),
    ]
end

function _ppbw_config(pt::PPBWorkloadPoint, scratch::String)
    outfile = joinpath(scratch, "$(pt.kind)_$(pt.case)_$(pt.mode).csv")
    return PPCConfig(
        profile = PPBW_PROFILE,
        worker = true,
        worker_case = pt.case,
        worker_mode = pt.mode,
        worker_threads = 1,
        worker_repeat = 1,
        worker_repeats = 1,
        worker_seed = PPBW_SEED,
        worker_mc_samples = pt.samples,
        worker_outfile = outfile,
        worker_parity = pt.kind === :parity,
        warmup = 0,
        # One worker: the workload must never spawn a Distributed pool. The
        # harness's own pool is provisioned only for process_workers >= 2, and
        # every mode's env caps the runtime's pool at this value
        # (SPACEAGORA_PERF_PROCS); the adaptive route withdraws below two.
        process_workers = 1,
        solver_mode = PPBConfig().solver_mode,
        parity_samples = 512,
    )
end

"""
    ppb_workload_run_point(pt, scratch)

Run one workload point the way a harness worker runs it. A pinned
process-route campaign cannot be run whole here (it would spawn a pool), so its
two halves are run separately: the coordinator's warm-up solve under the mode's
outer-split environment, and the pool worker's sample task. Adaptive campaign
modes run through the full worker path (with a one-worker cap, so the runner
stays in process) and also run the pool worker's sample task, which is what a
process-route sample executes on a worker.
"""
function ppb_workload_run_point(pt::PPBWorkloadPoint, scratch::String)
    cfg = _ppbw_config(pt, scratch)
    mode = ppc_mode_specs()[pt.mode]
    if pt.kind === :parity
        ppc_run_worker_parity(cfg)
    elseif mode.backend == "process" && pt.samples > 1
        withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=pt.samples)...) do
            ppc_solve_once(ppc_single_config(pt.case, cfg; seed=cfg.worker_seed - 1, mc_index=1), cfg)
        end
        ppc_process_sample_task(pt.case, cfg, pt.mode, 1, cfg.worker_seed)
    else
        ppc_run_worker_performance(cfg)
        if mode.backend == "auto" && pt.samples > 1
            ppc_process_sample_task(pt.case, cfg, pt.mode, 1, cfg.worker_seed)
        end
    end
    nprocs() == 1 || error("workload point $(pt.case)/$(pt.mode) started Distributed workers")
    return nothing
end

function ppb_workload_run(points=ppb_workload_points(); verbose::Bool=true)
    failed = PPBWorkloadPoint[]
    scratch = mktempdir(; prefix="spaceagora_paper_workload_")
    try
        withenv(_ppbw_state_env(scratch)...) do
            cd(scratch) do
                for (i, pt) in enumerate(points)
                    t0 = time()
                    try
                        ppb_workload_run_point(pt, scratch)
                        verbose && @info "paper workload point $(i)/$(length(points))" pt.phase pt.case pt.mode pt.kind seconds = round(time() - t0; digits = 1)
                    catch err
                        push!(failed, pt)
                        @warn "paper workload point failed; its first run will compile at run time" pt.case pt.mode pt.kind exception = (err, catch_backtrace())
                    end
                end
            end
        end
    finally
        rm(scratch; recursive=true, force=true)
    end
    return failed
end

@setup_workload begin
    points = ppb_workload_points()
    @compile_workload begin
        ppb_workload_run(points)
    end
    # Nothing the workload built may be serialized into the image. The harness
    # copy's own caches live in this module, so they would be; reset them.
    _PPC_ADAPTIVE_ROUTE_STATE[] = nothing
    _PPC_SAMPLE_FN_CACHE[] = nothing
    empty!(_PPC_GRAM_MODEL_CACHE)
    # SpaceAGORA's furnish set and planet caches belong to SpaceAGORA's image,
    # and a dependent package's precompile does not write another module's
    # state into its own image. The workload does construct SPICE-backed
    # planets, though, and this is the reset src/precompile_workload.jl
    # requires in that situation, so it is done here too rather than relying
    # on that rule.
    SpaceAGORA.SimulationModel.Planets._reset_furnished_kernels!()
    SpaceAGORA.SimulationCampaigns.reset_predictive_machine_constants!()
    SpaceAGORA.SimulationCampaigns.reset_campaign_route_state_persistence!()
end

end # module
