# Standalone single-point CLI for the outer-backend (threads vs. process)
# comparison at a FIXED LEO constellation size. Companion to
# leo_ensemble_outer_backend_scaling.jl (the sweep controller).
#
#   julia --project=. --threads=<T> leo_ensemble_outer_backend_scaling_worker.jl \
#       <n_sats> <outer_workers> <threads|process>
#
# outer_backend is one of:
#   threads -- run_constellation_ensemble(args; threads=outer_workers), the existing
#              fixed-Julia-thread-count API (matches leo_thread_scaling_by_mode.jl's
#              ensemble-route "standard" points). Needs --threads=<outer_workers> at
#              process startup.
#   process -- run_constellation_ensemble(args; threads=:auto, ...), routed to
#              SpaceAGORA.ParallelProcess's real OS-process backend via a route_tuning
#              that pins process_max_workers=outer_workers and forces the :process
#              route regardless of sample/mission-size thresholds, with a fresh
#              OuterRouteState() so no bandit exploration history carries over between
#              sweep points. --threads=1 is sufficient at process startup here -- the
#              process pool's own workers (each --threads=1) provide the parallelism,
#              not this coordinator process.
#
# Always mode="standard" (native GRAM point density): the workload class fixed by
# this session's GRAM-constructor-lock (data/GRAMSuite.jl, submodule) and
# GRAM/SPICE-lock-unification (src/simulation/runtime_services.jl) fixes -- previously
# this crashed outright under the ensemble route at outer_workers>=2. Compare the
# "threads" series here against leo_thread_scaling_by_mode.jl's own ensemble-route
# "standard" points (pre-fix: crash; post-fix: safe but lock-limited near ~1x
# regardless of worker count) for the full before/after story; "process" is the new
# route that sidesteps native GRAM's single process-wide lock entirely.

include(joinpath(@__DIR__, "leo_constellation_size_scaling_point.jl"))

length(ARGS) >= 3 || error(
    "Usage: julia --project=. --threads=<N> $(basename(@__FILE__)) <n_sats> <outer_workers> <threads|process>"
)
n_sats = parse(Int, ARGS[1])
outer_workers = parse(Int, ARGS[2])
outer_backend = ARGS[3]
outer_backend in ("threads", "process") || error(
    "outer_backend must be \"threads\" or \"process\", got \"$(outer_backend)\""
)

const MODE = "standard"

config_args = build_constellation_config(n_sats, MODE)
pairs = env_pairs_for(MODE, "ensemble")

function run_once()
    return withenv(pairs...) do
        if outer_backend == "threads"
            SpaceAGORA.run_constellation_ensemble(
                config_args; threads=outer_workers, return_solution=true, return_solver_metadata=true
            )
        else
            tuning = OuterRouteTuning(
                mc_process_min_samples=1,
                mc_process_min_mission_s=0.0,
                process_max_workers=outer_workers,
            )
            SpaceAGORA.run_constellation_ensemble(
                config_args;
                threads=:auto,
                route_state=OuterRouteState(),
                route_tuning=tuning,
                return_solution=true,
                return_solver_metadata=true
            )
        end
    end
end

warmup = run_once()
n_ok_warmup = length(warmup.successful)
println("warmup: $(n_ok_warmup)/$(n_sats) members succeeded")
flush(stdout)
n_ok_warmup == n_sats ||
    @warn "warmup did not report full success -- results below may not reflect a full propagation." n_sats outer_workers outer_backend

GC.gc()
t = @elapsed (result = run_once())
n_ok = length(result.successful)
println("repeat 1 ($(round(t; digits=4)) s): $(n_ok)/$(n_sats) members succeeded")
println(
    "median wall time (mode=$(MODE), route=ensemble, outer_backend=$(outer_backend), n_sat=$(n_sats), " *
    "outer_workers=$(outer_workers)): $(round(t; digits=4)) s"
)

# This process is a one-shot subprocess launched per sweep point by the controller --
# its own exit doesn't automatically clean up any SpaceAGORA.ParallelProcess pool
# workers it spawned via ensure_process_workers!. Left alone, each outer_backend
# ="process" point leaks its pool workers, and an orphaned pool worker can keep this
# process's own inherited stdout/stderr pipe open after this process exits, hanging
# the controller's `wait` on it. Shut the pool down explicitly before exiting.
outer_backend == "process" && SpaceAGORA.shutdown_process_pool!(SpaceAGORA.campaign_process_pool())
