# Where does the image's timed-repeat slowdown come from? A worker-shaped probe.
#
# Builds one paper point's configuration the way a harness worker does, solves it
# once, then times the right-hand side alone (the ODE function the solver calls,
# at a fixed state) in rounds, and prints where the native code of the RHS call
# chain lives: its address, and whether it was compiled in this process
# (time_compile > 0) or loaded from a package image (time_compile == 0).
#
# Launch it twice, once as a stock worker and once with the workload on the load
# path, one process at a time on a quiet machine:
#
#   julia --threads=1 --project=. image_layout_probe.jl p1 15
#   JULIA_LOAD_PATH="@:output/paper_workload/env:@v#.#:@stdlib" SPACEAGORA_PPC_WORKLOAD=1 \
#     julia --threads=1 --project=. image_layout_probe.jl p1 15
#
# Points: p1 (gravity_1sat_l50_vacuum_4150000s, serial), p3
# (independent_1sat_1hr, one pool-worker sample), p5 (mcgrid_16sat_8mc, one
# sample). The findings this produced are in README.md ("Why the timed repeats
# move").

const PPBL_REPO = normpath(joinpath(@__DIR__, "..", "..", "..", ".."))
const PPBL_PPC = joinpath(PPBL_REPO, "benchmarks", "studies", "parallelization_performance")
for f in ("cli.jl", "modes.jl", "cases.jl", "trajectory_parity.jl", "reporting.jl", "execution.jl")
    include(joinpath(PPBL_PPC, f))
end
using Statistics

const PPBL_POINTS = Dict(
    "p1" => (case="gravity_1sat_l50_vacuum_4150000s", mode="serial", iters=2000),
    "p3" => (case="independent_1sat_1hr", mode="outer_process", iters=20000),
    "p5" => (case="mcgrid_16sat_8mc", mode="outer_inner_static", iters=2000),
)

function ppbl_rhs_round(f, du, u, p, t, n)
    t0 = time_ns()
    for _ in 1:n
        f(du, u, p, t)
    end
    return (time_ns() - t0) / n
end

# The compiled specializations of the RHS call chain's functions for this
# point's parameter type, with the address of their native code and where that
# code came from.
function ppbl_code_addresses(ptype::Type)
    pstr = string(ptype)
    relevant(mi) = occursin(pstr, string(mi.specTypes)) ||
        (length(mi.specTypes.parameters) >= 2 && occursin(string(mi.specTypes.parameters[2]), pstr))
    E = SpaceAGORA.SimulationEngine
    PE = SpaceAGORA.SimulationModel.DynamicEffectors.PerturbationEffectors
    fns = [E.spacecraft_dynamics!, E._spacecraft_dynamics_dispatch!, E._evaluate_dynamic_effector,
           SpaceAGORA.SimulationModel.calcForceTorque, PE._harmonics_scalar_force_ii]
    rows = NamedTuple[]
    for fn in fns, m in methods(fn), mi in Base.specializations(m)
        (mi === nothing || !relevant(mi)) && continue
        ci = isdefined(mi, :cache) ? mi.cache : nothing
        while ci !== nothing
            sp = UInt(getfield(ci, :specptr))
            sp == 0 || push!(rows, (func=string(nameof(fn)), address=sp, from_image=ci.time_compile == 0))
            ci = isdefined(ci, :next) ? ci.next : nothing
        end
    end
    return rows
end

function ppbl_main(point::String, rounds::Int)
    spec = PPBL_POINTS[point]
    workload = any(p -> p.name == PPC_WORKLOAD_PACKAGE, keys(Base.loaded_modules))
    cfg = PPCConfig(profile="full", worker=true, worker_case=spec.case, worker_mode=spec.mode,
                    worker_threads=Threads.nthreads(), worker_seed=20260616, worker_mc_samples=1,
                    warmup=0, process_workers=1, solver_mode="auto_stiff", parity_samples=512)
    envp = ppc_mode_env_pairs(ppc_mode_specs()[spec.mode], cfg; outer_tasks=1)
    ptype = withenv(envp...) do
        sol = ppc_solve_once(ppc_single_config(spec.case, cfg; seed=20260616, mc_index=1), cfg).solution
        k = length(sol.t) ÷ 2 + 1
        u = copy(sol.u[k]); du = similar(u)
        f, p, t = sol.prob.f, sol.prob.p, sol.t[k]
        ppbl_rhs_round(f, du, u, p, t, 10)
        ns = [ppbl_rhs_round(f, du, u, p, t, spec.iters) for _ in 1:rounds]
        println("probe point=$(point) workload=$(workload) rhs_ns median=$(round(median(ns); digits=1)) " *
                "min=$(round(minimum(ns); digits=1)) max=$(round(maximum(ns); digits=1))")
        typeof(p)
    end
    rows = ppbl_code_addresses(ptype)
    for r in rows
        println("code func=$(r.func) address=0x$(string(r.address; base=16)) from_image=$(r.from_image)")
    end
    if !isempty(rows)
        lo, hi = extrema(r.address for r in rows)
        println("code span=$(round((hi - lo) / 2^20; digits=2)) MiB over $(length(rows)) specializations")
    end
end

ppbl_main(get(ARGS, 1, "p1"), parse(Int, get(ARGS, 2, "15")))
