# Does the campaign plan move a single bit of any trajectory?
#
# Runs one Monte Carlo campaign of a parallelization_performance case under
# several campaign plans -- the same samples, the same seeds, the same
# predictive-mode environment -- through the planner's own dispatch
# (`_predictive_dispatch`), and writes every sample's full state history (every
# saved time and every state component, raw Float64, in sample order) to one
# file per plan. Compare the files with `cmp`: a plan only decides which
# consumer runs which sample, so the files must be byte-identical.
#
# The repository is the ACTIVE PROJECT, not this file's location, so the same
# script can be run against another checkout of the code to compare commits:
#
#   julia --project=<checkout> --threads=8 \
#       benchmarks/studies/predictive_cost_terms/plan_parity.jl \
#       --case=mcgrid_8sat_16mc --plans=threads8,process2+3,none --out=<prefix>
#
# Plans: `threads<W>` (threads route, W tasks, one thread each), `process<W>+<L>`
# (W pool workers and L coordinator local slots), `none` (one consumer on the
# whole thread pool); a `b<B>` suffix (`threads2b4`, `process2+3b2`) gives
# each threads-route task or local slot an inner budget of B threads, which is
# how a campaign is dumped at budget 1 against a wider budget.

using Distributed

const PCP_REPO_ROOT = dirname(Base.active_project())
const PCP_PPC_DIR = joinpath(PCP_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(PCP_PPC_DIR, "cli.jl"))
include(joinpath(PCP_PPC_DIR, "modes.jl"))
include(joinpath(PCP_PPC_DIR, "cases.jl"))
include(joinpath(PCP_PPC_DIR, "trajectory_parity.jl"))
include(joinpath(PCP_PPC_DIR, "execution.jl"))

pcp_arg(key, default) = begin
    for a in ARGS
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const PCP_CASE = pcp_arg("case", "mcgrid_8sat_16mc")
const PCP_SAMPLES = parse(Int, pcp_arg("samples", "16"))
const PCP_PLANS = String.(split(pcp_arg("plans", "threads8,process2+3,none"), ","))
const PCP_OUT = pcp_arg("out", joinpath(tempdir(), "plan_parity"))
const PCP_SEED = parse(Int, pcp_arg("seed", "20260615"))

const SCamp = SpaceAGORA.SimulationCampaigns

function pcp_plan(spec::String, n::Int)
    spec == "none" && return SCamp._predictive_plan(:none, 1, 0, n, true, nothing)
    m = match(r"^threads(\d+)(?:b(\d+))?$", spec)
    if m !== nothing
        b = m[2] === nothing ? 0 : parse(Int, m[2])
        return SCamp._predictive_plan(:threads, parse(Int, m[1]), 0, n, b <= 1, nothing;
                                      inner_budget = b)
    end
    m = match(r"^process(\d+)\+(\d+)(?:b(\d+))?$", spec)
    if m !== nothing
        b = m[3] === nothing ? 0 : parse(Int, m[3])
        return SCamp._predictive_plan(:process, parse(Int, m[1]), parse(Int, m[2]), n,
                                      parse(Int, m[2]) == 0 && b <= 1, nothing; inner_budget = b)
    end
    error("unknown plan $(spec)")
end

function main()
    cfg = PPCConfig(profile="full", solver_mode="auto_stiff")
    mode = ppc_mode_specs()["predictive"]
    case_name = PCP_CASE
    jobs = [(i, PCP_SEED + i - 1) for i in 1:PCP_SAMPLES]
    # The harness's own split: a coordinator-side sample runs under the
    # coordinator's environment, a worker-side one under the worker's.
    sample_fn = job -> begin
        idx, seed = job
        run = () -> begin
            args = ppc_single_config(case_name, cfg; seed=seed, mc_index=idx)
            sol = ppc_solve_once(args, cfg).solution
            (t=Float64.(sol.t), u=[Float64[x for x in u] for u in sol.u], retcode=string(sol.retcode))
        end
        Distributed.myid() == 1 ? run() : withenv(run, ppc_mode_env_pairs(ppc_mode_specs()["predictive"], cfg)...)
    end
    pool_workers = maximum(s -> (m = match(r"^process(\d+)", s); m === nothing ? 0 : parse(Int, m[1])), PCP_PLANS)
    if pool_workers >= 2
        ids = ppc_ensure_process_workers!(pool_workers)
        SpaceAGORA.adopt_process_workers!(SpaceAGORA.campaign_process_pool(), ids)
    end
    tuning = SCamp._campaign_route_tuning()
    println("host=$(gethostname()) threads=$(Threads.nthreads()) case=$(case_name) samples=$(PCP_SAMPLES) " *
            "project=$(PCP_REPO_ROOT)")
    for spec in PCP_PLANS
        plan = pcp_plan(spec, PCP_SAMPLES)
        r = withenv(ppc_mode_env_pairs(mode, cfg; outer_tasks=1)...) do
            SCamp._predictive_dispatch(sample_fn, jobs, plan, tuning; fail_fast=false)
        end
        samples = sort(collect(r.samples); by=s -> s.index)
        all(s -> s.success, samples) || error("plan $(spec): a sample failed: " *
            join([sprint(showerror, s.error) for s in samples if !s.success], "; "))
        path = "$(PCP_OUT)_$(spec).bin"
        steps = 0
        open(path, "w") do io
            for s in samples
                write(io, s.value.t)
                for u in s.value.u
                    write(io, u)
                end
                steps += length(s.value.t)
            end
        end
        println("plan=$(spec) route=$(r.route) consumers=$(r.threads) local_slots=$(r.local_slots) " *
                "budget=$(SCamp._predictive_declared_budget(plan)) " *
                "samples=$(length(samples)) saved_steps=$(steps) bytes=$(filesize(path)) " *
                "retcodes=$(join(unique(s.value.retcode for s in samples), ",")) wall=$(round(r.elapsed_s; digits=3))s -> $(path)")
        flush(stdout)
    end
end

main()
