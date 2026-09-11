# Scenario 4 — Monte Carlo OF constellations: both parallelization axes at once.
#
# Claim it backs: given a FIXED thread budget on one machine, SpaceAGORA lets the
# user (or its own routing) partition that budget between outer sample-level
# workers and inner constellation-level threads — and the framework's
# oversubscription guards keep every split safe. The sweep holds total budget
# constant and moves the outer×inner boundary, answering "where should scarce
# cores go?" empirically: this composition of campaign-level and RHS-level
# parallelism under one resource governor is the capability the compared tools
# don't offer on a workstation.
#
#   julia --project=. benchmarks/studies/paper_scenarios/s4_hybrid_mc_constellation.jl
#
# Grid (env-overridable):
#   PS_MC_SAMPLES  samples per campaign                     (default 8)
#   PS_N_SATS      satellites per sample                    (default 16)
#   PS_THREADS     total thread budget T                    (default Sys.CPU_THREADS)
#   PS_MISSION_S   mission length                           (default 600 s)
#
# Splits: outer ∈ {1, 2, 4, ..., T} with inner = T ÷ outer (thread backend), plus a
# process-outer row (W = min(4, T) process workers, inner 1). Density is the
# lock-free GRAM surrogate so inner density threading is actually exercisable.

include(joinpath(@__DIR__, "common.jl"))

const SAMPLES = ps_env_int("PS_MC_SAMPLES", 8)
const N_SATS = ps_env_int("PS_N_SATS", 16)
const THREADS = ps_threads_default()
const MISSION_S = ps_env_float("PS_MISSION_S", 600.0)

function thread_splits(total::Int)::Vector{Tuple{Int, Int}}
    splits = Tuple{Int, Int}[]
    outer = 1
    while outer <= total
        push!(splits, (outer, max(1, fld(total, outer))))
        outer *= 2
    end
    return splits
end

function main()
    out_dir = ps_results_dir()
    splits = thread_splits(THREADS)
    println("S4 hybrid MC×constellation: samples=$(SAMPLES) n_sats=$(N_SATS) budget=$(THREADS) splits=$(splits) (+process)")
    rows = Dict{String, String}[]

    for (outer, inner) in splits
        label = "s4_threads_o$(outer)_i$(inner)"
        println("[run] $(label)")
        env = vcat(
            ps_mc_env(density="gram_surrogate", mission_s=MISSION_S, outer_backend="threads",
                      outer_workers=outer, inner_budget=inner),
            [
                "PS_WORKLOAD" => "mc",
                "PS_N_SATS" => string(N_SATS),
                "PS_MISSION_S" => string(MISSION_S),
                "PS_GRAVITY" => "l20",
                "PS_DENSITY" => "gram_surrogate",
                "PS_MC_SAMPLES" => string(SAMPLES),
                "PS_MC_BACKEND" => (outer == 1 ? "serial" : "threads"),
                "PS_MC_MEMBERS" => "0",
                "PS_OUTER_WORKERS" => string(outer),
                "PS_WARMUP" => string(ps_warmup()),
                "PS_REPEATS" => string(ps_repeats()),
            ],
        )
        result = ps_run_worker(label=label, julia_threads=THREADS, env=env, log_dir=out_dir)
        push!(rows, ps_row(vcat(ps_hostname_meta(), [
            "scenario" => "s4", "samples" => string(SAMPLES), "n_sats" => string(N_SATS),
            "backend" => "threads", "outer_workers" => string(outer), "inner_budget" => string(inner),
            "total_budget" => string(THREADS), "mission_s" => string(MISSION_S),
        ]), result))
        haskey(result, "median_s") && println("  -> median $(result["median_s"]) s")
    end

    proc_workers = min(4, THREADS)
    label = "s4_process_o$(proc_workers)_i1"
    println("[run] $(label)")
    env = vcat(
        ps_mc_env(density="gram_surrogate", mission_s=MISSION_S, outer_backend="process",
                  outer_workers=proc_workers, inner_budget=1),
        [
            "PS_WORKLOAD" => "mc",
            "PS_N_SATS" => string(N_SATS),
            "PS_MISSION_S" => string(MISSION_S),
            "PS_GRAVITY" => "l20",
            "PS_DENSITY" => "gram_surrogate",
            "PS_MC_SAMPLES" => string(SAMPLES),
            "PS_MC_BACKEND" => "process",
            "PS_MC_MEMBERS" => "0",
            "PS_OUTER_WORKERS" => string(proc_workers),
            "PS_WARMUP" => string(ps_warmup()),
            "PS_REPEATS" => string(ps_repeats()),
        ],
    )
    result = ps_run_worker(label=label, julia_threads=1, env=env, log_dir=out_dir)
    push!(rows, ps_row(vcat(ps_hostname_meta(), [
        "scenario" => "s4", "samples" => string(SAMPLES), "n_sats" => string(N_SATS),
        "backend" => "process", "outer_workers" => string(proc_workers), "inner_budget" => "1",
        "total_budget" => string(THREADS), "mission_s" => string(MISSION_S),
    ]), result))
    haskey(result, "median_s") && println("  -> median $(result["median_s"]) s")

    ps_write_csv(joinpath(out_dir, "s4_hybrid_mc_constellation.csv"), rows, [
        "hostname", "cpu_threads", "scenario", "samples", "n_sats", "backend", "outer_workers",
        "inner_budget", "total_budget", "mission_s",
        "ok", "median_s", "min_s", "max_s", "times_s", "maxrss_mb", "workers_rss_mb",
    ])
end

main()
