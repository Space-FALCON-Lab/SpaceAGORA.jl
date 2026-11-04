# Scenario 2 — High-fidelity constellation with native GRAM: the resource-constraint
# centerpiece.
#
# Claim it backs: when the atmosphere model is a non-thread-safe native library,
# the brute-force fix (one OS process per worker, each with its own 2–4 GB
# GRAM+SPICE image) is exactly what a RAM-constrained machine cannot afford — and
# SpaceAGORA's look-ahead density cache recovers most of the process-route speedup
# with ONE shared GRAM instance. Each point reports peak RSS (coordinator +, for
# the process mode, summed pool workers) alongside wall time, so the paper can show
# the speed-per-GB trade explicitly.
#
#   julia --project=. benchmarks/studies/paper_scenarios/s2_gram_atmosphere_modes.jl
#
# Grid (env-overridable):
#   PS_SIZES        constellation sizes                  (default 4,16,64)
#   PS_THREADS      parallel thread budget               (default Sys.CPU_THREADS)
#   PS_MISSION_S    mission length                       (default 600 s)
#   PS_PROC_WORKERS process-pool workers for the process mode (default min(8, threads))
#
# Modes per size:
#   serial_standard    --threads=1, direct native GRAM (freeze-per-step) — baseline
#   threads_standard   --threads=T, RHS threaded but density lock-bound — "safe, not fast"
#   threads_lookahead  --threads=T, vacuum-predicted look-ahead density cache
#   process_members    per-satellite propagation on W --threads=1 pool processes

include(joinpath(@__DIR__, "common.jl"))

const SIZES = [parse(Int, s) for s in split(ps_env_str("PS_SIZES", "4,16,64"), ",")]
const THREADS = ps_threads_default()
const MISSION_S = ps_env_float("PS_MISSION_S", 600.0)
const PROC_WORKERS = ps_env_int("PS_PROC_WORKERS", min(8, THREADS))

function point_env(mode::String, n::Int)::Tuple{Int, Vector{Pair{String, String}}}
    common = [
        "PS_N_SATS" => string(n),
        "PS_MISSION_S" => string(MISSION_S),
        "PS_GRAVITY" => "l20",
        "PS_WARMUP" => string(ps_warmup()),
        "PS_REPEATS" => string(ps_repeats()),
    ]
    if mode == "serial_standard"
        return 1, vcat(
            ps_constellation_env(rhs_mode="serial", density="gram_standard", mission_s=MISSION_S, inner_budget=1),
            common, ["PS_WORKLOAD" => "constellation", "PS_DENSITY" => "gram_standard"],
        )
    elseif mode == "threads_standard"
        return THREADS, vcat(
            ps_constellation_env(rhs_mode="flat", density="gram_standard", mission_s=MISSION_S, inner_budget=THREADS),
            common, ["PS_WORKLOAD" => "constellation", "PS_DENSITY" => "gram_standard"],
        )
    elseif mode == "threads_lookahead"
        return THREADS, vcat(
            ps_constellation_env(rhs_mode="flat", density="gram_lookahead", mission_s=MISSION_S, inner_budget=THREADS),
            common, ["PS_WORKLOAD" => "constellation", "PS_DENSITY" => "gram_lookahead"],
        )
    end
    # process_members: coordinator only dispatches (1 thread); each of the W pool
    # workers is an independent --threads=1 process with its own native GRAM.
    return 1, vcat(
        ps_mc_env(density="gram_standard", mission_s=MISSION_S, outer_backend="process",
                  outer_workers=PROC_WORKERS, inner_budget=1),
        common, [
            "PS_WORKLOAD" => "mc",
            "PS_DENSITY" => "gram_standard",
            "PS_MC_BACKEND" => "process",
            "PS_MC_MEMBERS" => "1",
            "PS_OUTER_WORKERS" => string(PROC_WORKERS),
        ],
    )
end

function main()
    out_dir = ps_results_dir()
    println("S2 GRAM atmosphere modes: sizes=$(SIZES) threads=$(THREADS) proc_workers=$(PROC_WORKERS) mission_s=$(MISSION_S)")
    rows = Dict{String, String}[]
    for n in SIZES, mode in ("serial_standard", "threads_standard", "threads_lookahead", "process_members")
        julia_threads, env = point_env(mode, n)
        label = "s2_n$(n)_$(mode)"
        println("[run] $(label) (--threads=$(julia_threads))")
        result = ps_run_worker(label=label, julia_threads=julia_threads, env=env, log_dir=out_dir)
        push!(rows, ps_row(vcat(ps_hostname_meta(), [
            "scenario" => "s2", "n_sats" => string(n), "mode" => mode,
            "julia_threads" => string(julia_threads), "proc_workers" => (mode == "process_members" ? string(PROC_WORKERS) : "0"),
            "mission_s" => string(MISSION_S),
        ]), result))
        haskey(result, "median_s") && println("  -> median $(result["median_s"]) s  rss $(get(result, "maxrss_mb", "?")) MB (+$(get(result, "workers_rss_mb", "0")) MB workers)")
    end
    ps_write_csv(joinpath(out_dir, "s2_gram_atmosphere_modes.csv"), rows, [
        "hostname", "cpu_threads", "scenario", "n_sats", "mode", "julia_threads", "proc_workers", "mission_s",
        "ok", "median_s", "min_s", "max_s", "times_s", "maxrss_mb", "workers_rss_mb",
    ])
end

main()
