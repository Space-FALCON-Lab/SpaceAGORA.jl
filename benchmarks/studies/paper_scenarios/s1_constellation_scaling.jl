# Scenario 1 — Constellation scaling baseline (Earth LEO, L20 harmonics, vacuum).
#
# Claim it backs: outer thread routing scales a coupled multi-satellite propagation
# near-linearly with constellation size on a single workstation, and the marginal
# memory cost per satellite is a state vector, not a process image. This is the
# clean-room curve every other scenario is interpreted against: no atmosphere, no
# native-library locks — any deviation from linear is parallelization overhead.
#
#   julia --project=. benchmarks/studies/paper_scenarios/s1_constellation_scaling.jl
#
# Grid (env-overridable):
#   PS_SIZES     comma list of constellation sizes   (default 1,4,16,64,256)
#   PS_THREADS   parallel thread budget              (default Sys.CPU_THREADS)
#   PS_MISSION_S mission length                      (default 1800 s)
#
# Modes per size: serial (--threads=1, serial RHS) and parallel (--threads=T,
# flat/batched RHS). Speedup(N) = serial(N) / parallel(N).

include(joinpath(@__DIR__, "common.jl"))

const SIZES = [parse(Int, s) for s in split(ps_env_str("PS_SIZES", "1,4,16,64,256"), ",")]
const THREADS = ps_threads_default()
const MISSION_S = ps_env_float("PS_MISSION_S", 1800.0)

function main()
    out_dir = ps_results_dir()
    println("S1 constellation scaling: sizes=$(SIZES) threads=$(THREADS) mission_s=$(MISSION_S)")
    rows = Dict{String, String}[]
    for n in SIZES, mode in ("serial", "parallel")
        julia_threads = mode == "serial" ? 1 : THREADS
        rhs_mode = mode == "serial" ? "serial" : "flat"
        label = "s1_n$(n)_$(mode)"
        println("[run] $(label) (--threads=$(julia_threads))")
        env = vcat(
            ps_constellation_env(rhs_mode=rhs_mode, density="none", mission_s=MISSION_S, inner_budget=julia_threads),
            [
                "PS_WORKLOAD" => "constellation",
                "PS_N_SATS" => string(n),
                "PS_MISSION_S" => string(MISSION_S),
                "PS_GRAVITY" => "l20",
                "PS_DENSITY" => "none",
                "PS_WARMUP" => string(ps_warmup()),
                "PS_REPEATS" => string(ps_repeats()),
            ],
        )
        result = ps_run_worker(label=label, julia_threads=julia_threads, env=env, log_dir=out_dir)
        push!(rows, ps_row(vcat(ps_hostname_meta(), [
            "scenario" => "s1", "n_sats" => string(n), "mode" => mode,
            "julia_threads" => string(julia_threads), "mission_s" => string(MISSION_S),
        ]), result))
        haskey(result, "median_s") && println("  -> median $(result["median_s"]) s")
    end
    ps_write_csv(joinpath(out_dir, "s1_constellation_scaling.csv"), rows, [
        "hostname", "cpu_threads", "scenario", "n_sats", "mode", "julia_threads", "mission_s",
        "ok", "median_s", "min_s", "max_s", "times_s", "maxrss_mb", "workers_rss_mb",
    ])
end

main()
