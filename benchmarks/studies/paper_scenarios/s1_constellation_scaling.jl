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
#   PS_SIZES     comma list of constellation sizes   (default 1,4,16,64,256,1024,2048,4096)
#   PS_THREADS   parallel thread budget              (default Sys.CPU_THREADS)
#   PS_MISSION_S mission length                      (default 1800 s)
#   PS_GRAVITY   gravity model: invsq | l20 | l50    (default l20)
#   PS_DENSITY   density model: none | gram_standard | gram_lookahead | gram_surrogate (default none)
#
# Modes per size: serial (--threads=1, serial RHS) and parallel (--threads=T,
# flat/batched RHS). Speedup(N) = serial(N) / parallel(N).
#
# PS_GRAVITY/PS_DENSITY exist to isolate per-satellite work cost as a variable:
# the default (l20, vacuum) is deliberately cheap per RHS call, so fixed
# per-call overhead (Polyester.@batch dispatch, the flat work-item queue
# rebuild -- see FABLE_FINDINGS.md D1) is a large fraction of wall time at
# small N, especially right where SPACEAGORA_RHS_BATCH_THREAD_THRESHOLD
# (default 16) flips batched RHS threading on. A heavier force model (e.g.
# PS_GRAVITY=l50, or a GRAM density mode) grows the numerator (real work per
# call) without moving that threshold, so it isolates whether the low-N/N=16
# efficiency dip is overhead-dominated cost or a fixed step in the routing
# policy. Non-default combinations write to a suffixed CSV/label set so they
# never clobber the (l20, none) paper-baseline results already on disk.

include(joinpath(@__DIR__, "common.jl"))

const SIZES = [parse(Int, s) for s in split(ps_env_str("PS_SIZES", "1,4,16,64,256,1024,2048,4096"), ",")]
const THREADS = ps_threads_default()
const MISSION_S = ps_env_float("PS_MISSION_S", 1800.0)
const GRAVITY = ps_env_str("PS_GRAVITY", "l20")
const DENSITY = ps_env_str("PS_DENSITY", "none")
const IS_BASELINE_WORKLOAD = (GRAVITY == "l20" && DENSITY == "none")
const VARIANT_SUFFIX = IS_BASELINE_WORKLOAD ? "" : "_$(GRAVITY)_$(DENSITY)"

function main()
    out_dir = ps_results_dir()
    println("S1 constellation scaling: sizes=$(SIZES) threads=$(THREADS) mission_s=$(MISSION_S) gravity=$(GRAVITY) density=$(DENSITY)")
    rows = Dict{String, String}[]
    for n in SIZES, mode in ("serial", "parallel")
        julia_threads = mode == "serial" ? 1 : THREADS
        rhs_mode = mode == "serial" ? "serial" : "flat"
        label = "s1_n$(n)_$(mode)$(VARIANT_SUFFIX)"
        println("[run] $(label) (--threads=$(julia_threads))")
        env = vcat(
            ps_constellation_env(rhs_mode=rhs_mode, density=DENSITY, mission_s=MISSION_S, inner_budget=julia_threads),
            [
                "PS_WORKLOAD" => "constellation",
                "PS_N_SATS" => string(n),
                "PS_MISSION_S" => string(MISSION_S),
                "PS_GRAVITY" => GRAVITY,
                "PS_DENSITY" => DENSITY,
                "PS_WARMUP" => string(ps_warmup()),
                "PS_REPEATS" => string(ps_repeats()),
            ],
        )
        result = ps_run_worker(label=label, julia_threads=julia_threads, env=env, log_dir=out_dir)
        push!(rows, ps_row(vcat(ps_hostname_meta(), [
            "scenario" => "s1", "n_sats" => string(n), "mode" => mode,
            "julia_threads" => string(julia_threads), "mission_s" => string(MISSION_S),
            "gravity" => GRAVITY, "density" => DENSITY,
        ]), result))
        haskey(result, "median_s") && println("  -> median $(result["median_s"]) s")
    end
    ps_write_csv(joinpath(out_dir, "s1_constellation_scaling$(VARIANT_SUFFIX).csv"), rows, [
        "hostname", "cpu_threads", "scenario", "n_sats", "mode", "julia_threads", "mission_s",
        "gravity", "density", "ok", "median_s", "min_s", "max_s", "times_s", "maxrss_mb", "workers_rss_mb",
    ])
end

main()
