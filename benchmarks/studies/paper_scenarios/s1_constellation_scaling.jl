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
#   PS_THREADS   parallel thread budget, comma list for a thread-count ladder
#                (default Sys.CPU_THREADS) -- e.g. PS_THREADS=1,2,4,8,12 sweeps
#                the same sizes across multiple thread budgets in one CSV write.
#                Serial runs once per size regardless of ladder length (it's
#                thread-count-independent); each ladder value adds one more
#                "parallel" row per size, distinguished by julia_threads.
#   PS_MISSION_S mission length                      (default 1800 s)
#   PS_GRAVITY   gravity model: invsq | l20 | l50    (default l20)
#   PS_DENSITY   density model: none | gram_standard | gram_lookahead | gram_surrogate (default none)
#   PS_RESULTS_SUFFIX  appended to the results/<hostname> dir (see common.jl) --
#                use for exploratory sweeps (e.g. a thread-count ladder) so they
#                land in their own pseudo-host directory instead of overwriting
#                the canonical per-host baseline CSV.
#
# Modes per size: serial (--threads=1, serial RHS) and parallel (--threads=T,
# flat/batched RHS) for each T in the PS_THREADS ladder. Speedup(N, T) =
# serial(N) / parallel(N, T).
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
const THREADS_LADDER = [parse(Int, s) for s in split(ps_env_str("PS_THREADS", string(Sys.CPU_THREADS)), ",")]
const MISSION_S = ps_env_float("PS_MISSION_S", 1800.0)
const GRAVITY = ps_env_str("PS_GRAVITY", "l20")
const DENSITY = ps_env_str("PS_DENSITY", "none")
const IS_BASELINE_WORKLOAD = (GRAVITY == "l20" && DENSITY == "none")
const VARIANT_SUFFIX = IS_BASELINE_WORKLOAD ? "" : "_$(GRAVITY)_$(DENSITY)"

function ps_run_point(n::Int, mode::String, julia_threads::Int, out_dir::String, rows::Vector{Dict{String, String}})
    rhs_mode = mode == "serial" ? "serial" : "flat"
    label = mode == "serial" ? "s1_n$(n)_serial$(VARIANT_SUFFIX)" : "s1_n$(n)_parallel_t$(julia_threads)$(VARIANT_SUFFIX)"
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

function main()
    out_dir = ps_results_dir()
    println("S1 constellation scaling: sizes=$(SIZES) threads=$(THREADS_LADDER) mission_s=$(MISSION_S) gravity=$(GRAVITY) density=$(DENSITY)")
    rows = Dict{String, String}[]
    for n in SIZES
        # Serial once per size -- thread-count-independent, so a multi-value
        # PS_THREADS ladder doesn't pay for it more than once.
        ps_run_point(n, "serial", 1, out_dir, rows)
        for t in THREADS_LADDER
            ps_run_point(n, "parallel", t, out_dir, rows)
        end
    end
    ps_write_csv(joinpath(out_dir, "s1_constellation_scaling$(VARIANT_SUFFIX).csv"), rows, [
        "hostname", "cpu_threads", "scenario", "n_sats", "mode", "julia_threads", "mission_s",
        "gravity", "density", "ok", "median_s", "min_s", "max_s", "times_s", "maxrss_mb", "workers_rss_mb",
    ])
end

main()
