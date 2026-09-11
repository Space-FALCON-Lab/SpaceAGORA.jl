# Scenario 5 — Routing-profile ladder R0–R5 across machines: the zero-tuning claim.
#
# Claim it backs: the utilization demonstrated by S1–S4 does not require expert
# hand-tuning. At a fixed thread budget, the adaptive profiles (R4/R5) match or
# closely approach the best hand-picked profile on every workload — run this same
# script on the M3 Pro MacBook, the Zen 5 desktop, and the Threadripper and the
# conclusion should hold on all three, turning "resource-constrained" from rhetoric
# into a cross-machine result.
#
#   julia --project=. benchmarks/studies/paper_scenarios/s5_routing_profile_ladder.jl
#
# Grid (env-overridable):
#   PS_THREADS   fixed thread budget for every profile (default Sys.CPU_THREADS)
#   PS_MISSION_S mission length                        (default 600 s)
#
# Workloads (three shapes that stress different routing axes):
#   const16_l20_vacuum    16 sat, L20 harmonics, vacuum      — outer-dominant
#   const16_gram_cache    16 sat, L20 + aero + look-ahead GRAM — outer + callback
#   single_l50_vacuum     1 sat, L50 harmonics, vacuum       — inner-only
#
# Profiles: R0 R1_a R2 R3 R4 R5 (R1_b omitted: process-outer routing applies to
# campaigns, not a single monolithic run_simulation). Applied via
# with_parallel_profile inside the worker; no manual RHS-mode forcing.

include(joinpath(@__DIR__, "common.jl"))

const THREADS = ps_threads_default()
const MISSION_S = ps_env_float("PS_MISSION_S", 600.0)
const PROFILES = ["R0", "R1_a", "R2", "R3", "R4", "R5"]

const WORKLOADS = [
    (name="const16_l20_vacuum", n_sats=16, gravity="l20", density="none"),
    (name="const16_gram_cache", n_sats=16, gravity="l20", density="gram_lookahead"),
    (name="single_l50_vacuum", n_sats=1, gravity="l50", density="none"),
]

function main()
    out_dir = ps_results_dir()
    println("S5 routing profile ladder: threads=$(THREADS) mission_s=$(MISSION_S) profiles=$(PROFILES)")
    rows = Dict{String, String}[]
    for wl in WORKLOADS, profile in PROFILES
        label = "s5_$(wl.name)_$(profile)"
        println("[run] $(label)")
        env = Pair{String, String}[
            "PS_WORKLOAD" => "constellation",
            "PS_N_SATS" => string(wl.n_sats),
            "PS_MISSION_S" => string(MISSION_S),
            "PS_GRAVITY" => wl.gravity,
            "PS_DENSITY" => wl.density,
            "PS_PROFILE" => profile,
            "PS_WARMUP" => string(ps_warmup()),
            "PS_REPEATS" => string(ps_repeats()),
        ]
        if wl.density == "gram_lookahead"
            # GRAM safety knobs only; routing itself comes from the profile.
            append!(env, [
                "SPACEAGORA_DENSITY_FREEZE_PER_STEP" => "0",
                "SPACEAGORA_VACUUM_GRAM_CACHE" => "1",
                "SPACEAGORA_VACUUM_GRAM_CACHE_NPOINTS" => "20",
                "SPACEAGORA_VACUUM_GRAM_CACHE_HORIZON_S" => string(MISSION_S + 500.0),
                "SPACEAGORA_VACUUM_GRAM_CACHE_DEVIATION_M" => "1e8",
            ])
        end
        result = ps_run_worker(label=label, julia_threads=THREADS, env=env, log_dir=out_dir)
        push!(rows, ps_row(vcat(ps_hostname_meta(), [
            "scenario" => "s5", "workload" => wl.name, "profile" => profile,
            "n_sats" => string(wl.n_sats), "julia_threads" => string(THREADS),
            "mission_s" => string(MISSION_S),
        ]), result))
        haskey(result, "median_s") && println("  -> median $(result["median_s"]) s")
    end
    ps_write_csv(joinpath(out_dir, "s5_routing_profile_ladder.csv"), rows, [
        "hostname", "cpu_threads", "scenario", "workload", "profile", "n_sats", "julia_threads",
        "mission_s", "ok", "median_s", "min_s", "max_s", "times_s", "maxrss_mb", "workers_rss_mb",
    ])
end

main()
