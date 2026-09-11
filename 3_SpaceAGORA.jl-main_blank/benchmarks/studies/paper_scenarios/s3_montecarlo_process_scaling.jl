# Scenario 3 — Monte Carlo throughput scaling across process workers (Earth GRAM).
#
# Claim it backs: for independent uncertainty samples, SpaceAGORA's process route
# delivers near-linear throughput up to the machine's core/RAM ceiling — and the
# engineering that makes it usable (locked GRAM construction, world-age fixes, GRAM
# model serialization, one-time per-worker JIT warmup outside the timed region) is
# what turns "embarrassingly parallel in principle" into measured speedup. Reported
# per-worker peak RSS quantifies the RAM cost per unit of process parallelism that
# Scenario 2's shared-instance route avoids.
#
#   julia --project=. benchmarks/studies/paper_scenarios/s3_montecarlo_process_scaling.jl
#
# Grid (env-overridable):
#   PS_MC_SAMPLES  samples per campaign            (default 32)
#   PS_WORKERS     comma worker ladder             (default 1,2,4,8 capped at PS_THREADS,
#                                                   extend to 16,32 on the big machines)
#   PS_MISSION_S   mission length per sample       (default 600 s)
#
# Each sample: 1 satellite, Earth L20 + aero + direct native GRAM (freeze-per-step),
# per-seed altitude/RAAN jitter. serial (in-process) anchors the ladder at W=1.

include(joinpath(@__DIR__, "common.jl"))

const SAMPLES = ps_env_int("PS_MC_SAMPLES", 32)
const THREADS = ps_threads_default()
const WORKER_LADDER = [w for w in
    [parse(Int, s) for s in split(ps_env_str("PS_WORKERS", "1,2,4,8"), ",")] if w <= THREADS]
const MISSION_S = ps_env_float("PS_MISSION_S", 600.0)

function main()
    out_dir = ps_results_dir()
    println("S3 Monte Carlo process scaling: samples=$(SAMPLES) workers=$(WORKER_LADDER) mission_s=$(MISSION_S)")
    rows = Dict{String, String}[]
    for w in WORKER_LADDER
        backend = w == 1 ? "serial" : "process"
        label = "s3_w$(w)_$(backend)"
        println("[run] $(label)")
        env = vcat(
            ps_mc_env(density="gram_standard", mission_s=MISSION_S, outer_backend=backend,
                      outer_workers=w, inner_budget=1),
            [
                "PS_WORKLOAD" => "mc",
                "PS_N_SATS" => "1",
                "PS_MISSION_S" => string(MISSION_S),
                "PS_GRAVITY" => "l20",
                "PS_DENSITY" => "gram_standard",
                "PS_MC_SAMPLES" => string(SAMPLES),
                "PS_MC_BACKEND" => backend,
                "PS_MC_MEMBERS" => "0",
                "PS_OUTER_WORKERS" => string(w),
                "PS_WARMUP" => string(ps_warmup()),
                "PS_REPEATS" => string(ps_repeats()),
            ],
        )
        result = ps_run_worker(label=label, julia_threads=1, env=env, log_dir=out_dir)
        push!(rows, ps_row(vcat(ps_hostname_meta(), [
            "scenario" => "s3", "samples" => string(SAMPLES), "workers" => string(w),
            "backend" => backend, "mission_s" => string(MISSION_S),
        ]), result))
        if haskey(result, "median_s")
            per_sample = parse(Float64, result["median_s"]) / SAMPLES
            println("  -> median $(result["median_s"]) s ($(round(per_sample; digits=3)) s/sample)")
        end
    end
    ps_write_csv(joinpath(out_dir, "s3_montecarlo_process_scaling.csv"), rows, [
        "hostname", "cpu_threads", "scenario", "samples", "workers", "backend", "mission_s",
        "ok", "median_s", "min_s", "max_s", "times_s", "maxrss_mb", "workers_rss_mb",
    ])
end

main()
