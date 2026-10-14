# Process-pool heap-growth verification: runs a P6p-shaped campaign
# (native-GRAM one-satellite Monte Carlo samples, process route) through the
# real harness dispatch path (ppc_run_sample_batch -> ppc_ensure_process_workers!)
# and prints markers an external `ps`-sampling loop can correlate against, so
# the caller (pool_rss_probe.sh) can report peak RSS per process and total,
# with and without the pool-worker heap-size-hint, back to back.
#
# Usage (normally invoked by pool_rss_probe.sh, not by hand), always under a
# hard memory cap:
#   systemd-run --user --scope -p MemoryMax=16G -p MemorySwapMax=0 -q -- \
#       julia --project=. --threads=1 benchmarks/studies/heap_contention/pool_rss_probe.jl \
#       --case=montecarlo_mars_gram_live --workers=4 --samples=800

const HC_STUDY_DIR = @__DIR__
const HC_REPO_ROOT = normpath(joinpath(HC_STUDY_DIR, "..", "..", ".."))
const HC_PPC_DIR = joinpath(HC_REPO_ROOT, "benchmarks", "studies", "parallelization_performance")

include(joinpath(HC_PPC_DIR, "cli.jl"))
include(joinpath(HC_PPC_DIR, "modes.jl"))
include(joinpath(HC_PPC_DIR, "cases.jl"))
include(joinpath(HC_PPC_DIR, "execution.jl"))
using Printf

hc_arg(args, key, default) = begin
    for a in args
        startswith(a, "--$key=") && return String(split(a, "=", limit=2)[2])
    end
    default
end

const HC_CASE    = hc_arg(ARGS, "case", "montecarlo_mars_gram_live")
const HC_WORKERS = parse(Int, hc_arg(ARGS, "workers", "4"))
const HC_SAMPLES = parse(Int, hc_arg(ARGS, "samples", "800"))
const HC_WARMUP  = parse(Int, hc_arg(ARGS, "warmup-samples", "2"))
const HC_PROFILE = hc_arg(ARGS, "profile", "smoke")

function main()
    cfg = PPCConfig(profile=HC_PROFILE, process_workers=HC_WORKERS, worker_seed=20260923)
    case_spec = ppc_case_catalog()[HC_CASE]
    mode = ppc_mode_specs()["outer_process"]

    @printf("host=%s case=%s workers=%d samples=%d pid=%d\n",
            gethostname(), HC_CASE, HC_WORKERS, HC_SAMPLES, getpid())
    println("MARKER pool_hint_env=", get(ENV, "SPACEAGORA_POOL_WORKER_HEAP_SIZE_HINT", "<unset>"))
    flush(stdout)

    println("MARKER WARMUP_START ", time())
    flush(stdout)
    ppc_run_sample_batch(case_spec, cfg, mode, HC_WARMUP)
    println("MARKER WARMUP_END ", time())
    flush(stdout)

    println("MARKER MAIN_BATCH_START ", time())
    flush(stdout)
    t0 = time()
    ppc_run_sample_batch(case_spec, cfg, mode, HC_SAMPLES)
    wall_s = time() - t0
    println("MARKER MAIN_BATCH_END ", time())
    @printf("MARKER wall_s=%.3f\n", wall_s)
    flush(stdout)
end

main()
