const _THREAD_SCALE_64_AERO_GRAM_REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

# Set SPACEAGORA_THREAD_SCALE_64_USE_AERO=0 to run the no-aero/no-GRAM variant
# (multi_64_nbody_srp: harmonics L20 + NBody Sun/Moon + SRP only).
# Defaults to 1 (full aero+GRAM case).
const _THREAD_SCALE_64_USE_AERO = let
    raw = lowercase(strip(get(ENV, "SPACEAGORA_THREAD_SCALE_64_USE_AERO", "1")))
    !(raw in ("0", "false", "no", "off"))
end

const _THREAD_SCALE_64_AERO_GRAM_CASE   = _THREAD_SCALE_64_USE_AERO ? "multi_64_aero_gram"  : "multi_64_nbody_srp"
const _THREAD_SCALE_64_AERO_GRAM_OUTDIR = joinpath(
    _THREAD_SCALE_64_AERO_GRAM_REPO_ROOT,
    "output",
    "performance",
    _THREAD_SCALE_64_USE_AERO ? "thread_scaling_64_aero_gram" : "thread_scaling_64_nbody_srp",
)

ENV["SPACEAGORA_THREAD_SCALE_CASE"] = get(
    ENV,
    "SPACEAGORA_THREAD_SCALE_CASE",
    _THREAD_SCALE_64_AERO_GRAM_CASE,
)
ENV["SPACEAGORA_THREAD_SCALE_OUTDIR"] = get(
    ENV,
    "SPACEAGORA_THREAD_SCALE_OUTDIR",
    _THREAD_SCALE_64_AERO_GRAM_OUTDIR,
)
ENV["SPACEAGORA_PERF_WARMUP_MISSION_SCALE"] = get(
    ENV,
    "SPACEAGORA_PERF_WARMUP_MISSION_SCALE",
    "0.1",
)

include(joinpath(@__DIR__, "performance_thread_scaling_1024.jl"))

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main_thread_scaling()
end
