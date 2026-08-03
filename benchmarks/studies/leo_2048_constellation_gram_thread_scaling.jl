# Thread-scaling sweep for the LEO 1024-satellite GRAM constellation benchmark
# (leo_2048_constellation_gram_scaling.jl -- N_SATS defaults to 1024 there, overridable
# via LEO_SCALING_N_SATS). Julia's thread count is fixed at process startup, so each
# point on the curve is a separate `julia --threads=N` subprocess running that script's
# "parallel" mode (warmup + 3 timed repeats each).
#
# NOT run automatically as part of any other task -- each thread count costs 4 full
# solves of a 1024-satellite, 10-minute-mission constellation with real GRAM (~25s/solve
# measured at 1 thread; the older ~490-570s/solve figure recorded elsewhere in this repo
# was for a 2048-satellite, 1-hour mission, a different scenario, not a like-for-like
# comparison), so the full default ladder still adds up across 7 thread counts x 4
# solves. Run explicitly when you have the time budget:
#
#   julia --project=. benchmarks/studies/leo_2048_constellation_gram_thread_scaling.jl
#
# Override the ladder with a comma-separated list of thread counts as the first arg:
#
#   julia --project=. benchmarks/studies/leo_2048_constellation_gram_thread_scaling.jl 1,2,4,8
#
# Uses "parallel" mode uniformly across every thread count (rather than switching to
# "serial" mode at threads=1) so the curve isolates thread-count effects on a single
# code path -- the flat/SIMD-batch route -- rather than conflating it with the
# separately-measured true-serial route (see leo_2048_constellation_gram_scaling.jl's
# own serial-vs-parallel comparison for that baseline instead).

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const WORKER_SCRIPT = joinpath(@__DIR__, "leo_2048_constellation_gram_scaling.jl")
const DEFAULT_THREAD_LADDER = [1, 2, 4, 8, 16, 32, 64]

function parse_thread_ladder(args::Vector{String})::Vector{Int}
    isempty(args) && return DEFAULT_THREAD_LADDER
    values = [parse(Int, strip(tok)) for tok in split(args[1], ",")]
    isempty(values) && throw(ArgumentError("Thread ladder must not be empty."))
    all(v -> v >= 1, values) || throw(ArgumentError("Thread ladder values must be >= 1, got $values."))
    maximum(values) <= 64 || throw(ArgumentError(
        "Thread ladder values must be <= 64 (this sweep is scoped to \"up to 64 threads\"); got max=$(maximum(values))."
    ))
    return values
end

struct ThreadScalingResult
    threads::Int
    median_s::Union{Nothing, Float64}
    output::String
end

function run_worker(threads::Int)::ThreadScalingResult
    julia_bin = Base.julia_cmd().exec[1]
    cmd = Cmd([julia_bin, "--project=$(REPO_ROOT)", "--threads=$(threads)", WORKER_SCRIPT, "parallel"])
    io = IOBuffer()
    run(pipeline(cmd; stdout=io, stderr=io); wait=true)
    output = String(take!(io))
    # Worker prints e.g. "median wall time (parallel, lookahead GRAM, 8 threads): 361.7 s"
    # -- two comma-separated fields before "N threads", not one, so `.*?` (not `[^,]+`)
    # is required to span both.
    m = match(r"median wall time \(.*?,\s*\d+\s*threads\):\s*([\d.]+)\s*s", output)
    median_s = m === nothing ? nothing : parse(Float64, m.captures[1])
    return ThreadScalingResult(threads, median_s, output)
end

function main()
    ladder = parse_thread_ladder(ARGS)
    println("Thread-scaling sweep for leo_2048_constellation_gram_scaling.jl (mode=parallel)")
    println("Thread ladder: $(ladder)")
    println("Estimated total wall time: tens of minutes or more (4 solves per thread count, ~25s/solve measured at 1 thread; parallel-mode scaling on this workload has historically been poor -- see project memory on the GRAM global lock -- so don't expect solves to get much cheaper at higher thread counts).")
    println()

    results = ThreadScalingResult[]
    for threads in ladder
        println("[run] threads=$(threads) ...")
        r = run_worker(threads)
        push!(results, r)
        if r.median_s === nothing
            println("  -> FAILED to parse median wall time from worker output; raw output follows:")
            println(r.output)
        else
            println("  -> median wall time: $(r.median_s) s")
        end
    end

    baseline = isempty(results) ? nothing : results[1].median_s
    println()
    println("Summary (mode=parallel, N_SATS=1024, mission=600s, baseline=threads=$(ladder[1])):")
    println(rpad("threads", 10), rpad("median_s", 12), rpad("speedup", 10), "efficiency")
    for r in results
        speedup = (baseline === nothing || r.median_s === nothing) ? NaN : baseline / r.median_s
        efficiency = isnan(speedup) ? NaN : speedup / r.threads
        println(
            rpad(string(r.threads), 10),
            rpad(r.median_s === nothing ? "FAILED" : string(round(r.median_s; digits=3)), 12),
            rpad(isnan(speedup) ? "-" : string(round(speedup; digits=3)), 10),
            isnan(efficiency) ? "-" : string(round(efficiency; digits=3))
        )
    end
end

main()
