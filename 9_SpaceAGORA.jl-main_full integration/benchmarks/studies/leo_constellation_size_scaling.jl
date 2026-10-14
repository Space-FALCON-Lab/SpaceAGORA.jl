# Constellation-size scaling sweep for the LEO GRAM constellation benchmark
# (leo_2048_constellation_gram_scaling.jl). Answers "how does wall time scale with the
# number of satellites?" holding orbit parameters, mission length, and thread count
# fixed at each mode -- as opposed to leo_2048_constellation_gram_thread_scaling.jl,
# which holds satellite count fixed and sweeps thread count instead.
#
# Julia's thread count is fixed at process startup, and "serial" mode is only a
# meaningful baseline at --threads=1 (see leo_2048_constellation_gram_scaling.jl's own
# header), so each (size, mode) pair is a separate `julia --threads=N` subprocess:
# --threads=1 for "serial", --threads=8 for "parallel". Constellation size is passed
# down via the LEO_SCALING_N_SATS env var, which the worker script reads to override
# its default N_SATS -- every other config knob (altitude, inclination, mission
# length, harmonics degree/order, tolerances) comes from the worker script's own
# `build_constellation_config`, so it's guaranteed identical to the leo constellation
# script's orbit parameters at every size.
#
# NOT run automatically as part of any other task -- the default ladder is 11 sizes x
# 2 modes = 22 full solves (warmup + 3 timed repeats each in the worker), so total wall
# time is on the order of the single largest-size runs times several (dominated by the
# 512/1024-satellite points; small sizes are fast). Run explicitly when you have the
# time budget:
#
#   julia --project=. benchmarks/studies/leo_constellation_size_scaling.jl
#
# Override the size ladder with a comma-separated list of satellite counts as the
# first arg, and/or restrict which mode(s) run via a second arg (serial|parallel|both,
# default both):
#
#   julia --project=. benchmarks/studies/leo_constellation_size_scaling.jl 1,4,16,64,256,1024
#   julia --project=. benchmarks/studies/leo_constellation_size_scaling.jl 1,4,16,64,256,1024 parallel

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const WORKER_SCRIPT = joinpath(@__DIR__, "leo_2048_constellation_gram_scaling.jl")
const DEFAULT_SIZE_LADDER = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
const SERIAL_THREADS = 1
const PARALLEL_THREADS = 8

function parse_size_ladder(raw::AbstractString)::Vector{Int}
    values = [parse(Int, strip(tok)) for tok in split(raw, ",")]
    isempty(values) && throw(ArgumentError("Satellite-count ladder must not be empty."))
    all(v -> v >= 1, values) || throw(ArgumentError("Satellite-count ladder values must be >= 1, got $values."))
    return values
end

function parse_modes(raw::AbstractString)::Vector{String}
    mode = lowercase(strip(raw))
    mode in ("serial", "parallel", "both") || throw(ArgumentError(
        "Mode must be one of: serial, parallel, both. Got \"$raw\"."
    ))
    return mode == "both" ? ["serial", "parallel"] : [mode]
end

struct SizeScalingResult
    n_sats::Int
    serial_median_s::Union{Nothing, Float64}
    parallel_median_s::Union{Nothing, Float64}
end

function run_worker(mode::String, n_sats::Int)::Union{Nothing, Float64}
    threads = mode == "serial" ? SERIAL_THREADS : PARALLEL_THREADS
    julia_bin = Base.julia_cmd().exec[1]
    cmd = Cmd([julia_bin, "--project=$(REPO_ROOT)", "--threads=$(threads)", WORKER_SCRIPT, mode])
    cmd = addenv(cmd, "LEO_SCALING_N_SATS" => string(n_sats))
    io = IOBuffer()
    run(pipeline(cmd; stdout=io, stderr=io); wait=true)
    output = String(take!(io))
    # Worker prints e.g. "median wall time (parallel, lookahead GRAM, 8 threads): 361.7 s"
    # -- two comma-separated fields before "N threads", not one, so `.*?` (not `[^,]+`)
    # is required to span both.
    m = match(r"median wall time \(.*?,\s*\d+\s*threads\):\s*([\d.]+)\s*s", output)
    if m === nothing
        println("  -> FAILED to parse median wall time from worker output; raw output follows:")
        println(output)
        return nothing
    end
    return parse(Float64, m.captures[1])
end

function main()
    ladder = length(ARGS) >= 1 ? parse_size_ladder(ARGS[1]) : DEFAULT_SIZE_LADDER
    modes = length(ARGS) >= 2 ? parse_modes(ARGS[2]) : parse_modes("both")

    println("Constellation-size scaling sweep for leo_2048_constellation_gram_scaling.jl")
    println("Satellite-count ladder: $(ladder)")
    println("Modes: $(modes) (serial @ $(SERIAL_THREADS) thread, parallel @ $(PARALLEL_THREADS) threads)")
    println()

    results = SizeScalingResult[]
    for n_sats in ladder
        serial_s = nothing
        parallel_s = nothing
        if "serial" in modes
            println("[run] n_sats=$(n_sats) mode=serial ...")
            serial_s = run_worker("serial", n_sats)
            println(serial_s === nothing ? "  -> FAILED" : "  -> median wall time: $(serial_s) s")
        end
        if "parallel" in modes
            println("[run] n_sats=$(n_sats) mode=parallel ...")
            parallel_s = run_worker("parallel", n_sats)
            println(parallel_s === nothing ? "  -> FAILED" : "  -> median wall time: $(parallel_s) s")
        end
        push!(results, SizeScalingResult(n_sats, serial_s, parallel_s))
    end

    println()
    println("Summary (orbit parameters/mission length fixed, see leo_2048_constellation_gram_scaling.jl's build_constellation_config):")
    header = rpad("n_sats", 10)
    "serial" in modes && (header *= rpad("serial_s", 12))
    "parallel" in modes && (header *= rpad("parallel_s", 12))
    ("serial" in modes && "parallel" in modes) && (header *= "speedup")
    println(header)

    for r in results
        line = rpad(string(r.n_sats), 10)
        if "serial" in modes
            line *= rpad(r.serial_median_s === nothing ? "FAILED" : string(round(r.serial_median_s; digits=3)), 12)
        end
        if "parallel" in modes
            line *= rpad(r.parallel_median_s === nothing ? "FAILED" : string(round(r.parallel_median_s; digits=3)), 12)
        end
        if "serial" in modes && "parallel" in modes
            speedup = (r.serial_median_s === nothing || r.parallel_median_s === nothing) ? NaN : r.serial_median_s / r.parallel_median_s
            line *= isnan(speedup) ? "-" : string(round(speedup; digits=3))
        end
        println(line)
    end
end

main()
