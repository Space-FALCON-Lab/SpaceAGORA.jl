const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

using Statistics
using StaticArrays

include(joinpath(REPO_ROOT, "src", "simulation_model", "SimulationModel.jl"))

const PP = SimulationModel.ParallelPolicy

Base.@kwdef struct MicrobenchConfig
    effectors::Int = 4096
    work_iters::Int = 16
    repeats::Int = 12
    warmup::Int = 3
    budget::Int = max(1, Threads.nthreads())
    allotments::Vector{Int} = [1, 2, 4, 8]
    csv_path::Union{Nothing, String} = nothing
end

@inline function _parse_positive_int(raw::AbstractString, name::String)::Int
    token = strip(String(raw))
    parsed = try
        parse(Int, token)
    catch
        throw(ArgumentError("$name must be an integer, got '$raw'"))
    end
    parsed > 0 || throw(ArgumentError("$name must be > 0, got '$parsed'"))
    return parsed
end

function _parse_allotments(raw::AbstractString)::Vector{Int}
    vals = Int[]
    for token in split(raw, ",")
        t = strip(token)
        isempty(t) && continue
        push!(vals, _parse_positive_int(t, "allotment"))
    end
    isempty(vals) && throw(ArgumentError("Allotments list cannot be empty."))
    return vals
end

function parse_microbench_cli()::MicrobenchConfig
    effectors = 4096
    work_iters = 16
    repeats = 12
    warmup = 3
    budget = max(1, Threads.nthreads())
    allotments = [1, 2, 4, 8]
    csv_path = nothing
    for arg in ARGS
        if startswith(arg, "--effectors=")
            effectors = _parse_positive_int(split(arg, "=", limit=2)[2], "--effectors")
        elseif startswith(arg, "--work-iters=")
            work_iters = _parse_positive_int(split(arg, "=", limit=2)[2], "--work-iters")
        elseif startswith(arg, "--repeats=")
            repeats = _parse_positive_int(split(arg, "=", limit=2)[2], "--repeats")
        elseif startswith(arg, "--warmup=")
            warmup = _parse_positive_int(split(arg, "=", limit=2)[2], "--warmup")
        elseif startswith(arg, "--budget=")
            budget = _parse_positive_int(split(arg, "=", limit=2)[2], "--budget")
        elseif startswith(arg, "--allotments=")
            allotments = _parse_allotments(split(arg, "=", limit=2)[2])
        elseif startswith(arg, "--csv=")
            path = strip(split(arg, "=", limit=2)[2])
            csv_path = isempty(path) ? nothing : abspath(path)
        else
            throw(ArgumentError(
                "Unknown argument '$arg'. Supported: --effectors=N --work-iters=N --repeats=N --warmup=N --budget=N --allotments=a,b,c --csv=/path/out.csv"
            ))
        end
    end
    return MicrobenchConfig(
        effectors=effectors,
        work_iters=work_iters,
        repeats=repeats,
        warmup=warmup,
        budget=budget,
        allotments=allotments,
        csv_path=csv_path
    )
end

@inline function _kernel(idx::Int, work_iters::Int)
    x = 1.0e-3 * idx
    y = 2.0e-3 * idx + 0.5
    z = 3.0e-3 * idx + 1.0
    @inbounds for k in 1:work_iters
        s = sin(x + 0.017 * k) + cos(y - 0.013 * k)
        x = muladd(s, 1.0e-4, x)
        y = muladd(s, 1.7e-4, y)
        z = muladd(s, 2.3e-4, z)
    end
    fx = x + y * 0.2
    fy = y + z * 0.3
    fz = z + x * 0.4
    tx = (y - z) * 0.7
    ty = (z - x) * 0.8
    tz = (x - y) * 0.9
    return fx, fy, fz, tx, ty, tz
end

function run_serial(num_items::Int, work_iters::Int)
    acc = MVector{6, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    @inbounds for idx in 1:num_items
        fx, fy, fz, tx, ty, tz = _kernel(idx, work_iters)
        acc[1] += fx
        acc[2] += fy
        acc[3] += fz
        acc[4] += tx
        acc[5] += ty
        acc[6] += tz
    end
    return SVector{6, Float64}(acc)
end

function run_atomic(num_items::Int, allotment::Int, work_iters::Int)
    fx = Threads.Atomic{Float64}(0.0)
    fy = Threads.Atomic{Float64}(0.0)
    fz = Threads.Atomic{Float64}(0.0)
    tx = Threads.Atomic{Float64}(0.0)
    ty = Threads.Atomic{Float64}(0.0)
    tz = Threads.Atomic{Float64}(0.0)
    PP.threaded_foreach(num_items, allotment) do idx
        vfx, vfy, vfz, vtx, vty, vtz = _kernel(idx, work_iters)
        Threads.atomic_add!(fx, vfx)
        Threads.atomic_add!(fy, vfy)
        Threads.atomic_add!(fz, vfz)
        Threads.atomic_add!(tx, vtx)
        Threads.atomic_add!(ty, vty)
        Threads.atomic_add!(tz, vtz)
    end
    return SVector{6, Float64}(fx[], fy[], fz[], tx[], ty[], tz[])
end

function run_threaded_reduce(num_items::Int, allotment::Int, work_iters::Int)
    reduced = PP.threaded_reduce(
        num_items,
        allotment,
        () -> MVector{6, Float64}(0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (local_acc, idx) -> begin
            fx, fy, fz, tx, ty, tz = _kernel(idx, work_iters)
            local_acc[1] += fx
            local_acc[2] += fy
            local_acc[3] += fz
            local_acc[4] += tx
            local_acc[5] += ty
            local_acc[6] += tz
            return nothing
        end,
        (dest, src) -> begin
            @inbounds for i in 1:6
                dest[i] += src[i]
            end
            return nothing
        end
    )
    return SVector{6, Float64}(reduced)
end

function _bench(f::Function, warmup::Int, repeats::Int)
    for _ in 1:warmup
        f()
    end
    samples_ms = Vector{Float64}(undef, repeats)
    out = f()
    for i in 1:repeats
        t0 = time_ns()
        out = f()
        samples_ms[i] = (time_ns() - t0) / 1.0e6
    end
    return (result=out, median_ms=median(samples_ms), mean_ms=mean(samples_ms), min_ms=minimum(samples_ms), max_ms=maximum(samples_ms))
end

@inline _max_abs_diff(a::SVector{6, Float64}, b::SVector{6, Float64}) = maximum(abs.(a .- b))

function _write_csv(path::String, rows)
    dir = dirname(path)
    isempty(dir) || mkpath(dir)
    open(path, "w") do io
        println(io, "allotment,atomic_median_ms,reduce_median_ms,reduce_speedup_vs_atomic,max_abs_diff")
        for row in rows
            println(
                io,
                string(
                    row.allotment, ",",
                    row.atomic_ms, ",",
                    row.reduce_ms, ",",
                    row.speedup, ",",
                    row.max_abs_diff
                )
            )
        end
    end
end

function main()
    cfg = parse_microbench_cli()
    budget = min(cfg.budget, max(1, Threads.nthreads()))
    allotments = sort(unique(max(1, a) for a in cfg.allotments))

    println("== Dynamic Effector Reduction Microbenchmark ==")
    println("Machine threads: ", Threads.nthreads())
    println("Inner budget: ", budget)
    println("Effectors: ", cfg.effectors)
    println("Work iters per effector: ", cfg.work_iters)
    println("Warmup/repeats: ", cfg.warmup, "/", cfg.repeats)
    println("Allotments: ", join(allotments, ", "))
    println("")

    withenv("SPACEAGORA_INNER_THREAD_BUDGET" => string(budget)) do
        serial_stats = _bench(() -> run_serial(cfg.effectors, cfg.work_iters), cfg.warmup, cfg.repeats)
        println("serial median: ", round(serial_stats.median_ms; digits=3), " ms")
        println("")

        rows = NamedTuple[]
        println("allotment | atomic_ms | reduce_ms | speedup | max_abs_diff")
        println("--------- | --------- | --------- | ------- | ------------")
        for allotment in allotments
            atomic_stats = _bench(() -> run_atomic(cfg.effectors, allotment, cfg.work_iters), cfg.warmup, cfg.repeats)
            reduce_stats = _bench(() -> run_threaded_reduce(cfg.effectors, allotment, cfg.work_iters), cfg.warmup, cfg.repeats)
            diff = _max_abs_diff(atomic_stats.result, reduce_stats.result)
            speedup = atomic_stats.median_ms / reduce_stats.median_ms
            push!(rows, (
                allotment=allotment,
                atomic_ms=round(atomic_stats.median_ms; digits=3),
                reduce_ms=round(reduce_stats.median_ms; digits=3),
                speedup=round(speedup; digits=3),
                max_abs_diff=diff
            ))
            println(
                lpad(string(allotment), 9), " | ",
                lpad(string(round(atomic_stats.median_ms; digits=3)), 9), " | ",
                lpad(string(round(reduce_stats.median_ms; digits=3)), 9), " | ",
                lpad(string(round(speedup; digits=3)), 7), " | ",
                diff
            )
        end

        println("")
        best_row = rows[argmax([r.speedup for r in rows])]
        println("best reduce speedup vs atomic: ", best_row.speedup, "x at allotment=", best_row.allotment)
        serial_best = serial_stats.median_ms / minimum(r.reduce_ms for r in rows)
        println("best reduce speedup vs serial: ", round(serial_best; digits=3), "x")

        if cfg.csv_path !== nothing
            _write_csv(cfg.csv_path, rows)
            println("wrote CSV: ", cfg.csv_path)
        end
    end
end

main()
