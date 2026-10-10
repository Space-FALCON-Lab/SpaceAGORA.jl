# Where one serial solve of the coupled case allocates.
#
#   julia --project=. benchmarks/studies/coupled_local_comparison/profile_allocations.jl \
#       [case=stack32_e6_actuated_saved] [profile=smoke] [sample_rate=0.001]
#
# Runs the harness's serial mode environment, warms up twice, then reports:
#   - total allocation and GC time of one solve, for the case and for the same
#     case with trajectory output off (the case name without `_saved`), which
#     separates saving from dynamics;
#   - a sampled allocation profile (Profile.Allocs) of one solve, aggregated by
#     the innermost SpaceAGORA source line on each sampled stack, by the
#     innermost line in any package, and by allocated type. Sampled bytes are
#     scaled by 1 / sample_rate, so they are estimates.
include(joinpath(@__DIR__, "..", "parallelization_performance.jl"))
using Profile, Printf

case = get(ARGS, 1, "stack32_e6_actuated_saved")
profile = get(ARGS, 2, "smoke")
rate = parse(Float64, get(ARGS, 3, "0.001"))
top = 25

cfg = parse_parallelization_performance_cli([profile, "--cases=$case", "--modes=serial", "--threads=1"])
src_root = joinpath(pkgdir(SpaceAGORA), "src") * "/"
short(file) = replace(file, src_root => "src/", r"^.*/packages/" => "")

function aggregate(allocs, key)
    acc = Dict{String, Tuple{Float64, Int}}()
    for a in allocs
        k = key(a)
        k === nothing && continue
        b, n = get(acc, k, (0.0, 0))
        acc[k] = (b + a.size, n + 1)
    end
    return sort!(collect(acc); by=x -> -x[2][1])
end

function innermost(a, keep)
    for sf in a.stacktrace
        sf.from_c && continue
        file = string(sf.file)
        keep(file) && return @sprintf("%s:%d %s", short(file), sf.line, sf.func)
    end
    return nothing
end

function report(title, rows, total_sampled)
    println("\n== ", title)
    for (k, (b, n)) in rows[1:min(top, end)]
        @printf("%6.1f%%  %10.1f MiB est  %7d samples  %s\n",
                100 * b / total_sampled, b / rate / 2^20, n, k)
    end
end

withenv(ppc_mode_env_pairs(ppc_mode_specs()["serial"], cfg; outer_tasks=1)...) do
    solve(name) = ppc_solve_once(ppc_single_config(name, cfg; seed=cfg.seed), cfg)
    names = endswith(case, "_saved") ? [case, replace(case, r"_saved$" => "")] : [case]
    println("== Totals for one solve after two warm-ups (profile=$profile, serial)")
    for name in names
        solve(name); solve(name)
        t = @timed solve(name)
        @printf("%-28s %8.3f s  %10.1f MiB  gc %6.3f s\n", name, t.time, t.bytes / 2^20, t.gctime)
    end

    Profile.Allocs.clear()
    Profile.Allocs.@profile sample_rate=rate solve(case)
    allocs = Profile.Allocs.fetch().allocs
    total = sum(a -> Float64(a.size), allocs; init=0.0)
    @printf("\n%d sampled allocations, %.1f MiB sampled, %.1f MiB estimated at sample_rate=%g\n",
            length(allocs), total / 2^20, total / rate / 2^20, rate)

    report("Innermost SpaceAGORA source line", aggregate(allocs, a -> innermost(a, f -> startswith(f, src_root))), total)
    report("Innermost line outside Base/Core (any package)",
           aggregate(allocs, a -> innermost(a, f -> !isempty(f) && !startswith(f, "./") &&
                                                    !occursin("/share/julia/base/", f) && !occursin("/julia/stdlib/", f))), total)
    report("Allocated type", aggregate(allocs, a -> first(string(a.type), 160)), total)
end
