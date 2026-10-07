# Summarize a run.sh output directory: timing, allocation and bit-identity.
#
#   julia --project=. benchmarks/studies/coupled_local_comparison/compare.jl <outdir>
#
# Writes <outdir>/summary.csv and prints the tables README.md describes.
using CSV, DataFrames, SHA, Statistics, Printf

out = abspath(ARGS[1])
rows = NamedTuple[]
for tag in sort(filter(d -> isdir(joinpath(out, d)), readdir(out))), point in sort(readdir(joinpath(out, tag)))
    dir = joinpath(out, tag, point)
    isfile(joinpath(dir, "done")) || continue
    raw_file = only(filter(f -> startswith(f, "parallelization_performance_raw_"), readdir(dir)))
    raw = CSV.read(joinpath(dir, raw_file), DataFrame)
    all(raw.success) || error("failed repeat in $dir")
    # Every timed solve of the point dumps its step times and final state (the
    # warm-up solves do not); one distinct digest means all of them agree bit for bit.
    states = joinpath(dir, "states")
    digests = unique(bytes2hex(sha256(read(joinpath(states, f)))) for f in readdir(states))
    plan = "rhs_plan_mode" in names(raw) ? join(unique(string.(raw.rhs_plan_mode)), "|") : ""
    push!(rows, (
        code=tag, case=raw.case[1], mode=raw.mode[1], threads=raw.thread_count[1],
        repeats=nrow(raw), median_s=median(raw.wall_time_s),
        min_s=minimum(raw.wall_time_s), max_s=maximum(raw.wall_time_s),
        alloc_mib=median(raw.sample_alloc_mb_sum), gc_s=median(raw.sample_gc_time_sum_s),
        solves_dumped=length(readdir(states)), distinct_states=length(digests),
        state_digest=length(digests) == 1 ? digests[1] : "MIXED", plan=plan,
        git_commit=string(raw.git_commit[1])[1:min(9, end)],
    ))
end
df = DataFrame(rows)
isempty(df) && error("no completed points under $out")

serial = Dict((r.code, r.case) => r.median_s for r in eachrow(df) if r.mode == "serial")
df.speedup_vs_serial = [get(serial, (r.code, r.case), NaN) / r.median_s for r in eachrow(df)]
pre = Dict((r.case, r.mode, r.threads) => r for r in eachrow(df) if r.code == "pre")
df.patch_speedup = [haskey(pre, (r.case, r.mode, r.threads)) ? pre[(r.case, r.mode, r.threads)].median_s / r.median_s : NaN
                    for r in eachrow(df)]
CSV.write(joinpath(out, "summary.csv"), df)

println("\n== Every point (median of timed repeats; speedup vs the same code's serial run)")
for r in eachrow(sort(df, [:case, :code, :threads, :mode]))
    @printf("%-28s %-4s %-18s t=%-2d %8.3f s [%7.3f, %7.3f]  x%5.2f  alloc %9.1f MiB  gc %6.3f s  states %d/%d  %s\n",
            r.case, r.code, r.mode, r.threads, r.median_s, r.min_s, r.max_s, r.speedup_vs_serial,
            r.alloc_mib, r.gc_s, r.distinct_states, r.solves_dumped, r.plan)
end

println("\n== Adaptive (predictive, R7) against the best static route at the same thread count")
static = ("inner_only", "rhs_serial", "rhs_satellite", "rhs_per_satellite", "rhs_flat")
for g in groupby(df, [:code, :case, :threads])
    s = filter(r -> r.mode in static, g)
    a = filter(r -> r.mode == "predictive", g)
    (isempty(s) || isempty(a)) && continue
    best = s[argmin(s.median_s), :]
    @printf("%-28s %-4s t=%-2d best static %-18s %8.3f s   adaptive %8.3f s   adaptive/best %5.2f\n",
            g.case[1], g.code[1], g.threads[1], best.mode, best.median_s, a.median_s[1], a.median_s[1] / best.median_s)
end

if any(==("pre"), df.code) && any(==("post"), df.code)
    println("\n== PR #222 effect (pre / post median, same route and thread count)")
    for r in eachrow(sort(filter(r -> r.code == "post" && !isnan(r.patch_speedup), df), [:case, :threads, :mode]))
        p = pre[(r.case, r.mode, r.threads)]
        @printf("%-28s %-18s t=%-2d %8.3f -> %8.3f s  x%5.2f   alloc %9.1f -> %9.1f MiB\n",
                r.case, r.mode, r.threads, p.median_s, r.median_s, r.patch_speedup, p.alloc_mib, r.alloc_mib)
    end
end

println("\n== Bit identity of step times and final state, per case and code")
for g in groupby(df, [:case, :code])
    @printf("%-28s %-4s points=%-3d all repeats identical within points: %-5s  one state across routes/threads: %s\n",
            g.case[1], g.code[1], nrow(g), all(==(1), g.distinct_states), length(unique(g.state_digest)) == 1)
end
for g in groupby(df, :case)
    codes = unique(g.code)
    length(codes) == 2 || continue
    same = Set(g.state_digest[g.code .== "pre"]) == Set(g.state_digest[g.code .== "post"])
    @printf("%-28s pre and post produce the same state set: %s\n", g.case[1], same)
end
