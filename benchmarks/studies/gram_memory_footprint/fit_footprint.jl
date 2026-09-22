# Fit an intercept and a per-spacecraft slope to one results CSV of the GRAM
# memory-footprint study, and print the table that
# `docs/architecture/gram_memory_footprint.md` and the routing constants in
# `src/parallel/routing/machine_topology.jl` are sourced from.
#
#   julia --project=. benchmarks/studies/gram_memory_footprint/fit_footprint.jl \
#       benchmarks/studies/gram_memory_footprint/results/gram_memory_footprint_<host>.csv
#
# Two quantities are fitted per path, both against `n_sats`:
#
#   maxrss_mb     the process's peak resident set: the conservative form, and
#                 the same column the paper_scenarios S1/S2 CSVs record, so the
#                 two studies are directly comparable. The fitted INTERCEPT is
#                 the package + SPICE + GRAM image plus one solve's fixed
#                 working set; the fitted SLOPE is the per-spacecraft term the
#                 routing constant charges. This is the fit the constants come
#                 from, because a worker that does not fit must not be started.
#   retained_mb   rss_retained_mb, live state after a full collection.
#                 Reported for contrast, never charged.
#
# `rss_base_mb` is printed per point but NOT subtracted: it is read right after
# a collection whose page release the kernel may or may not have completed, and
# it was seen to move by ~260 MB between two points of the same path. The peak
# is the stable reading; the intercept of its fit is the honest baseline.
#
# Ordinary least squares on (n_sats, y); R^2 is reported so a non-linear path is
# visible rather than averaged away.

using Printf
using Statistics

function read_csv(path::String)
    lines = readlines(path)
    isempty(lines) && error("empty CSV: $(path)")
    header = String.(split(lines[1], ","))
    rows = Dict{String, String}[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        fields = String.(split(line, ","))
        length(fields) == length(header) || continue
        push!(rows, Dict(zip(header, fields)))
    end
    return rows
end

num(row, key) = something(tryparse(Float64, get(row, key, "")), NaN)

"""OLS fit of `y = a + b*x`; returns (a, b, r2). `r2` is NaN with fewer than
three points or no variance in `x`."""
function ols(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    n >= 2 || return (NaN, NaN, NaN)
    x̄ = mean(x); ȳ = mean(y)
    sxx = sum((x .- x̄) .^ 2)
    sxx > 0 || return (NaN, NaN, NaN)
    b = sum((x .- x̄) .* (y .- ȳ)) / sxx
    a = ȳ - b * x̄
    n >= 3 || return (a, b, NaN)
    ss_res = sum((y .- (a .+ b .* x)) .^ 2)
    ss_tot = sum((y .- ȳ) .^ 2)
    r2 = ss_tot > 0 ? 1.0 - ss_res / ss_tot : NaN
    return (a, b, r2)
end

function main(paths::Vector{String})
    isempty(paths) && error("usage: fit_footprint.jl <results.csv> [...]")
    rows = reduce(vcat, read_csv.(paths))
    ok_rows = filter(r -> get(r, "ok", "false") == "true", rows)
    isempty(ok_rows) && error("no successful rows in $(paths)")

    by_path = Dict{String, Vector{Dict{String, String}}}()
    for r in ok_rows
        push!(get!(by_path, r["path"], Dict{String, String}[]), r)
    end

    println("host=$(ok_rows[1]["hostname"])  points=$(length(ok_rows))")
    println()
    @printf("%-18s %5s  %10s %10s %10s %10s %10s\n",
            "path", "N", "base_MB", "maxrss_MB", "peak-base", "retain_MB", "solve_s")
    for path in sort(collect(keys(by_path)))
        for r in sort(by_path[path]; by = r -> parse(Int, r["n_sats"]))
            base = num(r, "rss_base_mb")
            @printf("%-18s %5s  %10.1f %10.1f %10.1f %10.1f %10.1f\n",
                    path, r["n_sats"], base, num(r, "maxrss_mb"),
                    num(r, "maxrss_mb") - base, num(r, "rss_retained_mb"),
                    num(r, "solve_s"))
        end
    end

    println()
    @printf("%-18s %6s  %12s %12s %8s  %12s %12s %8s\n",
            "path", "points", "pk_int_MB", "pk_MB/sat", "R2", "rt_int_MB", "rt_MB/sat", "R2")
    for path in sort(collect(keys(by_path)))
        rs = by_path[path]
        x = Float64[parse(Float64, r["n_sats"]) for r in rs]
        wl = Float64[num(r, "maxrss_mb") for r in rs]
        rt = Float64[num(r, "rss_retained_mb") for r in rs]
        (a1, b1, r1) = ols(x, wl)
        (a2, b2, r2) = ols(x, rt)
        @printf("%-18s %6d  %12.1f %12.3f %8.4f  %12.1f %12.3f %8.4f\n",
                path, length(rs), a1, b1, r1, a2, b2, r2)
    end

    println()
    println("Package + SPICE + GRAM image (rss_base_mb), per path:")
    for path in sort(collect(keys(by_path)))
        b = Float64[num(r, "rss_base_mb") for r in by_path[path]]
        @printf("  %-18s min=%.1f MB  median=%.1f MB  max=%.1f MB  (n=%d)\n",
                path, minimum(b), median(b), maximum(b), length(b))
    end
    return nothing
end

main(String.(ARGS))
