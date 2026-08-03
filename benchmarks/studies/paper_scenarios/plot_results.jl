# Plots the paper_scenarios (S1-S5) suite for the JAIS parallelization paper.
#
# Reads whatever `results/<hostname>/sN_*.csv` files exist on disk -- no hostname
# or scenario is hardcoded, so this can be re-run at any point during the
# multi-machine data collection (partial per-machine coverage, machines added
# later, etc.) and it will pick up everything currently present. Per-host plots
# land in `results/<hostname>/plots/`; if two or more hosts have data, a
# cross-machine comparison set is also written to `results/plots_cross_machine/`.
#
# Usage:
#   julia --project=. benchmarks/studies/paper_scenarios/plot_results.jl [results_dir]
#
# results_dir defaults to benchmarks/studies/paper_scenarios/results.

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")  # headless GR before Plots loads a backend

using CSV
using DataFrames
using Plots
using Printf

const PSP_DIR = @__DIR__
const PSP_DEFAULT_RESULTS_ROOT = joinpath(PSP_DIR, "results")

# s1 is handled separately (psp_load_s1_variants) since it can have multiple
# gravity/density-variant CSVs per host; these four are always single-file.
const PSP_SCENARIO_FILES = Dict(
    "s2" => "s2_gram_atmosphere_modes.csv",
    "s3" => "s3_montecarlo_process_scaling.csv",
    "s4" => "s4_hybrid_mc_constellation.csv",
    "s5" => "s5_routing_profile_ladder.csv",
)

const PSP_S2_MODE_ORDER = ["serial_standard", "threads_standard", "threads_lookahead", "process_members"]
const PSP_S2_MODE_COLOR = Dict(
    "serial_standard" => :grey40, "threads_standard" => :darkorange,
    "threads_lookahead" => :steelblue, "process_members" => :firebrick,
)
const PSP_S2_MODE_LABEL = Dict(
    "serial_standard" => "Serial (native GRAM)",
    "threads_standard" => "Threaded (native GRAM)",
    "threads_lookahead" => "Threaded + look-ahead cache",
    "process_members" => "Process pool (1 GRAM/worker)",
)
const PSP_S2_MODE_MARKER = Dict(
    "threads_standard" => :diamond, "threads_lookahead" => :circle, "process_members" => :square,
)

const PSP_S5_PROFILE_ORDER = ["R0", "R1_a", "R2", "R3", "R4", "R5"]
const PSP_S5_PROFILE_COLOR = Dict(
    "R0" => :grey40, "R1_a" => :darkorange, "R2" => :seagreen,
    "R3" => :steelblue, "R4" => :purple, "R5" => :firebrick,
)
# Display forms for identifiers that appear verbatim in the source CSVs but
# read as internal code names rather than plot-ready labels.
const PSP_S5_PROFILE_DISPLAY = Dict("R1_a" => "R1a")
psp_profile_display(p::AbstractString)::String = get(PSP_S5_PROFILE_DISPLAY, p, p)

const PSP_S1_VARIANT_LABEL = Dict(
    "baseline" => "L20 gravity (baseline)",
    "l50_none" => "L50 gravity (heavier)",
)
psp_variant_label(v::AbstractString)::String = get(PSP_S1_VARIANT_LABEL, v, v)

const PSP_S5_WORKLOAD_LABEL = Dict(
    "const16_l20_vacuum" => "16-sat, L20 vacuum",
    "const16_gram_cache" => "16-sat, GRAM cache",
    "single_l50_vacuum" => "1-sat, L50",
)
psp_workload_label(w::AbstractString)::String = get(PSP_S5_WORKLOAD_LABEL, w, w)

const PSP_HOST_LABEL = Dict(
    "space-falcon-1" => "12-thread machine",
    "space-falcon-lab-TRX50-AERO-D" => "64-thread machine",
)
psp_host_label(h::AbstractString)::String = get(PSP_HOST_LABEL, h, h)

# ── Discovery / loading ─────────────────────────────────────────────────────

"""Hostname directories under `root` that contain at least one sN_*.csv file."""
function psp_discover_hosts(root::String)::Vector{String}
    isdir(root) || return String[]
    hosts = String[]
    for entry in sort(readdir(root))
        dir = joinpath(root, entry)
        isdir(dir) || continue
        any(f -> occursin(r"^s[1-5]_.*\.csv$", f), readdir(dir)) && push!(hosts, entry)
    end
    return hosts
end

"""Read one CSV, filtered to successful (ok=true) rows. Returns `nothing` if the
file doesn't exist, is empty, or has no successful rows -- callers treat that as
"not run (yet)", not an error."""
function psp_read_csv_ok(path::String)::Union{Nothing, DataFrame}
    isfile(path) || return nothing
    df = CSV.read(path, DataFrame; missingstring=["", "NA"], truestrings=["true"], falsestrings=["false"])
    nrow(df) == 0 && return nothing
    if "ok" in names(df)
        df = df[coalesce.(df.ok, false), :]
    end
    nrow(df) == 0 && return nothing
    return df
end

function psp_read_scenario(host_dir::String, key::String)::Union{Nothing, DataFrame}
    return psp_read_csv_ok(joinpath(host_dir, PSP_SCENARIO_FILES[key]))
end

# s1_constellation_scaling.jl writes a suffixed CSV per (gravity, density)
# combination (e.g. `s1_constellation_scaling_l50_none.csv`) so overhead-
# sensitivity variants never clobber the (l20, vacuum) paper baseline. Pick up
# every such file automatically and tag rows with a `variant` column (the
# un-suffixed baseline file gets "baseline") rather than hardcoding one name.
const PSP_S1_FILE_RE = r"^s1_constellation_scaling(?:_(.+))?\.csv$"

function psp_load_s1_variants(host_dir::String)::Union{Nothing, DataFrame}
    frames = DataFrame[]
    for file in sort(readdir(host_dir))
        m = match(PSP_S1_FILE_RE, file)
        m === nothing && continue
        variant = m.captures[1] === nothing ? "baseline" : String(m.captures[1])
        df = psp_read_csv_ok(joinpath(host_dir, file))
        df === nothing && continue
        df = copy(df)
        df[!, :variant] .= variant
        push!(frames, df)
    end
    isempty(frames) && return nothing
    return vcat(frames...; cols=:union)
end

function psp_load_host(host_dir::String)::Dict{String, DataFrame}
    data = Dict{String, DataFrame}()
    s1 = psp_load_s1_variants(host_dir)
    s1 === nothing || (data["s1"] = s1)
    for key in ("s2", "s3", "s4", "s5")
        df = psp_read_scenario(host_dir, key)
        df === nothing || (data[key] = df)
    end
    return data
end

# ── Shared plot style ───────────────────────────────────────────────────────

function psp_style(; kwargs...)
    base = (
        tickfont=Plots.font(9), guidefont=Plots.font(11), titlefont=Plots.font(12),
        legendfont=Plots.font(9), size=(760, 500), framestyle=:box,
        grid=true, gridalpha=0.25,
        left_margin=12Plots.PlotMeasures.mm, bottom_margin=10Plots.PlotMeasures.mm,
    )
    return merge(base, NamedTuple(kwargs))
end

function psp_safe(f::Function, label::String)::Vector{String}
    try
        return f()
    catch err
        @warn "Plot generation failed: $(label)" exception = (err, catch_backtrace())
        return String[]
    end
end

# ── S1: constellation scaling ───────────────────────────────────────────────

"""One variant's speedup curve as (n_sats, speedup, threads), or `nothing` if
the variant is missing serial or parallel rows."""
function psp_s1_speedup_curve(df::DataFrame)::Union{Nothing, NamedTuple}
    serial = sort(df[df.mode .== "serial", :], :n_sats)
    par = sort(df[df.mode .== "parallel", :], :n_sats)
    (nrow(serial) == 0 || nrow(par) == 0) && return nothing
    joined = innerjoin(
        rename(serial[:, [:n_sats, :median_s]], :median_s => :t_serial),
        rename(par[:, [:n_sats, :median_s]], :median_s => :t_parallel),
        on=:n_sats,
    )
    nrow(joined) == 0 && return nothing
    sort!(joined, :n_sats)
    joined.speedup = joined.t_serial ./ joined.t_parallel
    return (serial=serial, par=par, joined=joined, threads=par.julia_threads[1])
end

function psp_plot_s1_variant(df::DataFrame, outdir::String, host::String, variant::String)::Vector{String}
    paths = String[]
    curve = psp_s1_speedup_curve(df)
    curve === nothing && return paths
    (; serial, par, joined, threads) = curve
    suffix = variant == "baseline" ? "" : "_$(variant)"
    tag = variant == "baseline" ? "" : " [$(variant)]"

    p1 = Plots.plot(; psp_style(
        xlabel="Satellites (N)", ylabel="Median wall time (s)",
        xscale=:log10, yscale=:log10, legend=:topleft)...)
    Plots.plot!(p1, serial.n_sats, serial.median_s; label="Serial", marker=:circle, linewidth=2, color=:grey40)
    Plots.plot!(p1, par.n_sats, par.median_s; label="Parallel (T=$(threads))", marker=:circle, linewidth=2, color=:steelblue)
    path1 = joinpath(outdir, "s1_walltime$(suffix).pdf")
    Plots.savefig(p1, path1)
    push!(paths, path1)

    p2 = Plots.plot(; psp_style(
        xlabel="Satellites (N)", ylabel="Speedup (x)", xscale=:log10, legend=:topleft)...)
    Plots.plot!(p2, joined.n_sats, min.(joined.n_sats, threads);
        label="Ideal scaling (T=$(threads))", linestyle=:dash, color=:grey, linewidth=1)
    Plots.plot!(p2, joined.n_sats, joined.speedup; label="Measured", marker=:circle, linewidth=2, color=:steelblue)
    path2 = joinpath(outdir, "s1_speedup$(suffix).pdf")
    Plots.savefig(p2, path2)
    push!(paths, path2)
    return paths
end

"""Overlay every variant's speedup curve on one axis -- the overhead-
sensitivity comparison (e.g. cheap L20-vacuum vs. heavier L50 per-satellite
cost) isolating how much of the low-N/threshold-crossing inefficiency is
fixed dispatch overhead vs. genuine parallelization cost."""
function psp_plot_s1_overhead_sensitivity(df::DataFrame, outdir::String, host::String)::Vector{String}
    variants = sort(unique(df.variant))
    length(variants) < 2 && return String[]

    p = Plots.plot(; psp_style(
        xlabel="Satellites (N)", ylabel="Speedup vs. serial (x)", xscale=:log10, legend=:topleft)...)
    colors = Plots.palette(:tab10)
    drawn = false
    for (i, v) in enumerate(variants)
        curve = psp_s1_speedup_curve(df[df.variant .== v, :])
        curve === nothing && continue
        Plots.plot!(p, curve.joined.n_sats, curve.joined.speedup;
            label=psp_variant_label(v), marker=:circle, linewidth=2, color=colors[mod1(i, 10)])
        drawn = true
    end
    drawn || return String[]
    path = joinpath(outdir, "s1_overhead_sensitivity.pdf")
    Plots.savefig(p, path)
    return [path]
end

function psp_plot_s1(df::DataFrame, outdir::String, host::String)::Vector{String}
    paths = String[]
    variants = "variant" in names(df) ? sort(unique(df.variant)) : ["baseline"]
    for v in variants
        sub = "variant" in names(df) ? df[df.variant .== v, :] : df
        append!(paths, psp_plot_s1_variant(sub, outdir, host, v))
    end
    "variant" in names(df) && append!(paths, psp_plot_s1_overhead_sensitivity(df, outdir, host))
    return paths
end

# ── S2: GRAM atmosphere modes ───────────────────────────────────────────────

function psp_plot_s2(df::DataFrame, outdir::String, host::String)::Vector{String}
    paths = String[]
    modes = [m for m in PSP_S2_MODE_ORDER if m in unique(df.mode)]
    isempty(modes) && return paths

    p1 = Plots.plot(; psp_style(
        xlabel="Satellites (N)", ylabel="Median wall time (s)",
        xscale=:log10, yscale=:log10, legend=:topleft)...)
    for m in modes
        sub = sort(df[df.mode .== m, :], :n_sats)
        nrow(sub) == 0 && continue
        Plots.plot!(p1, sub.n_sats, sub.median_s;
            label=PSP_S2_MODE_LABEL[m], marker=:circle, linewidth=2, color=PSP_S2_MODE_COLOR[m])
    end
    path1 = joinpath(outdir, "s2_walltime.pdf")
    Plots.savefig(p1, path1)
    push!(paths, path1)

    base = sort(df[df.mode .== "serial_standard", [:n_sats, :median_s]], :n_sats)
    nrow(base) == 0 && return paths
    rename!(base, :median_s => :t_base)

    p2 = Plots.plot(; psp_style(
        xlabel="Satellites (N)", ylabel="Speedup (x)", xscale=:log10, legend=:topleft)...)
    Plots.plot!(p2, base.n_sats, fill(1.0, nrow(base));
        label="Serial baseline", linestyle=:dash, color=:grey, linewidth=1)
    for m in modes
        m == "serial_standard" && continue
        sub = sort(df[df.mode .== m, [:n_sats, :median_s]], :n_sats)
        nrow(sub) == 0 && continue
        joined = innerjoin(sub, base, on=:n_sats)
        nrow(joined) == 0 && continue
        joined.speedup = joined.t_base ./ joined.median_s
        Plots.plot!(p2, joined.n_sats, joined.speedup;
            label=PSP_S2_MODE_LABEL[m], marker=:circle, linewidth=2, color=PSP_S2_MODE_COLOR[m])
    end
    path2 = joinpath(outdir, "s2_speedup.pdf")
    Plots.savefig(p2, path2)
    push!(paths, path2)

    # Speed-per-GB centerpiece: does the look-ahead cache recover most of the
    # process-route speedup without paying the process-pool's per-GRAM-instance
    # RAM cost?
    p3 = Plots.plot(; psp_style(
        xlabel="Peak memory, coordinator + workers (GB)", ylabel="Speedup vs. serial (x)",
        legend=:topleft)...)
    any_points = false
    for m in modes
        m == "serial_standard" && continue
        sub = sort(df[df.mode .== m, [:n_sats, :median_s, :maxrss_mb, :workers_rss_mb]], :n_sats)
        nrow(sub) == 0 && continue
        joined = innerjoin(sub, base, on=:n_sats)
        nrow(joined) == 0 && continue
        joined.speedup = joined.t_base ./ joined.median_s
        joined.mem_gb = (joined.maxrss_mb .+ joined.workers_rss_mb) ./ 1024
        Plots.scatter!(p3, joined.mem_gb, joined.speedup;
            label=PSP_S2_MODE_LABEL[m], marker=get(PSP_S2_MODE_MARKER, m, :circle),
            markersize=7, color=PSP_S2_MODE_COLOR[m])
        for r in eachrow(joined)
            Plots.annotate!(p3, r.mem_gb, r.speedup,
                Plots.text("  N=$(r.n_sats)", 7, :left, color=PSP_S2_MODE_COLOR[m]))
        end
        any_points = true
    end
    if any_points
        # Left-aligned "  N=<n>" annotations extend to the right of their point;
        # widen the auto x-range so the rightmost label isn't clipped by the axis edge.
        xmin, xmax = Plots.xlims(p3)
        Plots.plot!(p3, xlims=(xmin, xmax + 0.18 * (xmax - xmin)))
        path3 = joinpath(outdir, "s2_speed_per_memory.pdf")
        Plots.savefig(p3, path3)
        push!(paths, path3)
    end
    return paths
end

# ── S3: Monte Carlo process scaling ─────────────────────────────────────────

function psp_plot_s3(df::DataFrame, outdir::String, host::String)::Vector{String}
    paths = String[]
    df = sort(df, :workers)
    nrow(df) == 0 && return paths
    samples = df.samples[1]
    t1_rows = df.median_s[df.workers .== 1]
    isempty(t1_rows) && return paths
    t1 = t1_rows[1]

    p1 = Plots.plot(; psp_style(
        xlabel="Process workers", ylabel="Median wall time (s)",
        xscale=:log10, yscale=:log10, legend=:topright)...)
    Plots.plot!(p1, df.workers, t1 ./ df.workers; label="Ideal linear scaling", linestyle=:dash, color=:grey, linewidth=1)
    Plots.plot!(p1, df.workers, df.median_s; label="Measured", marker=:circle, linewidth=2, color=:steelblue)
    path1 = joinpath(outdir, "s3_total_time.pdf")
    Plots.savefig(p1, path1)
    push!(paths, path1)

    p2 = Plots.plot(; psp_style(
        xlabel="Process workers", ylabel="Time per sample (s)", xscale=:log10, legend=:topright)...)
    Plots.plot!(p2, df.workers, df.median_s ./ samples;
        label="$(samples) samples per campaign", marker=:circle, linewidth=2, color=:steelblue)
    path2 = joinpath(outdir, "s3_per_sample_time.pdf")
    Plots.savefig(p2, path2)
    push!(paths, path2)

    p3 = Plots.plot(; psp_style(
        xlabel="Process workers", ylabel="Summed worker peak RSS (GB)", legend=:topleft)...)
    Plots.plot!(p3, df.workers, df.workers_rss_mb ./ 1024;
        label="Summed worker memory", marker=:circle, linewidth=2, color=:firebrick)
    path3 = joinpath(outdir, "s3_memory.pdf")
    Plots.savefig(p3, path3)
    push!(paths, path3)
    return paths
end

# ── S4: hybrid MC x constellation ───────────────────────────────────────────

function psp_plot_s4(df::DataFrame, outdir::String, host::String)::Vector{String}
    paths = String[]
    thr = sort(df[df.backend .== "threads", :], :outer_workers)
    nrow(thr) == 0 && return paths
    total_budget = thr.total_budget[1]

    p1 = Plots.plot(; psp_style(
        xlabel="Outer (MC) workers  [inner budget = $(total_budget) / outer]",
        ylabel="Median wall time (s)", xscale=:log10, yscale=:log10, legend=:topright,
        xticks=(thr.outer_workers, string.(thr.outer_workers)))...)
    Plots.plot!(p1, thr.outer_workers, thr.median_s;
        label="Threads (outer × inner = $(total_budget))", marker=:circle, linewidth=2, color=:steelblue)

    proc = df[df.backend .== "process", :]
    if nrow(proc) > 0
        Plots.scatter!(p1, proc.outer_workers, proc.median_s;
            label="Process route (outer=$(proc.outer_workers[1]))", marker=:diamond, markersize=8, color=:firebrick)
    end

    best_idx = argmin(thr.median_s)
    Plots.scatter!(p1, [thr.outer_workers[best_idx]], [thr.median_s[best_idx]];
        label="Best split (outer=$(thr.outer_workers[best_idx]), inner=$(thr.inner_budget[best_idx]))",
        marker=:star5, markersize=11, color=:gold)

    path1 = joinpath(outdir, "s4_partition.pdf")
    Plots.savefig(p1, path1)
    push!(paths, path1)
    return paths
end

# ── S5: routing profile ladder ──────────────────────────────────────────────

function psp_plot_s5(df::DataFrame, outdir::String, host::String)::Vector{String}
    paths = String[]
    workloads = sort(unique(df.workload))
    isempty(workloads) && return paths
    profiles = [p for p in PSP_S5_PROFILE_ORDER if p in unique(df.profile)]
    isempty(profiles) && return paths

    speedup = DataFrame(workload=String[], profile=String[], speedup=Float64[])
    for wl in workloads
        base_rows = df[(df.workload .== wl) .& (df.profile .== "R0"), :median_s]
        isempty(base_rows) && continue
        t0 = base_rows[1]
        for p in profiles
            r = df[(df.workload .== wl) .& (df.profile .== p), :median_s]
            isempty(r) && continue
            push!(speedup, (wl, p, t0 / r[1]))
        end
    end
    nrow(speedup) == 0 && return paths

    # Plots.jl (without StatsPlots) has no built-in grouped/dodged bar mode --
    # `bar_position` only accepts :stack/:overlay, which would sum profiles on
    # top of each other and misrepresent the comparison. Hand-roll dodge by
    # offsetting each profile's bars within its workload's group slot.
    n_wl, n_pr = length(workloads), length(profiles)
    group_width = 0.8
    bar_width = group_width / n_pr
    group_centers = collect(1:n_wl)

    p1 = Plots.plot(; psp_style(
        ylabel="Speedup vs. R0 (serial) (x)", legend=:outertopright,
        xticks=(group_centers, psp_workload_label.(workloads)), xrotation=15,
        tickfont=Plots.font(8), size=(860, 520),
        left_margin=14Plots.PlotMeasures.mm, bottom_margin=20Plots.PlotMeasures.mm,
        right_margin=4Plots.PlotMeasures.mm)...)
    for (j, pr) in enumerate(profiles)
        offset = (j - (n_pr + 1) / 2) * bar_width
        heights = [begin
            sub = speedup[(speedup.workload .== wl) .& (speedup.profile .== pr), :speedup]
            isempty(sub) ? 0.0 : sub[1]
        end for wl in workloads]
        Plots.bar!(p1, group_centers .+ offset, heights;
            bar_width=bar_width * 0.92, label=psp_profile_display(pr), color=PSP_S5_PROFILE_COLOR[pr], linecolor=:black)
    end
    path1 = joinpath(outdir, "s5_profile_speedup.pdf")
    Plots.savefig(p1, path1)
    push!(paths, path1)
    return paths
end

# ── Cross-machine comparison ────────────────────────────────────────────────

function psp_plot_cross_machine(all_data::Dict{String, Dict{String, DataFrame}}, outdir::String)::Vector{String}
    hosts = sort(collect(keys(all_data)))
    length(hosts) < 2 && return String[]
    mkpath(outdir)
    paths = String[]
    host_colors = Plots.palette(:tab10)

    paths_p1 = psp_safe("s1_speedup_by_machine") do
        p = Plots.plot(; psp_style(
            xlabel="Satellites (N)", ylabel="Speedup vs. serial (x)", xscale=:log10, legend=:topleft)...)
        drawn = false
        for (i, h) in enumerate(hosts)
            haskey(all_data[h], "s1") || continue
            d = all_data[h]["s1"]
            # Cross-machine comparison uses the paper-baseline (l20, vacuum)
            # variant only -- other hosts won't necessarily have run the same
            # overhead-sensitivity variants (e.g. l50), so mixing variants
            # here would join on n_sats across incompatible workloads.
            "variant" in names(d) && (d = d[d.variant .== "baseline", :])
            curve = psp_s1_speedup_curve(d)
            curve === nothing && continue
            Plots.plot!(p, curve.joined.n_sats, curve.joined.speedup;
                label=psp_host_label(h), marker=:circle, linewidth=2, color=host_colors[mod1(i, 10)])
            drawn = true
        end
        drawn || return String[]
        path = joinpath(outdir, "s1_speedup_by_machine.pdf")
        Plots.savefig(p, path)
        return [path]
    end
    append!(paths, paths_p1)

    paths_p2 = psp_safe("s3_per_sample_by_machine") do
        p = Plots.plot(; psp_style(
            xlabel="Process workers", ylabel="Time per sample (s)",
            xscale=:log10, yscale=:log10, legend=:topright)...)
        drawn = false
        for (i, h) in enumerate(hosts)
            haskey(all_data[h], "s3") || continue
            d = sort(all_data[h]["s3"], :workers)
            nrow(d) == 0 && continue
            Plots.plot!(p, d.workers, d.median_s ./ d.samples;
                label=psp_host_label(h), marker=:circle, linewidth=2, color=host_colors[mod1(i, 10)])
            drawn = true
        end
        drawn || return String[]
        path = joinpath(outdir, "s3_per_sample_by_machine.pdf")
        Plots.savefig(p, path)
        return [path]
    end
    append!(paths, paths_p2)

    paths_p3 = psp_safe("s5_profile_portability") do
        workload = "const16_l20_vacuum"
        p = Plots.plot(; psp_style(
            xlabel="Profile", ylabel="Speedup vs. R0 (x)", legend=:topleft)...)
        drawn = false
        for (i, h) in enumerate(hosts)
            haskey(all_data[h], "s5") || continue
            sub = all_data[h]["s5"]
            sub = sub[sub.workload .== workload, :]
            nrow(sub) == 0 && continue
            base = sub[sub.profile .== "R0", :median_s]
            isempty(base) && continue
            t0 = base[1]
            profiles = [pr for pr in PSP_S5_PROFILE_ORDER if pr in unique(sub.profile)]
            isempty(profiles) && continue
            sp = [t0 / sub[sub.profile .== pr, :median_s][1] for pr in profiles]
            Plots.plot!(p, psp_profile_display.(profiles), sp;
                label=psp_host_label(h), marker=:circle, linewidth=2, color=host_colors[mod1(i, 10)])
            drawn = true
        end
        drawn || return String[]
        path = joinpath(outdir, "s5_profile_portability.pdf")
        Plots.savefig(p, path)
        return [path]
    end
    append!(paths, paths_p3)

    return paths
end

# ── Report ───────────────────────────────────────────────────────────────────

function psp_write_host_report(path::String, host::String, data::Dict{String, DataFrame}, plot_paths::Vector{String})::Nothing
    open(path, "w") do io
        println(io, "# paper_scenarios plots -- $(host)")
        println(io)
        println(io, "Scenarios with data: $(join(sort(collect(keys(data))), ", "))")
        println(io, "Plots written: $(length(plot_paths))")
        println(io)
        for p in plot_paths
            println(io, "- $(basename(p))")
        end
    end
    return nothing
end

# ── Main ─────────────────────────────────────────────────────────────────────

function main(args::Vector{String}=ARGS)
    results_root = isempty(args) ? PSP_DEFAULT_RESULTS_ROOT : abspath(args[1])
    isdir(results_root) || throw(ArgumentError("Results directory does not exist: $(results_root)"))

    hosts = psp_discover_hosts(results_root)
    if isempty(hosts)
        println("[plot_results] No sN_*.csv files found under $(results_root).")
        return
    end
    println("[plot_results] hosts found: $(join(hosts, ", "))")

    all_data = Dict{String, Dict{String, DataFrame}}()
    for host in hosts
        host_dir = joinpath(results_root, host)
        data = psp_load_host(host_dir)
        all_data[host] = data
        if isempty(data)
            println("[plot_results] $(host): no successful rows in any scenario CSV, skipping.")
            continue
        end
        outdir = joinpath(host_dir, "plots")
        mkpath(outdir)

        plot_paths = String[]
        haskey(data, "s1") && append!(plot_paths, psp_safe(() -> psp_plot_s1(data["s1"], outdir, host), "$(host)/s1"))
        haskey(data, "s2") && append!(plot_paths, psp_safe(() -> psp_plot_s2(data["s2"], outdir, host), "$(host)/s2"))
        haskey(data, "s3") && append!(plot_paths, psp_safe(() -> psp_plot_s3(data["s3"], outdir, host), "$(host)/s3"))
        haskey(data, "s4") && append!(plot_paths, psp_safe(() -> psp_plot_s4(data["s4"], outdir, host), "$(host)/s4"))
        haskey(data, "s5") && append!(plot_paths, psp_safe(() -> psp_plot_s5(data["s5"], outdir, host), "$(host)/s5"))

        println("[plot_results] $(host): scenarios=$(join(sort(collect(keys(data))), ",")) plots=$(length(plot_paths)) -> $(outdir)")
        psp_write_host_report(joinpath(outdir, "README.md"), host, data, plot_paths)
    end

    cross_dir = joinpath(results_root, "plots_cross_machine")
    cross_paths = psp_plot_cross_machine(all_data, cross_dir)
    if !isempty(cross_paths)
        println("[plot_results] cross-machine: plots=$(length(cross_paths)) -> $(cross_dir)")
    end

    println("[plot_results] done.")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
