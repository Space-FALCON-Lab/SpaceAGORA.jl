# Build summary tables and plots for perturbation force magnitude / central-body gravity.
#
# Usage: julia --project=. benchmarks/studies/aerobraking_perturbation_mc/plot_perturbation_force_ratios.jl [run_dir]
#        Default: use the most recent timestamped run in output/aerobraking_perturbation_mc/.
# Options:
#   --planets mars,earth,venus,titan   (default: planets discovered in the run)
#   --out-dir PATH               directory for output plots/tables (default: <run_dir>/plots/)
#   --mass-scale SCALE           spacecraft mass scale for plotted slices (default: 1.0)
#   --inclination DEG            inclination held fixed for baseline/argp plots (default: 93)
#   --argp DEG                   argument of periapsis held fixed for baseline/inclination plots (default: 80)
#   --rebuild-summary            rebuild perturbation_force_ratio_summary.{csv,feather}
#   --metric peak|p95|p50|max_in_atmosphere
#   --plot-set all|time|heatmaps|rankings
#   --density-case nominal|low|high|all
#   --time-apo-panels N        max apoapsis columns for time/sweep panel grids (default: 6)
#   --summary-threads N          summary worker tasks (default: auto-threaded Julia)
#   --plot-threads N             plot worker tasks for independent PDF outputs (default: auto-threaded Julia)
#   --analytical-overlays none|basic|detailed|both

using Arrow
using CSV
using DataFrames
using Plots
using Printf
using Statistics

include(joinpath(@__DIR__, "analytical_perturbation_models.jl"))
using .AnalyticalPerturbationModels

gr()

const REPO_ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_OUTPUT_BASE = joinpath(REPO_ROOT, "output", "aerobraking_perturbation_mc")
const SUMMARY_BASENAME = "perturbation_force_ratio_summary"
const THREAD_REEXEC_ENV = "SPACEAGORA_AERO_PERTURB_THREAD_REEXEC"
const THREAD_REEXEC_DISABLE_ENV = "SPACEAGORA_AERO_PERTURB_NO_THREAD_REEXEC"
const AUTO_THREADS_ENV = "SPACEAGORA_AERO_PERTURB_AUTO_THREADS"
const PLOT_THREADS_ENV = "SPACEAGORA_AERO_PERTURB_PLOT_THREADS"

# Gravitational parameters (m^3/s^2) matching planets.jl.
const PLANET_GM = Dict(
    "mars"  => 4.282837285418775e13,
    "earth" => 3.98600436233e14,
    "venus" => 3.24858599e14,
    "titan" => 8.981e12,
)
const PLANET_RADIUS_M = Dict(
    "mars"  => 3.396190e6,
    "earth" => 6.3781366e6,
    "venus" => 6.0518e6,
    "titan" => 2.575e6,
)

const PERIAPSIS_REGIME_ORDER = ["shallow", "nominal", "deep"]
const PLANET_ORDER = ["mars", "earth", "venus", "titan"]
const METRIC_COLUMN = Dict(
    "peak" => :ratio_peak,
    "p95" => :ratio_p95,
    "p50" => :ratio_p50,
    "max_in_atmosphere" => :ratio_peak_in_atmosphere,
)
const GRID_LEFT_MARGIN = 12Plots.mm
const GRID_BOTTOM_MARGIN = 12Plots.mm
const GRID_TOP_MARGIN = 5Plots.mm
const DEFAULT_TIME_APO_PANELS = 6
const MAX_HEATMAP_XTICKS = 12
const ANALYTICAL_OVERLAY_MODES = ("none", "basic", "detailed", "both")
const PLOT_SAVE_LOCK = ReentrantLock()

# (dynamics_case, density_case, line_color, line_style)
const TRACE_STYLES = [
    ("j2",               "none",    "#1f77b4", :solid),
    ("harmonics_low",    "none",    "#2ca02c", :solid),
    ("srp",              "none",    "#ff7f0e", :solid),
    ("third_body_sun",   "none",    "#d62728", :solid),
    ("gram_aero",        "high",    "#7b2d8b", :solid),
    ("gram_aero",        "nominal", "#9e4fc4", :solid),
    ("gram_aero",        "low",     "#c99be0", :dot),
    ("full_environment", "high",    "#8c3a00", :solid),
    ("full_environment", "nominal", "#bf5e1a", :solid),
    ("full_environment", "low",     "#e0a070", :dot),
]

const ANALYTICAL_COMPARISON_STYLES = [
    (:j2, "j2", "none", "#1f77b4"),
    (:harmonics, "harmonics_low", "none", "#2ca02c"),
    (:srp, "srp", "none", "#ff7f0e"),
    (:third_body, "third_body_sun", "none", "#d62728"),
    (:drag, "gram_aero", "nominal", "#9e4fc4"),
]

function _truthy_env(name::String)::Bool
    lowercase(strip(get(ENV, name, ""))) in ("1", "true", "yes", "on")
end

function _ensure_multithreaded_cli!(args=ARGS)::Bool
    if abspath(PROGRAM_FILE) != @__FILE__
        return false
    end
    Threads.nthreads() > 1 && return false
    if _truthy_env(THREAD_REEXEC_ENV) || _truthy_env(THREAD_REEXEC_DISABLE_ENV)
        return false
    end

    requested = strip(get(ENV, AUTO_THREADS_ENV, "auto"))
    isempty(requested) && (requested = "auto")
    project = Base.active_project()
    parts = copy(Base.julia_cmd().exec)
    project !== nothing && push!(parts, "--project=$(project)")
    push!(parts, "--threads=$(requested)")
    push!(parts, abspath(PROGRAM_FILE))
    append!(parts, String.(args))

    println("Restarting with Julia --threads=$(requested) for multithreaded summary/plot generation...")
    flush(stdout)
    run(setenv(Cmd(parts), THREAD_REEXEC_ENV => "1"))
    exit(0)
    return true
end

Base.@kwdef struct PlotOptions
    run_dir::String
    planets::Vector{String}
    out_dir::String
    inferred_latest_run::Bool = false
    spacecraft_mass_scale::Float64 = 1.0
    inclination_deg::Float64 = 93.0
    argp_deg::Float64 = 80.0
    rebuild_summary::Bool = false
    metric::String = "peak"
    plot_set::String = "all"
    density_case::String = "nominal"
    time_apo_panels::Int = DEFAULT_TIME_APO_PANELS
    summary_threads::Int = Threads.nthreads()
    plot_threads::Int = Threads.nthreads()
    analytical_overlays::String = "both"
end

function _parse_case_dir(dirname::String)
    m = match(
        r"^case_\d+_(mars|earth|venus|titan)_(shallow|nominal|deep)_apo_(\d+)km_(.+)_density_(none|low|nominal|high)_ms(\d+p\d+)_inc(\d+)_aop(\d+)$",
        dirname,
    )
    m === nothing && return nothing
    return (;
        planet=m[1],
        periapsis=m[2],
        apoapsis_alt_km=parse(Float64, m[3]),
        dynamics=m[4],
        density=m[5],
        spacecraft_mass_scale=parse(Float64, replace(m[6], "p" => ".")),
        inclination_deg=parse(Float64, m[7]),
        argp_deg=parse(Float64, m[8]),
    )
end

function _discover_planets(run_dir::String)
    planets = Set{String}()
    for entry in readdir(run_dir; join=true)
        isdir(entry) || continue
        info = _parse_case_dir(basename(entry))
        info === nothing && continue
        push!(planets, info.planet)
    end
    return [planet for planet in PLANET_ORDER if planet in planets]
end

function _metric_column(metric::String)::Symbol
    haskey(METRIC_COLUMN, metric) || throw(ArgumentError("Unsupported --metric '$metric'. Use one of: $(join(keys(METRIC_COLUMN), ", "))."))
    return METRIC_COLUMN[metric]
end

function _summary_paths(out_dir::String)
    return (
        feather=joinpath(out_dir, SUMMARY_BASENAME * ".feather"),
        csv=joinpath(out_dir, SUMMARY_BASENAME * ".csv"),
    )
end

function _trajectory_path(case_dir::String)::String
    preferred = joinpath(case_dir, "trajectory_with_active_force.feather")
    isfile(preferred) && return preferred
    single_orbit = joinpath(case_dir, "trajectory_with_active_force_orbit_001.feather")
    isfile(single_orbit) && return single_orbit
    manifest = joinpath(case_dir, "orbit_chunk_manifest.feather")
    if isfile(manifest)
        chunks = DataFrame(Arrow.Table(manifest))
        if nrow(chunks) > 0 && :path in propertynames(chunks)
            sort!(chunks, :orbit)
            first_path = String(chunks.path[1])
            isfile(first_path) && return first_path
            rel_path = joinpath(case_dir, basename(first_path))
            isfile(rel_path) && return rel_path
        end
    end
    return ""
end

function _all_case_entries(run_dir::String)
    rows = NamedTuple[]
    for entry in readdir(run_dir; join=true)
        isdir(entry) || continue
        info = _parse_case_dir(basename(entry))
        info === nothing && continue
        path = _trajectory_path(entry)
        isempty(path) && continue
        push!(rows, merge(info, (; case_dir=entry, trajectory_path=path)))
    end
    sort!(rows, by=r -> (r.planet, r.periapsis, r.apoapsis_alt_km, r.spacecraft_mass_scale, r.inclination_deg, r.argp_deg, r.dynamics, r.density))
    return rows
end

function _density_allowed(dynamics::AbstractString, density::AbstractString, density_case::AbstractString)::Bool
    density_case == "all" && return true
    dynamics in ("gram_aero", "full_environment") && return density == density_case
    return density == "none"
end

function _trace_styles_for_density(density_case::String)
    return [style for style in TRACE_STYLES if _density_allowed(style[1], style[2], density_case)]
end

function _filter_summary(df::DataFrame, opts::PlotOptions)::DataFrame
    isempty(df) && return df
    mask = map(eachrow(df)) do row
        String(row.planet) in opts.planets &&
        abs(Float64(row.spacecraft_mass_scale) - opts.spacecraft_mass_scale) <= 0.01 &&
        abs(Float64(row.inclination_deg) - opts.inclination_deg) <= 0.5 &&
        abs(Float64(row.argp_deg) - opts.argp_deg) <= 0.5 &&
        _density_allowed(String(row.dynamics_case), String(row.density_case), opts.density_case)
    end
    return df[mask, :]
end

function _has_columns(tbl, names_needed::Vector{Symbol})::Bool
    available = Set(Symbol.(propertynames(tbl)))
    return all(name -> name in available, names_needed)
end

function _central_gravity_force(mass, gm::Float64, r2)
    return mass .* gm ./ r2
end

function _force_ratio(fmag, mass, gm::Float64, r2)
    grav = _central_gravity_force(mass, gm, r2)
    ratio = fmag ./ grav
    return ifelse.((fmag .== 0.0) .| .!isfinite.(ratio), NaN, ratio)
end

function _info_for_analytical(info)
    periapsis = hasproperty(info, :periapsis) ? getproperty(info, :periapsis) : getproperty(info, :periapsis_regime)
    return (
        planet=String(info.planet),
        periapsis=String(periapsis),
        apoapsis_alt_km=Float64(info.apoapsis_alt_km),
        spacecraft_mass_scale=Float64(info.spacecraft_mass_scale),
    )
end

function _ratio_for_case(info, tbl, gm::Float64)
    dynamics = hasproperty(info, :dynamics) ? getproperty(info, :dynamics) : getproperty(info, :dynamics_case)
    if String(dynamics) == "harmonics_low"
        return simulated_model_ratio_series(_info_for_analytical(info), tbl, :harmonics, collect(eachindex(tbl.time)))
    end
    r2 = tbl.sc1_pos_1 .^ 2 .+ tbl.sc1_pos_2 .^ 2 .+ tbl.sc1_pos_3 .^ 2
    return collect(_force_ratio(tbl.active_perturbation_force_mag, tbl.sc1_mass, gm, r2))
end

function _periapsis_altitude_m(planet::String, regime::String)::Float64
    if planet == "mars"
        return Dict("shallow" => 150e3, "nominal" => 125e3, "deep" => 110e3)[regime]
    elseif planet == "venus"
        return Dict("shallow" => 180e3, "nominal" => 150e3, "deep" => 135e3)[regime]
    elseif planet == "titan"
        return Dict("shallow" => 900e3, "nominal" => 720e3, "deep" => 650e3)[regime]
    end
    return Dict("shallow" => 180e3, "nominal" => 145e3, "deep" => 120e3)[regime]
end

function _orbit_period_s(info, gm::Float64)::Float64
    analytical_info = _info_for_analytical(info)
    planet = analytical_info.planet
    rp_m = PLANET_RADIUS_M[planet] + _periapsis_altitude_m(planet, analytical_info.periapsis)
    ra_m = PLANET_RADIUS_M[planet] + analytical_info.apoapsis_alt_km * 1e3
    a_m = 0.5 * (rp_m + ra_m)
    return 2π * sqrt(a_m^3 / gm)
end

function _time_over_period(info, tbl, gm::Float64)::Vector{Float64}
    period_s = _orbit_period_s(info, gm)
    period_s > 0.0 && isfinite(period_s) || return fill(NaN, length(tbl.time))
    return collect(Float64.(tbl.time)) ./ period_s
end

function _load_ratio(info, gm::Float64)
    feather_path = String(info.trajectory_path)
    tbl = Arrow.Table(feather_path)
    needed = [:time, :sc1_pos_1, :sc1_pos_2, :sc1_pos_3, :sc1_mass, :active_perturbation_force_mag]
    if !_has_columns(tbl, needed)
        missing_names = [String(name) for name in needed if !(name in Symbol.(propertynames(tbl)))]
        error("missing column(s): $(join(missing_names, ", "))")
    end

    ratio = _ratio_for_case(info, tbl, gm)
    return _time_over_period(info, tbl, gm), ratio
end

function _finite_max(values)::Float64
    finite_values = filter(isfinite, values)
    isempty(finite_values) && return NaN
    return maximum(finite_values)
end

function _safe_quantile(values, p::Float64)::Float64
    finite_values = filter(isfinite, values)
    isempty(finite_values) && return NaN
    return quantile(finite_values, p)
end

function _first_peak_index(values, peak::Float64)::Int
    isfinite(peak) || return 0
    for i in eachindex(values)
        isfinite(values[i]) && values[i] == peak && return Int(i)
    end
    return 0
end

function _in_atmosphere_vector(tbl, n::Int)
    if :in_atmosphere in Symbol.(propertynames(tbl))
        return Bool.(collect(tbl.in_atmosphere))
    end
    return fill(false, n)
end

function _column_or_default(tbl, name::Symbol, idx::Int, default)
    name in Symbol.(propertynames(tbl)) || return default
    return getproperty(tbl, name)[idx]
end

function _summarize_case(info)::NamedTuple
    gm = PLANET_GM[info.planet]
    tbl = Arrow.Table(info.trajectory_path)
    needed = [:time, :sc1_pos_1, :sc1_pos_2, :sc1_pos_3, :sc1_mass, :active_perturbation_force_mag]
    if !_has_columns(tbl, needed)
        missing_names = [String(name) for name in needed if !(name in Symbol.(propertynames(tbl)))]
        error("missing column(s): $(join(missing_names, ", "))")
    end

    ratio = _ratio_for_case(info, tbl, gm)
    finite_count = count(isfinite, ratio)
    peak = _finite_max(ratio)
    peak_idx = _first_peak_index(ratio, peak)
    in_atmosphere = _in_atmosphere_vector(tbl, length(ratio))
    ratio_in = [ratio[i] for i in eachindex(ratio) if in_atmosphere[i]]
    ratio_out = [ratio[i] for i in eachindex(ratio) if !in_atmosphere[i]]

    return (
        planet=info.planet,
        periapsis_regime=info.periapsis,
        apoapsis_alt_km=info.apoapsis_alt_km,
        dynamics_case=info.dynamics,
        density_case=info.density,
        spacecraft_mass_scale=info.spacecraft_mass_scale,
        inclination_deg=info.inclination_deg,
        argp_deg=info.argp_deg,
        ratio_peak=peak,
        ratio_p95=_safe_quantile(ratio, 0.95),
        ratio_p50=_safe_quantile(ratio, 0.50),
        ratio_peak_in_atmosphere=_finite_max(ratio_in),
        ratio_peak_out_of_atmosphere=_finite_max(ratio_out),
        time_of_ratio_peak_s=peak_idx == 0 ? NaN : Float64(tbl.time[peak_idx]),
        altitude_at_ratio_peak_m=peak_idx == 0 ? NaN : Float64(_column_or_default(tbl, :altitude_m, peak_idx, NaN)),
        orbit_of_ratio_peak=peak_idx == 0 ? missing : _column_or_default(tbl, :orbit, peak_idx, missing),
        finite_ratio_sample_count=finite_count,
        trajectory_path=info.trajectory_path,
    )
end

function _summary_worker_count(requested_threads::Int, case_count::Int)::Int
    requested_threads > 0 || throw(ArgumentError("--summary-threads must be > 0; got $requested_threads."))
    available = Threads.nthreads()
    if requested_threads > available
        @warn "Requested summary threads exceed available Julia threads; capping" requested_threads available
    end
    return max(1, min(requested_threads, available, case_count))
end

function _plot_worker_count(requested_threads::Int, job_count::Int)::Int
    requested_threads > 0 || throw(ArgumentError("--plot-threads must be > 0; got $requested_threads."))
    job_count == 0 && return 0
    available = Threads.nthreads()
    if requested_threads > available
        @warn "Requested plot threads exceed available Julia threads; capping" requested_threads available
    end
    return max(1, min(requested_threads, available, job_count))
end

function _summarize_cases_threaded(cases, requested_threads::Int)
    worker_count = _summary_worker_count(requested_threads, length(cases))
    println("Summary threads: $worker_count / $(Threads.nthreads()) available")
    rows = Vector{Any}(nothing, length(cases))
    next_case = Threads.Atomic{Int}(1)
    completed = Threads.Atomic{Int}(0)
    progress_lock = ReentrantLock()

    tasks = map(1:worker_count) do _
        Threads.@spawn begin
            while true
                idx = Threads.atomic_add!(next_case, 1)
                idx > length(cases) && break
                info = cases[idx]
                try
                    rows[idx] = _summarize_case(info)
                catch err
                    @warn "Skipping unreadable trajectory force history" path=info.trajectory_path exception=(err, catch_backtrace())
                finally
                    done = Threads.atomic_add!(completed, 1) + 1
                    if done == 1 || done % 100 == 0 || done == length(cases)
                        lock(progress_lock) do
                            @printf("  summarized %d/%d\n", done, length(cases))
                            flush(stdout)
                        end
                    end
                end
            end
        end
    end
    foreach(fetch, tasks)
    return [row for row in rows if row !== nothing]
end

function _run_plot_jobs_threaded(jobs::Vector{<:NamedTuple}, requested_threads::Int)
    worker_count = _plot_worker_count(requested_threads, length(jobs))
    worker_count == 0 && return nothing
    println("Plot threads   : $worker_count / $(Threads.nthreads()) available for $(length(jobs)) job(s)")
    next_job = Threads.Atomic{Int}(1)
    completed = Threads.Atomic{Int}(0)
    print_lock = ReentrantLock()
    failures = Vector{Any}(nothing, length(jobs))

    tasks = map(1:worker_count) do _
        Threads.@spawn begin
            while true
                idx = Threads.atomic_add!(next_job, 1)
                idx > length(jobs) && break
                job = jobs[idx]
                try
                    lock(print_lock) do
                        println("Plotting $(job.label)...")
                        flush(stdout)
                    end
                    job.fn()
                catch err
                    failures[idx] = (label=job.label, err=err, bt=catch_backtrace())
                    @warn "Plot job failed" label=job.label exception=(err, catch_backtrace())
                finally
                    done = Threads.atomic_add!(completed, 1) + 1
                    lock(print_lock) do
                        @printf("  plotted %d/%d\n", done, length(jobs))
                        flush(stdout)
                    end
                end
            end
        end
    end
    foreach(fetch, tasks)

    failed = [failure for failure in failures if failure !== nothing]
    isempty(failed) || error("$(length(failed)) plot job(s) failed; see warnings above.")
    return nothing
end

function _summary_cache_matches_cases(df::DataFrame, cases)::Bool
    nrow(df) == length(cases) || return false
    :trajectory_path in propertynames(df) || return false
    cached_paths = Set(abspath.(String.(df.trajectory_path)))
    case_paths = Set(abspath.(String[case.trajectory_path for case in cases]))
    return cached_paths == case_paths
end

function build_force_ratio_summary(run_dir::String, out_dir::String; rebuild::Bool=false, summary_threads::Int=Threads.nthreads())::DataFrame
    paths = _summary_paths(out_dir)
    cases = _all_case_entries(run_dir)
    isempty(cases) && error("No trajectory_with_active_force*.feather files found in $run_dir")

    if !rebuild && isfile(paths.feather)
        cached_df = DataFrame(Arrow.Table(paths.feather))
        if _summary_cache_matches_cases(cached_df, cases)
            println("Using cached summary: $(paths.feather)")
            return cached_df
        end
        println("Cached summary is stale for current run ($(nrow(cached_df)) cached row(s), $(length(cases)) discovered trajectory file(s)); rebuilding.")
    end

    println("Building force-ratio summary from $(length(cases)) trajectory file(s)...")
    rows = _summarize_cases_threaded(cases, summary_threads)
    df = DataFrame(rows)
    mkpath(out_dir)
    Arrow.write(paths.feather, df)
    CSV.write(paths.csv, df)
    println("Summary saved: $(paths.feather)")
    println("Summary saved: $(paths.csv)")
    return df
end

function _present_traces(df::DataFrame, density_case::String)
    isempty(df) && return Tuple{String, String, String, Symbol}[]
    present = Set((String(row.dynamics_case), String(row.density_case)) for row in eachrow(df))
    return [style for style in _trace_styles_for_density(density_case) if (style[1], style[2]) in present]
end

function _trace_label(planet::String, dynamics::String, density::String)
    density_label = if density == "high"
        " rho high"
    elseif density == "nominal"
        " rho nominal"
    elseif density == "low"
        " rho low"
    else
        ""
    end

    base = if dynamics == "j2"
        "J2"
    elseif dynamics == "harmonics_low"
        planet == "titan" ? "Harmonics 5x5" : "Harmonics 20x20"
    elseif dynamics == "srp"
        "SRP"
    elseif dynamics == "third_body_sun"
        planet == "earth" ? "3rd Body (Sun+Moon)" : planet == "titan" ? "3rd Body (Saturn+Sun)" : "3rd Body (Sun)"
    elseif dynamics == "gram_aero"
        "Aero"
    elseif dynamics == "full_environment"
        if planet == "earth"
            "Full Env (Sun+Moon+20x20+SRP+Aero)"
        elseif planet == "titan"
            "Full Env (Saturn+Sun+5x5+SRP+Aero)"
        else
            "Full Env (Sun+20x20+SRP+Aero)"
        end
    else
        dynamics
    end
    return base * density_label
end

function _trace_label_for_rows(df::DataFrame, dynamics::String, density::String)
    planets = Set(String.(df.planet))
    density_label = if density == "high"
        " rho high"
    elseif density == "nominal"
        " rho nominal"
    elseif density == "low"
        " rho low"
    else
        ""
    end
    if dynamics == "harmonics_low" && "titan" in planets && length(planets) > 1
        return "Harmonics (20x20; Titan 5x5)" * density_label
    elseif dynamics == "third_body_sun" && "titan" in planets && length(planets) > 1
        return "3rd Body (planet-specific)" * density_label
    elseif dynamics == "full_environment" && "titan" in planets && length(planets) > 1
        return "Full Env (planet-specific)" * density_label
    end
    planet = isempty(planets) ? "earth" : first(planets)
    return _trace_label(planet, dynamics, density)
end

function _downsample(v::AbstractVector, n::Int=2000)
    length(v) <= n && return v
    step = max(1, div(length(v), n))
    return v[1:step:end]
end

function _downsample_indices(n::Int, max_count::Int=2000)
    n <= max_count && return collect(1:n)
    step = max(1, div(n, max_count))
    return collect(1:step:n)
end

function _analytical_model_key(dynamics::AbstractString)
    dynamics == "j2" && return :j2
    dynamics == "harmonics_low" && return :harmonics
    dynamics == "srp" && return :srp
    dynamics == "third_body_sun" && return :third_body
    dynamics == "gram_aero" && return :drag
    dynamics == "full_environment" && return :full_environment
    return nothing
end

function _analytical_info(row)
    return (
        planet=String(row.planet),
        periapsis=String(row.periapsis_regime),
        apoapsis_alt_km=Float64(row.apoapsis_alt_km),
        spacecraft_mass_scale=Float64(row.spacecraft_mass_scale),
    )
end

function _analytical_series_for_row(row, idxs)
    tbl = Arrow.Table(String(row.trajectory_path))
    return analytical_series(_analytical_info(row), tbl, idxs)
end

function _plot_analytical_overlay!(sp, row, x_periods, idxs, model_key, color, mode::String)
    mode == "none" && return false
    model_key === nothing && return false
    series = try
        _analytical_series_for_row(row, idxs)
    catch err
        @warn "Skipping analytical overlay" path=String(row.trajectory_path) exception=(err, catch_backtrace())
        return false
    end

    xs = x_periods[idxs]
    plotted = false
    if mode in ("basic", "both")
        ys = series.basic[model_key]
        finite = isfinite.(ys)
        if any(finite)
            plot!(sp, xs[finite], ys[finite];
                color=color, linestyle=:dashdot, linewidth=0.8, alpha=0.85, label=nothing)
            plotted = true
        end
    end
    if mode in ("detailed", "both")
        ys = series.detailed[model_key]
        finite = isfinite.(ys)
        if any(finite)
            plot!(sp, xs[finite], ys[finite];
                color=color, linestyle=:dash, linewidth=0.8, alpha=0.85, label=nothing)
            plotted = true
        end
    end
    return plotted
end

function _analytical_comparison_styles(opts::PlotOptions)
    styles = NamedTuple[]
    for (model_key, dynamics, default_density, color) in ANALYTICAL_COMPARISON_STYLES
        if model_key == :drag
            push!(styles, (model_key=model_key, dynamics=dynamics, density="nominal", color=color, label_suffix=" rho nominal"))
        else
            push!(styles, (model_key=model_key, dynamics=dynamics, density=default_density, color=color, label_suffix=""))
        end
    end
    return styles
end

function _validate_pdf(path::String)
    isfile(path) || error("PDF was not created: $path")
    filesize(path) > 0 || error("PDF is empty: $path")

    pdfinfo = Sys.which("pdfinfo")
    pdfinfo === nothing && return true

    output = try
        read(`$pdfinfo $path`, String)
    catch err
        error("PDF validation failed for $path: $err")
    end
    m = match(r"Pages:\s+(\d+)", output)
    m === nothing && error("PDF validation failed for $path: page count not found")
    parse(Int, m[1]) > 0 || error("PDF validation failed for $path: zero pages")
    return true
end

function _save_pdf(fig, out_path::String)
    mkpath(dirname(out_path))
    raw_path = tempname(dirname(out_path)) * ".raw.pdf"
    normalized_path = tempname(dirname(out_path)) * ".pdf"
    try
        lock(PLOT_SAVE_LOCK) do
            savefig(fig, raw_path)
        end
        _validate_pdf(raw_path)

        gs = Sys.which("gs")
        if gs !== nothing
            run(`$gs -q -dBATCH -dNOPAUSE -sDEVICE=pdfwrite -sOutputFile=$normalized_path $raw_path`)
            _validate_pdf(normalized_path)
            mv(normalized_path, out_path; force=true)
        else
            mv(raw_path, out_path; force=true)
            raw_path = ""
        end
        _validate_pdf(out_path)
    finally
        !isempty(raw_path) && isfile(raw_path) && rm(raw_path; force=true)
        isfile(normalized_path) && rm(normalized_path; force=true)
    end
    return out_path
end

function _ordered_periapsis(values)
    present = Set(String.(values))
    ordered = [regime for regime in PERIAPSIS_REGIME_ORDER if regime in present]
    extras = sort([regime for regime in present if !(regime in ordered)])
    return vcat(ordered, extras)
end

function _representative_values(values::AbstractVector{<:Real}, max_count::Int)
    vals = sort(unique(Float64.(values)))
    max_count > 0 || throw(ArgumentError("--time-apo-panels must be > 0; got $max_count."))
    length(vals) <= max_count && return vals
    max_count == 1 && return [first(vals)]
    idxs = unique(round.(Int, range(1, length(vals); length=max_count)))
    return vals[idxs]
end

function _values_label(values)::String
    return join(string.(round.(Int, values)), ", ")
end

function _xticks_for_dense_grid(xs::AbstractVector{<:Real}; max_ticks::Int=MAX_HEATMAP_XTICKS)
    vals = sort(unique(Float64.(xs)))
    ticks = _representative_values(vals, min(max_ticks, length(vals)))
    return (ticks, string.(round.(Int, ticks)))
end

function _metric_label(metric::String)
    return metric == "peak" ? "Peak force / gravity" :
           metric == "p95" ? "P95 force / gravity" :
           metric == "p50" ? "Median force / gravity" :
           "Peak in-atmosphere force / gravity"
end

function plot_planet_time_histories(run_dir::String, planet::String, summary_df::DataFrame, opts::PlotOptions)
    planet_df = summary_df[summary_df.planet .== planet, :]
    if isempty(planet_df)
        println("  No force-ratio summary rows for $planet with the selected filters; skipping time histories.")
        return nothing
    end

    trace_styles = _present_traces(planet_df, opts.density_case)
    isempty(trace_styles) && return nothing
    periapsis_regimes = _ordered_periapsis(planet_df.periapsis_regime)
    all_apo_alts = sort(unique(Float64.(planet_df.apoapsis_alt_km)))
    apo_alts = _representative_values(all_apo_alts, opts.time_apo_panels)
    if length(apo_alts) < length(all_apo_alts)
        println("  Dense apoapsis grid: plotting representative time-history columns for $(uppercasefirst(planet)): $(_values_label(apo_alts)) km (all $(length(all_apo_alts)) values remain in summaries/heatmaps).")
    end
    n_apo = length(apo_alts)
    n_peri = length(periapsis_regimes)
    subplots = Plots.Plot[]

    for (peri_idx, periapsis) in enumerate(periapsis_regimes)
        for (apo_idx, apo_alt_km) in enumerate(apo_alts)
            subtitle = "$(uppercasefirst(periapsis)) / $(round(Int, apo_alt_km)) km apo"
            sp = plot(
                title=subtitle,
                titlefontsize=6,
                xlabel=peri_idx == n_peri ? "Time / orbital period" : "",
                ylabel=apo_idx == 1 ? "Force / gravity" : "",
                yscale=:log10,
                ylims=(1e-9, 1.0),
                legend=false,
                framestyle=:box,
                grid=true,
                gridalpha=0.25,
                guidefontsize=7,
                tickfontsize=5,
            )

            any_data = false
            for (dynamics, density, color, lstyle) in trace_styles
                rows = planet_df[
                    (planet_df.periapsis_regime .== periapsis) .&
                    (planet_df.apoapsis_alt_km .== apo_alt_km) .&
                    (planet_df.dynamics_case .== dynamics) .&
                    (planet_df.density_case .== density),
                    :,
                ]
                isempty(rows) && continue
                x_periods, ratio = try
                    _load_ratio(rows[1, :], PLANET_GM[planet])
                catch err
                    path = String(rows.trajectory_path[1])
                    @warn "Skipping unreadable trajectory force history" path exception=(err, catch_backtrace())
                    continue
                end
                idxs = _downsample_indices(length(x_periods))
                _plot_analytical_overlay!(sp, rows[1, :], x_periods, idxs, _analytical_model_key(dynamics), color, opts.analytical_overlays)
                plot!(sp, x_periods[idxs], ratio[idxs];
                    color=color, linestyle=lstyle, linewidth=0.55, alpha=0.85, label=nothing)
                any_data = true
            end
            !any_data && annotate!(sp, 0.5, 5e-5, text("no data", :center, 6, :gray))
            push!(subplots, sp)
        end
    end

    leg = plot(
        framestyle=:none, grid=false, legend=:inside, legendfontsize=7,
        background_color_inside=:transparent, background_color_outside=:transparent,
        xlims=(0, 1), ylims=(0, 1),
    )
    for (dynamics, density, color, lstyle) in trace_styles
        label = _trace_label(planet, dynamics, density)
        plot!(leg, [NaN], [NaN]; color=color, linestyle=lstyle, linewidth=1.4, label="$(label) simulation")
        if _analytical_model_key(dynamics) !== nothing
            if opts.analytical_overlays in ("basic", "both")
                plot!(leg, [NaN], [NaN]; color=color, linestyle=:dashdot, linewidth=1.2, label="$(label) analytical basic")
            end
            if opts.analytical_overlays in ("detailed", "both")
                plot!(leg, [NaN], [NaN]; color=color, linestyle=:dash, linewidth=1.2, label="$(label) analytical detailed")
            end
        end
    end

    layout = @eval Plots.@layout [Plots.grid($n_peri, $n_apo) a{0.24w}]
    title_suffix = length(apo_alts) < length(all_apo_alts) ? " (representative apoapsis cuts)" : ""
    fig = plot(
        subplots..., leg;
        layout=layout,
        size=(max(900, 260 * n_apo), 300 * n_peri),
        plot_title="$(uppercasefirst(planet)) - Perturbation Force / Central-Body Gravity$(title_suffix)",
        plot_titlefontsize=11,
        left_margin=GRID_LEFT_MARGIN,
        bottom_margin=GRID_BOTTOM_MARGIN,
        top_margin=GRID_TOP_MARGIN,
    )

    ms_tag = replace(@sprintf("ms%.2f", opts.spacecraft_mass_scale), "." => "p")
    inc_tag = @sprintf("inc%03.0f", opts.inclination_deg)
    aop_tag = @sprintf("aop%03.0f", opts.argp_deg)
    density_tag = opts.density_case == "all" ? "all_density" : "density_$(opts.density_case)"
    out_path = joinpath(opts.out_dir, "perturbation_force_ratios_$(planet)_$(ms_tag)_$(inc_tag)_$(aop_tag)_$(density_tag).pdf")
    _save_pdf(fig, out_path)
    println("  Saved: $out_path")
    return out_path
end

function _analytical_comparison_source_rows(planet_df::DataFrame, periapsis, apo_alt_km, opts::PlotOptions)
    rows = NamedTuple[]
    for style in _analytical_comparison_styles(opts)
        local_rows = planet_df[
            (planet_df.periapsis_regime .== periapsis) .&
            (planet_df.apoapsis_alt_km .== apo_alt_km) .&
            (planet_df.dynamics_case .== style.dynamics) .&
            (planet_df.density_case .== style.density),
            :,
        ]
        isempty(local_rows) && continue
        push!(rows, (model_key=style.model_key, row=local_rows[1, :], color=style.color, label_suffix=style.label_suffix))
    end
    return rows
end

function plot_analytical_basic_vs_detailed_time_histories(planet::String, summary_df::DataFrame, opts::PlotOptions)
    opts.analytical_overlays == "none" && return nothing
    planet_df = summary_df[summary_df.planet .== planet, :]
    if isempty(planet_df)
        println("  No analytical comparison rows for $planet with the selected filters; skipping.")
        return nothing
    end

    labels = analytical_labels()
    periapsis_regimes = _ordered_periapsis(planet_df.periapsis_regime)
    all_apo_alts = sort(unique(Float64.(planet_df.apoapsis_alt_km)))
    apo_alts = _representative_values(all_apo_alts, opts.time_apo_panels)
    n_apo = length(apo_alts)
    n_peri = length(periapsis_regimes)
    subplots = Plots.Plot[]

    for (peri_idx, periapsis) in enumerate(periapsis_regimes)
        for (apo_idx, apo_alt_km) in enumerate(apo_alts)
            subtitle = "$(uppercasefirst(String(periapsis))) / $(round(Int, apo_alt_km)) km apo"
            sp = plot(
                title=subtitle,
                titlefontsize=6,
                xlabel=peri_idx == n_peri ? "Time / orbital period" : "",
                ylabel=apo_idx == 1 ? "Force / gravity or Pi" : "",
                yscale=:log10,
                ylims=(1e-9, 1.0),
                legend=false,
                framestyle=:box,
                grid=true,
                gridalpha=0.25,
                guidefontsize=7,
                tickfontsize=5,
            )

            any_data = false
            for item in _analytical_comparison_source_rows(planet_df, periapsis, apo_alt_km, opts)
                row = item.row
                tbl = Arrow.Table(String(row.trajectory_path))
                n = length(tbl.time)
                idxs = _downsample_indices(n)
                x_periods = _time_over_period(row, tbl, PLANET_GM[String(row.planet)])
                series = try
                    analytical_series(_analytical_info(row), tbl, idxs)
                catch err
                    @warn "Skipping analytical basic-vs-detailed comparison" path=String(row.trajectory_path) exception=(err, catch_backtrace())
                    continue
                end
                xs = x_periods[idxs]
                simulation = simulated_model_ratio_series(_analytical_info(row), tbl, item.model_key, idxs)
                basic = series.basic[item.model_key]
                detailed = series.detailed[item.model_key]
                finite_simulation = isfinite.(simulation)
                finite_basic = isfinite.(basic)
                finite_detailed = isfinite.(detailed)
                if any(finite_basic)
                    plot!(sp, xs[finite_basic], basic[finite_basic];
                        color=item.color, linestyle=:dashdot, linewidth=0.8, alpha=0.85, label=nothing)
                    any_data = true
                end
                if any(finite_detailed)
                    plot!(sp, xs[finite_detailed], detailed[finite_detailed];
                        color=item.color, linestyle=:dash, linewidth=1.2, alpha=0.95, label=nothing)
                    any_data = true
                end
                if any(finite_simulation)
                    plot!(sp, xs[finite_simulation], simulation[finite_simulation];
                        color=item.color, linestyle=:solid, linewidth=0.55, alpha=0.85, label=nothing)
                    any_data = true
                end
            end
            !any_data && annotate!(sp, 0.5, 5e-5, text("no data", :center, 6, :gray))
            push!(subplots, sp)
        end
    end

    leg = plot(
        framestyle=:none, grid=false, legend=:inside, legendfontsize=6,
        background_color_inside=:transparent, background_color_outside=:transparent,
        xlims=(0, 1), ylims=(0, 1),
    )
    for style in _analytical_comparison_styles(opts)
        model_label = labels[style.model_key] * style.label_suffix
        plot!(leg, [NaN], [NaN]; color=style.color, linestyle=:solid, linewidth=0.8, label="Simulated $(model_label)")
        plot!(leg, [NaN], [NaN]; color=style.color, linestyle=:dashdot, linewidth=1.2, label="Basic $(model_label)")
        plot!(leg, [NaN], [NaN]; color=style.color, linestyle=:dash, linewidth=1.6, label="Detailed $(model_label)")
    end

    layout = @eval Plots.@layout [Plots.grid($n_peri, $n_apo) a{0.30w}]
    title_suffix = length(apo_alts) < length(all_apo_alts) ? " (representative apoapsis cuts)" : ""
    fig = plot(
        subplots..., leg;
        layout=layout,
        size=(max(900, 260 * n_apo), 300 * n_peri),
        plot_title="$(uppercasefirst(planet)) - Simulated vs Analytical Perturbation Parameters$(title_suffix)",
        plot_titlefontsize=11,
        left_margin=GRID_LEFT_MARGIN,
        bottom_margin=GRID_BOTTOM_MARGIN,
        top_margin=GRID_TOP_MARGIN,
    )

    ms_tag = replace(@sprintf("ms%.2f", opts.spacecraft_mass_scale), "." => "p")
    inc_tag = @sprintf("inc%03.0f", opts.inclination_deg)
    aop_tag = @sprintf("aop%03.0f", opts.argp_deg)
    density_tag = "density_nominal"
    out_path = joinpath(opts.out_dir, "analytical_basic_vs_detailed_$(planet)_$(ms_tag)_$(inc_tag)_$(aop_tag)_$(density_tag).pdf")
    _save_pdf(fig, out_path)
    println("  Saved: $out_path")
    return out_path
end

function _heatmap_matrix(df::DataFrame, metric_col::Symbol, apo_alts, periapsis_regimes)
    z = fill(NaN, length(periapsis_regimes), length(apo_alts))
    for (iy, periapsis) in enumerate(periapsis_regimes)
        for (ix, apo) in enumerate(apo_alts)
            rows = df[(df.periapsis_regime .== periapsis) .& (df.apoapsis_alt_km .== apo), :]
            isempty(rows) && continue
            vals = filter(isfinite, Float64.(rows[!, metric_col]))
            isempty(vals) && continue
            z[iy, ix] = log10(maximum(vals))
        end
    end
    return z
end

function plot_regime_heatmaps(summary_df::DataFrame, opts::PlotOptions)
    metric_col = _metric_column(opts.metric)
    trace_styles = _present_traces(summary_df, opts.density_case)
    paths = String[]
    for (dynamics, density, _, _) in trace_styles
        local_df = summary_df[(summary_df.dynamics_case .== dynamics) .& (summary_df.density_case .== density), :]
        isempty(local_df) && continue
        subplots = Plots.Plot[]
        for planet in opts.planets
            planet_df = local_df[local_df.planet .== planet, :]
            isempty(planet_df) && continue
            apo_alts = sort(unique(Float64.(planet_df.apoapsis_alt_km)))
            periapsis_regimes = _ordered_periapsis(planet_df.periapsis_regime)
            z = _heatmap_matrix(planet_df, metric_col, apo_alts, periapsis_regimes)
            xticks = _xticks_for_dense_grid(apo_alts)
            sp = heatmap(
                apo_alts,
                periapsis_regimes,
                z;
                title=uppercasefirst(planet),
                titlefontsize=8,
                xlabel="Apoapsis altitude (km)",
                ylabel="Regime",
                color=:viridis,
                clims=(-9, 0),
                colorbar_title="log10 ratio",
                framestyle=:box,
                guidefontsize=7,
                tickfontsize=6,
                xticks=xticks,
                xrotation=35,
            )
            push!(subplots, sp)
        end
        isempty(subplots) && continue
        fig = plot(
            subplots...;
            layout=(length(subplots), 1),
            size=(1000, max(320, 280 * length(subplots))),
            plot_title="$(_trace_label_for_rows(local_df, dynamics, density)) - $(_metric_label(opts.metric))",
            plot_titlefontsize=11,
            left_margin=12Plots.mm,
            bottom_margin=16Plots.mm,
        )
        density_tag = density == "none" ? "none" : density
        out_path = joinpath(opts.out_dir, "perturbation_force_ratio_heatmap_$(opts.metric)_$(dynamics)_$(density_tag).pdf")
        _save_pdf(fig, out_path)
        println("  Saved: $out_path")
        push!(paths, out_path)
    end
    return paths
end

function plot_cross_planet_comparison(summary_df::DataFrame, opts::PlotOptions)
    metric_col = _metric_column(opts.metric)
    trace_styles = _present_traces(summary_df, opts.density_case)
    periapsis_regimes = _ordered_periapsis(summary_df.periapsis_regime)
    apo_xticks = _xticks_for_dense_grid(summary_df.apoapsis_alt_km)
    planet_colors = Dict("mars" => "#b4492c", "earth" => "#2b6cb0", "venus" => "#c18f25", "titan" => "#6a4c93")
    subplots = Plots.Plot[]

    for (idx, periapsis) in enumerate(periapsis_regimes)
        sp = plot(
            title=uppercasefirst(periapsis),
            titlefontsize=8,
            xlabel="Apoapsis altitude (km)",
            ylabel=idx == 1 ? _metric_label(opts.metric) : "",
            yscale=:log10,
            ylims=(1e-9, 1.0),
            legend=idx == 1 ? :outertopright : false,
            framestyle=:box,
            grid=true,
            gridalpha=0.25,
            guidefontsize=8,
            tickfontsize=7,
            xticks=apo_xticks,
            xrotation=35,
        )
        for planet in opts.planets
            planet_df = summary_df[(summary_df.planet .== planet) .& (summary_df.periapsis_regime .== periapsis), :]
            isempty(planet_df) && continue
            for (dynamics, density, _, lstyle) in trace_styles
                local_df = planet_df[(planet_df.dynamics_case .== dynamics) .& (planet_df.density_case .== density), :]
                isempty(local_df) && continue
                sort!(local_df, :apoapsis_alt_km)
                xs = Float64.(local_df.apoapsis_alt_km)
                ys = Float64.(local_df[!, metric_col])
                finite = isfinite.(ys)
                any(finite) || continue
                label = "$(uppercasefirst(planet)) $(_trace_label(planet, dynamics, density))"
                plot!(sp, xs[finite], ys[finite];
                    color=get(planet_colors, planet, "#333333"),
                    linestyle=lstyle,
                    marker=:circle,
                    markersize=2.5,
                    linewidth=1.0,
                    alpha=0.85,
                    label=label)
            end
        end
        push!(subplots, sp)
    end

    isempty(subplots) && return nothing
    fig = plot(
        subplots...;
        layout=(length(subplots), 1),
        size=(1300, max(420, 320 * length(subplots))),
        plot_title="Cross-Planet Perturbation Force Ratio - $(_metric_label(opts.metric))",
        plot_titlefontsize=12,
        left_margin=14Plots.mm,
        bottom_margin=10Plots.mm,
    )
    out_path = joinpath(opts.out_dir, "perturbation_force_ratio_cross_planet_$(opts.metric).pdf")
    _save_pdf(fig, out_path)
    println("  Saved: $out_path")
    return out_path
end

function write_rankings(summary_df::DataFrame, opts::PlotOptions; top_n::Int=10)
    metric_col = _metric_column(opts.metric)
    rows = NamedTuple[]
    grouped = groupby(summary_df, [:planet, :periapsis_regime])
    for group in grouped
        local_df = group[isfinite.(Float64.(group[!, metric_col])), :]
        isempty(local_df) && continue
        sort!(local_df, metric_col; rev=true)
        for (rank, row) in enumerate(eachrow(first(local_df, min(top_n, nrow(local_df)))))
            push!(rows, (
                planet=String(row.planet),
                periapsis_regime=String(row.periapsis_regime),
                rank=rank,
                metric=opts.metric,
                metric_value=Float64(row[metric_col]),
                apoapsis_alt_km=Float64(row.apoapsis_alt_km),
                dynamics_case=String(row.dynamics_case),
                density_case=String(row.density_case),
                spacecraft_mass_scale=Float64(row.spacecraft_mass_scale),
                inclination_deg=Float64(row.inclination_deg),
                argp_deg=Float64(row.argp_deg),
                time_of_ratio_peak_s=Float64(row.time_of_ratio_peak_s),
                altitude_at_ratio_peak_m=Float64(row.altitude_at_ratio_peak_m),
            ))
        end
    end
    ranking_df = DataFrame(rows)
    out_csv = joinpath(opts.out_dir, "perturbation_force_ratio_rankings_$(opts.metric).csv")
    CSV.write(out_csv, ranking_df)
    println("  Saved: $out_csv")
    return ranking_df, out_csv
end

function plot_rankings(ranking_df::DataFrame, opts::PlotOptions)
    isempty(ranking_df) && return nothing
    label(row) = "$(row.planet) $(row.periapsis_regime) / $(round(Int, row.apoapsis_alt_km)) km / $(row.dynamics_case) $(row.density_case)"
    top = first(sort(ranking_df, :metric_value; rev=true), min(25, nrow(ranking_df)))
    labels = label.(eachrow(top))
    values = Float64.(top.metric_value)
    fig = bar(
        reverse(labels),
        reverse(values);
        orientation=:h,
        xscale=:log10,
        xlims=(1e-9, 1.0),
        legend=false,
        xlabel=_metric_label(opts.metric),
        ylabel="",
        title="Top Perturbation Force Ratios",
        titlefontsize=12,
        guidefontsize=8,
        tickfontsize=6,
        size=(1300, 900),
        left_margin=58Plots.mm,
        bottom_margin=12Plots.mm,
        color="#4c78a8",
    )
    out_path = joinpath(opts.out_dir, "perturbation_force_ratio_rankings_$(opts.metric).pdf")
    _save_pdf(fig, out_path)
    println("  Saved: $out_path")
    return out_path
end

function plot_orbit_parameter_effects(summary_df::DataFrame, planet::String, parameter::Symbol, opts::PlotOptions)
    parameter in (:inclination, :argp) || throw(ArgumentError("Unsupported sweep parameter: $parameter"))
    planet_df = summary_df[
        (summary_df.planet .== planet) .&
        (abs.(summary_df.spacecraft_mass_scale .- opts.spacecraft_mass_scale) .<= 0.01) .&
        map((d, rho) -> _density_allowed(String(d), String(rho), opts.density_case), summary_df.dynamics_case, summary_df.density_case),
        :,
    ]
    parameter_label = parameter == :inclination ? "Inclination" : "Argument of Periapsis"
    parameter_tag = parameter == :inclination ? "inclination" : "argp"
    fixed_label = parameter == :inclination ? @sprintf("argp=%.0f deg", opts.argp_deg) : @sprintf("inc=%.0f deg", opts.inclination_deg)

    if isempty(planet_df)
        println("  No $parameter_label sweep rows for $planet with the selected filters; skipping.")
        return nothing
    end

    trace_styles = _present_traces(planet_df, opts.density_case)
    isempty(trace_styles) && return nothing
    all_apo_alts = sort(unique(Float64.(planet_df.apoapsis_alt_km)))
    apo_alts = _representative_values(all_apo_alts, opts.time_apo_panels)
    if length(apo_alts) < length(all_apo_alts)
        println("  Dense apoapsis grid: plotting representative $parameter_label sweep columns for $(uppercasefirst(planet)): $(_values_label(apo_alts)) km.")
    end
    periapsis_regimes = _ordered_periapsis(planet_df.periapsis_regime)
    parameter_values = sort(unique(Float64.(parameter == :inclination ? planet_df.inclination_deg : planet_df.argp_deg)))
    n_apo = length(apo_alts)
    n_peri = length(periapsis_regimes)
    subplots = Plots.Plot[]

    for (peri_idx, periapsis) in enumerate(periapsis_regimes)
        for (apo_idx, apo_alt_km) in enumerate(apo_alts)
            subtitle = "$(uppercasefirst(periapsis)) / $(round(Int, apo_alt_km)) km apo"
            sp = plot(
                title=subtitle,
                titlefontsize=6,
                xlabel=peri_idx == n_peri ? "$parameter_label (deg)" : "",
                ylabel=apo_idx == 1 ? "Peak force / gravity" : "",
                yscale=:log10,
                ylims=(1e-9, 1.0),
                legend=false,
                framestyle=:box,
                grid=true,
                gridalpha=0.25,
                guidefontsize=7,
                tickfontsize=5,
            )

            any_data = false
            for (dynamics, density, color, lstyle) in trace_styles
                xs = Float64[]
                ys = Float64[]
                for value in parameter_values
                    mask = if parameter == :inclination
                        (planet_df.periapsis_regime .== periapsis) .&
                        (planet_df.apoapsis_alt_km .== apo_alt_km) .&
                        (planet_df.dynamics_case .== dynamics) .&
                        (planet_df.density_case .== density) .&
                        (planet_df.inclination_deg .== value) .&
                        (abs.(planet_df.argp_deg .- opts.argp_deg) .<= 0.5)
                    else
                        (planet_df.periapsis_regime .== periapsis) .&
                        (planet_df.apoapsis_alt_km .== apo_alt_km) .&
                        (planet_df.dynamics_case .== dynamics) .&
                        (planet_df.density_case .== density) .&
                        (planet_df.argp_deg .== value) .&
                        (abs.(planet_df.inclination_deg .- opts.inclination_deg) .<= 0.5)
                    end
                    rows = planet_df[mask, :]
                    isempty(rows) && continue
                    peak_ratio = Float64(rows.ratio_peak[1])
                    isfinite(peak_ratio) || continue
                    push!(xs, value)
                    push!(ys, peak_ratio)
                end
                isempty(xs) && continue
                plot!(sp, xs, ys;
                    color=color, linestyle=lstyle, marker=:circle, markersize=2.2,
                    linewidth=1.0, alpha=0.9, label=nothing)
                any_data = true
            end
            !any_data && annotate!(sp, 0.5, 5e-5, text("no data", :center, 6, :gray))
            push!(subplots, sp)
        end
    end

    leg = plot(
        framestyle=:none, grid=false, legend=:inside, legendfontsize=7,
        background_color_inside=:transparent, background_color_outside=:transparent,
        xlims=(0, 1), ylims=(0, 1),
    )
    for (dynamics, density, color, lstyle) in trace_styles
        plot!(leg, [NaN], [NaN]; color=color, linestyle=lstyle, marker=:circle,
            markersize=2.2, linewidth=1.4, label=_trace_label(planet, dynamics, density))
    end

    layout = @eval Plots.@layout [Plots.grid($n_peri, $n_apo) a{0.14w}]
    title_suffix = length(apo_alts) < length(all_apo_alts) ? " (representative apoapsis cuts)" : ""
    fig = plot(
        subplots..., leg;
        layout=layout,
        size=(max(900, 260 * n_apo), 300 * n_peri),
        plot_title="$(uppercasefirst(planet)) - Peak Perturbation Force / Gravity vs $parameter_label ($fixed_label)$(title_suffix)",
        plot_titlefontsize=11,
        left_margin=GRID_LEFT_MARGIN,
        bottom_margin=GRID_BOTTOM_MARGIN,
        top_margin=GRID_TOP_MARGIN,
    )

    ms_tag = replace(@sprintf("ms%.2f", opts.spacecraft_mass_scale), "." => "p")
    fixed_tag = parameter == :inclination ? @sprintf("aop%03.0f", opts.argp_deg) : @sprintf("inc%03.0f", opts.inclination_deg)
    density_tag = opts.density_case == "all" ? "all_density" : "density_$(opts.density_case)"
    out_path = joinpath(opts.out_dir, "perturbation_force_ratios_vs_$(parameter_tag)_$(planet)_$(ms_tag)_$(fixed_tag)_$(density_tag).pdf")
    _save_pdf(fig, out_path)
    println("  Saved: $out_path")
    return out_path
end

function _most_recent_run(base::String)
    isdir(base) || error("Output directory not found: $base")
    entries = filter(
        d -> isdir(joinpath(base, d)) && occursin(r"^\d{8}_\d{6}$", d),
        readdir(base),
    )
    isempty(entries) && error("No timestamped run directories found in $base")
    return joinpath(base, last(sort(entries)))
end

function _parse_args(args)::PlotOptions
    args = collect(args)
    planets = String[]
    out_dir = nothing
    positional = String[]
    spacecraft_mass_scale = 1.0
    inclination_deg = 93.0
    argp_deg = 80.0
    rebuild_summary = false
    metric = "peak"
    plot_set = "all"
    density_case = "nominal"
    time_apo_panels = DEFAULT_TIME_APO_PANELS
    summary_threads = parse(Int, get(ENV, "SPACEAGORA_AERO_PERTURB_SUMMARY_THREADS", string(Threads.nthreads())))
    plot_threads = parse(Int, get(ENV, PLOT_THREADS_ENV, string(Threads.nthreads())))
    analytical_overlays = "both"

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--planets"
            i += 1; planets = split(args[i], ",") .|> strip .|> String
        elseif a == "--out-dir"
            i += 1; out_dir = abspath(args[i])
        elseif a == "--mass-scale"
            i += 1; spacecraft_mass_scale = parse(Float64, args[i])
        elseif a == "--inclination"
            i += 1; inclination_deg = parse(Float64, args[i])
        elseif a == "--argp"
            i += 1; argp_deg = parse(Float64, args[i])
        elseif a == "--rebuild-summary"
            rebuild_summary = true
        elseif a == "--metric"
            i += 1; metric = String(strip(args[i]))
        elseif a == "--plot-set"
            i += 1; plot_set = String(strip(args[i]))
        elseif a == "--density-case"
            i += 1; density_case = String(strip(args[i]))
        elseif a == "--time-apo-panels"
            i += 1; time_apo_panels = parse(Int, args[i])
        elseif a == "--summary-threads"
            i += 1; summary_threads = parse(Int, args[i])
        elseif a == "--plot-threads"
            i += 1; plot_threads = parse(Int, args[i])
        elseif a == "--analytical-overlays"
            i += 1; analytical_overlays = String(strip(args[i]))
        elseif !startswith(a, "-")
            push!(positional, a)
        else
            throw(ArgumentError("Unknown option: $a"))
        end
        i += 1
    end

    inferred_latest_run = isempty(positional)
    run_dir = inferred_latest_run ? _most_recent_run(DEFAULT_OUTPUT_BASE) : abspath(positional[1])
    isempty(planets) && (planets = _discover_planets(run_dir))
    isempty(planets) && (planets = copy(PLANET_ORDER))
    out_dir === nothing && (out_dir = joinpath(run_dir, "plots"))
    _metric_column(metric)
    plot_set in ("all", "time", "heatmaps", "rankings") || throw(ArgumentError("Unsupported --plot-set '$plot_set'."))
    density_case in ("nominal", "low", "high", "all") || throw(ArgumentError("Unsupported --density-case '$density_case'."))
    analytical_overlays in ANALYTICAL_OVERLAY_MODES || throw(ArgumentError("Unsupported --analytical-overlays '$analytical_overlays'. Use one of: $(join(ANALYTICAL_OVERLAY_MODES, ", "))."))
    time_apo_panels > 0 || throw(ArgumentError("--time-apo-panels must be > 0; got $time_apo_panels."))
    summary_threads > 0 || throw(ArgumentError("--summary-threads must be > 0; got $summary_threads."))
    plot_threads > 0 || throw(ArgumentError("--plot-threads must be > 0; got $plot_threads."))
    unknown_planets = setdiff(planets, PLANET_ORDER)
    isempty(unknown_planets) || throw(ArgumentError("Unsupported planet(s): $(join(unknown_planets, ", "))."))

    return PlotOptions(
        run_dir=run_dir,
        planets=planets,
        out_dir=out_dir,
        inferred_latest_run=inferred_latest_run,
        spacecraft_mass_scale=spacecraft_mass_scale,
        inclination_deg=inclination_deg,
        argp_deg=argp_deg,
        rebuild_summary=rebuild_summary,
        metric=metric,
        plot_set=plot_set,
        density_case=density_case,
        time_apo_panels=time_apo_panels,
        summary_threads=summary_threads,
        plot_threads=plot_threads,
        analytical_overlays=analytical_overlays,
    )
end

function main(args=ARGS)
    opts = _parse_args(args)
    isdir(opts.run_dir) || error("Not a directory: $(opts.run_dir)")
    mkpath(opts.out_dir)
    println("Run directory : $(opts.run_dir)")
    opts.inferred_latest_run && println("Dataset       : most recent timestamped run under $(DEFAULT_OUTPUT_BASE)")
    println("Output        : $(opts.out_dir)")
    println("Planets       : $(join(opts.planets, ", "))")
    println("Filter        : mass_scale=$(opts.spacecraft_mass_scale)  inclination=$(opts.inclination_deg) deg  argp=$(opts.argp_deg) deg  density_case=$(opts.density_case)")
    println("Metric        : $(opts.metric)")
    println("Plot set      : $(opts.plot_set)")
    println("Time apo cuts : $(opts.time_apo_panels)")
    println("Summary tasks : $(opts.summary_threads)")
    println("Plot tasks    : $(opts.plot_threads)")
    println("Analytical    : $(opts.analytical_overlays)")

    summary_df_all = build_force_ratio_summary(opts.run_dir, opts.out_dir; rebuild=opts.rebuild_summary, summary_threads=opts.summary_threads)
    summary_df = _filter_summary(summary_df_all, opts)
    filtered_csv = joinpath(opts.out_dir, "perturbation_force_ratio_summary_filtered.csv")
    CSV.write(filtered_csv, summary_df)
    println("Filtered summary saved: $filtered_csv ($(nrow(summary_df)) row(s))")

    if opts.plot_set in ("all", "time")
        jobs = NamedTuple[]
        for planet in opts.planets
            planet_all = summary_df_all[summary_df_all.planet .== planet, :]
            push!(jobs, (label="time histories for $planet", fn=let planet=planet
                () -> plot_planet_time_histories(opts.run_dir, planet, summary_df, opts)
            end))
            push!(jobs, (label="analytical comparison for $planet", fn=let planet=planet
                () -> plot_analytical_basic_vs_detailed_time_histories(planet, summary_df, opts)
            end))
            push!(jobs, (label="inclination sweep for $planet", fn=let planet=planet, planet_all=copy(planet_all)
                () -> plot_orbit_parameter_effects(planet_all, planet, :inclination, opts)
            end))
            push!(jobs, (label="argument-of-periapsis sweep for $planet", fn=let planet=planet, planet_all=copy(planet_all)
                () -> plot_orbit_parameter_effects(planet_all, planet, :argp, opts)
            end))
        end
        _run_plot_jobs_threaded(jobs, opts.plot_threads)
    end

    if opts.plot_set in ("all", "heatmaps")
        jobs = NamedTuple[
            (label="regime heatmaps", fn=() -> plot_regime_heatmaps(summary_df, opts)),
            (label="cross-planet comparison", fn=() -> plot_cross_planet_comparison(summary_df, opts)),
        ]
        _run_plot_jobs_threaded(jobs, opts.plot_threads)
    end

    if opts.plot_set in ("all", "rankings")
        println("Writing and plotting rankings...")
        flush(stdout)
        ranking_df, _ = write_rankings(summary_df, opts)
        plot_rankings(ranking_df, opts)
    end

    println("Done.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    _ensure_multithreaded_cli!()
    main()
end
