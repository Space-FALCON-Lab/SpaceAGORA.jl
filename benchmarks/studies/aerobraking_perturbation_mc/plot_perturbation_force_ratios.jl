# For each available planet, plots a grid of perturbation force magnitude / central-body gravity.
# Rows = periapsis regimes found in the run, columns = apoapsis altitudes found in the run.
# Each subplot traces the non-two_body dynamics/density cases found for that
# (planet, periapsis, apoapsis) slot.
# Reads trajectory_with_active_force.feather files written by study.jl.
#
# Usage: julia --project=. benchmarks/studies/aerobraking_perturbation_mc/plot_perturbation_force_ratios.jl <run_dir>
#        Omit <run_dir> to use the most recent run in output/aerobraking_perturbation_mc/.
# Options:
#   --planets mars,earth,venus   (default: planets discovered in the run)
#   --out-dir PATH               directory for output PDFs (default: <run_dir>/plots/)
#   --mass-scale SCALE           spacecraft mass scale for plotted slices (default: 1.0)
#   --inclination DEG            inclination held fixed for baseline/argp plots (default: 93)
#   --argp DEG                   argument of periapsis held fixed for baseline/inclination plots (default: 80)

using Arrow
using DataFrames
using Plots
using Printf

gr()

const REPO_ROOT = abspath(joinpath(@__DIR__, "..", "..", ".."))
const DEFAULT_OUTPUT_BASE = joinpath(REPO_ROOT, "output", "aerobraking_perturbation_mc")

# Gravitational parameters (m³/s²) matching planets.jl
const PLANET_GM = Dict(
    "mars"  => 4.282837285418775e13,
    "earth" => 3.98600436233e14,
    "venus" => 3.24858599e14,
)

const PERIAPSIS_REGIME_ORDER = ["shallow", "nominal", "deep"]
const PLANET_ORDER = ["mars", "earth", "venus"]
const GRID_LEFT_MARGIN = 12Plots.mm
const GRID_BOTTOM_MARGIN = 12Plots.mm
const GRID_TOP_MARGIN = 5Plots.mm

# (dynamics_case, density_case, line_color, line_style)
const TRACE_STYLES = [
    ("j2",               "none",    "#1f77b4", :solid),
    ("harmonics_low",    "none",    "#2ca02c", :solid),
    ("srp",              "none",    "#ff7f0e", :solid),
    ("third_body_sun",   "none",    "#d62728", :solid),
    ("gram_aero",        "high",    "#7b2d8b", :solid),
    ("gram_aero",        "nominal", "#9e4fc4", :dash),
    ("gram_aero",        "low",     "#c99be0", :dot),
    ("full_environment", "high",    "#8c3a00", :solid),
    ("full_environment", "nominal", "#bf5e1a", :dash),
    ("full_environment", "low",     "#e0a070", :dot),
]

function _parse_case_dir(dirname::String)
    m = match(
        r"^case_\d+_(mars|earth|venus)_(shallow|nominal|deep)_apo_(\d+)km_(.+)_density_(none|low|nominal|high)_ms(\d+p\d+)_inc(\d+)_aop(\d+)$",
        dirname,
    )
    m === nothing && return nothing
    return (;
        planet=m[1], periapsis=m[2], apoapsis_alt_km=parse(Float64, m[3]),
        dynamics=m[4], density=m[5],
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

function _build_index(run_dir::String, planet::String;
        spacecraft_mass_scale::Float64=1.0,
        inclination_deg::Float64=93.0,
        argp_deg::Float64=80.0)
    index = Dict{Tuple{String, Float64, String, String}, String}()
    apo_alts = Set{Float64}()
    periapsis_regimes = Set{String}()
    for entry in readdir(run_dir; join=true)
        isdir(entry) || continue
        info = _parse_case_dir(basename(entry))
        info === nothing && continue
        info.planet != planet && continue
        abs(info.spacecraft_mass_scale - spacecraft_mass_scale) > 0.01 && continue
        abs(info.inclination_deg - inclination_deg) > 0.5 && continue
        abs(info.argp_deg - argp_deg) > 0.5 && continue
        p = joinpath(entry, "trajectory_with_active_force.feather")
        isfile(p) || continue
        index[(info.periapsis, info.apoapsis_alt_km, info.dynamics, info.density)] = p
        push!(apo_alts, info.apoapsis_alt_km)
        push!(periapsis_regimes, info.periapsis)
    end
    ordered_periapsis_regimes = [regime for regime in PERIAPSIS_REGIME_ORDER if regime in periapsis_regimes]
    return index, ordered_periapsis_regimes, sort!(collect(apo_alts))
end

function _build_parameter_sweep_index(run_dir::String, planet::String, parameter::Symbol;
        spacecraft_mass_scale::Float64=1.0,
        inclination_deg::Float64=93.0,
        argp_deg::Float64=80.0)
    parameter in (:inclination, :argp) || throw(ArgumentError("Unsupported sweep parameter: $parameter"))

    index = Dict{Tuple{String, Float64, String, String, Float64}, String}()
    apo_alts = Set{Float64}()
    periapsis_regimes = Set{String}()
    parameter_values = Set{Float64}()

    for entry in readdir(run_dir; join=true)
        isdir(entry) || continue
        info = _parse_case_dir(basename(entry))
        info === nothing && continue
        info.planet != planet && continue
        abs(info.spacecraft_mass_scale - spacecraft_mass_scale) > 0.01 && continue

        value = if parameter == :inclination
            abs(info.argp_deg - argp_deg) > 0.5 && continue
            info.inclination_deg
        else
            abs(info.inclination_deg - inclination_deg) > 0.5 && continue
            info.argp_deg
        end

        p = joinpath(entry, "trajectory_with_active_force.feather")
        isfile(p) || continue
        index[(info.periapsis, info.apoapsis_alt_km, info.dynamics, info.density, value)] = p
        push!(apo_alts, info.apoapsis_alt_km)
        push!(periapsis_regimes, info.periapsis)
        push!(parameter_values, value)
    end

    ordered_periapsis_regimes = [regime for regime in PERIAPSIS_REGIME_ORDER if regime in periapsis_regimes]
    return index, ordered_periapsis_regimes, sort!(collect(apo_alts)), sort!(collect(parameter_values))
end

function _present_traces(index::Dict{Tuple{String, Float64, String, String}, String})
    present = Set((dynamics, density) for (_, _, dynamics, density) in keys(index))
    return [style for style in TRACE_STYLES if (style[1], style[2]) in present]
end

function _present_traces(index::Dict{Tuple{String, Float64, String, String, Float64}, String})
    present = Set((dynamics, density) for (_, _, dynamics, density, _) in keys(index))
    return [style for style in TRACE_STYLES if (style[1], style[2]) in present]
end

function _trace_label(planet::String, dynamics::String, density::String)
    density_label = if density == "high"
        " rho +25%"
    elseif density == "nominal"
        " rho nominal"
    elseif density == "low"
        " rho -25%"
    else
        ""
    end

    base = if dynamics == "j2"
        "J2"
    elseif dynamics == "harmonics_low"
        "Harmonics 20x20"
    elseif dynamics == "srp"
        "SRP"
    elseif dynamics == "third_body_sun"
        planet == "earth" ? "3rd Body (Sun+Moon)" : "3rd Body (Sun)"
    elseif dynamics == "gram_aero"
        "Aero"
    elseif dynamics == "full_environment"
        planet == "earth" ? "Full Env (Sun+Moon+20x20+SRP+Aero)" : "Full Env (Sun+20x20+SRP+Aero)"
    else
        dynamics
    end
    return base * density_label
end

function _has_columns(tbl, names_needed::Vector{Symbol})::Bool
    available = Set(Symbol.(propertynames(tbl)))
    return all(name -> name in available, names_needed)
end

function _load_ratio(feather_path::String, gm::Float64)
    tbl = Arrow.Table(feather_path)
    needed = [:time, :sc1_pos_1, :sc1_pos_2, :sc1_pos_3, :sc1_mass, :active_perturbation_force_mag]
    if !_has_columns(tbl, needed)
        missing_names = [String(name) for name in needed if !(name in Symbol.(propertynames(tbl)))]
        error("missing column(s): $(join(missing_names, ", "))")
    end
    pos1 = tbl.sc1_pos_1
    pos2 = tbl.sc1_pos_2
    pos3 = tbl.sc1_pos_3
    mass = tbl.sc1_mass
    fmag = tbl.active_perturbation_force_mag
    t_hr = collect(tbl.time) ./ 3600.0

    r2 = pos1 .^ 2 .+ pos2 .^ 2 .+ pos3 .^ 2
    grav = mass .* gm ./ r2
    ratio = fmag ./ grav
    # zero-force timesteps (outside atmosphere) → NaN so log-scale segments aren't drawn
    ratio = ifelse.(fmag .== 0.0, NaN, ratio)
    return t_hr, collect(ratio)
end

function _load_peak_ratio(feather_path::String, gm::Float64)
    _, ratio = _load_ratio(feather_path, gm)
    finite_ratio = filter(isfinite, ratio)
    isempty(finite_ratio) && return NaN
    return maximum(finite_ratio)
end

function _downsample(v::AbstractVector, n::Int=2000)
    length(v) <= n && return v
    step = max(1, div(length(v), n))
    return v[1:step:end]
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
        savefig(fig, raw_path)
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

function plot_planet(run_dir::String, planet::String; output_dir::String=run_dir,
        spacecraft_mass_scale::Float64=1.0, inclination_deg::Float64=93.0, argp_deg::Float64=80.0)
    gm = PLANET_GM[planet]
    index, periapsis_regimes, apo_alts = _build_index(run_dir, planet;
        spacecraft_mass_scale=spacecraft_mass_scale,
        inclination_deg=inclination_deg,
        argp_deg=argp_deg)

    if isempty(index)
        println("  No trajectory_with_active_force.feather files for $planet with the selected filters; skipping.")
        return nothing
    end

    trace_styles = _present_traces(index)
    if isempty(trace_styles)
        println("  No plotted perturbation cases for $planet with the selected filters; skipping.")
        return nothing
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
                xlabel=peri_idx == n_peri ? "Time (hr)" : "",
                ylabel=apo_idx == 1 ? "Force / gravity" : "",
                yscale=:log10,
                legend=false,
                ylims=(1e-9, 1.0),
                framestyle=:box,
                grid=true,
                gridalpha=0.25,
                guidefontsize=7,
                tickfontsize=5,
            )

            any_data = false
            for (dynamics, density, color, lstyle) in trace_styles
                path = get(index, (periapsis, apo_alt_km, dynamics, density), nothing)
                path === nothing && continue
                t_hr, ratio = try
                    _load_ratio(path, gm)
                catch err
                    @warn "Skipping unreadable trajectory force history" path exception=(err, catch_backtrace())
                    continue
                end
                plot!(sp, _downsample(t_hr), _downsample(ratio);
                    color=color, linestyle=lstyle, linewidth=0.8, alpha=0.9, label=nothing)
                any_data = true
            end
            !any_data && annotate!(sp, 0.5, 5e-5, text("no data", :center, 6, :gray))

            push!(subplots, sp)
        end
    end

    # Shared legend panel
    leg = plot(
        framestyle=:none, grid=false, legend=:inside, legendfontsize=7,
        background_color_inside=:transparent, background_color_outside=:transparent,
        xlims=(0, 1), ylims=(0, 1),
    )
    for (dynamics, density, color, lstyle) in trace_styles
        label = _trace_label(planet, dynamics, density)
        plot!(leg, [NaN], [NaN]; color=color, linestyle=lstyle, linewidth=1.4, label=label)
    end

    layout = @eval Plots.@layout [Plots.grid($n_peri, $n_apo) a{0.14w}]
    fig = plot(
        subplots..., leg;
        layout=layout,
        size=(max(900, 260 * n_apo), 300 * n_peri),
        plot_title="$(uppercasefirst(planet))  —  Perturbation Force / Central-Body Gravity",
        plot_titlefontsize=11,
        left_margin=GRID_LEFT_MARGIN,
        bottom_margin=GRID_BOTTOM_MARGIN,
        top_margin=GRID_TOP_MARGIN,
    )

    ms_tag = replace(@sprintf("ms%.2f", spacecraft_mass_scale), "." => "p")
    inc_tag = @sprintf("inc%03.0f", inclination_deg)
    aop_tag = @sprintf("aop%03.0f", argp_deg)
    out_path = joinpath(output_dir, "perturbation_force_ratios_$(planet)_$(ms_tag)_$(inc_tag)_$(aop_tag).pdf")
    _save_pdf(fig, out_path)
    println("  Saved: $out_path")
    return out_path
end

function plot_orbit_parameter_effects(run_dir::String, planet::String, parameter::Symbol; output_dir::String=run_dir,
        spacecraft_mass_scale::Float64=1.0, inclination_deg::Float64=93.0, argp_deg::Float64=80.0)
    gm = PLANET_GM[planet]
    index, periapsis_regimes, apo_alts, parameter_values = _build_parameter_sweep_index(run_dir, planet, parameter;
        spacecraft_mass_scale=spacecraft_mass_scale,
        inclination_deg=inclination_deg,
        argp_deg=argp_deg)

    parameter_label = parameter == :inclination ? "Inclination" : "Argument of Periapsis"
    parameter_tag = parameter == :inclination ? "inclination" : "argp"
    fixed_label = parameter == :inclination ? @sprintf("argp=%.0f°", argp_deg) : @sprintf("inc=%.0f°", inclination_deg)

    if isempty(index)
        println("  No $parameter_label sweep force histories for $planet with the selected filters; skipping.")
        return nothing
    end

    trace_styles = _present_traces(index)
    if isempty(trace_styles)
        println("  No plotted $parameter_label sweep perturbation cases for $planet with the selected filters; skipping.")
        return nothing
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
                xlabel=peri_idx == n_peri ? "$parameter_label (deg)" : "",
                ylabel=apo_idx == 1 ? "Peak force / gravity" : "",
                yscale=:log10,
                legend=false,
                ylims=(1e-9, 1.0),
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
                    path = get(index, (periapsis, apo_alt_km, dynamics, density, value), nothing)
                    path === nothing && continue
                    peak_ratio = try
                        _load_peak_ratio(path, gm)
                    catch err
                        @warn "Skipping unreadable trajectory force history" path exception=(err, catch_backtrace())
                        continue
                    end
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
        label = _trace_label(planet, dynamics, density)
        plot!(leg, [NaN], [NaN]; color=color, linestyle=lstyle, marker=:circle,
            markersize=2.2, linewidth=1.4, label=label)
    end

    layout = @eval Plots.@layout [Plots.grid($n_peri, $n_apo) a{0.14w}]
    fig = plot(
        subplots..., leg;
        layout=layout,
        size=(max(900, 260 * n_apo), 300 * n_peri),
        plot_title="$(uppercasefirst(planet))  —  Peak Perturbation Force / Gravity vs $parameter_label ($fixed_label)",
        plot_titlefontsize=11,
        left_margin=GRID_LEFT_MARGIN,
        bottom_margin=GRID_BOTTOM_MARGIN,
        top_margin=GRID_TOP_MARGIN,
    )

    ms_tag = replace(@sprintf("ms%.2f", spacecraft_mass_scale), "." => "p")
    fixed_tag = parameter == :inclination ? @sprintf("aop%03.0f", argp_deg) : @sprintf("inc%03.0f", inclination_deg)
    out_path = joinpath(output_dir, "perturbation_force_ratios_vs_$(parameter_tag)_$(planet)_$(ms_tag)_$(fixed_tag).pdf")
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

function _parse_args(args)
    args = collect(args)
    planets = String[]
    out_dir = nothing
    positional = String[]
    spacecraft_mass_scale = 1.0
    inclination_deg = 93.0
    argp_deg = 80.0

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--planets"
            i += 1
            planets = split(args[i], ",") .|> strip .|> String
        elseif a == "--out-dir"
            i += 1
            out_dir = abspath(args[i])
        elseif a == "--mass-scale"
            i += 1; spacecraft_mass_scale = parse(Float64, args[i])
        elseif a == "--inclination"
            i += 1; inclination_deg = parse(Float64, args[i])
        elseif a == "--argp"
            i += 1; argp_deg = parse(Float64, args[i])
        elseif !startswith(a, "-")
            push!(positional, a)
        end
        i += 1
    end
    run_dir = isempty(positional) ? _most_recent_run(DEFAULT_OUTPUT_BASE) : abspath(positional[1])
    isempty(planets) && (planets = _discover_planets(run_dir))
    isempty(planets) && (planets = copy(PLANET_ORDER))
    out_dir === nothing && (out_dir = joinpath(run_dir, "plots"))
    return run_dir, planets, out_dir, spacecraft_mass_scale, inclination_deg, argp_deg
end

function main(args=ARGS)
    run_dir, planets, out_dir, mass_scale, inclination_deg, argp_deg = _parse_args(args)
    isdir(run_dir) || error("Not a directory: $run_dir")
    mkpath(out_dir)
    println("Run directory : $run_dir")
    println("Output        : $out_dir")
    println("Planets       : $(join(planets, ", "))")
    println("Filter        : mass_scale=$(mass_scale)  inclination=$(inclination_deg)°  argp=$(argp_deg)°")
    for planet in planets
        println("Plotting $planet...")
        flush(stdout)
        plot_planet(run_dir, planet; output_dir=out_dir,
            spacecraft_mass_scale=mass_scale, inclination_deg=inclination_deg, argp_deg=argp_deg)
        plot_orbit_parameter_effects(run_dir, planet, :inclination; output_dir=out_dir,
            spacecraft_mass_scale=mass_scale, inclination_deg=inclination_deg, argp_deg=argp_deg)
        plot_orbit_parameter_effects(run_dir, planet, :argp; output_dir=out_dir,
            spacecraft_mass_scale=mass_scale, inclination_deg=inclination_deg, argp_deg=argp_deg)
    end
    println("Done.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
