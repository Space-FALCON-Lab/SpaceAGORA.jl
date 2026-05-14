# For each planet, plots a 3×3 figure of perturbation force magnitude / central-body gravity.
# Rows = periapsis regimes (shallow / nominal / deep), columns = apoapsis regimes (low / medium / high).
# Each subplot traces all non-two_body dynamics cases found for that (planet, periapsis, apoapsis) slot.
# Reads trajectory_with_active_force.feather files; run convert_to_feather.jl first.
#
# Usage: julia --project=. benchmarks/studies/aerobraking_perturbation_mc/plot_perturbation_force_ratios.jl <run_dir>
#        Omit <run_dir> to use the most recent run in output/aerobraking_perturbation_mc/.
# Options:
#   --planets mars,earth,venus   (default: all three)
#   --out-dir PATH               directory for output PNGs (default: <run_dir>/plots/)

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

const PERIAPSIS_REGIMES = ["shallow", "nominal", "deep"]
const APOAPSIS_REGIMES  = ["low", "medium", "high"]

# (dynamics_case, density_case, line_color, line_style, legend_label)
const TRACE_STYLES = [
    ("j2",               "none",    "#1f77b4", :solid, "J2"),
    ("harmonics_low",    "none",    "#2ca02c", :solid, "Harmonics 4×4"),
    ("srp",              "none",    "#ff7f0e", :solid, "SRP"),
    ("third_body_sun",   "none",    "#d62728", :solid, "3rd Body (Sun)"),
    ("gram_aero",        "high",    "#7b2d8b", :solid, "Aero ρ+25%"),
    ("gram_aero",        "nominal", "#9e4fc4", :dash,  "Aero ρ nominal"),
    ("gram_aero",        "low",     "#c99be0", :dot,   "Aero ρ−25%"),
    ("full_environment", "high",    "#8c3a00", :solid, "Full Env ρ+25%"),
    ("full_environment", "nominal", "#bf5e1a", :dash,  "Full Env ρ nominal"),
    ("full_environment", "low",     "#e0a070", :dot,   "Full Env ρ−25%"),
]

function _parse_case_dir(dirname::String)
    m = match(
        r"^case_\d+_(mars|earth|venus)_(shallow|nominal|deep)_apo_(low|medium|high)_(.+)_density_(none|low|nominal|high)$",
        dirname,
    )
    m === nothing && return nothing
    return (; planet=m[1], periapsis=m[2], apoapsis=m[3], dynamics=m[4], density=m[5])
end

function _build_index(run_dir::String, planet::String)
    index = Dict{NTuple{4, String}, String}()
    for entry in readdir(run_dir; join=true)
        isdir(entry) || continue
        info = _parse_case_dir(basename(entry))
        info === nothing && continue
        info.planet != planet && continue
        p = joinpath(entry, "trajectory_with_active_force.feather")
        isfile(p) || continue
        index[(info.periapsis, info.apoapsis, info.dynamics, info.density)] = p
    end
    return index
end

function _load_ratio(feather_path::String, gm::Float64)
    tbl = Arrow.Table(feather_path)
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

function _downsample(v::AbstractVector, n::Int=2000)
    length(v) <= n && return v
    step = max(1, div(length(v), n))
    return v[1:step:end]
end

function plot_planet(run_dir::String, planet::String; output_dir::String=run_dir)
    gm = PLANET_GM[planet]
    index = _build_index(run_dir, planet)

    if isempty(index)
        println("  No feather files for $planet — run convert_to_feather.jl first.")
        return nothing
    end

    subplots = Plots.Plot[]

    for periapsis in PERIAPSIS_REGIMES
        for apoapsis in APOAPSIS_REGIMES
            subtitle = "$(uppercasefirst(periapsis)) periapsis / $(uppercasefirst(apoapsis)) apoapsis"
            sp = plot(
                title=subtitle,
                titlefontsize=7,
                xlabel="",
                ylabel="",
                yscale=:log10,
                legend=false,
                ylims=(1e-9, 1.0),
                framestyle=:box,
                grid=true,
                gridalpha=0.25,
                tickfontsize=6,
            )

            any_data = false
            for (dynamics, density, color, lstyle, _) in TRACE_STYLES
                path = get(index, (periapsis, apoapsis, dynamics, density), nothing)
                path === nothing && continue
                t_hr, ratio = _load_ratio(path, gm)
                plot!(sp, _downsample(t_hr), _downsample(ratio);
                    color=color, linestyle=lstyle, linewidth=0.8, alpha=0.9, label=nothing)
                any_data = true
            end
            !any_data && annotate!(sp, 0.5, 5e-5, text("no data", :center, 7, :gray))

            push!(subplots, sp)
        end
    end

    # Shared legend panel
    leg = plot(
        framestyle=:none, grid=false, legend=:inside, legendfontsize=7,
        background_color_inside=:transparent, background_color_outside=:transparent,
        xlims=(0, 1), ylims=(0, 1),
    )
    for (_, _, color, lstyle, label) in TRACE_STYLES
        plot!(leg, [NaN], [NaN]; color=color, linestyle=lstyle, linewidth=1.4, label=label)
    end

    fig = plot(
        subplots..., leg;
        layout=@layout([grid(3, 3) a{0.16w}]),
        size=(1500, 900),
        plot_title="$(uppercasefirst(planet))  —  Perturbation Force / Central-Body Gravity",
        plot_titlefontsize=11,
        left_margin=8Plots.mm,
        bottom_margin=6Plots.mm,
        top_margin=4Plots.mm,
    )

    # Shared axis labels via annotation on the overall figure (use inner margins instead)
    out_path = joinpath(output_dir, "perturbation_force_ratios_$(planet).pdf")
    savefig(fig, out_path)
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
    planets = ["mars", "earth", "venus"]
    out_dir = nothing
    positional = String[]

    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--planets"
            i += 1
            planets = split(args[i], ",") .|> strip .|> String
        elseif a == "--out-dir"
            i += 1
            out_dir = abspath(args[i])
        elseif !startswith(a, "-")
            push!(positional, a)
        end
        i += 1
    end
    run_dir = isempty(positional) ? _most_recent_run(DEFAULT_OUTPUT_BASE) : abspath(positional[1])
    out_dir === nothing && (out_dir = joinpath(run_dir, "plots"))
    return run_dir, planets, out_dir
end

function main(args=ARGS)
    run_dir, planets, out_dir = _parse_args(args)
    isdir(run_dir) || error("Not a directory: $run_dir")
    mkpath(out_dir)
    println("Run directory : $run_dir")
    println("Output        : $out_dir")
    println("Planets       : $(join(planets, ", "))")
    for planet in planets
        println("Plotting $planet...")
        flush(stdout)
        plot_planet(run_dir, planet; output_dir=out_dir)
    end
    println("Done.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
