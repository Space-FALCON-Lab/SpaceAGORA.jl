# Generate IPPW-poster-oriented post-processing products from an existing
# aerobraking perturbation Monte Carlo run.
#
# Usage:
#   julia --project=. benchmarks/studies/aerobraking_perturbation_mc/generate_ippw_poster_plots.jl
#   julia --project=. benchmarks/studies/aerobraking_perturbation_mc/generate_ippw_poster_plots.jl output/aerobraking_perturbation_mc/20260519_121929
#   julia --project=. benchmarks/studies/aerobraking_perturbation_mc/generate_ippw_poster_plots.jl base_run supplement_run
#   julia --project=. benchmarks/studies/aerobraking_perturbation_mc/generate_ippw_poster_plots.jl base_run nominal_inclination_argp_run shallow_deep_inclination_argp_run periapsis_altitude_run
#
# The analytical contour here is an additive prediction from the isolated
# perturbation cases in the run, not the final closed-form secular Xi model.

using Arrow
using CSV
using DataFrames
using LinearAlgebra
using Plots
using Printf
using Statistics

include(joinpath(@__DIR__, "analytical_perturbation_models.jl"))
using .AnalyticalPerturbationModels

const PLANETS = ["mars", "venus", "earth", "titan"]
const PLANET_LABEL = Dict("mars" => "Mars", "venus" => "Venus", "earth" => "Earth", "titan" => "Titan")
const REGIME_LABEL = Dict("shallow" => "Shallow", "nominal" => "Nominal", "deep" => "Deep")
const POSTER_DENSITY = "nominal"
const DEFAULT_MASS_SCALE = 1.0
const DEFAULT_INCLINATION_DEG = 93.0
const DEFAULT_ARGP_DEG = 80.0
const PHASE_ARGP_DEGS = [0.0, 90.0, 180.0, 270.0]
const COMBINED_PHASE_ARGP_DEGS = [0.0, 180.0]
const PHASE_INCLINATIONS_DEG = sort(unique(vcat(
    collect(range(15.0, 165.0; step=7.5)),
    [DEFAULT_INCLINATION_DEG],
)))
const PHASE_APOAPSIS_ALTITUDES_KM = Dict(
    "mars" => collect(range(1_000.0, 30_000.0; length=21)),
    "venus" => collect(range(1_000.0, 40_000.0; length=21)),
    "earth" => collect(range(1_000.0, 60_000.0; length=21)),
    "titan" => collect(range(1_000.0, 40_000.0; length=21)),
)
const POSTER_PERIAPSIS = Dict("mars" => "nominal", "venus" => "nominal", "earth" => "nominal", "titan" => "nominal")
const POSTER_APOAPSIS_KM = Dict("mars" => 5000.0, "venus" => 10000.0, "earth" => 36000.0, "titan" => 10000.0)
const PERIAPSIS_ALTITUDE_KM = Dict(
    "mars" => Dict("shallow" => 150.0, "nominal" => 125.0, "deep" => 110.0),
    "venus" => Dict("shallow" => 180.0, "nominal" => 150.0, "deep" => 135.0),
    "earth" => Dict("shallow" => 180.0, "nominal" => 145.0, "deep" => 120.0),
    "titan" => Dict("shallow" => 900.0, "nominal" => 720.0, "deep" => 650.0),
)
const STABLE_EPSILON_KM = 0.1
const ZLK_SECONDS_PER_YEAR = 365.25 * 86400.0

gr()

function latest_run_dir()
    root = joinpath(REPO_ROOT, "output", "aerobraking_perturbation_mc")
    isdir(root) || throw(ArgumentError("No output directory found at $root"))
    runs = filter(name -> isdir(joinpath(root, name)), readdir(root))
    isempty(runs) && throw(ArgumentError("No timestamped runs found under $root"))
    return joinpath(root, sort(runs)[end])
end

function ensure_out_dir(run_dir::String)
    out_dir = joinpath(run_dir, "poster")
    mkpath(out_dir)
    return out_dir
end

function read_table(path::String)::DataFrame
    isfile(path) || throw(ArgumentError("Missing required table: $path"))
    return DataFrame(Arrow.Table(path))
end

function read_results(run_dirs::Vector{String})::DataFrame
    dfs = [read_table(joinpath(run_dir, "results.feather")) for run_dir in run_dirs]
    isempty(dfs) && throw(ArgumentError("At least one run directory is required."))
    return vcat(dfs...; cols=:union)
end

function finite_values(x)
    vals = Float64[]
    for value in skipmissing(x)
        f = Float64(value)
        isfinite(f) && push!(vals, f)
    end
    return vals
end

function representative_apo(rows::DataFrame, planet::String)::Float64
    available = sort(unique(Float64.(rows.apoapsis_alt_km)))
    target = POSTER_APOAPSIS_KM[planet]
    isempty(available) && return target
    return available[argmin(abs.(available .- target))]
end

function row_key(row)
    return (
        String(row.planet),
        String(row.periapsis_regime),
        Float64(row.apoapsis_alt_km),
        Float64(row.spacecraft_mass_scale),
        Float64(row.inclination_deg),
        Float64(row.argp_deg),
    )
end

function two_body_delta_lookup(results::DataFrame)
    tb = results[results.dynamics_case .== "two_body", :]
    lookup = Dict{Tuple{String, String, Float64, Float64, Float64, Float64}, Float64}()
    for row in eachrow(tb)
        lookup[row_key(row)] = Float64(row.delta_rp_alt_m)
    end
    return lookup
end

function perturbation_delta(row, baseline_lookup)
    value = Float64(row.delta_rp_alt_m)
    base = get(baseline_lookup, row_key(row), NaN)
    return isfinite(base) ? value - base : value
end

function additive_prediction_table(results::DataFrame)
    baseline = two_body_delta_lookup(results)
    wanted = Set(["j2", "harmonics_low", "third_body_sun", "gram_aero"])
    rows = NamedTuple[]
    grouped = groupby(results, [:planet, :periapsis_regime, :apoapsis_alt_km, :spacecraft_mass_scale, :inclination_deg, :argp_deg])
    for group in grouped
        parts = Dict{String, Float64}()
        initial_e = Float64(first(group.initial_e))
        initial_a_m = Float64(first(group.initial_a_m))
        requested_rp_alt_m = hasproperty(group, :requested_rp_alt_m) ? Float64(first(group.requested_rp_alt_m)) : NaN
        full_nominal = NaN
        for row in eachrow(group)
            dyn = String(row.dynamics_case)
            density = String(row.density_case)
            if dyn in wanted && (density == "none" || density == POSTER_DENSITY)
                parts[dyn] = perturbation_delta(row, baseline)
            end
            if dyn == "full_environment" && density == POSTER_DENSITY
                full_nominal = Float64(row.delta_rp_alt_m)
            end
        end
        isfinite(full_nominal) || continue
        has_additive_parts = all(has_key -> haskey(parts, has_key), wanted)
        pred_grav_m = has_additive_parts ? parts["j2"] + parts["harmonics_low"] + parts["third_body_sun"] : NaN
        pred_drag_m = has_additive_parts ? parts["gram_aero"] : NaN
        pred_total_m = has_additive_parts ? pred_grav_m + pred_drag_m : NaN
        push!(rows, (
            planet=String(first(group.planet)),
            periapsis_regime=String(first(group.periapsis_regime)),
            apoapsis_alt_km=Float64(first(group.apoapsis_alt_km)),
            spacecraft_mass_scale=Float64(first(group.spacecraft_mass_scale)),
            inclination_deg=Float64(first(group.inclination_deg)),
            argp_deg=Float64(first(group.argp_deg)),
            initial_e=initial_e,
            initial_a_m=initial_a_m,
            requested_rp_alt_m=requested_rp_alt_m,
            j2_delta_rp_m=get(parts, "j2", NaN),
            harmonics_delta_rp_m=get(parts, "harmonics_low", NaN),
            third_body_delta_rp_m=get(parts, "third_body_sun", NaN),
            drag_delta_rp_m=pred_drag_m,
            additive_gravity_delta_rp_m=pred_grav_m,
            additive_total_delta_rp_m=pred_total_m,
            full_environment_delta_rp_m=full_nominal,
            additive_error_m=isfinite(full_nominal) ? pred_total_m - full_nominal : NaN,
        ))
    end
    return DataFrame(rows)
end

function case_trajectory(row)::Union{Arrow.Table, Nothing}
    path = hasproperty(row, :trajectory_with_active_force_feather) ? String(row.trajectory_with_active_force_feather) : ""
    if isempty(path) || !isfile(path)
        case_dir = hasproperty(row, :case_dir) ? String(row.case_dir) : ""
        path = joinpath(case_dir, "trajectory_with_active_force.feather")
    end
    isfile(path) || return nothing
    return Arrow.Table(path)
end

function max_density(tbl)
    hasproperty(tbl, :density_kg_m3) || return NaN
    vals = finite_values(tbl.density_kg_m3)
    return isempty(vals) ? NaN : maximum(vals)
end

function inferred_beta(tbl, mass_scale::Float64)
    ref_area = AnalyticalPerturbationModels._spacecraft_reference_area_m2()
    mass = AnalyticalPerturbationModels._nominal_mass_kg(mass_scale)
    return AnalyticalPerturbationModels._drag_beta_from_history(tbl, mass, ref_area)
end

function pi_values_for_reference(row)
    tbl = case_trajectory(row)
    tbl === nothing && return nothing
    info = (
        planet=String(row.planet),
        periapsis=String(row.periapsis_regime),
        apoapsis_alt_km=Float64(row.apoapsis_alt_km),
        spacecraft_mass_scale=Float64(row.spacecraft_mass_scale),
    )
    ctx = AnalyticalPerturbationModels._context(info, tbl)
    rho_p = max_density(tbl)
    beta = inferred_beta(tbl, Float64(row.spacecraft_mass_scale))
    pi_j2 = AnalyticalPerturbationModels._basic_j2(ctx)
    pi_h = AnalyticalPerturbationModels._basic_harmonics(ctx)
    pi_3b = AnalyticalPerturbationModels._basic_third_body(ctx)
    pi_d = isfinite(rho_p) && beta > 0 ? rho_p * ctx.rp_m / beta : NaN
    return (
        planet=String(row.planet),
        periapsis_regime=String(row.periapsis_regime),
        apoapsis_alt_km=Float64(row.apoapsis_alt_km),
        initial_e=Float64(row.initial_e),
        initial_a_km=Float64(row.initial_a_m) / 1e3,
        pi_j2=pi_j2,
        pi_h=pi_h,
        pi_3b=pi_3b,
        pi_d=pi_d,
        rho_p_kg_m3=rho_p,
        beta_kg_m2=beta,
    )
end

function build_pi_table(results::DataFrame, out_dir::String)
    rows = NamedTuple[]
    for planet in PLANETS
        peri = POSTER_PERIAPSIS[planet]
        planet_rows = results[
            (results.planet .== planet) .&
            (results.periapsis_regime .== peri) .&
            (results.dynamics_case .== "full_environment") .&
            (results.density_case .== POSTER_DENSITY) .&
            (results.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
            (results.inclination_deg .== DEFAULT_INCLINATION_DEG) .&
            (results.argp_deg .== DEFAULT_ARGP_DEG),
            :,
        ]
        isempty(planet_rows) && continue
        apo = representative_apo(planet_rows, planet)
        selected = planet_rows[planet_rows.apoapsis_alt_km .== apo, :]
        isempty(selected) && continue
        values = pi_values_for_reference(first(eachrow(selected)))
        values === nothing || push!(rows, values)
    end
    df = DataFrame(rows)
    CSV.write(joinpath(out_dir, "pi_parameter_reference_values.csv"), df)
    plot_pi_table(df, out_dir)
    return df
end

function plot_pi_table(df::DataFrame, out_dir::String)
    isempty(df) && return nothing
    params = [
        ("Pi_J2", :pi_j2),
        ("Pi_H", :pi_h),
        ("Pi_3b", :pi_3b),
        ("Pi_D", :pi_d),
    ]
    subplots = Plots.Plot[]
    colors = ["#b4492c", "#c18f25", "#2b6cb0", "#6a4c93"]
    for (label, col) in params
        vals = [max(Float64(v), 1e-18) for v in df[!, col]]
        sp = bar(
            [PLANET_LABEL[String(planet)] for planet in df.planet],
            vals;
            yscale=:log10,
            title=label,
            legend=false,
            color=colors[1:nrow(df)],
            ylabel="value",
            framestyle=:box,
            grid=true,
            tickfontsize=7,
            guidefontsize=8,
            titlefontsize=10,
        )
        push!(subplots, sp)
    end
    fig = plot(subplots...; layout=(2, 2), size=(1000, 700), plot_title="Reference Dimensionless Parameters")
    out = joinpath(out_dir, "pi_parameter_reference_values.pdf")
    savefig(fig, out)
    return out
end

function pivot_grid(df::DataFrame, xcol::Symbol, ycol::Symbol, zcol::Symbol)
    xs = sort(unique(Float64.(df[!, xcol])))
    ys = sort(unique(Float64.(df[!, ycol])))
    z = fill(NaN, length(ys), length(xs))
    for row in eachrow(df)
        ix = findfirst(==(Float64(row[xcol])), xs)
        iy = findfirst(==(Float64(row[ycol])), ys)
        ix === nothing || iy === nothing || (z[iy, ix] = Float64(row[zcol]) / 1e3)
    end
    return xs, ys, z
end

function classify_delta_hp_km(value_km::Float64; epsilon_km::Float64=STABLE_EPSILON_KM)::Int
    !isfinite(value_km) && return 0
    value_km > epsilon_km && return 1
    value_km < -epsilon_km && return -1
    return 0
end

function regime_name(value_km::Float64; epsilon_km::Float64=STABLE_EPSILON_KM)::String
    class = classify_delta_hp_km(value_km; epsilon_km)
    class > 0 && return "raising"
    class < 0 && return "lowering"
    return "stable"
end

function regime_grid(z_km; epsilon_km::Float64=STABLE_EPSILON_KM)
    out = fill(NaN, size(z_km))
    for idx in eachindex(z_km)
        value = z_km[idx]
        isfinite(value) || continue
        out[idx] = classify_delta_hp_km(value; epsilon_km)
    end
    return out
end

function periapsis_sort_key(value)::Float64
    label = String(value)
    order = Dict("shallow" => 1.0, "nominal" => 2.0, "deep" => 3.0)
    haskey(order, label) && return order[label]
    m = match(r"^p(\d+)$", label)
    m === nothing || return 100.0 + parse(Float64, m.captures[1])
    return 1_000.0
end

function periapsis_display_label(value)::String
    label = String(value)
    return get(REGIME_LABEL, label, label)
end

function periapsis_altitude_km_for_row(row, planet::String)::Float64
    if hasproperty(row, :requested_rp_alt_m)
        alt_m = Float64(row.requested_rp_alt_m)
        isfinite(alt_m) && alt_m > 0.0 && return alt_m / 1e3
    end
    return PERIAPSIS_ALTITUDE_KM[planet][String(row.periapsis_regime)]
end

function altitude_tick_labels_km(values)
    vals = sort(unique(Float64.(values)))
    labels = [
        abs(v - round(v)) < 1e-6 ? @sprintf("%.0f", v) : @sprintf("%.1f", v)
        for v in vals
    ]
    return vals, labels
end

function altitude_category_panel(df::DataFrame, xcol::Symbol, title::String, xlabel::String)
    sp = plot(
        title=title,
        xlabel=xlabel,
        ylabel="periapsis alt. (km)",
        legend=false,
        framestyle=:none,
        grid=true,
        titlefontsize=7,
        guidefontsize=7,
        tickfontsize=5,
    )
    isempty(df) && return sp
    alts, labels = altitude_tick_labels_km(df.periapsis_alt_km)
    alt_index = Dict(alt => idx for (idx, alt) in enumerate(alts))
    xs = sort(unique(Float64.(df[!, xcol])))
    colors = Dict(-1 => "#c84630", 0 => "#edf6e5", 1 => "#2b6cb0")
    xvals = Float64[]
    yvals = Float64[]
    cvals = String[]
    for row in eachrow(df)
        alt = Float64(row.periapsis_alt_km)
        isfinite(alt) || continue
        push!(xvals, Float64(row[xcol]))
        push!(yvals, Float64(alt_index[alt]))
        class = classify_delta_hp_km(Float64(row.full_environment_delta_rp_m) / 1e3)
        push!(cvals, colors[class])
    end
    scatter!(
        sp,
        xvals,
        yvals;
        marker=:rect,
        markersize=7.5,
        markerstrokewidth=0.35,
        markerstrokecolor=:black,
        color=cvals,
        alpha=0.96,
    )
    xticks!(sp, xs)
    yticks!(sp, (collect(1:length(alts)), labels))
    ylims!(sp, 0.35, length(alts) + 0.65)
    xlims!(sp, minimum(xs) - 8.0, maximum(xs) + 8.0)
    annotate!(sp, 0.03, 0.94, text(@sprintf("%d/%d cells", length(xvals), length(xs) * length(alts)), :left, 5, :black, :relative))
    return sp
end

function add_zero_contour!(sp, xs, ys, z; label="")
    finite = filter(isfinite, vec(z))
    isempty(finite) && return sp
    minimum(finite) <= 0.0 <= maximum(finite) || return sp
    contour!(sp, xs, ys, z; levels=[0.0], color=:black, linewidth=2.0, linestyle=:dash, label=label)
    return sp
end

function slice_label(planet::String, df::DataFrame, xcol::Symbol, ycol::Symbol)
    isempty(df) && return ""
    peri = POSTER_PERIAPSIS[planet]
    fixed = String[]
    push!(fixed, "peri=$(REGIME_LABEL[peri])")
    push!(fixed, @sprintf("mass=%.2g", DEFAULT_MASS_SCALE))
    if ycol != :inclination_deg
        push!(fixed, @sprintf("i=%.0f deg", DEFAULT_INCLINATION_DEG))
    end
    if ycol != :apoapsis_alt_km && length(unique(df.apoapsis_alt_km)) == 1
        push!(fixed, @sprintf("apo=%.0f km", first(df.apoapsis_alt_km)))
    end
    if ycol != :initial_e && length(unique(round.(df.initial_e; digits=3))) == 1
        push!(fixed, @sprintf("e=%.3f", first(df.initial_e)))
    end
    return join(fixed, ", ")
end

function regime_slice(additive::DataFrame, planet::String)
    peri = POSTER_PERIAPSIS[planet]
    base = additive[
        (additive.planet .== planet) .&
        (additive.periapsis_regime .== peri) .&
        (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE),
        :,
    ]
    if planet in ("mars", "titan")
        return base[base.inclination_deg .== DEFAULT_INCLINATION_DEG, :], :argp_deg, :initial_e, "omega (deg)", "eccentricity"
    elseif planet == "venus"
        apo = representative_apo(base, planet)
        return base[(base.apoapsis_alt_km .== apo) .& (base.argp_deg .!= NaN), :], :argp_deg, :inclination_deg, "omega (deg)", "inclination (deg)"
    else
        return base[base.inclination_deg .== DEFAULT_INCLINATION_DEG, :], :argp_deg, :apoapsis_alt_km, "omega (deg)", "apoapsis altitude (km)"
    end
end

function plot_regime_maps(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    regime_colors = cgrad(["#c84630", "#edf6e5", "#2b6cb0"], [0.0, 0.5, 1.0])
    for planet in PLANETS
        df, xcol, ycol, xlabel, ylabel = regime_slice(additive, planet)
        if isempty(df)
            push!(subplots, plot(title=PLANET_LABEL[planet], framestyle=:none))
            continue
        end
        xs, ys, z_full = pivot_grid(df, xcol, ycol, :full_environment_delta_rp_m)
        _, _, z_pred = pivot_grid(df, xcol, ycol, :additive_total_delta_rp_m)
        z_regime = regime_grid(z_full)
        title = "$(PLANET_LABEL[planet])\n$(slice_label(planet, df, xcol, ycol))"
        sp = heatmap(
            xs,
            ys,
            z_regime;
            title=title,
            xlabel=xlabel,
            ylabel=ylabel,
            color=regime_colors,
            clims=(-1.5, 1.5),
            colorbar=false,
            framestyle=:box,
            guidefontsize=8,
            tickfontsize=7,
            titlefontsize=8,
        )
        add_zero_contour!(sp, xs, ys, z_full; label="numerical Delta hp=0")
        add_zero_contour!(sp, xs, ys, z_pred; label="additive prediction")
        scatter!(sp, Float64.(df[!, xcol]), Float64.(df[!, ycol]); marker=:circle, markersize=2.0, color=:black, alpha=0.45, label="")
        push!(subplots, sp)
    end
    legend = plot(framestyle=:none, grid=false, xlims=(0, 1), ylims=(0, 1), legend=:inside, background_color_inside=:transparent)
    scatter!(legend, [NaN], [NaN]; color="#2b6cb0", marker=:rect, markersize=9, label="raising")
    scatter!(legend, [NaN], [NaN]; color="#edf6e5", marker=:rect, markersize=9, label="stable")
    scatter!(legend, [NaN], [NaN]; color="#c84630", marker=:rect, markersize=9, label="lowering")
    plot!(legend, [NaN], [NaN]; color=:black, linewidth=2.0, linestyle=:dash, label="zero contour")
    fig = plot(
        subplots..., legend;
        layout=@layout([a b; c d; e{0.08h}]),
        size=(1200, 1050),
        plot_title=@sprintf("Numerical Per-Orbit Delta hp Regime Maps (stable: |Delta hp| <= %.2f km/orbit)", STABLE_EPSILON_KM),
        plot_titlefontsize=13,
        left_margin=10Plots.mm,
        bottom_margin=10Plots.mm,
    )
    out = joinpath(out_dir, "regime_maps_delta_hp_2x2.pdf")
    savefig(fig, out)
    return out
end

function plot_delta_hp_vs_argp(additive::DataFrame, out_dir::String)
    colors = Dict("mars" => "#b4492c", "venus" => "#c18f25", "earth" => "#2b6cb0", "titan" => "#6a4c93")
    subplots = Plots.Plot[]
    for planet in PLANETS
        peri = POSTER_PERIAPSIS[planet]
        base = additive[
            (additive.planet .== planet) .&
            (additive.periapsis_regime .== peri) .&
            (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
            (additive.inclination_deg .== DEFAULT_INCLINATION_DEG),
            :,
        ]
        isempty(base) && continue
        apo = representative_apo(base, planet)
        df = sort(base[base.apoapsis_alt_km .== apo, :], :argp_deg)
        isempty(df) && continue
        ys = Float64.(df.full_environment_delta_rp_m) ./ 1e3
        ypred = Float64.(df.additive_total_delta_rp_m) ./ 1e3
        finite = filter(isfinite, vcat(ys, ypred))
        lim = isempty(finite) ? 1.0 : max(STABLE_EPSILON_KM * 3.0, quantile(abs.(finite), 0.96) * 1.15)
        sp = plot(
            xlabel="omega (deg)",
            ylabel="Delta hp (km/orbit)",
            title=@sprintf("%s: %s, apo=%.0f km, i=%.0f deg", PLANET_LABEL[planet], REGIME_LABEL[peri], apo, DEFAULT_INCLINATION_DEG),
            framestyle=:box,
            grid=true,
            legend=:topright,
            size=(550, 380),
            ylims=(-lim, lim),
            titlefontsize=8,
            guidefontsize=8,
            tickfontsize=7,
        )
        hspan!(sp, [STABLE_EPSILON_KM, lim]; color="#dceeff", alpha=0.32, label="raising")
        hspan!(sp, [-lim, -STABLE_EPSILON_KM]; color="#f8ded8", alpha=0.32, label="lowering")
        hspan!(sp, [-STABLE_EPSILON_KM, STABLE_EPSILON_KM]; color="#edf6e5", alpha=0.70, label="stable")
        hline!(sp, [0.0]; color=:black, linestyle=:dot, linewidth=1.2, label="")
        plot!(sp, df.argp_deg, ys; color=colors[planet], marker=:circle, linewidth=2.2, label="numerical")
        plot!(sp, df.argp_deg, ypred; color=:black, linestyle=:dash, linewidth=1.4, label="additive")
        push!(subplots, sp)
    end
    sp = plot(
        subplots...;
        layout=(2, 2),
        size=(1200, 820),
        plot_title=@sprintf("Reference Delta hp(omega) Slices (stable band: +/- %.2f km/orbit)", STABLE_EPSILON_KM),
        plot_titlefontsize=13,
        left_margin=9Plots.mm,
        bottom_margin=9Plots.mm,
    )
    out = joinpath(out_dir, "delta_hp_vs_argp_reference.pdf")
    savefig(sp, out)
    return out
end

function write_reference_corridor_summary(additive::DataFrame, out_dir::String)
    rows = NamedTuple[]
    for planet in PLANETS
        peri = POSTER_PERIAPSIS[planet]
        base = additive[
            (additive.planet .== planet) .&
            (additive.periapsis_regime .== peri) .&
            (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
            (additive.inclination_deg .== DEFAULT_INCLINATION_DEG),
            :,
        ]
        isempty(base) && continue
        apo = representative_apo(base, planet)
        df = sort(base[base.apoapsis_alt_km .== apo, :], :argp_deg)
        isempty(df) && continue
        classes = [regime_name(Float64(v) / 1e3) for v in df.full_environment_delta_rp_m]
        n = length(classes)
        stable_omega = [@sprintf("%.0f", Float64(row.argp_deg)) for (idx, row) in enumerate(eachrow(df)) if classes[idx] == "stable"]
        raising_omega = [@sprintf("%.0f", Float64(row.argp_deg)) for (idx, row) in enumerate(eachrow(df)) if classes[idx] == "raising"]
        lowering_omega = [@sprintf("%.0f", Float64(row.argp_deg)) for (idx, row) in enumerate(eachrow(df)) if classes[idx] == "lowering"]
        push!(rows, (
            planet=planet,
            periapsis_regime=peri,
            apoapsis_alt_km=apo,
            inclination_deg=DEFAULT_INCLINATION_DEG,
            stable_epsilon_km=STABLE_EPSILON_KM,
            raising_fraction=count(==("raising"), classes) / n,
            stable_fraction=count(==("stable"), classes) / n,
            lowering_fraction=count(==("lowering"), classes) / n,
            stable_omega_deg=join(stable_omega, ";"),
            raising_omega_deg=join(raising_omega, ";"),
            lowering_omega_deg=join(lowering_omega, ";"),
            min_delta_hp_km=minimum(Float64.(df.full_environment_delta_rp_m)) / 1e3,
            max_delta_hp_km=maximum(Float64.(df.full_environment_delta_rp_m)) / 1e3,
        ))
    end
    df = DataFrame(rows)
    CSV.write(joinpath(out_dir, "reference_corridor_summary.csv"), df)
    return df
end

function regime_fraction_table(additive::DataFrame; epsilon_km::Float64=STABLE_EPSILON_KM)
    rows = NamedTuple[]
    for group in groupby(additive, [:planet, :periapsis_regime])
        vals = finite_values(group.full_environment_delta_rp_m) ./ 1e3
        isempty(vals) && continue
        n = length(vals)
        classes = [classify_delta_hp_km(v; epsilon_km) for v in vals]
        push!(rows, (
            planet=String(first(group.planet)),
            periapsis_regime=String(first(group.periapsis_regime)),
            n=n,
            raising_fraction=count(==(1), classes) / n,
            stable_fraction=count(==(0), classes) / n,
            lowering_fraction=count(==(-1), classes) / n,
            median_delta_hp_km=median(vals),
            p95_abs_delta_hp_km=quantile(abs.(vals), 0.95),
        ))
    end
    return DataFrame(rows)
end

function plot_regime_fraction_summary(additive::DataFrame, out_dir::String)
    df = regime_fraction_table(additive)
    CSV.write(joinpath(out_dir, "full_survey_regime_fractions.csv"), df)
    subplots = Plots.Plot[]
    colors = Dict("lowering" => "#c84630", "stable" => "#edf6e5", "raising" => "#2b6cb0")
    for planet in PLANETS
        local_df = df[df.planet .== planet, :]
        isempty(local_df) && continue
        sort!(local_df, :periapsis_regime; by=periapsis_sort_key)
        labels = [periapsis_display_label(x) for x in local_df.periapsis_regime]
        sp = plot(
            title=PLANET_LABEL[planet],
            ylabel="fraction of full sweep",
            xlims=(0.4, length(labels) + 0.6),
            ylims=(0, 1),
            legend=planet == "mars" ? :outertopright : false,
            framestyle=:box,
            grid=true,
            xticks=(1:length(labels), labels),
            titlefontsize=9,
            guidefontsize=8,
            tickfontsize=7,
        )
        for (idx, row) in enumerate(eachrow(local_df))
            bottom = 0.0
            segments = (
                ("lowering", Float64(row.lowering_fraction)),
                ("stable", Float64(row.stable_fraction)),
                ("raising", Float64(row.raising_fraction)),
            )
            for (name, fraction) in segments
                top = bottom + fraction
                shape = Shape(
                    [idx - 0.32, idx + 0.32, idx + 0.32, idx - 0.32],
                    [bottom, bottom, top, top],
                )
                plot!(sp, shape; color=colors[name], linecolor=:black, linewidth=0.3,
                    label=(planet == "mars" && idx == 1) ? name : "")
                bottom = top
            end
        end
        push!(subplots, sp)
    end
    fig = plot(
        subplots...;
        layout=(2, 2),
        size=(1100, 800),
        plot_title=@sprintf("Full 5D Survey Regime Fractions (stable: |Delta hp| <= %.2f km/orbit)", STABLE_EPSILON_KM),
        plot_titlefontsize=13,
        left_margin=10Plots.mm,
        bottom_margin=9Plots.mm,
    )
    out = joinpath(out_dir, "full_survey_regime_fractions.pdf")
    savefig(fig, out)
    return out
end

function atlas_panel(df::DataFrame, xcol::Symbol, ycol::Symbol, title::String, xlabel::String, ylabel::String)
    colors = Dict(-1 => "#c84630", 0 => "#edf6e5", 1 => "#2b6cb0")
    sp = plot(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        legend=false,
        framestyle=:box,
        grid=true,
        titlefontsize=7,
        guidefontsize=7,
        tickfontsize=5,
    )
    isempty(df) && return sp
    xs, ys, z_full = pivot_grid(df, xcol, ycol, :full_environment_delta_rp_m)
    filled = count(isfinite, z_full)
    total = length(z_full)
    dense_enough = filled >= 0.7 * total && length(xs) >= 3 && length(ys) >= 3
    if dense_enough
        heatmap!(sp, xs, ys, regime_grid(z_full); color=cgrad(["#c84630", "#edf6e5", "#2b6cb0"], [0.0, 0.5, 1.0]),
            clims=(-1.5, 1.5), colorbar=false)
        add_zero_contour!(sp, xs, ys, z_full)
    end
    classes = classify_delta_hp_km.(Float64.(df.full_environment_delta_rp_m) ./ 1e3)
    scatter!(
        sp,
        Float64.(df[!, xcol]),
        Float64.(df[!, ycol]);
        marker=:circle,
        markersize=dense_enough ? 1.7 : 4.0,
        markerstrokewidth=0.25,
        markerstrokecolor=:black,
        color=[colors[c] for c in classes],
        alpha=dense_enough ? 0.65 : 0.95,
    )
    annotate!(sp, 0.03, 0.94, text(@sprintf("%d/%d cells", filled, total), :left, 5, :black, :relative))
    if !dense_enough
        annotate!(sp, 0.03, 0.84, text("sparse", :left, 5, :black, :relative))
    end
    return sp
end

function plot_omega_e_atlas(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        for peri in ["shallow", "nominal", "deep"]
            df = additive[
                (additive.planet .== planet) .&
                (additive.periapsis_regime .== peri) .&
                (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
                (additive.inclination_deg .== DEFAULT_INCLINATION_DEG),
                :,
            ]
            title = @sprintf("%s %s, i=%.0f deg", PLANET_LABEL[planet], REGIME_LABEL[peri], DEFAULT_INCLINATION_DEG)
            push!(subplots, atlas_panel(df, :argp_deg, :initial_e, title, "omega (deg)", "e"))
        end
    end
    fig = plot(
        subplots...;
        layout=(4, 3),
        size=(1350, 1200),
        plot_title="Slice Atlas: omega vs eccentricity, fixed inclination",
        plot_titlefontsize=13,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    out = joinpath(out_dir, "slice_atlas_omega_vs_e.pdf")
    savefig(fig, out)
    return out
end

function plot_omega_inclination_reference_atlas(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        for peri in ["shallow", "nominal", "deep"]
            base = additive[
                (additive.planet .== planet) .&
                (additive.periapsis_regime .== peri) .&
                (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE),
                :,
            ]
            isempty(base) && continue
            apo = representative_apo(base, planet)
            df = base[base.apoapsis_alt_km .== apo, :]
            title = @sprintf("%s %s, apo=%.0f km", PLANET_LABEL[planet], REGIME_LABEL[peri], apo)
            push!(subplots, atlas_panel(df, :argp_deg, :inclination_deg, title, "omega (deg)", "inclination (deg)"))
        end
    end
    fig = plot(
        subplots...;
        layout=(4, 3),
        size=(1350, 1200),
        plot_title="Slice Atlas: omega vs inclination, representative apoapsis",
        plot_titlefontsize=13,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    out = joinpath(out_dir, "slice_atlas_omega_vs_inclination_reference_apo.pdf")
    savefig(fig, out)
    return out
end

function plot_omega_apo_atlas(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        for peri in ["shallow", "nominal", "deep"]
            df = additive[
                (additive.planet .== planet) .&
                (additive.periapsis_regime .== peri) .&
                (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
                (additive.inclination_deg .== DEFAULT_INCLINATION_DEG),
                :,
            ]
            title = @sprintf("%s %s, i=%.0f deg", PLANET_LABEL[planet], REGIME_LABEL[peri], DEFAULT_INCLINATION_DEG)
            push!(subplots, atlas_panel(df, :argp_deg, :apoapsis_alt_km, title, "omega (deg)", "apoapsis (km)"))
        end
    end
    fig = plot(
        subplots...;
        layout=(4, 3),
        size=(1350, 1200),
        plot_title="Slice Atlas: omega vs apoapsis altitude, fixed inclination",
        plot_titlefontsize=13,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    out = joinpath(out_dir, "slice_atlas_omega_vs_apoapsis.pdf")
    savefig(fig, out)
    return out
end

function plot_omega_periapsis_atlas(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        base = additive[
            (additive.planet .== planet) .&
            (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
            (additive.inclination_deg .== DEFAULT_INCLINATION_DEG),
            :,
        ]
        isempty(base) && continue
        apo = representative_apo(base, planet)
        df = base[base.apoapsis_alt_km .== apo, :]
        mapped = copy(df)
        mapped[!, :periapsis_alt_km] = [periapsis_altitude_km_for_row(row, planet) for row in eachrow(mapped)]
        title = @sprintf("%s, apo=%.0f km, i=%.0f deg", PLANET_LABEL[planet], apo, DEFAULT_INCLINATION_DEG)
        push!(subplots, atlas_panel(mapped, :argp_deg, :periapsis_alt_km, title, "omega (deg)", "periapsis alt. (km)"))
    end
    fig = plot(
        subplots...;
        layout=(2, 2),
        size=(1100, 800),
        plot_title="Slice Atlas: omega vs periapsis altitude, representative apoapsis",
        plot_titlefontsize=13,
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    out = joinpath(out_dir, "slice_atlas_omega_vs_periapsis_regime.pdf")
    savefig(fig, out)
    return out
end

function plot_inclination_e_atlas(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        for peri in ["shallow", "nominal", "deep"]
            df = additive[
                (additive.planet .== planet) .&
                (additive.periapsis_regime .== peri) .&
                (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
                (additive.argp_deg .== DEFAULT_ARGP_DEG),
                :,
            ]
            title = @sprintf("%s %s, omega=%.0f deg", PLANET_LABEL[planet], REGIME_LABEL[peri], DEFAULT_ARGP_DEG)
            push!(subplots, atlas_panel(df, :inclination_deg, :initial_e, title, "inclination (deg)", "e"))
        end
    end
    fig = plot(
        subplots...;
        layout=(4, 3),
        size=(1350, 1200),
        plot_title="Slice Atlas: inclination vs eccentricity, fixed omega",
        plot_titlefontsize=13,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    out = joinpath(out_dir, "slice_atlas_inclination_vs_e.pdf")
    savefig(fig, out)
    return out
end

function plot_inclination_apo_atlas(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        for peri in ["shallow", "nominal", "deep"]
            df = additive[
                (additive.planet .== planet) .&
                (additive.periapsis_regime .== peri) .&
                (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
                (additive.argp_deg .== DEFAULT_ARGP_DEG),
                :,
            ]
            title = @sprintf("%s %s, omega=%.0f deg", PLANET_LABEL[planet], REGIME_LABEL[peri], DEFAULT_ARGP_DEG)
            push!(subplots, atlas_panel(df, :inclination_deg, :apoapsis_alt_km, title, "inclination (deg)", "apoapsis (km)"))
        end
    end
    fig = plot(
        subplots...;
        layout=(4, 3),
        size=(1350, 1200),
        plot_title="Slice Atlas: inclination vs apoapsis altitude, fixed omega",
        plot_titlefontsize=13,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    out = joinpath(out_dir, "slice_atlas_inclination_vs_apoapsis.pdf")
    savefig(fig, out)
    return out
end

function plot_inclination_omega_reference_atlas(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        for peri in ["shallow", "nominal", "deep"]
            base = additive[
                (additive.planet .== planet) .&
                (additive.periapsis_regime .== peri) .&
                (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE),
                :,
            ]
            isempty(base) && continue
            apo = representative_apo(base, planet)
            df = base[base.apoapsis_alt_km .== apo, :]
            title = @sprintf("%s %s, apo=%.0f km", PLANET_LABEL[planet], REGIME_LABEL[peri], apo)
            push!(subplots, atlas_panel(df, :inclination_deg, :argp_deg, title, "inclination (deg)", "omega (deg)"))
        end
    end
    fig = plot(
        subplots...;
        layout=(4, 3),
        size=(1350, 1200),
        plot_title="Slice Atlas: inclination vs omega, representative apoapsis",
        plot_titlefontsize=13,
        left_margin=7Plots.mm,
        bottom_margin=7Plots.mm,
    )
    out = joinpath(out_dir, "slice_atlas_inclination_vs_omega_reference_apo.pdf")
    savefig(fig, out)
    return out
end

function plot_inclination_periapsis_atlas(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        base = additive[
            (additive.planet .== planet) .&
            (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
            (additive.argp_deg .== DEFAULT_ARGP_DEG),
            :,
        ]
        isempty(base) && continue
        apo = representative_apo(base, planet)
        df = base[base.apoapsis_alt_km .== apo, :]
        mapped = copy(df)
        mapped[!, :periapsis_alt_km] = [periapsis_altitude_km_for_row(row, planet) for row in eachrow(mapped)]
        title = @sprintf("%s, apo=%.0f km, omega=%.0f deg", PLANET_LABEL[planet], apo, DEFAULT_ARGP_DEG)
        push!(subplots, atlas_panel(mapped, :inclination_deg, :periapsis_alt_km, title, "inclination (deg)", "periapsis alt. (km)"))
    end
    fig = plot(
        subplots...;
        layout=(2, 2),
        size=(1100, 800),
        plot_title="Slice Atlas: inclination vs periapsis altitude, representative apoapsis",
        plot_titlefontsize=13,
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    out = joinpath(out_dir, "slice_atlas_inclination_vs_periapsis_regime.pdf")
    savefig(fig, out)
    return out
end

function plot_slice_atlas(additive::DataFrame, out_dir::String)
    paths = String[]
    push!(paths, plot_omega_e_atlas(additive, out_dir))
    push!(paths, plot_omega_apo_atlas(additive, out_dir))
    push!(paths, plot_omega_inclination_reference_atlas(additive, out_dir))
    push!(paths, plot_omega_periapsis_atlas(additive, out_dir))
    push!(paths, plot_inclination_e_atlas(additive, out_dir))
    push!(paths, plot_inclination_apo_atlas(additive, out_dir))
    push!(paths, plot_inclination_omega_reference_atlas(additive, out_dir))
    push!(paths, plot_inclination_periapsis_atlas(additive, out_dir))
    return paths
end

function priority_panel_title(planet::String, subtitle::String)
    return "$(PLANET_LABEL[planet])\n$(subtitle)"
end

function plot_priority_omega_vs_apoapsis(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        peri = POSTER_PERIAPSIS[planet]
        df = additive[
            (additive.planet .== planet) .&
            (additive.periapsis_regime .== peri) .&
            (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
            (additive.inclination_deg .== DEFAULT_INCLINATION_DEG),
            :,
        ]
        title = priority_panel_title(planet, @sprintf("%s, i=%.0f deg", REGIME_LABEL[peri], DEFAULT_INCLINATION_DEG))
        push!(subplots, atlas_panel(df, :argp_deg, :apoapsis_alt_km, title, "omega (deg)", "apoapsis (km)"))
    end
    fig = plot(
        subplots...;
        layout=(2, 2),
        size=(1150, 850),
        plot_title="Priority Slice: omega vs apoapsis altitude",
        plot_titlefontsize=13,
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    out = joinpath(out_dir, "priority_slice_omega_vs_apoapsis_nominal_i093.pdf")
    savefig(fig, out)
    return out
end

function plot_priority_omega_vs_eccentricity(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        peri = POSTER_PERIAPSIS[planet]
        df = additive[
            (additive.planet .== planet) .&
            (additive.periapsis_regime .== peri) .&
            (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
            (additive.inclination_deg .== DEFAULT_INCLINATION_DEG),
            :,
        ]
        title = priority_panel_title(planet, @sprintf("%s, i=%.0f deg", REGIME_LABEL[peri], DEFAULT_INCLINATION_DEG))
        push!(subplots, atlas_panel(df, :argp_deg, :initial_e, title, "omega (deg)", "eccentricity"))
    end
    fig = plot(
        subplots...;
        layout=(2, 2),
        size=(1150, 850),
        plot_title="Priority Slice: omega vs eccentricity",
        plot_titlefontsize=13,
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    out = joinpath(out_dir, "priority_slice_omega_vs_eccentricity_nominal_i093.pdf")
    savefig(fig, out)
    return out
end

function plot_priority_inclination_vs_apoapsis(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        peri = POSTER_PERIAPSIS[planet]
        df = additive[
            (additive.planet .== planet) .&
            (additive.periapsis_regime .== peri) .&
            (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
            (additive.argp_deg .== DEFAULT_ARGP_DEG),
            :,
        ]
        title = priority_panel_title(planet, @sprintf("%s, omega=%.0f deg", REGIME_LABEL[peri], DEFAULT_ARGP_DEG))
        push!(subplots, atlas_panel(df, :inclination_deg, :apoapsis_alt_km, title, "inclination (deg)", "apoapsis (km)"))
    end
    fig = plot(
        subplots...;
        layout=(2, 2),
        size=(1150, 850),
        plot_title="Priority Slice: inclination vs apoapsis altitude",
        plot_titlefontsize=13,
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    out = joinpath(out_dir, "priority_slice_inclination_vs_apoapsis_nominal_aop080.pdf")
    savefig(fig, out)
    return out
end

function plot_priority_omega_vs_periapsis_altitude(additive::DataFrame, out_dir::String)
    subplots = Plots.Plot[]
    for planet in PLANETS
        base = additive[
            (additive.planet .== planet) .&
            (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
            (additive.inclination_deg .== DEFAULT_INCLINATION_DEG),
            :,
        ]
        isempty(base) && continue
        apo = representative_apo(base, planet)
        df = copy(base[base.apoapsis_alt_km .== apo, :])
        df[!, :periapsis_alt_km] = [
            periapsis_altitude_km_for_row(row, planet)
            for row in eachrow(df)
        ]
        title = priority_panel_title(planet, @sprintf("apo=%.0f km, i=%.0f deg", apo, DEFAULT_INCLINATION_DEG))
        push!(subplots, atlas_panel(df, :argp_deg, :periapsis_alt_km, title, "omega (deg)", "periapsis alt. (km)"))
    end
    fig = plot(
        subplots...;
        layout=(2, 2),
        size=(1150, 850),
        plot_title="Priority Slice: omega vs periapsis altitude",
        plot_titlefontsize=13,
        left_margin=8Plots.mm,
        bottom_margin=8Plots.mm,
    )
    out = joinpath(out_dir, "priority_slice_omega_vs_periapsis_altitude_reference_apo_i093.pdf")
    savefig(fig, out)
    return out
end

function plot_priority_slices(additive::DataFrame, out_dir::String)
    paths = String[]
    push!(paths, plot_priority_omega_vs_apoapsis(additive, out_dir))
    push!(paths, plot_priority_omega_vs_eccentricity(additive, out_dir))
    push!(paths, plot_priority_inclination_vs_apoapsis(additive, out_dir))
    push!(paths, plot_priority_omega_vs_periapsis_altitude(additive, out_dir))
    return paths
end

function phase_grid_coverage(additive::DataFrame, out_dir::String)
    rows = NamedTuple[]
    for planet in PLANETS
        peri = POSTER_PERIAPSIS[planet]
        base = phase_grid_subset(additive[
            (additive.planet .== planet) .&
            (additive.periapsis_regime .== peri) .&
            (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE),
            :,
        ], planet)
        for argp in PHASE_ARGP_DEGS
            df = base[base.argp_deg .== argp, :]
            xs = isempty(df) ? Float64[] : sort(unique(Float64.(df.apoapsis_alt_km)))
            ys = isempty(df) ? Float64[] : sort(unique(Float64.(df.inclination_deg)))
            total = length(xs) * length(ys)
            push!(rows, (
                planet=planet,
                periapsis_regime=peri,
                argp_deg=argp,
                rows=nrow(df),
                apoapsis_levels=length(xs),
                inclination_levels=length(ys),
                nominal_grid_cells=total,
                coverage_fraction=total > 0 ? nrow(df) / total : 0.0,
            ))
        end
    end
    df = DataFrame(rows)
    CSV.write(joinpath(out_dir, "phase_inclination_apoapsis_coverage.csv"), df)
    return df
end

function plot_phase_inclination_apoapsis(additive::DataFrame, out_dir::String)
    paths = String[]
    combined_subplots = Plots.Plot[]
    phase_grid_coverage(additive, out_dir)
    for (planet_idx, planet) in enumerate(PLANETS)
        peri = POSTER_PERIAPSIS[planet]
        base = phase_grid_subset(additive[
            (additive.planet .== planet) .&
            (additive.periapsis_regime .== peri) .&
            (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE),
            :,
        ], planet)
        isempty(base) && continue
        subplots = Plots.Plot[]
        for (argp_idx, argp) in enumerate(PHASE_ARGP_DEGS)
            df = base[base.argp_deg .== argp, :]
            title = @sprintf("AoP=%.0f deg", argp)
            ylabel = argp_idx == 1 ? "$(PLANET_LABEL[planet])\nInclination (deg)" : ""
            combined_xlabel = planet_idx == length(PLANETS) ? "Apoapsis (km)" : ""
            panel = phase_box_panel(df, :apoapsis_alt_km, :inclination_deg, title, "Apoapsis (km)", ylabel)
            push!(subplots, panel)
            if argp in COMBINED_PHASE_ARGP_DEGS
                combined_title = planet_idx == 1 ? title : ""
                combined_panel = phase_box_panel(df, :apoapsis_alt_km, :inclination_deg, combined_title, combined_xlabel, ylabel)
                push!(combined_subplots, combined_panel)
            end
        end
        fig = plot(
            subplots...;
            layout=(1, 4),
            size=(1700, 420),
            left_margin=16Plots.mm,
            right_margin=8Plots.mm,
            top_margin=8Plots.mm,
            bottom_margin=18Plots.mm,
        )
        out = joinpath(out_dir, @sprintf("phase_dependence_inclination_vs_apoapsis_%s_%s_full_environment.pdf", planet, peri))
        savefig(fig, out)
        push!(paths, out)
    end
    if !isempty(combined_subplots)
        fig = plot(
            phase_regime_legend_panel(), combined_subplots...;
            layout=@layout([l{0.08h}; a b; c d; e f; g h]),
            size=(1200, 2600),
            left_margin=18Plots.mm,
            right_margin=10Plots.mm,
            top_margin=8Plots.mm,
            bottom_margin=20Plots.mm,
        )
        out = joinpath(out_dir, "phase_dependence_inclination_vs_apoapsis_all_planets_full_environment.pdf")
        savefig(fig, out)
        push!(paths, out)
    end
    return paths
end

function phase_regime_legend_panel()
    sp = plot(
        xlims=(0, 1),
        ylims=(0, 1),
        framestyle=:none,
        grid=false,
        legend=false,
        xticks=false,
        yticks=false,
    )
    legend_box = Shape(
        [0.08, 0.98, 0.98, 0.08],
        [0.12, 0.12, 0.90, 0.90],
    )
    plot!(sp, legend_box; color=:white, fillalpha=0.0, linecolor=:black, linewidth=0.7)
    entries = [
        (0.22, "#c84630", "Lowering"),
        (0.50, "#edf6e5", "Stable"),
        (0.78, "#2b6cb0", "Raising"),
    ]
    for (x, color, label) in entries
        box = Shape(
            [x - 0.045, x + 0.015, x + 0.015, x - 0.045],
            [0.36, 0.36, 0.66, 0.66],
        )
        plot!(sp, box; color=color, linecolor=:black, linewidth=0.5)
        annotate!(sp, x + 0.035, 0.51, text(label, :left, 18, :black))
    end
    return sp
end

function phase_value_in_grid(value, grid; atol::Float64=1e-6)::Bool
    f = Float64(value)
    return any(abs(f - Float64(g)) <= atol for g in grid)
end

function phase_grid_subset(df::DataFrame, planet::String)::DataFrame
    isempty(df) && return df
    apo_grid = PHASE_APOAPSIS_ALTITUDES_KM[planet]
    keep = [
        phase_value_in_grid(row.apoapsis_alt_km, apo_grid) &&
        phase_value_in_grid(row.inclination_deg, PHASE_INCLINATIONS_DEG) &&
        phase_value_in_grid(row.argp_deg, PHASE_ARGP_DEGS)
        for row in eachrow(df)
    ]
    return df[keep, :]
end

function phase_sparse_ticks(values; max_ticks::Int=6)
    vals = sort(unique(Float64.(values)))
    isempty(vals) && return (Float64[], String[])
    nticks = min(max_ticks, length(vals))
    idxs = unique(round.(Int, range(1, length(vals); length=nticks)))
    ticks = vals[idxs]
    labels = [
        abs(v - round(v)) < 1e-6 ? @sprintf("%.0f", v) : @sprintf("%.1f", v)
        for v in ticks
    ]
    return ticks, labels
end

function phase_box_panel(df::DataFrame, xcol::Symbol, ycol::Symbol, title::String, xlabel::String, ylabel::String)
    sp = plot(
        title=title,
        xlabel=xlabel,
        ylabel=ylabel,
        legend=false,
        framestyle=:box,
        grid=true,
        titlefontsize=18,
        guidefontsize=18,
        tickfontsize=18,
        xrotation=45,
    )
    isempty(df) && return sp

    xs, ys, z_full = pivot_grid(df, xcol, ycol, :full_environment_delta_rp_m)
    heatmap!(
        sp,
        xs,
        ys,
        regime_grid(z_full);
        color=cgrad(["#c84630", "#edf6e5", "#2b6cb0"], [0.0, 0.5, 1.0]),
        clims=(-1.5, 1.5),
        colorbar=false,
    )
    xticks!(sp, phase_sparse_ticks(xs; max_ticks=5))
    yticks!(sp, phase_sparse_ticks(ys; max_ticks=6))
    return sp
end

function plot_equilibrium_terms(additive::DataFrame, out_dir::String)
    df = additive[
        (additive.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
        (additive.inclination_deg .== DEFAULT_INCLINATION_DEG),
        :,
    ]
    colors = Dict("mars" => "#b4492c", "venus" => "#c18f25", "earth" => "#2b6cb0", "titan" => "#6a4c93")
    sp = plot(
        xlabel="drag contribution, Delta hp drag (km/orbit)",
        ylabel="gravity contribution, Delta hp grav (km/orbit)",
        title="Additive Equilibrium Diagram from Isolated Perturbation Cases",
        framestyle=:box,
        grid=true,
        legend=:outertopright,
        size=(900, 720),
    )
    vals = vcat(finite_values(df.drag_delta_rp_m), finite_values(df.additive_gravity_delta_rp_m)) ./ 1e3
    lim = isempty(vals) ? 1.0 : max(1.0, quantile(abs.(vals), 0.98))
    plot!(sp, [-lim, lim], [lim, -lim]; color=:black, linestyle=:dash, linewidth=2, label="Delta hp total=0")
    for planet in PLANETS
        local_df = df[df.planet .== planet, :]
        isempty(local_df) && continue
        scatter!(sp, Float64.(local_df.drag_delta_rp_m) ./ 1e3, Float64.(local_df.additive_gravity_delta_rp_m) ./ 1e3;
            markerstrokewidth=0, markersize=3, alpha=0.55, color=colors[planet], label=PLANET_LABEL[planet])
    end
    xlims!(sp, -lim, lim)
    ylims!(sp, -lim, lim)
    out = joinpath(out_dir, "equilibrium_drag_vs_gravity_additive.pdf")
    savefig(sp, out)
    return out
end

function write_short_period_offset(additive::DataFrame, out_dir::String)
    rows = NamedTuple[]
    for group in groupby(additive, [:planet, :periapsis_regime])
        errs = finite_values(group.additive_error_m)
        full = finite_values(group.full_environment_delta_rp_m)
        isempty(errs) && continue
        push!(rows, (
            planet=String(first(group.planet)),
            periapsis_regime=String(first(group.periapsis_regime)),
            n=length(errs),
            mean_error_km=mean(errs) / 1e3,
            median_error_km=median(errs) / 1e3,
            rms_error_km=sqrt(mean(abs2, errs)) / 1e3,
            p95_abs_error_km=quantile(abs.(errs), 0.95) / 1e3,
            p95_abs_full_delta_km=isempty(full) ? NaN : quantile(abs.(full), 0.95) / 1e3,
        ))
    end
    df = DataFrame(rows)
    CSV.write(joinpath(out_dir, "short_period_additive_offset_summary.csv"), df)

    sp = scatter(
        Float64.(additive.full_environment_delta_rp_m) ./ 1e3,
        Float64.(additive.additive_total_delta_rp_m) ./ 1e3;
        group=additive.planet,
        xlabel="full-environment Delta hp (km/orbit)",
        ylabel="additive prediction Delta hp (km/orbit)",
        title="One-Orbit Additive Prediction Check",
        framestyle=:box,
        grid=true,
        markersize=2.5,
        markerstrokewidth=0,
        alpha=0.55,
        size=(760, 700),
    )
    vals = vcat(finite_values(additive.full_environment_delta_rp_m), finite_values(additive.additive_total_delta_rp_m)) ./ 1e3
    lim = isempty(vals) ? 1.0 : max(1.0, quantile(abs.(vals), 0.98))
    plot!(sp, [-lim, lim], [-lim, lim]; color=:black, linestyle=:dash, label="1:1")
    xlims!(sp, -lim, lim)
    ylims!(sp, -lim, lim)
    savefig(sp, joinpath(out_dir, "short_period_additive_prediction_check.pdf"))
    return df
end

function titan_zlk_estimate(results::DataFrame, out_dir::String)
    rows = NamedTuple[]
    planet = AnalyticalPerturbationModels._planet_cache("titan").planet
    nbody = AnalyticalPerturbationModels._planet_cache("titan").nbody_model
    saturn_idx = findfirst(==("Saturn"), String.(nbody.body_names))
    saturn_idx === nothing && return DataFrame()
    mu_saturn = Float64(nbody.body_mus[saturn_idx])
    mu_titan = AnalyticalPerturbationModels._planet_mu(planet)
    base = results[
        (results.planet .== "titan") .&
        map(x -> String(x) in ("shallow", "nominal", "deep"), results.periapsis_regime) .&
        (results.dynamics_case .== "full_environment") .&
        (results.density_case .== POSTER_DENSITY) .&
        (results.spacecraft_mass_scale .== DEFAULT_MASS_SCALE) .&
        (results.inclination_deg .== DEFAULT_INCLINATION_DEG) .&
        (results.argp_deg .== DEFAULT_ARGP_DEG),
        :,
    ]
    for row in eachrow(base)
        tbl = case_trajectory(row)
        tbl === nothing && continue
        info = (
            planet="titan",
            periapsis=String(row.periapsis_regime),
            apoapsis_alt_km=Float64(row.apoapsis_alt_km),
            spacecraft_mass_scale=Float64(row.spacecraft_mass_scale),
        )
        ctx = AnalyticalPerturbationModels._context(info, tbl)
        positions = AnalyticalPerturbationModels._third_body_positions(ctx, 0.0)
        a_out = norm(positions[saturn_idx])
        period_inner = 2 * pi * sqrt(ctx.semi_major_m^3 / mu_titan)
        period_outer = 2 * pi * sqrt(a_out^3 / (mu_saturn + mu_titan))
        e_outer = 0.0
        t_zlk = (mu_saturn + mu_titan) / mu_saturn * period_outer^2 / period_inner * (1.0 - e_outer^2)^(1.5)
        push!(rows, (
            periapsis_regime=String(row.periapsis_regime),
            apoapsis_alt_km=Float64(row.apoapsis_alt_km),
            initial_e=Float64(row.initial_e),
            inner_period_days=period_inner / 86400.0,
            saturn_titan_period_days=period_outer / 86400.0,
            titan_saturn_distance_km=a_out / 1e3,
            t_zlk_years=t_zlk / ZLK_SECONDS_PER_YEAR,
        ))
    end
    df = DataFrame(rows)
    CSV.write(joinpath(out_dir, "titan_zlk_timescale_estimate.csv"), df)
    return df
end

function write_readme(out_dir::String)
    text = """
    # IPPW poster generated products

    Generated by `benchmarks/studies/aerobraking_perturbation_mc/generate_ippw_poster_plots.jl`.

    Key files:

    - `pi_parameter_reference_values.csv` / `.pdf`: reference Pi values for the poster table.
    - `additive_perturbation_prediction.csv`: isolated-case additive decomposition used for overlays.
    - `regime_maps_delta_hp_2x2.pdf`: body-specific raising/stable/lowering regime maps with zero-contour overlays.
    - `delta_hp_vs_argp_reference.pdf`: reference Delta h_p(omega) curves with stable bands.
    - `full_survey_regime_fractions.csv` / `.pdf`: collapsed full-sweep raising/stable/lowering fractions.
    - `reference_corridor_summary.csv`: readable omega ranges for each reference Delta h_p slice.
    - `slice_atlas_omega_vs_*.pdf`: visual diagnostics for selecting clean omega-horizontal slices.
    - `slice_atlas_inclination_vs_*.pdf`: visual diagnostics for selecting clean inclination-horizontal slices.
    - `priority_slice_*.pdf`: cleaner poster-candidate slices for omega/apoapsis,
      omega/eccentricity, inclination/apoapsis, and omega/periapsis-altitude.
    - `phase_dependence_inclination_vs_apoapsis_all_planets_full_environment.pdf`:
      4x2 phase-dependence regime map across planets and selected AoP values.
    - `equilibrium_drag_vs_gravity_additive.pdf`: drag-vs-gravity equilibrium schematic populated from data.
    - `short_period_additive_offset_summary.csv` / `.pdf`: one-orbit additive-prediction offset check.
    - `titan_zlk_timescale_estimate.csv`: Titan ZLK period estimate from the reference geometry.

    Stable means `|Delta h_p| <= $(STABLE_EPSILON_KM) km/orbit`.
    The body-specific maps intentionally use different y-axes so each panel shows the most
    informative 2D slice through the 5D survey; the full-sweep fraction plot summarizes
    the collapsed 5D coverage.

    Note: the plotted additive zero contour is derived from isolated numerical perturbation cases
    and is not the final closed-form secular Xi equilibrium locus.
    """
    open(joinpath(out_dir, "README.md"), "w") do io
        write(io, text)
    end
end

function main(args=ARGS)
    run_dirs = isempty(args) ? [latest_run_dir()] : abspath.(collect(args))
    out_dir = ensure_out_dir(run_dirs[1])
    println("[ippw-poster] run_dirs = $(join(run_dirs, ", "))")
    println("[ippw-poster] out_dir = $out_dir")

    results = read_results(run_dirs)
    additive = additive_prediction_table(results)
    CSV.write(joinpath(out_dir, "additive_perturbation_prediction.csv"), additive)

    pi_df = build_pi_table(results, out_dir)
    CSV.write(joinpath(out_dir, "pi_parameter_reference_values.csv"), pi_df)

    plot_regime_maps(additive, out_dir)
    plot_delta_hp_vs_argp(additive, out_dir)
    write_reference_corridor_summary(additive, out_dir)
    plot_regime_fraction_summary(additive, out_dir)
    plot_slice_atlas(additive, out_dir)
    plot_priority_slices(additive, out_dir)
    plot_phase_inclination_apoapsis(additive, out_dir)
    plot_equilibrium_terms(additive, out_dir)
    write_short_period_offset(additive, out_dir)
    titan_zlk_estimate(results, out_dir)
    write_readme(out_dir)

    println("[ippw-poster] wrote poster products to $out_dir")
    return out_dir
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
