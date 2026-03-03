const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const BUILDER_SCRIPT = joinpath(REPO_ROOT, "GRAMSuite.jl/GRAM Suite 2.0", "simulation", "GRAM", "build_offline_static_grids.jl")

using Statistics
using Printf
using Dates
using CSV
using DataFrames

# Reuse GRAM initialization/config + adaptive altitude helpers.
include(BUILDER_SCRIPT)

const PLANET_ALT_CAP_KM = Dict(
    "jupiter" => 1000.0
)

function parse_cli_sweep(args::Vector{String})
    opts = Dict{String, String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || error("Unsupported argument '$arg'. Use --key=value.")
        body = arg[3:end]
        if occursin("=", body)
            k, v = split(body, "=", limit=2)
            opts[k] = v
        else
            k = body
            if i == length(args) || startswith(args[i + 1], "--")
                opts[k] = "true"
            else
                i += 1
                opts[k] = args[i]
            end
        end
        i += 1
    end
    return opts
end

@inline parse_bool_sweep(x::AbstractString) = lowercase(strip(String(x))) in ("1", "true", "yes", "on")

function parse_float_list(raw::AbstractString)
    vals = Float64[]
    for t in split(String(raw), ",")
        s = strip(t)
        isempty(s) && continue
        push!(vals, parse(Float64, s))
    end
    isempty(vals) && error("Expected non-empty comma-separated float list, got '$raw'.")
    return vals
end

@inline function pctl(x::Vector{Float64}, p::Float64)
    isempty(x) && return 0.0
    xs = sort(x)
    idx = clamp(Int(ceil(clamp(p, 0.0, 1.0) * length(xs))), 1, length(xs))
    return xs[idx]
end

function make_cfg(opts::Dict{String, String})
    return (
        nalt=parse_int(get(opts, "nalt", "81")),
        nlat=parse_int(get(opts, "nlat", "73")),
        nlon=parse_int(get(opts, "nlon", "145")),
        alt_min_km=parse_float(get(opts, "alt-min-km", "0.0")),
        alt_max_km=parse_float(get(opts, "profile-max-alt-km", "600.0")),
        elapsed_time_s=parse_float(get(opts, "elapsed-time-s", "0.0")),
        year=parse_int(get(opts, "year", "2001")),
        month=parse_int(get(opts, "month", "1")),
        day=parse_int(get(opts, "day", "6")),
        hour=parse_int(get(opts, "hour", "19")),
        minute=parse_int(get(opts, "minute", "0")),
        second=parse_float(get(opts, "second", "32.0")),
        seed=parse_int(get(opts, "seed", "1001"))
    )
end

function build_rho_profile(root::String, planet::String, cfg; alt_step_km::Float64, lat_step_deg::Float64, lon_step_deg::Float64)
    alt_km = linear_axis_nodes(cfg.alt_min_km, cfg.alt_max_km, alt_step_km; periodic=false)
    lat_deg = linear_axis_nodes(-90.0, 90.0, lat_step_deg; periodic=false)
    lon_deg = linear_axis_nodes(0.0, 360.0, lon_step_deg; periodic=true)

    @printf("[%s] Building rho profile (%d x %d x %d)\n", planet, length(alt_km), length(lat_deg), length(lon_deg))
    mean_logrho = Vector{Float64}(undef, length(alt_km))
    p95_rho = Vector{Float64}(undef, length(alt_km))
    max_rho = Vector{Float64}(undef, length(alt_km))

    atmos = atmosphere_for_planet(root, planet)
    try
        configure_common!(atmos, cfg)
        configure_planet_specific!(planet, atmos)

        @inbounds for ia in eachindex(alt_km)
            if ia == 1 || ia == length(alt_km) || ia % max(1, fld(length(alt_km), 10)) == 0
                @printf("[%s]   profile altitude slice %d / %d\n", planet, ia, length(alt_km))
            end
            h = alt_km[ia]
            samples = Float64[]
            sizehint!(samples, length(lat_deg) * length(lon_deg))
            for lat in lat_deg
                for lon in lon_deg
                    st = sample_state!(atmos, h, lat, lon, cfg.elapsed_time_s, planet)
                    push!(samples, max(st.rho, 1e-30))
                end
            end
            mean_logrho[ia] = mean(log.(samples))
            p95_rho[ia] = pctl(samples, 0.95)
            max_rho[ia] = maximum(samples)
        end
    finally
        close!(atmos)
    end

    return (alt_km=alt_km, mean_logrho=mean_logrho, p95_rho=p95_rho, max_rho=max_rho)
end

function evaluate_adaptive_combo(ref_alt_km::Vector{Float64}, ref_mean_logrho::Vector{Float64}; alt_max_km::Float64, min_step_km::Float64, max_step_km::Float64)
    idx = findall(a -> a <= alt_max_km + 1e-9, ref_alt_km)
    length(idx) >= 2 || error("Not enough reference points up to alt_max_km=$alt_max_km")

    alt = ref_alt_km[idx]
    logrho = ref_mean_logrho[idx]

    nodes = adaptive_altitude_nodes(alt, logrho; min_step_km=min_step_km, max_step_km=max_step_km)
    node_logrho = [_interp1_clamped(alt, logrho, h) for h in nodes]
    interp_logrho = [_interp1_clamped(nodes, node_logrho, h) for h in alt]

    rho_ref = exp.(logrho)
    rho_hat = exp.(interp_logrho)
    abs_err = abs.(rho_hat .- rho_ref)
    rel_err = abs_err ./ max.(rho_ref, 1e-30)

    steps = diff(nodes)
    return (
        alt_max_km=alt_max_km,
        max_step_km=max_step_km,
        nalt_nodes=length(nodes),
        alt_step_min_km=isempty(steps) ? 0.0 : minimum(steps),
        alt_step_max_km=isempty(steps) ? 0.0 : maximum(steps),
        rho_abs_max=maximum(abs_err),
        rho_abs_p95=pctl(abs_err, 0.95),
        rho_rel_max=maximum(rel_err),
        rho_rel_p95=pctl(rel_err, 0.95),
        rho_rel_p99=pctl(rel_err, 0.99)
    )
end

function pick_recommendation(df::DataFrame; suggested_alt_km::Float64, rel_p95_limit::Float64)
    if nrow(df) == 0
        return nothing
    end
    ok = filter(r -> r.rho_rel_p95 <= rel_p95_limit, eachrow(df))
    if !isempty(ok)
        ok_df = DataFrame(ok)
        sort!(ok_df, [:alt_max_km, :max_step_km], rev=[false, true])
        candidates = filter(r -> r.alt_max_km >= suggested_alt_km, eachrow(ok_df))
        return isempty(candidates) ? ok_df[1, :] : DataFrame(candidates)[1, :]
    end
    fallback = deepcopy(df)
    sort!(fallback, [:rho_rel_p95, :nalt_nodes, :max_step_km], rev=[false, false, true])
    return fallback[1, :]
end

function run_sweep()
    opts = parse_cli_sweep(copy(ARGS))

    planet = lowercase(strip(get(opts, "planet", "mars")))
    planet in SUPPORTED_PLANETS || error("Unsupported --planet=$planet. Supported: $(SUPPORTED_PLANETS)")
    planet_alt_cap_km = haskey(opts, "planet-cap-km") ? parse_float(opts["planet-cap-km"]) : get(PLANET_ALT_CAP_KM, planet, Inf)

    alt_max_list_raw = parse_float_list(get(opts, "alt-max-list-km", "120,200,300,450,600"))
    max_step_list = parse_float_list(get(opts, "max-step-list-km", "4,6,8,10,12,16,20,24,30,40,50,60,80,100,125,150,175,200,225,250"))
    min_step_km = parse_float(get(opts, "min-step-km", "0.5"))

    profile_alt_step_km = parse_float(get(opts, "profile-alt-step-km", "2.0"))
    profile_lat_step_deg = parse_float(get(opts, "profile-lat-step-deg", "15.0"))
    profile_lon_step_deg = parse_float(get(opts, "profile-lon-step-deg", "15.0"))

    rho_floor_kgm3 = parse_float(get(opts, "rho-floor-kgm3", "1e-10"))
    rel_p95_limit = parse_float(get(opts, "rho-rel-p95-limit", "0.10"))

    out_dir = get(opts, "out-dir", joinpath(REPO_ROOT, "output", "gram_adaptive_sweep"))
    mkpath(out_dir)

    cfg = make_cfg(opts)
    if !haskey(opts, "alt-min-km") && planet == "earth"
        cfg = merge(cfg, (alt_min_km=5.0,))
    end
    if isfinite(planet_alt_cap_km)
        cfg = merge(cfg, (alt_max_km=min(cfg.alt_max_km, planet_alt_cap_km),))
    end

    alt_max_list = isfinite(planet_alt_cap_km) ? min.(alt_max_list_raw, planet_alt_cap_km) : alt_max_list_raw
    alt_max_list = sort(unique(filter(a -> a > cfg.alt_min_km + 1e-9, alt_max_list)))
    isempty(alt_max_list) && error("No valid --alt-max-list-km values after applying constraints.")

    max_requested_alt = maximum(alt_max_list)
    cfg_alt_max_km = max(cfg.alt_max_km, max_requested_alt)
    if isfinite(planet_alt_cap_km)
        cfg_alt_max_km = min(cfg_alt_max_km, planet_alt_cap_km)
    end
    cfg = merge(cfg, (alt_max_km=cfg_alt_max_km,))

    libext = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")
    libpath = get(opts, "lib", joinpath(GRAM_ROOT, "Build", "lib", "libGRAM.$libext"))
    spice_path = get(opts, "spice", joinpath(GRAM_ROOT, "SPICE"))

    earth_merra2_override = haskey(opts, "earth-merra2-path") ? normpath(opts["earth-merra2-path"]) : nothing
    earth_merra2_path = resolve_earth_merra2_path(GRAM_ROOT, cfg.month; override=earth_merra2_override)
    cfg_earth_merra2_path[] = earth_merra2_path

    set_library!(libpath)
    initialize!(spice_path)

    profile = build_rho_profile(
        GRAM_ROOT,
        planet,
        cfg;
        alt_step_km=profile_alt_step_km,
        lat_step_deg=profile_lat_step_deg,
        lon_step_deg=profile_lon_step_deg
    )

    idx_suggest = findlast(profile.p95_rho .>= rho_floor_kgm3)
    suggested_alt_raw_km = idx_suggest === nothing ? cfg.alt_min_km : profile.alt_km[idx_suggest]
    suggested_alt_km = isfinite(planet_alt_cap_km) ? min(suggested_alt_raw_km, planet_alt_cap_km) : suggested_alt_raw_km

    rows = NamedTuple[]
    for alt_max in sort(unique(alt_max_list))
        max_allowed_step_km = max(min_step_km, alt_max)
        max_step_candidates = sort(unique(min.(max_step_list, max_allowed_step_km)))
        for max_step in max_step_candidates
            max_step >= min_step_km || continue
            metrics = evaluate_adaptive_combo(
                profile.alt_km,
                profile.mean_logrho;
                alt_max_km=alt_max,
                min_step_km=min_step_km,
                max_step_km=max_step
            )
            push!(rows, merge(metrics, (
                planet=planet,
                rho_floor_kgm3=rho_floor_kgm3,
                rho_rel_p95_limit=rel_p95_limit,
                planet_alt_cap_km=planet_alt_cap_km,
                suggested_alt_raw_km=suggested_alt_raw_km,
                suggested_alt_km=suggested_alt_km,
                pass_rel_p95=metrics.rho_rel_p95 <= rel_p95_limit
            )))
        end
    end

    df = DataFrame(rows)
    sort!(df, [:alt_max_km, :max_step_km])

    rec = pick_recommendation(df; suggested_alt_km=suggested_alt_km, rel_p95_limit=rel_p95_limit)

    profile_df = DataFrame(
        planet=fill(planet, length(profile.alt_km)),
        alt_km=profile.alt_km,
        mean_logrho=profile.mean_logrho,
        rho_p95_kgm3=profile.p95_rho,
        rho_max_kgm3=profile.max_rho
    )

    run_tag = Dates.format(now(), "yyyymmdd_HHMMSS")
    summary_csv = joinpath(out_dir, "$(planet)_adaptive_altitude_rho_sweep_$(run_tag).csv")
    profile_csv = joinpath(out_dir, "$(planet)_rho_profile_$(run_tag).csv")
    CSV.write(summary_csv, df)
    CSV.write(profile_csv, profile_df)

    println("\nAdaptive rho sweep summary ($planet):")
    show(df, allrows=true, allcols=true)

    println("\n\nDensity profile suggestion:")
    println(@sprintf("  rho_floor = %.3e kg/m^3", rho_floor_kgm3))
    if isfinite(planet_alt_cap_km)
        println(@sprintf("  planet_alt_cap_km = %.2f", planet_alt_cap_km))
        println(@sprintf("  suggested_alt_raw_km = %.2f", suggested_alt_raw_km))
    end
    println(@sprintf("  suggested_alt_km = %.2f", suggested_alt_km))

    if rec !== nothing
        println("\nRecommended combo:")
        println(@sprintf("  alt_max_km=%.2f max_step_km=%.2f min_step_km=%.2f", rec.alt_max_km, rec.max_step_km, min_step_km))
        println(@sprintf("  nalt_nodes=%d rho_rel_p95=%.4f rho_rel_max=%.4f", rec.nalt_nodes, rec.rho_rel_p95, rec.rho_rel_max))
    end

    println("\nSaved:")
    println("  $summary_csv")
    println("  $profile_csv")

    return (summary=df, profile=profile_df, recommendation=rec, summary_csv=summary_csv, profile_csv=profile_csv)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_sweep()
end
