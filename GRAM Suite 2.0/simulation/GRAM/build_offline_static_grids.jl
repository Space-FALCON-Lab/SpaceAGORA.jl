#!/usr/bin/env julia

using Dates
using Serialization
using Printf
using Statistics

const GRAM_ROOT = normpath(dirname(dirname(@__DIR__)))
include(joinpath(GRAM_ROOT, "Julia", "GRAM.jl"))
using .GRAM

const SUPPORTED_PLANETS = ("earth", "mars", "venus", "jupiter", "uranus", "neptune", "titan")
const cfg_earth_merra2_path = Ref{Union{Nothing, String}}(nothing)

const PLANET_BODIES = Dict(
    "earth" => BODY_EARTH,
    "mars" => BODY_MARS,
    "venus" => BODY_VENUS,
    "jupiter" => BODY_JUPITER,
    "uranus" => BODY_URANUS,
    "neptune" => BODY_NEPTUNE,
    "titan" => BODY_TITAN
)

# Tuned defaults for smaller/faster surrogates while preserving low-alt fidelity.
const DEFAULT_SURROGATE_BANDS = Dict(
    "earth" => "0:100:1,100:200:2,200:300:4",
    "mars" => "0:100:1,100:200:2,200:300:4",
    "venus" => "0:100:1,100:200:2,200:300:4",
    "titan" => "0:100:1,100:200:2,200:300:4",
    "jupiter" => "0:200:2,200:600:5,600:2000:20,2000:10000:100",
    "uranus" => "0:200:2,200:600:5,600:2000:20,2000:10000:100",
    "neptune" => "0:200:2,200:600:5,600:2000:20,2000:10000:100"
)

function parse_cli(args::Vector{String})
    opts = Dict{String, String}()
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || error("Unsupported argument '$arg'. Use --key=value.")
        body = arg[3:end]
        if occursin("=", body)
            key, value = split(body, "=", limit=2)
            opts[key] = value
        else
            key = body
            if i == length(args) || startswith(args[i + 1], "--")
                opts[key] = "true"
            else
                i += 1
                opts[key] = args[i]
            end
        end
        i += 1
    end
    return opts
end

parse_bool(s::AbstractString) = lowercase(strip(String(s))) in ("1", "true", "yes", "on")
parse_int(s::AbstractString) = parse(Int, strip(String(s)))
parse_float(s::AbstractString) = parse(Float64, strip(String(s)))

function parse_planets(raw::AbstractString)
    txt = lowercase(strip(String(raw)))
    if isempty(txt) || txt in ("all", "*")
        return collect(SUPPORTED_PLANETS)
    end
    planets = String[]
    for token in split(txt, ",")
        p = lowercase(strip(token))
        p in SUPPORTED_PLANETS || error("Unsupported planet '$p'. Supported: $(SUPPORTED_PLANETS)")
        p in planets || push!(planets, p)
    end
    isempty(planets) && error("No valid planets parsed from '$raw'.")
    return planets
end

function parse_surrogate_source(raw::AbstractString)::Symbol
    x = lowercase(strip(String(raw)))
    if x in ("grid", "from-grid")
        return :grid
    elseif x in ("direct", "gram", "from-gram")
        return :direct
    end
    error("Unsupported --surrogate-source='$raw'. Use one of: grid, direct.")
end

function parse_surrogate_precision(raw::AbstractString)::DataType
    x = lowercase(strip(String(raw)))
    if x in ("float32", "f32", "single")
        return Float32
    elseif x in ("float64", "f64", "double")
        return Float64
    end
    error("Unsupported --surrogate-precision='$raw'. Use one of: float32, float64.")
end

function parse_altitude_bands(spec::AbstractString)
    txt = strip(String(spec))
    isempty(txt) && error("Altitude band spec cannot be empty.")
    bands = NTuple{3, Float64}[]
    for token in split(txt, ",")
        t = strip(token)
        isempty(t) && continue
        parts = split(t, ":")
        length(parts) == 3 || error("Invalid altitude band '$t'. Expected a:b:step (km).")
        a = parse_float(parts[1])
        b = parse_float(parts[2])
        step = parse_float(parts[3])
        step > 0.0 || error("Altitude band step must be > 0, got $step in '$t'.")
        b > a || error("Altitude band upper bound must be > lower bound in '$t'.")
        push!(bands, (a, b, step))
    end
    isempty(bands) && error("No valid altitude bands parsed from '$spec'.")
    return bands
end

function clip_altitude_bands(bands::Vector{NTuple{3, Float64}}, alt_min_km::Float64, alt_max_km::Float64)
    clipped = NTuple{3, Float64}[]
    for (a, b, step) in bands
        lo = max(a, alt_min_km)
        hi = min(b, alt_max_km)
        if hi > lo
            push!(clipped, (lo, hi, step))
        end
    end
    isempty(clipped) && push!(clipped, (alt_min_km, alt_max_km, max(1e-3, alt_max_km - alt_min_km)))
    sort!(clipped, by=x -> x[1])
    return clipped
end

function altitude_nodes_from_bands(
    bands::Vector{NTuple{3, Float64}},
    alt_min_km::Float64,
    alt_max_km::Float64
)
    clipped = clip_altitude_bands(bands, alt_min_km, alt_max_km)
    nodes = Float64[alt_min_km]
    for (a, b, step) in clipped
        x = max(a, alt_min_km)
        if abs(x - nodes[end]) > 1e-9
            push!(nodes, x)
        end
        while x + step < b - 1e-9
            x += step
            abs(x - nodes[end]) > 1e-9 && push!(nodes, x)
        end
        abs(b - nodes[end]) > 1e-9 && push!(nodes, b)
    end
    abs(alt_max_km - nodes[end]) > 1e-9 && push!(nodes, alt_max_km)

    sort!(nodes)
    dedup = Float64[]
    for x in nodes
        if isempty(dedup) || abs(x - dedup[end]) > 1e-9
            push!(dedup, x)
        end
    end
    return dedup
end

function linear_axis_nodes(min_v::Float64, max_v::Float64, step_v::Float64; periodic::Bool=false)
    step_v > 0.0 || error("Axis step must be > 0, got $step_v.")
    nodes = Float64[]
    x = min_v
    while x <= max_v + 1e-9
        push!(nodes, x)
        x += step_v
    end
    if periodic
        if !isempty(nodes) && abs(nodes[end] - max_v) <= 1e-9
            pop!(nodes)
        end
        isempty(nodes) && push!(nodes, min_v)
    else
        if isempty(nodes) || abs(nodes[end] - max_v) > 1e-9
            push!(nodes, max_v)
        end
    end
    return nodes
end

function bands_to_string(bands::Vector{NTuple{3, Float64}})
    parts = String[]
    for (a, b, s) in bands
        push!(parts, @sprintf("%.6g:%.6g:%.6g", a, b, s))
    end
    return join(parts, ",")
end

@inline _clamp01(x::Float64) = clamp(x, 0.0, 1.0)

@inline function _safe_quantile(v::Vector{Float64}, p::Float64)::Float64
    isempty(v) && return 0.0
    return quantile(v, _clamp01(p))
end

function bands_from_nodes(nodes::Vector{Float64}; rel_tol::Float64=1e-6)
    length(nodes) >= 2 || error("Need at least two altitude nodes to form bands.")
    b = NTuple{3, Float64}[]
    i0 = 1
    step_ref = nodes[2] - nodes[1]
    for i in 2:(length(nodes) - 1)
        step_i = nodes[i + 1] - nodes[i]
        same_step = abs(step_i - step_ref) <= rel_tol * max(1.0, abs(step_ref))
        if !same_step
            push!(b, (nodes[i0], nodes[i], step_ref))
            i0 = i
            step_ref = step_i
        end
    end
    push!(b, (nodes[i0], nodes[end], step_ref))
    return b
end

@inline function _interp1_clamped(xnodes::Vector{Float64}, ynodes::Vector{Float64}, x::Float64)::Float64
    n = length(xnodes)
    n == length(ynodes) || error("x/ y profile lengths mismatch.")
    n >= 2 || error("Need at least two profile nodes for interpolation.")
    xq = clamp(x, xnodes[1], xnodes[end])
    i0 = clamp(searchsortedlast(xnodes, xq), 1, n - 1)
    i1 = i0 + 1
    x0 = xnodes[i0]
    x1 = xnodes[i1]
    if x1 == x0
        return ynodes[i0]
    end
    w = clamp((xq - x0) / (x1 - x0), 0.0, 1.0)
    return ynodes[i0] + w * (ynodes[i1] - ynodes[i0])
end

function density_profile_from_grid_payload(grid_payload::Dict{String, Any})
    get(grid_payload, "status", "error") == "ok" || error("Grid payload is not usable (status=$(get(grid_payload, "status", "unknown"))).")
    grid = grid_payload["grid"]
    fields = grid_payload["fields"]
    alt_km = Vector{Float64}(grid["alt_km"])
    rho_grid = fields["density_kgm3"]

    nalt = length(alt_km)
    nalt >= 2 || error("Grid payload requires at least two altitude nodes.")
    logrho = Vector{Float64}(undef, nalt)
    @inbounds for ia in 1:nalt
        slice = rho_grid[ia, :, :]
        s = 0.0
        n = length(slice)
        for v in slice
            s += log(max(Float64(v), 1e-30))
        end
        logrho[ia] = s / max(1, n)
    end
    return alt_km, logrho
end

function density_profile_direct(
    root::String,
    planet::String,
    cfg;
    pilot_alt_step_km::Float64,
    pilot_lat_step_deg::Float64,
    pilot_lon_step_deg::Float64
)
    alt_km = linear_axis_nodes(cfg.alt_min_km, cfg.alt_max_km, pilot_alt_step_km; periodic=false)
    lat_deg = linear_axis_nodes(-90.0, 90.0, pilot_lat_step_deg; periodic=false)
    lon_deg = linear_axis_nodes(0.0, 360.0, pilot_lon_step_deg; periodic=true)

    @printf("[%s] Building adaptive pilot profile (%d x %d x %d)\n", planet, length(alt_km), length(lat_deg), length(lon_deg))
    logrho = Vector{Float64}(undef, length(alt_km))

    atmos = atmosphere_for_planet(root, planet)
    try
        configure_common!(atmos, cfg)
        configure_planet_specific!(planet, atmos)

        @inbounds for ia in eachindex(alt_km)
            if ia == 1 || ia == length(alt_km) || ia % max(1, fld(length(alt_km), 8)) == 0
                @printf("[%s]   adaptive pilot altitude slice %d / %d\n", planet, ia, length(alt_km))
            end
            h = alt_km[ia]
            s = 0.0
            n = 0
            for lat in lat_deg
                for lon in lon_deg
                    st = sample_state!(atmos, h, lat, lon, cfg.elapsed_time_s, planet)
                    s += log(max(st.rho, 1e-30))
                    n += 1
                end
            end
            logrho[ia] = s / max(1, n)
        end
    finally
        close!(atmos)
    end

    return alt_km, logrho
end

function adaptive_altitude_nodes(
    alt_profile_km::Vector{Float64},
    logrho_profile::Vector{Float64};
    min_step_km::Float64,
    max_step_km::Float64
)
    n = length(alt_profile_km)
    n == length(logrho_profile) || error("Adaptive profile length mismatch.")
    n >= 2 || error("Adaptive profile requires at least two points.")
    min_step_km > 0.0 || error("Adaptive min step must be > 0.")
    max_step_km >= min_step_km || error("Adaptive max step must be >= min step.")

    grads = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        if i == 1
            dh = max(1e-9, alt_profile_km[2] - alt_profile_km[1])
            grads[i] = abs((logrho_profile[2] - logrho_profile[1]) / dh)
        elseif i == n
            dh = max(1e-9, alt_profile_km[n] - alt_profile_km[n - 1])
            grads[i] = abs((logrho_profile[n] - logrho_profile[n - 1]) / dh)
        else
            dh = max(1e-9, alt_profile_km[i + 1] - alt_profile_km[i - 1])
            grads[i] = abs((logrho_profile[i + 1] - logrho_profile[i - 1]) / dh)
        end
    end

    rho_lo = _safe_quantile(logrho_profile, 0.15)
    rho_hi = _safe_quantile(logrho_profile, 0.85)
    grad_lo = _safe_quantile(grads, 0.15)
    grad_hi = _safe_quantile(grads, 0.85)
    (rho_hi - rho_lo) < 1e-12 && (rho_hi = rho_lo + 1.0)
    (grad_hi - grad_lo) < 1e-12 && (grad_hi = grad_lo + 1.0)

    score = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        s_rho = _clamp01((logrho_profile[i] - rho_lo) / (rho_hi - rho_lo))
        s_grad = _clamp01((grads[i] - grad_lo) / (grad_hi - grad_lo))
        score[i] = max(s_rho, s_grad)
    end

    alt_min_km = alt_profile_km[1]
    alt_max_km = alt_profile_km[end]
    nodes = Float64[alt_min_km]
    while nodes[end] < alt_max_km - 1e-9
        h = nodes[end]
        s = _interp1_clamped(alt_profile_km, score, h)
        dz = max_step_km - s * (max_step_km - min_step_km)
        dz = clamp(dz, min_step_km, max_step_km)
        nxt = min(h + dz, alt_max_km)
        if nxt <= h + 1e-10
            nxt = min(h + min_step_km, alt_max_km)
        end
        push!(nodes, nxt)
    end

    dedup = Float64[]
    for x in nodes
        if isempty(dedup) || abs(x - dedup[end]) > 1e-9
            push!(dedup, x)
        end
    end
    return dedup
end

function resolve_surrogate_bands(
    opts::Dict{String, String},
    planet::String,
    alt_min_km::Float64,
    alt_max_km::Float64
)
    key_specific = "surrogate-alt-bands-$planet"
    spec = if haskey(opts, key_specific)
        opts[key_specific]
    elseif haskey(opts, "surrogate-alt-bands")
        opts["surrogate-alt-bands"]
    else
        get(DEFAULT_SURROGATE_BANDS, planet, @sprintf("%.6g:%.6g:1", alt_min_km, alt_max_km))
    end
    parsed = parse_altitude_bands(spec)
    return clip_altitude_bands(parsed, alt_min_km, alt_max_km)
end

function data_path_for_planet(root::String, planet::String)
    if planet == "earth"
        return joinpath(root, "Earth", "data")
    elseif planet == "mars"
        return joinpath(root, "Mars", "data")
    end
    return nothing
end

function atmosphere_for_planet(root::String, planet::String)
    body = PLANET_BODIES[planet]
    data_path = data_path_for_planet(root, planet)
    if data_path === nothing
        return create_atmosphere(body)
    end
    return create_atmosphere(body; data_path=data_path)
end

function _find_merra2_file(path::AbstractString, month::Integer)
    file = joinpath(path, "All Mean", @sprintf("MERRA2All_%02d.bin", month))
    return isfile(file) ? file : nothing
end

function resolve_earth_merra2_path(root::String, month::Integer; override::Union{Nothing, String}=nothing)
    root_clean = normpath(root)
    if !isempty(root_clean) && (endswith(root_clean, "/") || endswith(root_clean, "\\"))
        root_clean = root_clean[1:end-1]
    end
    candidates = String[]
    if override !== nothing
        push!(candidates, normpath(override))
    end
    if haskey(ENV, "EARTH_MERRA2_PATH") && !isempty(strip(ENV["EARTH_MERRA2_PATH"]))
        push!(candidates, normpath(ENV["EARTH_MERRA2_PATH"]))
    end

    root_parent = dirname(root_clean)
    root_grandparent = dirname(root_parent)
    append!(candidates, [
        joinpath(root_clean, "Earth", "data", "MERRA2data"),
        joinpath(root_parent, "GRAM_Data", "Earth", "data", "MERRA2data"),
        joinpath(root_grandparent, "SpaceAGORA.jl-rollback_struct", "GRAM_Data", "Earth", "data", "MERRA2data"),
        joinpath(root_grandparent, "SpaceAGORA.jl-rollback_struct_init", "GRAM_Data", "Earth", "data", "MERRA2data")
    ])

    for cand in candidates
        isdir(cand) || continue
        _find_merra2_file(cand, month) !== nothing && return cand
    end

    for cand in candidates
        isdir(cand) || continue
        all_mean = joinpath(cand, "All Mean")
        isdir(all_mean) || continue
        any(occursin(r"MERRA2All_\d\d\.bin$", f) for f in readdir(all_mean; join=true)) && return cand
    end

    return nothing
end

function configure_common!(atmos::Atmosphere, cfg)
    set_start_time!(
        atmos;
        year=cfg.year,
        month=cfg.month,
        day=cfg.day,
        hour=cfg.hour,
        minute=cfg.minute,
        seconds=cfg.second,
        scale=1,
        frame=1
    )
    set_seed!(atmos, cfg.seed)
    set_perturbation_action!(atmos, false)
    set_perturbation_scales!(
        atmos;
        density_scale=0.0,
        ew_wind_scale=0.0,
        ns_wind_scale=0.0,
        vertical_wind_scale=0.0
    )
    return nothing
end

function configure_planet_specific!(planet::String, atmos::Atmosphere)
    if planet == "earth"
        if cfg_earth_merra2_path[] === nothing
            error("Earth MERRA2 data not found. Set --earth-merra2-path=<.../MERRA2data> or EARTH_MERRA2_PATH.")
        end
        earth_set_merra2_path!(atmos, cfg_earth_merra2_path[]::String)
    end
    if planet == "mars"
        if isdefined(GRAM, :set_mgcm_dust_levels!)
            set_mgcm_dust_levels!(atmos; constant_level=0.0, min_level=0.0, max_level=0.0)
        end
        if isdefined(GRAM, :set_dust_storm!)
            set_dust_storm!(
                atmos;
                longitude_sun=0.0,
                duration=0.0,
                intensity=0.0,
                max_radius=0.0,
                latitude=0.0,
                longitude=0.0
            )
        end
        if isdefined(GRAM, :set_dust_density!)
            set_dust_density!(atmos; nu=0.0, diameter=0.0, density=0.0)
        end
    end
    return nothing
end

@inline function sample_state!(
    atmos::Atmosphere,
    h_km::Float64,
    lat_deg::Float64,
    lon_deg::Float64,
    elapsed_time_s::Float64,
    planet::String
)
    set_position!(
        atmos;
        height=h_km,
        latitude=lat_deg,
        longitude=lon_deg,
        elapsed_time=elapsed_time_s
    )
    err = update!(atmos)
    err == 0 || error("GRAM update failed (planet=$planet, alt_km=$h_km, lat_deg=$lat_deg, lon_deg=$lon_deg): $(get_error_message())")

    dyn = get_dynamics_state(atmos)
    wnd = get_winds_state(atmos)
    return (
        rho=Float64(dyn.density),
        temp=Float64(dyn.temperature),
        wind_ew=Float64(wnd.ewWind),
        wind_ns=Float64(wnd.nsWind),
        wind_up=Float64(wnd.verticalWind)
    )
end

function build_planet_grid(root::String, planet::String, cfg)
    @printf("[%s] Building full grid (%d x %d x %d)\n", planet, cfg.nalt, cfg.nlat, cfg.nlon)
    atmos = atmosphere_for_planet(root, planet)
    try
        configure_common!(atmos, cfg)
        configure_planet_specific!(planet, atmos)

        alt_km = collect(range(cfg.alt_min_km, cfg.alt_max_km, length=cfg.nalt))
        lat_deg = collect(range(-90.0, 90.0, length=cfg.nlat))
        lon_deg = collect(range(0.0, 360.0, length=cfg.nlon + 1))[1:end-1]
        dims = (cfg.nalt, cfg.nlat, cfg.nlon)

        rho = Array{Float64, 3}(undef, dims...)
        temp = Array{Float64, 3}(undef, dims...)
        wind_ew = Array{Float64, 3}(undef, dims...)
        wind_ns = Array{Float64, 3}(undef, dims...)
        wind_up = Array{Float64, 3}(undef, dims...)

        t0_ns = time_ns()
        @inbounds for ia in eachindex(alt_km)
            if ia == 1 || ia == length(alt_km) || ia % max(1, fld(length(alt_km), 10)) == 0
                @printf("[%s]   full-grid altitude slice %d / %d\n", planet, ia, length(alt_km))
            end
            h = alt_km[ia]
            for ilat in eachindex(lat_deg)
                lat = lat_deg[ilat]
                for ilon in eachindex(lon_deg)
                    lon = lon_deg[ilon]
                    s = sample_state!(atmos, h, lat, lon, cfg.elapsed_time_s, planet)
                    rho[ia, ilat, ilon] = s.rho
                    temp[ia, ilat, ilon] = s.temp
                    wind_ew[ia, ilat, ilon] = s.wind_ew
                    wind_ns[ia, ilat, ilon] = s.wind_ns
                    wind_up[ia, ilat, ilon] = s.wind_up
                end
            end
        end

        elapsed_s = (time_ns() - t0_ns) * 1e-9
        @printf("[%s] Full grid done in %.3f s\n", planet, elapsed_s)

        return Dict{String, Any}(
            "status" => "ok",
            "format" => "spaceagora_gram_static_grid_v1",
            "type" => "full_grid",
            "planet" => planet,
            "created_at_utc" => string(now(UTC)),
            "elapsed_s" => elapsed_s,
            "wind_mode" => "on_mean",
            "monte_carlo" => "off",
            "dust" => "off",
            "elapsed_time_s" => cfg.elapsed_time_s,
            "seed" => cfg.seed,
            "grid" => Dict(
                "alt_km" => alt_km,
                "lat_deg" => lat_deg,
                "lon_deg" => lon_deg
            ),
            "fields" => Dict(
                "density_kgm3" => rho,
                "temperature_K" => temp,
                "wind_ew_ms" => wind_ew,
                "wind_ns_ms" => wind_ns,
                "wind_up_ms" => wind_up
            )
        )
    finally
        close!(atmos)
    end
end

@inline function _axis_segment(nodes::Vector{Float64}, x::Float64)
    n = length(nodes)
    n >= 2 || error("Axis must have at least two nodes.")
    xq = clamp(x, nodes[1], nodes[end])
    i0 = clamp(searchsortedlast(nodes, xq), 1, n - 1)
    i1 = i0 + 1
    x0 = nodes[i0]
    x1 = nodes[i1]
    w = x1 == x0 ? 0.0 : (xq - x0) / (x1 - x0)
    return i0, i1, clamp(w, 0.0, 1.0)
end

@inline function _lon_segment(nodes::Vector{Float64}, lon_deg::Float64)
    n = length(nodes)
    n >= 2 || error("Longitude axis must have at least two nodes.")
    period = 360.0
    xq = mod(lon_deg, period)
    xq < 0 && (xq += period)
    i0 = searchsortedlast(nodes, xq)
    i0 = i0 == 0 ? n : clamp(i0, 1, n)
    i1 = i0 == n ? 1 : i0 + 1
    x0 = nodes[i0]
    x1 = i1 == 1 ? nodes[1] + period : nodes[i1]
    xq_adj = i1 == 1 && xq < x0 ? xq + period : xq
    w = x1 == x0 ? 0.0 : (xq_adj - x0) / (x1 - x0)
    return i0, i1, clamp(w, 0.0, 1.0)
end

@inline _lerp(a::Float64, b::Float64, w::Float64) = a + w * (b - a)

@inline function _trilerp(
    c000::Float64, c100::Float64, c010::Float64, c110::Float64,
    c001::Float64, c101::Float64, c011::Float64, c111::Float64,
    wa::Float64, wb::Float64, wc::Float64
)
    c00 = _lerp(c000, c100, wa)
    c10 = _lerp(c010, c110, wa)
    c01 = _lerp(c001, c101, wa)
    c11 = _lerp(c011, c111, wa)
    c0 = _lerp(c00, c10, wb)
    c1 = _lerp(c01, c11, wb)
    return _lerp(c0, c1, wc)
end

function lookup_payload_state(payload::Dict{String, Any}, alt_km::Float64, lat_deg::Float64, lon_deg::Float64)
    get(payload, "status", "error") == "ok" || error("Input payload is not usable (status=$(get(payload, "status", "unknown"))).")
    grid = payload["grid"]
    fields = payload["fields"]
    alt_nodes = grid["alt_km"]::Vector{Float64}
    lat_nodes = grid["lat_deg"]::Vector{Float64}
    lon_nodes = grid["lon_deg"]::Vector{Float64}

    ia0, ia1, wa = _axis_segment(alt_nodes, alt_km)
    ilat0, ilat1, wb = _axis_segment(lat_nodes, lat_deg)
    ilon0, ilon1, wc = _lon_segment(lon_nodes, lon_deg)

    rho_grid = fields["density_kgm3"]::Array{Float64, 3}
    temp_grid = fields["temperature_K"]::Array{Float64, 3}
    w_ew_grid = fields["wind_ew_ms"]::Array{Float64, 3}
    w_ns_grid = fields["wind_ns_ms"]::Array{Float64, 3}
    w_up_grid = fields["wind_up_ms"]::Array{Float64, 3}

    rho = _trilerp(
        rho_grid[ia0, ilat0, ilon0], rho_grid[ia1, ilat0, ilon0],
        rho_grid[ia0, ilat1, ilon0], rho_grid[ia1, ilat1, ilon0],
        rho_grid[ia0, ilat0, ilon1], rho_grid[ia1, ilat0, ilon1],
        rho_grid[ia0, ilat1, ilon1], rho_grid[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    temp = _trilerp(
        temp_grid[ia0, ilat0, ilon0], temp_grid[ia1, ilat0, ilon0],
        temp_grid[ia0, ilat1, ilon0], temp_grid[ia1, ilat1, ilon0],
        temp_grid[ia0, ilat0, ilon1], temp_grid[ia1, ilat0, ilon1],
        temp_grid[ia0, ilat1, ilon1], temp_grid[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    w_ew = _trilerp(
        w_ew_grid[ia0, ilat0, ilon0], w_ew_grid[ia1, ilat0, ilon0],
        w_ew_grid[ia0, ilat1, ilon0], w_ew_grid[ia1, ilat1, ilon0],
        w_ew_grid[ia0, ilat0, ilon1], w_ew_grid[ia1, ilat0, ilon1],
        w_ew_grid[ia0, ilat1, ilon1], w_ew_grid[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    w_ns = _trilerp(
        w_ns_grid[ia0, ilat0, ilon0], w_ns_grid[ia1, ilat0, ilon0],
        w_ns_grid[ia0, ilat1, ilon0], w_ns_grid[ia1, ilat1, ilon0],
        w_ns_grid[ia0, ilat0, ilon1], w_ns_grid[ia1, ilat0, ilon1],
        w_ns_grid[ia0, ilat1, ilon1], w_ns_grid[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    w_up = _trilerp(
        w_up_grid[ia0, ilat0, ilon0], w_up_grid[ia1, ilat0, ilon0],
        w_up_grid[ia0, ilat1, ilon0], w_up_grid[ia1, ilat1, ilon0],
        w_up_grid[ia0, ilat0, ilon1], w_up_grid[ia1, ilat0, ilon1],
        w_up_grid[ia0, ilat1, ilon1], w_up_grid[ia1, ilat1, ilon1],
        wa, wb, wc
    )

    return rho, temp, w_ew, w_ns, w_up
end

function build_surrogate_from_grid_payload(
    grid_payload::Dict{String, Any},
    planet::String,
    cfg,
    s_cfg
)
    Tnum = s_cfg.precision
    alt_km = s_cfg.alt_nodes === nothing ? altitude_nodes_from_bands(s_cfg.alt_bands, cfg.alt_min_km, cfg.alt_max_km) : Vector{Float64}(s_cfg.alt_nodes)
    lat_deg = linear_axis_nodes(-90.0, 90.0, s_cfg.lat_step_deg; periodic=false)
    lon_deg = linear_axis_nodes(0.0, 360.0, s_cfg.lon_step_deg; periodic=true)
    dims = (length(alt_km), length(lat_deg), length(lon_deg))

    @printf("[%s] Building surrogate from grid (%d x %d x %d)\n", planet, dims[1], dims[2], dims[3])

    rho = Array{Tnum, 3}(undef, dims...)
    temp = Array{Tnum, 3}(undef, dims...)
    wind_ew = Array{Tnum, 3}(undef, dims...)
    wind_ns = Array{Tnum, 3}(undef, dims...)
    wind_up = Array{Tnum, 3}(undef, dims...)

    t0_ns = time_ns()
    @inbounds for ia in eachindex(alt_km)
        if ia == 1 || ia == length(alt_km) || ia % max(1, fld(length(alt_km), 10)) == 0
            @printf("[%s]   surrogate altitude slice %d / %d\n", planet, ia, length(alt_km))
        end
        h = alt_km[ia]
        for ilat in eachindex(lat_deg)
            lat = lat_deg[ilat]
            for ilon in eachindex(lon_deg)
                lon = lon_deg[ilon]
                r, t, we, wn, wu = lookup_payload_state(grid_payload, h, lat, lon)
                rho[ia, ilat, ilon] = Tnum(r)
                temp[ia, ilat, ilon] = Tnum(t)
                wind_ew[ia, ilat, ilon] = Tnum(we)
                wind_ns[ia, ilat, ilon] = Tnum(wn)
                wind_up[ia, ilat, ilon] = Tnum(wu)
            end
        end
    end
    elapsed_s = (time_ns() - t0_ns) * 1e-9
    @printf("[%s] Surrogate-from-grid done in %.3f s\n", planet, elapsed_s)
    alt_steps = length(alt_km) >= 2 ? diff(alt_km) : Float64[]

    return Dict{String, Any}(
        "status" => "ok",
        "format" => "spaceagora_gram_surrogate_trilinear_v1",
        "type" => "surrogate_trilinear",
        "source" => "grid",
        "planet" => planet,
        "created_at_utc" => string(now(UTC)),
        "elapsed_s" => elapsed_s,
        "wind_mode" => "on_mean",
        "monte_carlo" => "off",
        "dust" => "off",
        "elapsed_time_s" => cfg.elapsed_time_s,
        "seed" => cfg.seed,
        "surrogate" => Dict(
            "precision" => string(Tnum),
            "lat_step_deg" => s_cfg.lat_step_deg,
            "lon_step_deg" => s_cfg.lon_step_deg,
            "alt_mode" => s_cfg.alt_mode,
            "alt_bands_km" => bands_to_string(s_cfg.alt_bands),
            "alt_nodes_count" => length(alt_km),
            "alt_step_min_km" => isempty(alt_steps) ? 0.0 : minimum(alt_steps),
            "alt_step_max_km" => isempty(alt_steps) ? 0.0 : maximum(alt_steps)
        ),
        "grid" => Dict(
            "alt_km" => alt_km,
            "lat_deg" => lat_deg,
            "lon_deg" => lon_deg
        ),
        "fields" => Dict(
            "density_kgm3" => rho,
            "temperature_K" => temp,
            "wind_ew_ms" => wind_ew,
            "wind_ns_ms" => wind_ns,
            "wind_up_ms" => wind_up
        )
    )
end

function build_surrogate_direct(root::String, planet::String, cfg, s_cfg)
    Tnum = s_cfg.precision
    alt_km = s_cfg.alt_nodes === nothing ? altitude_nodes_from_bands(s_cfg.alt_bands, cfg.alt_min_km, cfg.alt_max_km) : Vector{Float64}(s_cfg.alt_nodes)
    lat_deg = linear_axis_nodes(-90.0, 90.0, s_cfg.lat_step_deg; periodic=false)
    lon_deg = linear_axis_nodes(0.0, 360.0, s_cfg.lon_step_deg; periodic=true)
    dims = (length(alt_km), length(lat_deg), length(lon_deg))

    @printf("[%s] Building surrogate direct from GRAM (%d x %d x %d)\n", planet, dims[1], dims[2], dims[3])

    rho = Array{Tnum, 3}(undef, dims...)
    temp = Array{Tnum, 3}(undef, dims...)
    wind_ew = Array{Tnum, 3}(undef, dims...)
    wind_ns = Array{Tnum, 3}(undef, dims...)
    wind_up = Array{Tnum, 3}(undef, dims...)

    atmos = atmosphere_for_planet(root, planet)
    try
        configure_common!(atmos, cfg)
        configure_planet_specific!(planet, atmos)
        t0_ns = time_ns()
        @inbounds for ia in eachindex(alt_km)
            if ia == 1 || ia == length(alt_km) || ia % max(1, fld(length(alt_km), 10)) == 0
                @printf("[%s]   surrogate-direct altitude slice %d / %d\n", planet, ia, length(alt_km))
            end
            h = alt_km[ia]
            for ilat in eachindex(lat_deg)
                lat = lat_deg[ilat]
                for ilon in eachindex(lon_deg)
                    lon = lon_deg[ilon]
                    s = sample_state!(atmos, h, lat, lon, cfg.elapsed_time_s, planet)
                    rho[ia, ilat, ilon] = Tnum(s.rho)
                    temp[ia, ilat, ilon] = Tnum(s.temp)
                    wind_ew[ia, ilat, ilon] = Tnum(s.wind_ew)
                    wind_ns[ia, ilat, ilon] = Tnum(s.wind_ns)
                    wind_up[ia, ilat, ilon] = Tnum(s.wind_up)
                end
            end
        end
        elapsed_s = (time_ns() - t0_ns) * 1e-9
        @printf("[%s] Surrogate-direct done in %.3f s\n", planet, elapsed_s)
        alt_steps = length(alt_km) >= 2 ? diff(alt_km) : Float64[]

        return Dict{String, Any}(
            "status" => "ok",
            "format" => "spaceagora_gram_surrogate_trilinear_v1",
            "type" => "surrogate_trilinear",
            "source" => "direct",
            "planet" => planet,
            "created_at_utc" => string(now(UTC)),
            "elapsed_s" => elapsed_s,
            "wind_mode" => "on_mean",
            "monte_carlo" => "off",
            "dust" => "off",
            "elapsed_time_s" => cfg.elapsed_time_s,
            "seed" => cfg.seed,
            "surrogate" => Dict(
                "precision" => string(Tnum),
                "lat_step_deg" => s_cfg.lat_step_deg,
                "lon_step_deg" => s_cfg.lon_step_deg,
                "alt_mode" => s_cfg.alt_mode,
                "alt_bands_km" => bands_to_string(s_cfg.alt_bands),
                "alt_nodes_count" => length(alt_km),
                "alt_step_min_km" => isempty(alt_steps) ? 0.0 : minimum(alt_steps),
                "alt_step_max_km" => isempty(alt_steps) ? 0.0 : maximum(alt_steps)
            ),
            "grid" => Dict(
                "alt_km" => alt_km,
                "lat_deg" => lat_deg,
                "lon_deg" => lon_deg
            ),
            "fields" => Dict(
                "density_kgm3" => rho,
                "temperature_K" => temp,
                "wind_ew_ms" => wind_ew,
                "wind_ns_ms" => wind_ns,
                "wind_up_ms" => wind_up
            )
        )
    finally
        close!(atmos)
    end
end

function load_payload(path::String)
    open(path, "r") do io
        payload = deserialize(io)
        payload isa Dict || error("Expected Dict payload in '$path'.")
        return Dict{String, Any}(payload)
    end
end

function save_payload(path::String, payload::Dict{String, Any})
    mkpath(dirname(path))
    open(path, "w") do io
        serialize(io, payload)
    end
    return path
end

toml_escape(s::AbstractString) = replace(String(s), "\\" => "\\\\", "\"" => "\\\"", "\n" => " ")

function save_summary(path::String, summary_rows::Vector{Dict{String, Any}}, cfg, opts_meta::Dict{String, Any}, planets)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "created_at_utc = \"$(now(UTC))\"")
        println(io, "format = \"spaceagora_gram_offline_products_v2\"")
        println(io, "wind = \"on\"")
        println(io, "monte_carlo = \"off\"")
        println(io, "dust = \"off\"")
        println(io, "elapsed_time_s = $(cfg.elapsed_time_s)")
        println(io, "full_grid_nalt = $(cfg.nalt)")
        println(io, "full_grid_nlat = $(cfg.nlat)")
        println(io, "full_grid_nlon = $(cfg.nlon)")
        println(io, "alt_min_km = $(cfg.alt_min_km)")
        println(io, "alt_max_km = $(cfg.alt_max_km)")
        println(io, "planets = [\"$(join(planets, "\", \""))\"]")
        println(io, "build_grid = $(opts_meta["build_grid"])")
        println(io, "build_surrogate = $(opts_meta["build_surrogate"])")
        println(io, "surrogate_only = $(opts_meta["surrogate_only"])")
        println(io, "surrogate_source = \"$(opts_meta["surrogate_source"])\"")
        println(io, "surrogate_precision = \"$(opts_meta["surrogate_precision"])\"")
        println(io, "surrogate_lat_step_deg = $(opts_meta["surrogate_lat_step_deg"])")
        println(io, "surrogate_lon_step_deg = $(opts_meta["surrogate_lon_step_deg"])")
        if haskey(opts_meta, "surrogate_adaptive_alt")
            println(io, "surrogate_adaptive_alt = $(opts_meta["surrogate_adaptive_alt"])")
            println(io, "surrogate_adaptive_min_step_km = $(opts_meta["surrogate_adaptive_min_step_km"])")
            println(io, "surrogate_adaptive_max_step_km = $(opts_meta["surrogate_adaptive_max_step_km"])")
            println(io, "surrogate_adaptive_pilot_alt_step_km = $(opts_meta["surrogate_adaptive_pilot_alt_step_km"])")
            println(io, "surrogate_adaptive_pilot_lat_step_deg = $(opts_meta["surrogate_adaptive_pilot_lat_step_deg"])")
            println(io, "surrogate_adaptive_pilot_lon_step_deg = $(opts_meta["surrogate_adaptive_pilot_lon_step_deg"])")
        end
        println(io)

        for row in summary_rows
            println(io, "[[planet]]")
            println(io, "name = \"$(toml_escape(row["planet"]))\"")
            println(io, "status = \"$(toml_escape(row["status"]))\"")
            println(io, "grid_status = \"$(toml_escape(get(row, "grid_status", "skipped")))\"")
            println(io, "surrogate_status = \"$(toml_escape(get(row, "surrogate_status", "skipped")))\"")
            if haskey(row, "grid_file")
                println(io, "grid_file = \"$(toml_escape(row["grid_file"]))\"")
            end
            if haskey(row, "grid_elapsed_s")
                println(io, "grid_elapsed_s = $(row["grid_elapsed_s"])")
            end
            if haskey(row, "surrogate_file")
                println(io, "surrogate_file = \"$(toml_escape(row["surrogate_file"]))\"")
            end
            if haskey(row, "surrogate_elapsed_s")
                println(io, "surrogate_elapsed_s = $(row["surrogate_elapsed_s"])")
            end
            if haskey(row, "surrogate_source")
                println(io, "surrogate_source = \"$(toml_escape(row["surrogate_source"]))\"")
            end
            if haskey(row, "surrogate_alt_bands_km")
                println(io, "surrogate_alt_bands_km = \"$(toml_escape(row["surrogate_alt_bands_km"]))\"")
            end
            if haskey(row, "surrogate_alt_mode")
                println(io, "surrogate_alt_mode = \"$(toml_escape(row["surrogate_alt_mode"]))\"")
            end
            if haskey(row, "surrogate_alt_nodes")
                println(io, "surrogate_alt_nodes = $(row["surrogate_alt_nodes"])")
            end
            if haskey(row, "surrogate_alt_step_min_km")
                println(io, "surrogate_alt_step_min_km = $(row["surrogate_alt_step_min_km"])")
            end
            if haskey(row, "surrogate_alt_step_max_km")
                println(io, "surrogate_alt_step_max_km = $(row["surrogate_alt_step_max_km"])")
            end
            if haskey(row, "surrogate_alt_profile_source")
                println(io, "surrogate_alt_profile_source = \"$(toml_escape(row["surrogate_alt_profile_source"]))\"")
            end
            if haskey(row, "error")
                println(io, "error = \"$(toml_escape(row["error"]))\"")
            end
            println(io)
        end
    end
    return path
end

function main()
    opts = parse_cli(copy(ARGS))

    planets = parse_planets(get(opts, "planets", "all"))
    out_dir = get(opts, "out-dir", joinpath(GRAM_ROOT, "simulation", "GRAM", "static_grids"))
    strict = parse_bool(get(opts, "strict", "false"))

    surrogate_only = parse_bool(get(opts, "surrogate-only", "false"))
    build_grid = parse_bool(get(opts, "build-grid", "true"))
    build_surrogate = parse_bool(get(opts, "build-surrogate", "false"))
    if surrogate_only
        build_grid = false
        build_surrogate = true
        haskey(opts, "surrogate-source") || (opts["surrogate-source"] = "direct")
    end

    cfg = (
        nalt=parse_int(get(opts, "nalt", "81")),
        nlat=parse_int(get(opts, "nlat", "73")),
        nlon=parse_int(get(opts, "nlon", "145")),
        alt_min_km=parse_float(get(opts, "alt-min-km", "0.0")),
        alt_max_km=parse_float(get(opts, "alt-max-km", "2000.0")),
        elapsed_time_s=parse_float(get(opts, "elapsed-time-s", "0.0")),
        year=parse_int(get(opts, "year", "2001")),
        month=parse_int(get(opts, "month", "11")),
        day=parse_int(get(opts, "day", "6")),
        hour=parse_int(get(opts, "hour", "19")),
        minute=parse_int(get(opts, "minute", "0")),
        second=parse_float(get(opts, "second", "32.0")),
        seed=parse_int(get(opts, "seed", "1001"))
    )

    s_source = parse_surrogate_source(get(opts, "surrogate-source", "grid"))
    s_precision = parse_surrogate_precision(get(opts, "surrogate-precision", "float32"))
    s_lat_step_deg = parse_float(get(opts, "surrogate-lat-step-deg", "1.75"))
    s_lon_step_deg = parse_float(get(opts, "surrogate-lon-step-deg", "1.75"))
    s_adaptive_alt = parse_bool(get(opts, "surrogate-adaptive-alt", "false"))
    s_adaptive_min_step_km = parse_float(get(opts, "surrogate-adaptive-min-step-km", "0.5"))
    s_adaptive_max_step_km = parse_float(get(opts, "surrogate-adaptive-max-step-km", "6.0"))
    s_adaptive_pilot_alt_step_km = parse_float(get(opts, "surrogate-adaptive-pilot-alt-step-km", "2.0"))
    s_adaptive_pilot_lat_step_deg = parse_float(get(opts, "surrogate-adaptive-pilot-lat-step-deg", "10.0"))
    s_adaptive_pilot_lon_step_deg = parse_float(get(opts, "surrogate-adaptive-pilot-lon-step-deg", "10.0"))
    s_dir = get(opts, "surrogate-dir", joinpath(out_dir, "surrogates"))

    earth_merra2_override = haskey(opts, "earth-merra2-path") ? normpath(opts["earth-merra2-path"]) : nothing
    earth_merra2_path = resolve_earth_merra2_path(GRAM_ROOT, cfg.month; override=earth_merra2_override)
    if earth_merra2_path === nothing
        @warn "Earth MERRA2 dataset not found. Earth products will fail unless --earth-merra2-path (or EARTH_MERRA2_PATH) is provided."
    else
        @info "Using Earth MERRA2 path." path=earth_merra2_path
    end
    cfg_earth_merra2_path[] = earth_merra2_path

    cfg.nalt >= 2 || error("--nalt must be >= 2")
    cfg.nlat >= 2 || error("--nlat must be >= 2")
    cfg.nlon >= 2 || error("--nlon must be >= 2")
    cfg.alt_max_km > cfg.alt_min_km || error("--alt-max-km must be greater than --alt-min-km")
    s_lat_step_deg > 0.0 || error("--surrogate-lat-step-deg must be > 0")
    s_lon_step_deg > 0.0 || error("--surrogate-lon-step-deg must be > 0")
    s_adaptive_min_step_km > 0.0 || error("--surrogate-adaptive-min-step-km must be > 0")
    s_adaptive_max_step_km >= s_adaptive_min_step_km || error("--surrogate-adaptive-max-step-km must be >= --surrogate-adaptive-min-step-km")
    s_adaptive_pilot_alt_step_km > 0.0 || error("--surrogate-adaptive-pilot-alt-step-km must be > 0")
    s_adaptive_pilot_lat_step_deg > 0.0 || error("--surrogate-adaptive-pilot-lat-step-deg must be > 0")
    s_adaptive_pilot_lon_step_deg > 0.0 || error("--surrogate-adaptive-pilot-lon-step-deg must be > 0")

    if !build_grid && !build_surrogate
        @warn "Both --build-grid=false and --build-surrogate=false. Nothing to do."
        return nothing
    end

    libext = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")
    libpath = get(opts, "lib", joinpath(GRAM_ROOT, "Build", "lib", "libGRAM.$libext"))
    spice_path = get(opts, "spice", joinpath(GRAM_ROOT, "SPICE"))
    summary_file = get(opts, "summary-file", "grid_build_summary.toml")

    mkpath(out_dir)
    build_surrogate && mkpath(s_dir)

    set_library!(libpath)
    initialize!(spice_path)

    summary_rows = Dict{String, Any}[]

    for planet in planets
        row = Dict{String, Any}(
            "planet" => planet,
            "status" => "ok",
            "grid_status" => build_grid ? "pending" : "skipped",
            "surrogate_status" => build_surrogate ? "pending" : "skipped"
        )
        grid_payload = nothing
        grid_file_path = joinpath(out_dir, "$(planet)_grid.jls")

        if build_grid
            try
                grid_payload = build_planet_grid(GRAM_ROOT, planet, cfg)
                save_payload(grid_file_path, grid_payload)
                row["grid_status"] = "ok"
                row["grid_file"] = basename(grid_file_path)
                row["grid_elapsed_s"] = grid_payload["elapsed_s"]
            catch err
                err_msg = sprint(showerror, err)
                row["status"] = "error"
                row["grid_status"] = "error"
                row["error"] = err_msg
                fail_payload = Dict{String, Any}(
                    "status" => "error",
                    "format" => "spaceagora_gram_static_grid_v1",
                    "type" => "full_grid",
                    "planet" => planet,
                    "created_at_utc" => string(now(UTC)),
                    "error" => err_msg
                )
                save_payload(grid_file_path, fail_payload)
                @warn "Failed to build full grid." planet error=err_msg
                strict && rethrow(err)
            end
        end

        if build_surrogate
            surrogate_file_path = joinpath(s_dir, "$(planet)_surrogate.jls")
            try
                surrogate_payload = nothing
                source_used = s_source
                src = nothing
                src_valid = false

                if s_source == :grid
                    src = grid_payload
                    if !(src isa Dict{String, Any} && get(src, "status", "error") == "ok")
                        if isfile(grid_file_path)
                            loaded = load_payload(grid_file_path)
                            if get(loaded, "status", "error") == "ok"
                                src = loaded
                            end
                        end
                    end
                    src_valid = src isa Dict{String, Any} && get(src, "status", "error") == "ok"
                end

                alt_mode = "bands"
                alt_profile_source = "none"
                if s_adaptive_alt
                    profile_alt = Float64[]
                    profile_logrho = Float64[]
                    if src_valid
                        profile_alt, profile_logrho = density_profile_from_grid_payload(src)
                        alt_profile_source = "grid"
                    else
                        profile_alt, profile_logrho = density_profile_direct(
                            GRAM_ROOT,
                            planet,
                            cfg;
                            pilot_alt_step_km=s_adaptive_pilot_alt_step_km,
                            pilot_lat_step_deg=s_adaptive_pilot_lat_step_deg,
                            pilot_lon_step_deg=s_adaptive_pilot_lon_step_deg
                        )
                        alt_profile_source = "direct"
                    end
                    s_alt_nodes = adaptive_altitude_nodes(
                        profile_alt,
                        profile_logrho;
                        min_step_km=s_adaptive_min_step_km,
                        max_step_km=s_adaptive_max_step_km
                    )
                    s_bands = bands_from_nodes(s_alt_nodes)
                    alt_mode = "adaptive_profile"
                else
                    s_bands = resolve_surrogate_bands(opts, planet, cfg.alt_min_km, cfg.alt_max_km)
                    s_alt_nodes = altitude_nodes_from_bands(s_bands, cfg.alt_min_km, cfg.alt_max_km)
                end

                s_cfg = (
                    alt_bands=s_bands,
                    alt_nodes=s_alt_nodes,
                    alt_mode=alt_mode,
                    lat_step_deg=s_lat_step_deg,
                    lon_step_deg=s_lon_step_deg,
                    precision=s_precision
                )

                if s_source == :grid
                    if src_valid
                        surrogate_payload = build_surrogate_from_grid_payload(src, planet, cfg, s_cfg)
                    else
                        @warn "Grid source unavailable for surrogate-from-grid; falling back to direct surrogate build." planet
                        surrogate_payload = build_surrogate_direct(GRAM_ROOT, planet, cfg, s_cfg)
                        source_used = :direct
                    end
                else
                    surrogate_payload = build_surrogate_direct(GRAM_ROOT, planet, cfg, s_cfg)
                end

                save_payload(surrogate_file_path, surrogate_payload)
                row["surrogate_status"] = "ok"
                row["surrogate_file"] = basename(surrogate_file_path)
                row["surrogate_elapsed_s"] = surrogate_payload["elapsed_s"]
                row["surrogate_source"] = String(source_used)
                row["surrogate_alt_bands_km"] = bands_to_string(s_bands)
                row["surrogate_alt_mode"] = alt_mode
                row["surrogate_alt_nodes"] = length(s_alt_nodes)
                if length(s_alt_nodes) >= 2
                    steps = diff(s_alt_nodes)
                    row["surrogate_alt_step_min_km"] = minimum(steps)
                    row["surrogate_alt_step_max_km"] = maximum(steps)
                end
                row["surrogate_alt_profile_source"] = alt_profile_source
            catch err
                err_msg = sprint(showerror, err)
                row["status"] = "error"
                row["surrogate_status"] = "error"
                row["error"] = haskey(row, "error") ? "$(row["error"]) | surrogate: $err_msg" : err_msg
                fail_payload = Dict{String, Any}(
                    "status" => "error",
                    "format" => "spaceagora_gram_surrogate_trilinear_v1",
                    "type" => "surrogate_trilinear",
                    "planet" => planet,
                    "created_at_utc" => string(now(UTC)),
                    "error" => err_msg
                )
                save_payload(surrogate_file_path, fail_payload)
                @warn "Failed to build surrogate product." planet error=err_msg
                strict && rethrow(err)
            end
        end

        push!(summary_rows, row)
    end

    opts_meta = Dict{String, Any}(
        "build_grid" => build_grid,
        "build_surrogate" => build_surrogate,
        "surrogate_only" => surrogate_only,
        "surrogate_source" => String(s_source),
        "surrogate_precision" => string(s_precision),
        "surrogate_lat_step_deg" => s_lat_step_deg,
        "surrogate_lon_step_deg" => s_lon_step_deg,
        "surrogate_adaptive_alt" => s_adaptive_alt,
        "surrogate_adaptive_min_step_km" => s_adaptive_min_step_km,
        "surrogate_adaptive_max_step_km" => s_adaptive_max_step_km,
        "surrogate_adaptive_pilot_alt_step_km" => s_adaptive_pilot_alt_step_km,
        "surrogate_adaptive_pilot_lat_step_deg" => s_adaptive_pilot_lat_step_deg,
        "surrogate_adaptive_pilot_lon_step_deg" => s_adaptive_pilot_lon_step_deg
    )

    summary_path = joinpath(out_dir, summary_file)
    save_summary(summary_path, summary_rows, cfg, opts_meta, planets)

    println("Wrote summary file: $summary_path")
    for row in summary_rows
        grid_status = get(row, "grid_status", "skipped")
        surr_status = get(row, "surrogate_status", "skipped")
        grid_file = get(row, "grid_file", "-")
        surr_file = get(row, "surrogate_file", "-")
        @printf("  %-8s status=%-5s grid=%-8s (%s) surrogate=%-8s (%s)\n",
            row["planet"], row["status"], grid_status, grid_file, surr_status, surr_file)
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
