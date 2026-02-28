# include("../utils/Reference_system.jl")

# using .SimulationModel
using SatelliteToolbox
using StaticArrays
using LinearAlgebra
using Serialization
using ..SimulationModel: GRAM_LOCK
# export NoAtmosphereModel, ExponentialAtmosphereModel, GRAMAtmosphereModel, PolynomialFitAtmosphereModel, NRLMSISE00AtmosphereModel
# export getDensity
# using .AbstractTypes: AbstractPlanet, AbstractDensityModel

struct NoAtmosphereModel <: AbstractDensityModel
    # No additional fields needed for a simple no-atmosphere model
end

struct ExponentialAtmosphereModel <: AbstractDensityModel
    ρ_ref::Float64 # Reference density at reference height (kg/m³)
    h_ref::Float64 # Reference height (m)
    H::Float64     # Scale height (m)
end

struct GRAMAtmosphereModel{G, GA} <: AbstractDensityModel
    gram::G
    gram_atmosphere::GA
    gram_root::String
    gram_data_root::String
    spice_root::String
    planet_name::String
    initial_time::InitialTime
end

struct GRAMAtmosphereModelSurrugate{M} <: AbstractDensityModel
    base_model::M
    surrogate_file::String
    point_fallback_below_m::Union{Nothing, Float64}
end

const GRAMAtmosphereModelSurrogate = GRAMAtmosphereModelSurrugate

struct PolynomialFitAtmosphereModel <: AbstractDensityModel
    polyfit_coeffs::Vector{Float64} # Coefficients for the polynomial fit (ordered from highest degree to lowest)
end

struct NRLMSISE00AtmosphereModel <: AbstractDensityModel
    # No additional fields needed for NRLMSISE-00, but we can add any necessary parameters or configurations here if needed in the future
end

const _GRAM_WRAPPER = Ref{Any}(nothing)
const _GRAM_WRAPPER_FILE = Ref{String}("")
const _GRAM_SEED_WARNING_EMITTED = Ref(false)
const _GRAM_WIND_WARNING_EMITTED = Ref(false)
const _GRAM_LOCK_OFF_WARNING_EMITTED = Ref(false)
const _GRAM_STATIC_GRID_LOGGED = Ref(false)
const _GRAM_STATIC_GRID_CACHE = Dict{Any, Any}()
const _GRAM_STATIC_GRID_LOCK = ReentrantLock()
const _GRAM_STATIC_GRID_PREBUILD_IN_PROGRESS = Ref(false)
const _GRAM_SUPPORTED_PLANETS = ("earth", "mars", "venus", "titan", "jupiter", "uranus", "neptune")
const _GRAM_PLANET_DIR_NAMES = Dict{String, String}(
    "earth" => "Earth",
    "mars" => "Mars",
    "venus" => "Venus",
    "titan" => "Titan",
    "jupiter" => "Jupiter",
    "uranus" => "Uranus",
    "neptune" => "Neptune"
)
const _GRAM_OFFLINE_SURROGATE_CACHE = Dict{Any, Any}()
const _GRAM_OFFLINE_SURROGATE_LOCK = ReentrantLock()
const _GRAM_OFFLINE_SURROGATE_WARNED = Dict{String, Bool}()
const _GRAM_OFFLINE_SURROGATE_WARNED_LOCK = ReentrantLock()
const _GRAM_OFFLINE_SURROGATE_POINT_FALLBACK_BELOW_M_DEFAULT = Dict{String, Float64}(
    "titan" => 260_000.0
)

struct GRAMStaticGridKey
    planet_name::String
    alt_min_m::Float64
    alt_max_m::Float64
    n_alt::Int
    n_lat::Int
    n_lon::Int
    elapsed_time_s::Float64
    include_wind::Bool
end

struct GRAMStaticGrid
    key::GRAMStaticGridKey
    alt_nodes::Vector{Float64}
    lat_nodes::Vector{Float64}
    lon_nodes::Vector{Float64}
    rho::Array{Float64, 3}
    T::Array{Float64, 3}
    wind_e::Array{Float64, 3}
    wind_n::Array{Float64, 3}
    wind_u::Array{Float64, 3}
end

struct GRAMOfflineSurrogate
    planet_name::String
    source_file::String
    alt_nodes_m::Vector{Float64}
    lat_nodes_rad::Vector{Float64}
    lon_nodes_rad::Vector{Float64}
    rho::Array{Float64, 3}
    T::Array{Float64, 3}
    wind_e::Array{Float64, 3}
    wind_n::Array{Float64, 3}
    wind_u::Array{Float64, 3}
end

@inline function _parse_bool_env(name::String, default::Bool)::Bool
    raw = lowercase(strip(get(ENV, name, default ? "1" : "0")))
    if raw in ("1", "true", "yes", "on")
        return true
    elseif raw in ("0", "false", "no", "off")
        return false
    end
    throw(ArgumentError("Invalid $name='$raw'. Use one of: 1/0, true/false, yes/no, on/off."))
end

@inline function _parse_float_env(name::String, default::Float64)::Float64
    raw = strip(get(ENV, name, string(default)))
    parsed = try
        parse(Float64, raw)
    catch
        throw(ArgumentError("$name must be a floating-point value, got '$raw'"))
    end
    return parsed
end

@inline function _parse_int_env(name::String, default::Int)::Int
    raw = strip(get(ENV, name, string(default)))
    parsed = try
        parse(Int, raw)
    catch
        throw(ArgumentError("$name must be an integer value, got '$raw'"))
    end
    return parsed
end

@inline function _gram_warn_once(key::String, msg::String; kwargs...)
    emit = false
    lock(_GRAM_OFFLINE_SURROGATE_WARNED_LOCK) do
        if !get(_GRAM_OFFLINE_SURROGATE_WARNED, key, false)
            _GRAM_OFFLINE_SURROGATE_WARNED[key] = true
            emit = true
        end
    end
    emit && @warn msg kwargs...
    return nothing
end

@inline function Base.getproperty(model::GRAMAtmosphereModelSurrugate, name::Symbol)
    if name === :base_model || name === :surrogate_file || name === :point_fallback_below_m
        return getfield(model, name)
    end
    return getproperty(getfield(model, :base_model), name)
end

@inline function Base.propertynames(model::GRAMAtmosphereModelSurrugate, private::Bool=false)
    wrapped = propertynames(getfield(model, :base_model), private)
    return (:base_model, :surrogate_file, :point_fallback_below_m, wrapped...)
end

@inline _gram_static_grid_enabled() = _parse_bool_env("SPACEAGORA_GRAM_STATIC_GRID", false)
@inline _gram_static_grid_prebuild_all_planets_enabled() = _parse_bool_env("SPACEAGORA_GRAM_STATIC_GRID_PREBUILD_ALL_PLANETS", false)
@inline _gram_static_grid_prebuild_strict() = _parse_bool_env("SPACEAGORA_GRAM_STATIC_GRID_PREBUILD_STRICT", false)

@inline _spaceagora_repo_root() = normpath(joinpath(@__DIR__, "..", ".."))
@inline _gram_lib_extension() = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")
@inline function _gram_use_global_lock()::Bool
    mode = lowercase(strip(get(ENV, "SPACEAGORA_GRAM_GLOBAL_LOCK", "on")))
    if mode in ("on", "true", "1", "yes")
        return true
    elseif mode in ("off", "false", "0", "no")
        if Threads.nthreads() > 1
            if !_GRAM_LOCK_OFF_WARNING_EMITTED[]
                _GRAM_LOCK_OFF_WARNING_EMITTED[] = true
                @warn "SPACEAGORA_GRAM_GLOBAL_LOCK=off is unsafe with Threads.nthreads()>1. Keeping GRAM global lock enabled."
            end
            return true
        end
        return false
    end
    throw(ArgumentError("Unsupported SPACEAGORA_GRAM_GLOBAL_LOCK='$mode'. Use one of: on, off."))
end

@inline function _gram_normalize_planet_name(planet::AbstractString)::String
    key = lowercase(strip(String(planet)))
    key in _GRAM_SUPPORTED_PLANETS || throw(ArgumentError("Unsupported planet '$planet'. Supported planets: $(_GRAM_SUPPORTED_PLANETS)"))
    return key
end

@inline function _gram_planet_dir_name(planet::String)::String
    name = get(_GRAM_PLANET_DIR_NAMES, planet, "")
    isempty(name) && throw(ArgumentError("Unsupported GRAM planet '$planet' for surrogate path resolution."))
    return name
end

@inline function _gram_default_surrogate_file(planet::String)::String
    planet_dir = _gram_planet_dir_name(planet)
    preferred = joinpath(
        _spaceagora_repo_root(),
        "GRAM Suite 2.0",
        "simulation",
        "GRAM",
        "static_grids",
        "$(planet)_surrogate.jls"
    )
    if isfile(preferred)
        return preferred
    end
    # Compatibility fallback for per-planet copies.
    return joinpath(_spaceagora_repo_root(), "GRAM Suite 2.0", planet_dir, "$(planet)_surrogate.jls")
end

@inline function _gram_offline_surrogate_file(model::GRAMAtmosphereModelSurrugate)::String
    return model.surrogate_file
end

@inline function _gram_offline_surrogate_point_fallback_below_m(model::GRAMAtmosphereModel)::Union{Nothing, Float64}
    return get(_GRAM_OFFLINE_SURROGATE_POINT_FALLBACK_BELOW_M_DEFAULT, lowercase(strip(model.planet_name)), nothing)
end

@inline function _gram_offline_surrogate_point_fallback_below_m(model::GRAMAtmosphereModelSurrugate)::Union{Nothing, Float64}
    return model.point_fallback_below_m
end

@inline function _gram_parse_static_grid_planets(raw::AbstractString)::Vector{String}
    v = lowercase(strip(String(raw)))
    if isempty(v) || v in ("all", "*")
        return collect(_GRAM_SUPPORTED_PLANETS)
    end
    planets = String[]
    for token in split(v, ",")
        name = _gram_normalize_planet_name(token)
        name in planets || push!(planets, name)
    end
    isempty(planets) && throw(ArgumentError("No valid planets provided for static-grid prebuild."))
    return planets
end

@inline function _gram_static_grid_planets_from_env()::Vector{String}
    raw = get(ENV, "SPACEAGORA_GRAM_STATIC_GRID_PLANETS", "all")
    return _gram_parse_static_grid_planets(raw)
end

@inline function _gram_static_grid_key(model::GRAMAtmosphereModel, wind::Bool)::GRAMStaticGridKey
    alt_min_default = 0.0
    alt_max_default = 2_000_000.0
    alt_min_m = _parse_float_env("SPACEAGORA_GRAM_STATIC_GRID_ALT_MIN_M", alt_min_default)
    alt_max_m = _parse_float_env("SPACEAGORA_GRAM_STATIC_GRID_ALT_MAX_M", alt_max_default)
    if alt_max_m <= alt_min_m
        throw(ArgumentError("SPACEAGORA_GRAM_STATIC_GRID_ALT_MAX_M must be > ALT_MIN_M; got $alt_max_m <= $alt_min_m"))
    end

    n_alt = max(3, _parse_int_env("SPACEAGORA_GRAM_STATIC_GRID_NALT", 31))
    n_lat = max(3, _parse_int_env("SPACEAGORA_GRAM_STATIC_GRID_NLAT", 37))
    n_lon = max(8, _parse_int_env("SPACEAGORA_GRAM_STATIC_GRID_NLON", 73))
    elapsed_time_s = _parse_float_env("SPACEAGORA_GRAM_STATIC_GRID_ELAPSED_TIME_S", 0.0)
    include_wind = _parse_bool_env("SPACEAGORA_GRAM_STATIC_GRID_WIND", wind)

    return GRAMStaticGridKey(
        lowercase(strip(model.planet_name)),
        alt_min_m,
        alt_max_m,
        n_alt,
        n_lat,
        n_lon,
        elapsed_time_s,
        include_wind
    )
end

@inline function _gram_lerp(a::Float64, b::Float64, w::Float64)::Float64
    return a + w * (b - a)
end

@inline function _gram_trilerp(c000::Float64, c100::Float64, c010::Float64, c110::Float64, c001::Float64, c101::Float64, c011::Float64, c111::Float64, wa::Float64, wb::Float64, wc::Float64)::Float64
    c00 = _gram_lerp(c000, c100, wa)
    c10 = _gram_lerp(c010, c110, wa)
    c01 = _gram_lerp(c001, c101, wa)
    c11 = _gram_lerp(c011, c111, wa)
    c0 = _gram_lerp(c00, c10, wb)
    c1 = _gram_lerp(c01, c11, wb)
    return _gram_lerp(c0, c1, wc)
end

@inline function _gram_axis_segment(nodes::Vector{Float64}, x::Float64)::Tuple{Int, Int, Float64}
    n = length(nodes)
    n >= 2 || throw(ArgumentError("Grid axis must have at least 2 nodes."))
    xq = clamp(x, nodes[1], nodes[end])
    i0 = clamp(searchsortedlast(nodes, xq), 1, n - 1)
    i1 = i0 + 1
    x0 = nodes[i0]
    x1 = nodes[i1]
    w = x1 == x0 ? 0.0 : (xq - x0) / (x1 - x0)
    return i0, i1, clamp(w, 0.0, 1.0)
end

@inline function _gram_lon_segment(lon_nodes::Vector{Float64}, lon::Float64)::Tuple{Int, Int, Float64}
    n = length(lon_nodes)
    n >= 2 || throw(ArgumentError("Longitude grid must have at least 2 nodes."))
    period = 2pi
    lonq = mod(lon, period)
    lonq < 0 && (lonq += period)
    i0 = searchsortedlast(lon_nodes, lonq)
    i0 = i0 == 0 ? n : clamp(i0, 1, n)
    i1 = i0 == n ? 1 : i0 + 1
    lon0 = lon_nodes[i0]
    lon1 = i1 == 1 ? lon_nodes[1] + period : lon_nodes[i1]
    lonq_adj = i1 == 1 && lonq < lon0 ? lonq + period : lonq
    w = lon1 == lon0 ? 0.0 : (lonq_adj - lon0) / (lon1 - lon0)
    return i0, i1, clamp(w, 0.0, 1.0)
end

@inline function _gram_axis_segment_checked(nodes::Vector{Float64}, x::Float64)::Union{Nothing, Tuple{Int, Int, Float64}}
    n = length(nodes)
    n >= 2 || throw(ArgumentError("Grid axis must have at least 2 nodes."))
    tol = 1e-12 * max(abs(nodes[1]), abs(nodes[end]), 1.0)
    if x < nodes[1] - tol || x > nodes[end] + tol
        return nothing
    end
    xq = clamp(x, nodes[1], nodes[end])
    i0 = clamp(searchsortedlast(nodes, xq), 1, n - 1)
    i1 = i0 + 1
    x0 = nodes[i0]
    x1 = nodes[i1]
    w = x1 == x0 ? 0.0 : (xq - x0) / (x1 - x0)
    return i0, i1, clamp(w, 0.0, 1.0)
end

function _gram_static_grid_build(model::GRAMAtmosphereModel, key::GRAMStaticGridKey)::GRAMStaticGrid
    alt_nodes = collect(range(key.alt_min_m, key.alt_max_m, length=key.n_alt))
    lat_nodes = collect(range(-0.5pi, 0.5pi, length=key.n_lat))
    lon_nodes = collect(range(0.0, 2pi, length=key.n_lon + 1))[1:end-1]

    dims = (key.n_alt, key.n_lat, key.n_lon)
    rho = Array{Float64, 3}(undef, dims...)
    T = Array{Float64, 3}(undef, dims...)
    wind_e = Array{Float64, 3}(undef, dims...)
    wind_n = Array{Float64, 3}(undef, dims...)
    wind_u = Array{Float64, 3}(undef, dims...)

    t0_ns = time_ns()
    if !_GRAM_STATIC_GRID_LOGGED[]
        _GRAM_STATIC_GRID_LOGGED[] = true
        @warn "GRAM static grid interpolation enabled. Initial grid build can be expensive but avoids time-dependent GRAM calls afterwards."
    end

    @inbounds for ia in eachindex(alt_nodes)
        h = alt_nodes[ia]
        for ilat in eachindex(lat_nodes)
            lat = lat_nodes[ilat]
            for ilon in eachindex(lon_nodes)
                lon = lon_nodes[ilon]
                ρ, Ti, wv = if _gram_use_global_lock()
                    lock(GRAM_LOCK) do
                        _gram_density_state(model, h, lat, lon, key.elapsed_time_s, key.include_wind)
                    end
                else
                    _gram_density_state(model, h, lat, lon, key.elapsed_time_s, key.include_wind)
                end
                rho[ia, ilat, ilon] = ρ
                T[ia, ilat, ilon] = Ti
                wind_e[ia, ilat, ilon] = wv[1]
                wind_n[ia, ilat, ilon] = wv[2]
                wind_u[ia, ilat, ilon] = wv[3]
            end
        end
    end

    @info "Built GRAM static grid." planet=key.planet_name n_alt=key.n_alt n_lat=key.n_lat n_lon=key.n_lon elapsed_s=(time_ns() - t0_ns) * 1e-9

    return GRAMStaticGrid(key, alt_nodes, lat_nodes, lon_nodes, rho, T, wind_e, wind_n, wind_u)
end

function _gram_static_grid_get_or_build!(model::GRAMAtmosphereModel, wind::Bool)::GRAMStaticGrid
    key = _gram_static_grid_key(model, wind)
    lock(_GRAM_STATIC_GRID_LOCK) do
        if haskey(_GRAM_STATIC_GRID_CACHE, key)
            return _GRAM_STATIC_GRID_CACHE[key]::GRAMStaticGrid
        end
        grid = _gram_static_grid_build(model, key)
        _GRAM_STATIC_GRID_CACHE[key] = grid
        return grid
    end
end

function clear_gram_static_grid_cache!()
    lock(_GRAM_STATIC_GRID_LOCK) do
        empty!(_GRAM_STATIC_GRID_CACHE)
        _GRAM_STATIC_GRID_LOGGED[] = false
    end
    return nothing
end

function _gram_model_for_planet(base_model::GRAMAtmosphereModel, planet::String)::GRAMAtmosphereModel
    if planet == base_model.planet_name
        return base_model
    end
    return GRAMAtmosphereModel(
        gram_root_directory=base_model.gram_root,
        gram_data_directory=base_model.gram_data_root,
        spice_directory=base_model.spice_root,
        planet_name=planet,
        initial_time=base_model.initial_time
    )
end

function precompute_gram_static_grids!(
    base_model::GRAMAtmosphereModel;
    planets::Union{Nothing, AbstractVector{<:AbstractString}}=nothing,
    wind::Bool=true
)
    planet_list = planets === nothing ?
        _gram_static_grid_planets_from_env() :
        _gram_parse_static_grid_planets(join(planets, ","))

    if _GRAM_STATIC_GRID_PREBUILD_IN_PROGRESS[]
        return nothing
    end

    _GRAM_STATIC_GRID_PREBUILD_IN_PROGRESS[] = true
    try
        @info "Precomputing GRAM static grids." planets=planet_list
        strict = _gram_static_grid_prebuild_strict()
        failed = String[]
        for planet in planet_list
            try
                model = _gram_model_for_planet(base_model, planet)
                _gram_static_grid_get_or_build!(model, wind)
            catch err
                if strict
                    rethrow(err)
                end
                push!(failed, planet)
                @warn "Skipping GRAM static-grid prebuild for planet." planet error=sprint(showerror, err)
            end
        end
        isempty(failed) || @warn "GRAM static-grid prebuild completed with skipped planets." skipped=failed
    finally
        _GRAM_STATIC_GRID_PREBUILD_IN_PROGRESS[] = false
    end
    return nothing
end

@inline function _gram_static_grid_eval(grid::GRAMStaticGrid, h::Float64, lat::Float64, lon::Float64)::Tuple{Float64, Float64, SVector{3, Float64}}
    ia0, ia1, wa = _gram_axis_segment(grid.alt_nodes, h)
    ilat0, ilat1, wb = _gram_axis_segment(grid.lat_nodes, lat)
    ilon0, ilon1, wc = _gram_lon_segment(grid.lon_nodes, lon)

    ρ = _gram_trilerp(
        grid.rho[ia0, ilat0, ilon0], grid.rho[ia1, ilat0, ilon0],
        grid.rho[ia0, ilat1, ilon0], grid.rho[ia1, ilat1, ilon0],
        grid.rho[ia0, ilat0, ilon1], grid.rho[ia1, ilat0, ilon1],
        grid.rho[ia0, ilat1, ilon1], grid.rho[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    Ti = _gram_trilerp(
        grid.T[ia0, ilat0, ilon0], grid.T[ia1, ilat0, ilon0],
        grid.T[ia0, ilat1, ilon0], grid.T[ia1, ilat1, ilon0],
        grid.T[ia0, ilat0, ilon1], grid.T[ia1, ilat0, ilon1],
        grid.T[ia0, ilat1, ilon1], grid.T[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    wE = _gram_trilerp(
        grid.wind_e[ia0, ilat0, ilon0], grid.wind_e[ia1, ilat0, ilon0],
        grid.wind_e[ia0, ilat1, ilon0], grid.wind_e[ia1, ilat1, ilon0],
        grid.wind_e[ia0, ilat0, ilon1], grid.wind_e[ia1, ilat0, ilon1],
        grid.wind_e[ia0, ilat1, ilon1], grid.wind_e[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    wN = _gram_trilerp(
        grid.wind_n[ia0, ilat0, ilon0], grid.wind_n[ia1, ilat0, ilon0],
        grid.wind_n[ia0, ilat1, ilon0], grid.wind_n[ia1, ilat1, ilon0],
        grid.wind_n[ia0, ilat0, ilon1], grid.wind_n[ia1, ilat0, ilon1],
        grid.wind_n[ia0, ilat1, ilon1], grid.wind_n[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    wU = _gram_trilerp(
        grid.wind_u[ia0, ilat0, ilon0], grid.wind_u[ia1, ilat0, ilon0],
        grid.wind_u[ia0, ilat1, ilon0], grid.wind_u[ia1, ilat1, ilon0],
        grid.wind_u[ia0, ilat0, ilon1], grid.wind_u[ia1, ilat0, ilon1],
        grid.wind_u[ia0, ilat1, ilon1], grid.wind_u[ia1, ilat1, ilon1],
        wa, wb, wc
    )

    return ρ, Ti, SVector{3, Float64}(wE, wN, wU)
end

@inline function _gram_offline_surrogate_key(model::GRAMAtmosphereModelSurrugate)::NamedTuple{(:planet, :file), Tuple{String, String}}
    return (planet=lowercase(strip(model.planet_name)), file=_gram_offline_surrogate_file(model))
end

function _gram_load_offline_surrogate(file::String, planet::String)::GRAMOfflineSurrogate
    payload = open(file, "r") do io
        deserialize(io)
    end
    payload isa Dict || throw(ArgumentError("Expected Dict payload in '$file'."))
    dict = Dict{String, Any}(payload)

    get(dict, "status", "error") == "ok" || throw(ArgumentError("Payload status is not ok in '$file'."))
    get(dict, "type", "") == "surrogate_trilinear" || throw(ArgumentError("Unsupported payload type in '$file'."))

    payload_planet = lowercase(strip(String(get(dict, "planet", ""))))
    payload_planet == planet || throw(ArgumentError("Payload planet '$payload_planet' does not match expected '$planet' in '$file'."))

    grid = dict["grid"]
    fields = dict["fields"]
    alt_nodes_m = Float64.(Vector{Float64}(grid["alt_km"])) .* 1e3
    lat_nodes_rad = deg2rad.(Float64.(Vector{Float64}(grid["lat_deg"])))
    lon_nodes_rad = deg2rad.(Float64.(Vector{Float64}(grid["lon_deg"])))

    rho = Float64.(fields["density_kgm3"])
    T = Float64.(fields["temperature_K"])
    wind_e = Float64.(fields["wind_ew_ms"])
    wind_n = Float64.(fields["wind_ns_ms"])
    wind_u = Float64.(fields["wind_up_ms"])

    dims = size(rho)
    size(T) == dims || throw(ArgumentError("Temperature grid dimensions do not match density in '$file'."))
    size(wind_e) == dims || throw(ArgumentError("EW wind grid dimensions do not match density in '$file'."))
    size(wind_n) == dims || throw(ArgumentError("NS wind grid dimensions do not match density in '$file'."))
    size(wind_u) == dims || throw(ArgumentError("Vertical wind grid dimensions do not match density in '$file'."))
    length(alt_nodes_m) == dims[1] || throw(ArgumentError("Altitude axis size does not match fields in '$file'."))
    length(lat_nodes_rad) == dims[2] || throw(ArgumentError("Latitude axis size does not match fields in '$file'."))
    length(lon_nodes_rad) == dims[3] || throw(ArgumentError("Longitude axis size does not match fields in '$file'."))

    return GRAMOfflineSurrogate(payload_planet, file, alt_nodes_m, lat_nodes_rad, lon_nodes_rad, rho, T, wind_e, wind_n, wind_u)
end

function _gram_offline_surrogate_get_or_load!(model::GRAMAtmosphereModelSurrugate)::GRAMOfflineSurrogate
    key = _gram_offline_surrogate_key(model)
    if !isfile(key.file)
        throw(ArgumentError(
            "GRAM surrogate payload not found for planet='$(key.planet)': $(key.file). " *
            "Create/copy the file or pass `surrogate_file=...` to GRAMAtmosphereModelSurrugate."
        ))
    end

    lock(_GRAM_OFFLINE_SURROGATE_LOCK) do
        if haskey(_GRAM_OFFLINE_SURROGATE_CACHE, key)
            return _GRAM_OFFLINE_SURROGATE_CACHE[key]::GRAMOfflineSurrogate
        end
        surrogate = _gram_load_offline_surrogate(key.file, key.planet)
        _GRAM_OFFLINE_SURROGATE_CACHE[key] = surrogate
        @info "Loaded GRAM surrogate payload." planet=surrogate.planet_name file=surrogate.source_file
        return surrogate
    end
end

function _gram_offline_surrogate_eval(
    surrogate::GRAMOfflineSurrogate,
    h::Float64,
    lat::Float64,
    lon::Float64
)::Union{Nothing, Tuple{Float64, Float64, SVector{3, Float64}}}
    ia = _gram_axis_segment_checked(surrogate.alt_nodes_m, h)
    ia === nothing && return nothing
    ilat = _gram_axis_segment_checked(surrogate.lat_nodes_rad, lat)
    ilat === nothing && return nothing
    ilon0, ilon1, wc = _gram_lon_segment(surrogate.lon_nodes_rad, lon)

    ia0, ia1, wa = ia
    ilat0, ilat1, wb = ilat

    ρ = _gram_trilerp(
        surrogate.rho[ia0, ilat0, ilon0], surrogate.rho[ia1, ilat0, ilon0],
        surrogate.rho[ia0, ilat1, ilon0], surrogate.rho[ia1, ilat1, ilon0],
        surrogate.rho[ia0, ilat0, ilon1], surrogate.rho[ia1, ilat0, ilon1],
        surrogate.rho[ia0, ilat1, ilon1], surrogate.rho[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    Ti = _gram_trilerp(
        surrogate.T[ia0, ilat0, ilon0], surrogate.T[ia1, ilat0, ilon0],
        surrogate.T[ia0, ilat1, ilon0], surrogate.T[ia1, ilat1, ilon0],
        surrogate.T[ia0, ilat0, ilon1], surrogate.T[ia1, ilat0, ilon1],
        surrogate.T[ia0, ilat1, ilon1], surrogate.T[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    wE = _gram_trilerp(
        surrogate.wind_e[ia0, ilat0, ilon0], surrogate.wind_e[ia1, ilat0, ilon0],
        surrogate.wind_e[ia0, ilat1, ilon0], surrogate.wind_e[ia1, ilat1, ilon0],
        surrogate.wind_e[ia0, ilat0, ilon1], surrogate.wind_e[ia1, ilat0, ilon1],
        surrogate.wind_e[ia0, ilat1, ilon1], surrogate.wind_e[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    wN = _gram_trilerp(
        surrogate.wind_n[ia0, ilat0, ilon0], surrogate.wind_n[ia1, ilat0, ilon0],
        surrogate.wind_n[ia0, ilat1, ilon0], surrogate.wind_n[ia1, ilat1, ilon0],
        surrogate.wind_n[ia0, ilat0, ilon1], surrogate.wind_n[ia1, ilat0, ilon1],
        surrogate.wind_n[ia0, ilat1, ilon1], surrogate.wind_n[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    wU = _gram_trilerp(
        surrogate.wind_u[ia0, ilat0, ilon0], surrogate.wind_u[ia1, ilat0, ilon0],
        surrogate.wind_u[ia0, ilat1, ilon0], surrogate.wind_u[ia1, ilat1, ilon0],
        surrogate.wind_u[ia0, ilat0, ilon1], surrogate.wind_u[ia1, ilat0, ilon1],
        surrogate.wind_u[ia0, ilat1, ilon1], surrogate.wind_u[ia1, ilat1, ilon1],
        wa, wb, wc
    )
    return ρ, Ti, SVector{3, Float64}(wE, wN, wU)
end

function clear_gram_offline_surrogate_cache!()
    lock(_GRAM_OFFLINE_SURROGATE_LOCK) do
        empty!(_GRAM_OFFLINE_SURROGATE_CACHE)
    end
    lock(_GRAM_OFFLINE_SURROGATE_WARNED_LOCK) do
        empty!(_GRAM_OFFLINE_SURROGATE_WARNED)
    end
    return nothing
end

function _as_abspath(path::AbstractString)::String
    if isempty(path)
        return ""
    end
    return isabspath(path) ? normpath(path) : normpath(joinpath(_spaceagora_repo_root(), path))
end

@inline function _is_gram_root(path::AbstractString)::Bool
    return isdir(joinpath(path, "Build")) && isdir(joinpath(path, "Julia"))
end

function _resolve_gram_root(gram_root_directory::String, gram_directory::String)::String
    candidates = String[]

    !isempty(gram_root_directory) && push!(candidates, gram_root_directory)
    !isempty(get(ENV, "GRAM_ROOT", "")) && push!(candidates, ENV["GRAM_ROOT"])

    if !isempty(gram_directory)
        push!(candidates, gram_directory)
        push!(candidates, joinpath(gram_directory, ".."))
        push!(candidates, joinpath(gram_directory, "..", ".."))
    end

    push!(candidates, "GRAM Suite 2.0")

    for candidate in candidates
        candidate_abs = _as_abspath(candidate)
        _is_gram_root(candidate_abs) && return candidate_abs
    end

    throw(ArgumentError(
        "Unable to locate GRAM Suite root. Set `gram_root_directory` or `ENV[\"GRAM_ROOT\"]` to a path containing `Build/` and `Julia/`."
    ))
end

function _load_gram_wrapper!(gram_root::String)::Tuple{Any, Bool}
    wrapper_file = normpath(joinpath(gram_root, "Julia", "GRAM.jl"))
    isfile(wrapper_file) || throw(ArgumentError("GRAM Julia wrapper not found at: $wrapper_file"))

    loaded_now = false
    if _GRAM_WRAPPER[] === nothing
        Base.include(@__MODULE__, wrapper_file)
        _GRAM_WRAPPER[] = Base.invokelatest(getfield, @__MODULE__, :GRAM)
        _GRAM_WRAPPER_FILE[] = wrapper_file
        loaded_now = true
    elseif _GRAM_WRAPPER_FILE[] != wrapper_file
        throw(ArgumentError(
            "GRAM wrapper already loaded from $(_GRAM_WRAPPER_FILE[]). Restart Julia to load a different wrapper path."
        ))
    end

    return _GRAM_WRAPPER[], loaded_now
end

function _resolve_spice_directory(gram_root::String, spice_directory::String, gram_data_directory::String)::String
    candidates = String[]
    !isempty(spice_directory) && push!(candidates, spice_directory)
    push!(candidates, joinpath(gram_root, "SPICE"))
    !isempty(gram_data_directory) && push!(candidates, joinpath(gram_data_directory, "SPICE"))
    push!(candidates, "GRAM Suite 2.0/SPICE")
    push!(candidates, "GRAM_Data/SPICE")

    for candidate in candidates
        candidate_abs = _as_abspath(candidate)
        isdir(candidate_abs) && return candidate_abs
    end

    throw(ArgumentError("Unable to find a valid SPICE directory for GRAM initialization."))
end

function _resolve_gram_data_root(gram_root::String, gram_data_directory::String)::String
    candidates = String[]
    push!(candidates, gram_root)
    !isempty(gram_data_directory) && push!(candidates, gram_data_directory)
    push!(candidates, "GRAM Suite 2.0")
    push!(candidates, "GRAM_Data")

    for candidate in candidates
        candidate_abs = _as_abspath(candidate)
        if isdir(joinpath(candidate_abs, "Earth", "data")) || isdir(joinpath(candidate_abs, "Mars", "data"))
            return candidate_abs
        end
    end

    throw(ArgumentError("Unable to find GRAM data directories (expected `<root>/Earth/data` or `<root>/Mars/data`)."))
end

function _gram_body(gram, planet_name::String)::Int
    key = lowercase(strip(planet_name))
    key == "earth" && return gram.BODY_EARTH
    key == "mars" && return gram.BODY_MARS
    key == "venus" && return gram.BODY_VENUS
    key == "titan" && return gram.BODY_TITAN
    key == "jupiter" && return gram.BODY_JUPITER
    key == "uranus" && return gram.BODY_URANUS
    key == "neptune" && return gram.BODY_NEPTUNE
    throw(ArgumentError("Unsupported GRAM planet_name '$planet_name'."))
end

function _gram_data_path(gram_data_root::String, planet_name::String)::Union{Nothing, String}
    key = lowercase(strip(planet_name))
    if key == "earth" || key == "mars"
        planet_dir = uppercasefirst(key)
        data_path = joinpath(gram_data_root, planet_dir, "data")
        isdir(data_path) || throw(ArgumentError("GRAM data path not found: $data_path"))
        return data_path
    end
    return nothing
end

function GRAMAtmosphereModel(;
    gram_directory::String="GRAM Suite 2.0",
    gram_data_directory::String="GRAM Suite 2.0",
    gram_root_directory::String="",
    gram_library_path::String="",
    spice_directory::String="",
    planet_name::String="earth",
    seed::Int=1001,
    initial_time::InitialTime=InitialTime(),
    gram_min_relative_step_size::Union{Nothing, Real}=nothing,
    gram_perturbation_scales::Union{Nothing, Real, NTuple{4, Real}}=nothing,
    mars_map_year::Union{Nothing, Integer}=nothing,
    mars_mgcm_dust_levels::Union{Nothing, NTuple{3, Real}}=nothing,
    mars_dust_storm::Union{Nothing, NTuple{6, Real}}=nothing,
    mars_f107::Union{Nothing, Real}=nothing,
    mars_wind_scales::Union{Nothing, NTuple{2, Real}}=nothing,
    mars_mola_heights::Union{Nothing, Bool}=nothing,
    mars_min_max::Union{Nothing, Integer}=nothing
)
    gram_root = _resolve_gram_root(gram_root_directory, gram_directory)
    gram, loaded_now = _load_gram_wrapper!(gram_root)
    if loaded_now
        return Base.invokelatest(
            GRAMAtmosphereModel;
            gram_directory=gram_directory,
            gram_data_directory=gram_data_directory,
            gram_root_directory=gram_root_directory,
            gram_library_path=gram_library_path,
            spice_directory=spice_directory,
            planet_name=planet_name,
            seed=seed,
            initial_time=initial_time,
            gram_min_relative_step_size=gram_min_relative_step_size,
            gram_perturbation_scales=gram_perturbation_scales,
            mars_map_year=mars_map_year,
            mars_mgcm_dust_levels=mars_mgcm_dust_levels,
            mars_dust_storm=mars_dust_storm,
            mars_f107=mars_f107,
            mars_wind_scales=mars_wind_scales,
            mars_mola_heights=mars_mola_heights,
            mars_min_max=mars_min_max
        )
    end

    if !isempty(gram_library_path)
        gram.set_library!(_as_abspath(gram_library_path))
    else
        local_lib = joinpath(gram_root, "Build", "lib", "libGRAM.$(_gram_lib_extension())")
        if isfile(local_lib)
            gram.set_library!(local_lib)
        end
    end

    spice_root = _resolve_spice_directory(gram_root, spice_directory, gram_data_directory)
    gram.initialize!(spice_root)

    gram_data_root = _resolve_gram_data_root(gram_root, gram_data_directory)
    planet_key = lowercase(strip(planet_name))
    body = _gram_body(gram, planet_key)
    data_path = _gram_data_path(gram_data_root, planet_key)
    gram_atmosphere = data_path === nothing ? gram.create_atmosphere(body) : gram.create_atmosphere(body; data_path=data_path)

    gram.set_start_time!(
        gram_atmosphere;
        year=Int(initial_time.year),
        month=Int(initial_time.month),
        day=Int(initial_time.day),
        hour=Int(initial_time.hour),
        minute=Int(initial_time.minute),
        seconds=Float64(initial_time.second),
        scale=1,
        frame=1
    )

    if isdefined(gram, :set_seed!)
        gram.set_seed!(gram_atmosphere, seed)
    elseif seed != 1001 && !_GRAM_SEED_WARNING_EMITTED[]
        _GRAM_SEED_WARNING_EMITTED[] = true
        @warn "GRAMAtmosphereModel seed is ignored by the current Julia GRAM wrapper."
    end

    if gram_min_relative_step_size !== nothing && isdefined(gram, :set_min_relative_step_size!)
        gram.set_min_relative_step_size!(gram_atmosphere, Float64(gram_min_relative_step_size))
    end

    if gram_perturbation_scales !== nothing && isdefined(gram, :set_perturbation_scales!)
        if gram_perturbation_scales isa Real
            scale = Float64(gram_perturbation_scales)
            gram.set_perturbation_scales!(
                gram_atmosphere;
                density_scale=scale,
                ew_wind_scale=scale,
                ns_wind_scale=scale,
                vertical_wind_scale=scale
            )
        else
            s = gram_perturbation_scales
            gram.set_perturbation_scales!(
                gram_atmosphere;
                density_scale=Float64(s[1]),
                ew_wind_scale=Float64(s[2]),
                ns_wind_scale=Float64(s[3]),
                vertical_wind_scale=Float64(s[4])
            )
        end
    end

    if planet_key == "mars"
        if mars_map_year !== nothing && isdefined(gram, :set_map_year!)
            gram.set_map_year!(gram_atmosphere, Int(mars_map_year))
        end
        if mars_mgcm_dust_levels !== nothing && isdefined(gram, :set_mgcm_dust_levels!)
            gram.set_mgcm_dust_levels!(
                gram_atmosphere;
                constant_level=Float64(mars_mgcm_dust_levels[1]),
                min_level=Float64(mars_mgcm_dust_levels[2]),
                max_level=Float64(mars_mgcm_dust_levels[3])
            )
        end
        if mars_dust_storm !== nothing && isdefined(gram, :set_dust_storm!)
            ds = mars_dust_storm
            gram.set_dust_storm!(
                gram_atmosphere;
                longitude_sun=Float64(ds[1]),
                duration=Float64(ds[2]),
                intensity=Float64(ds[3]),
                max_radius=Float64(ds[4]),
                latitude=Float64(ds[5]),
                longitude=Float64(ds[6])
            )
        end
        if mars_f107 !== nothing && isdefined(gram, :set_f107!)
            gram.set_f107!(gram_atmosphere, Float64(mars_f107))
        end
        if mars_wind_scales !== nothing && isdefined(gram, :set_wind_scales!)
            gram.set_wind_scales!(
                gram_atmosphere;
                mean_winds=Float64(mars_wind_scales[1]),
                boundary_layer_winds=Float64(mars_wind_scales[2])
            )
        end
        if mars_mola_heights !== nothing && isdefined(gram, :set_mola_heights!)
            gram.set_mola_heights!(gram_atmosphere, mars_mola_heights)
        end
        if mars_min_max !== nothing && isdefined(gram, :set_min_max!)
            gram.set_min_max!(gram_atmosphere, Int(mars_min_max))
        end
    end

    model = GRAMAtmosphereModel(
        gram,
        gram_atmosphere,
        gram_root,
        gram_data_root,
        spice_root,
        planet_key,
        initial_time
    )
    if _gram_static_grid_enabled() && _gram_static_grid_prebuild_all_planets_enabled() && !_GRAM_STATIC_GRID_PREBUILD_IN_PROGRESS[]
        wind_enabled = _parse_bool_env("SPACEAGORA_GRAM_STATIC_GRID_WIND", true)
        precompute_gram_static_grids!(model; wind=wind_enabled)
    end
    return model
end

function GRAMAtmosphereModelSurrugate(;
    surrogate_file::String="",
    point_fallback_below_m::Union{Nothing, Real}=nothing,
    kwargs...
)
    base_model = GRAMAtmosphereModel(; kwargs...)

    file = isempty(strip(surrogate_file)) ?
        _gram_default_surrogate_file(base_model.planet_name) :
        _as_abspath(surrogate_file)
    point_fallback = if point_fallback_below_m === nothing
        _gram_offline_surrogate_point_fallback_below_m(base_model)
    else
        value = Float64(point_fallback_below_m)
        value >= 0.0 || throw(ArgumentError("point_fallback_below_m must be >= 0.0 m, got $value"))
        value
    end

    return GRAMAtmosphereModelSurrugate(base_model, file, point_fallback)
end

function Base.deepcopy_internal(model::GRAMAtmosphereModel, stackdict::IdDict)
    if haskey(stackdict, model)
        return stackdict[model]
    end

    copied = GRAMAtmosphereModel(
        gram_root_directory=model.gram_root,
        gram_data_directory=model.gram_data_root,
        spice_directory=model.spice_root,
        planet_name=model.planet_name,
        initial_time=model.initial_time
    )
    stackdict[model] = copied
    return copied
end

function Base.deepcopy_internal(model::GRAMAtmosphereModelSurrugate, stackdict::IdDict)
    if haskey(stackdict, model)
        return stackdict[model]
    end

    copied = GRAMAtmosphereModelSurrugate(
        deepcopy(model.base_model),
        model.surrogate_file,
        model.point_fallback_below_m
    )
    stackdict[model] = copied
    return copied
end
function interp(a, b, x)
    """

    """
    
    # check delta == diff b and a
    if (abs(b-a) > 20.0)
        if b <= 360.0 && b >= 350.0
            b = 360.0 - b
        elseif a <= 360.0 && a >= 350.0
            a = 360.0 - a
        end
    end

    value = x * (b - a) + a

    return value
end

function temperature_linear(h, p)
    """

    """

    # into atmosphere
    # if config.cnf.drag_state == true
    #     T = p.T
    # else
    #     T = p.T
    # end

    return p.T_ref
end

function getDensity(model::NoAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    # Return zero density and a default temperature for the no-atmosphere model
    ρ = 0.0
    T = p.args.environment_model.planet.T_ref
    wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    return ρ, T, wind_vec
end

"""
    Get the density using the Exponential atmosphere model.

    Parameters
    ----------
    model : ExponentialAtmosphereModel
        Model of the atmosphere.
    h : Float64
        Height in meters.
    lat : Float64
        Latitude in radians.
    lon : Float64
        Longitude in radians.
    el_time : Float64
        Elapsed time since the start of the simulation in seconds.
    wind : Bool
        Whether to include wind effects in the density calculation.

    Returns
    -------
    ρ : Float64
        Density in kg/m³.
    T : Float64
        Temperature in Kelvin.
    wind : SVector{3, Float64}
        Wind vector in m/s.
"""
function getDensity(model::ExponentialAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}
    """

    """

    ρ = p.ρ_ref * exp.((p.h_ref .- h)/p.H)

    T = temperature_linear(h, p)

    wind_vec = SVector{3, Float64}(0, 0, 0)

    return ρ, T, wind_vec
end

function getDensity(model::PolynomialFitAtmosphereModel, h::Union{Float64, Vector{Float64}}, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params

    if typeof(h) != Float64 || length(h) != 1
        # polyfit = model.polyfit_coeffs
        power = zeros(length(model.polyfit_coeffs),length(h))
        # Convert height from meters to kilometers
        h *= 1e-3
        # Calculate the polynomial value at height h
        for i=1:length(polyfit)
            power[i,:] = (h).^(length(polyfit)-i)
        end

        # Calculate the exponent term of the density using the polynomial coefficients
        exponent = zeros(length(h))
        @inbounds for j=eachindex(h)
            exponent[j] = sum(polyfit .* power[:,j])
        end

        # Calculate the density
        ρ = exp.(exponent)
        T = p.args.environment_model.planet.T_ref

        wind_vec = SVector{3, Float64}(0, 0, 0)

        return ρ, T, wind_vec
    else
        # polyfit = model.polyfit_coeffs
        power = MVector{length(model.polyfit_coeffs)}(zeros(length(model.polyfit_coeffs)))
        # Convert height from meters to kilometers
        h *= 1e-3
        # Calculate the polynomial value at height h
        @inbounds for i=eachindex(model.polyfit_coeffs)
            power[i] = h^(length(model.polyfit_coeffs)-i)
        end
        # Calculate the exponent term of the density using the polynomial coefficients
        exponent = sum(model.polyfit_coeffs .* power)
        # Calculate the density
        ρ = exp(exponent)
        T = p.args.environment_model.planet.T_ref

        wind_vec = SVector{3, Float64}(0, 0, 0)

        return ρ, T, wind_vec
    end
end

@inline function _gram_density_state(
    model::GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    # Resolve GRAM callables in latest world to avoid Julia 1.12 world-age binding warnings under hot reload.
    set_position! = Base.invokelatest(getfield, model.gram, Symbol("set_position!"))
    Base.invokelatest(
        set_position!,
        model.gram_atmosphere;
        height=h * 1e-3,
        latitude=rad2deg(lat),
        longitude=rad2deg(lon),
        elapsed_time=el_time
    )

    update! = Base.invokelatest(getfield, model.gram, Symbol("update!"))
    err = Base.invokelatest(update!, model.gram_atmosphere)
    if err != 0
        get_error_message = Base.invokelatest(getfield, model.gram, :get_error_message)
        throw(ErrorException("GRAM update failed (code=$err): $(Base.invokelatest(get_error_message))"))
    end

    get_dynamics_state = Base.invokelatest(getfield, model.gram, :get_dynamics_state)
    atmos = Base.invokelatest(get_dynamics_state, model.gram_atmosphere)
    rho_local = Float64(atmos.density)
    T_local = Float64(atmos.temperature)
    if isdefined(model.gram, :get_winds_state)
        get_winds_state = Base.invokelatest(getfield, model.gram, :get_winds_state)
        winds = Base.invokelatest(get_winds_state, model.gram_atmosphere)
        wind_vec_local = if wind
            SVector{3, Float64}(
                Float64(winds.perturbedEWWind),
                Float64(winds.perturbedNSWind),
                Float64(winds.perturbedVerticalWind)
            )
        else
            SVector{3, Float64}(
                Float64(winds.ewWind),
                Float64(winds.nsWind),
                Float64(winds.verticalWind)
            )
        end
        return rho_local, T_local, wind_vec_local
    end

    if !_GRAM_WIND_WARNING_EMITTED[]
        _GRAM_WIND_WARNING_EMITTED[] = true
        @warn "Loaded GRAM wrapper does not expose wind state. Returning zero wind."
    end
    return rho_local, T_local, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline function _gram_point_density(
    model::GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return if _gram_use_global_lock()
        lock(GRAM_LOCK) do
            _gram_density_state(model, h, lat, lon, el_time, wind)
        end
    else
        _gram_density_state(model, h, lat, lon, el_time, wind)
    end
end

@inline function _gram_point_density(
    model::GRAMAtmosphereModelSurrugate,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return _gram_point_density(model.base_model, h, lat, lon, el_time, wind)
end

function getDensity(model::GRAMAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    # cnf = params.cnf
    # println("In GRAM density model")
    # Check atmospheric interface conditions
    # planet_radius = p.args.environment_model.planet.Rp_e
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0
    if h > 2000.0e3
        rho = 0.0
        T = p.args.environment_model.planet.T_ref
        wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif _gram_static_grid_enabled()
        # Static piecewise-linear interpolation on a precomputed (alt, lat, lon) grid.
        # This mode intentionally ignores elapsed time to trade fidelity for throughput.
        grid = _gram_static_grid_get_or_build!(model, wind)
        grid_alt_max_m = last(grid.alt_nodes)
        if h > grid_alt_max_m
            rho = 0.0
            T = p.args.environment_model.planet.T_ref
            wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
            _gram_warn_once(
                "above_ceiling_zero_rho_static:$(grid.key.planet_name)",
                "State altitude is above GRAM static-grid ceiling. Returning rho=0 for this sample.",
                planet=grid.key.planet_name,
                h_m=h,
                grid_alt_max_m=grid_alt_max_m
            )
        else
            rho, T, wind_vec = _gram_static_grid_eval(grid, h, lat, lon)
        end
    elseif !drag_state && !p.args.mission_configuration.keplerian
        rho, T, wind_vec = density_polyfit(h, p)
    else
        rho, T, wind_vec = _gram_point_density(model, h, lat, lon, el_time, wind)
    end
    # println("el_time: ", el_time, "h: ", h, " lat: ", rad2deg(lat), " lon: ", rad2deg(lon), " rho: ", rho, " T: ", T, " wind_vec: ", wind_vec)
    # println("rho: ", rho)
    return rho, T, wind_vec
end

function getDensity(model::GRAMAtmosphereModelSurrugate, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0
    if h > 2000.0e3
        return 0.0, p.args.environment_model.planet.T_ref, SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif !drag_state && !p.args.mission_configuration.keplerian
        return density_polyfit(h, p)
    end

    surrogate = _gram_offline_surrogate_get_or_load!(model)
    point_fallback_below_m = _gram_offline_surrogate_point_fallback_below_m(model)
    if point_fallback_below_m !== nothing && h <= point_fallback_below_m
        _gram_warn_once(
            "point_fallback_below:$(surrogate.planet_name)",
            "State altitude is below configured surrogate point-fallback altitude. Using point-to-point GRAM for this sample.",
            planet=surrogate.planet_name,
            h_m=h,
            fallback_below_m=point_fallback_below_m
        )
        return _gram_point_density(model, h, lat, lon, el_time, wind)
    end

    s_eval = _gram_offline_surrogate_eval(surrogate, h, lat, lon)
    if s_eval !== nothing
        return s_eval
    end

    grid_alt_max_m = last(surrogate.alt_nodes_m)
    if h > grid_alt_max_m
        _gram_warn_once(
            "above_ceiling_zero_rho:$(surrogate.planet_name)",
            "State altitude is above GRAM surrogate ceiling. Returning rho=0 for this sample.",
            planet=surrogate.planet_name,
            h_m=h,
            grid_alt_max_m=grid_alt_max_m
        )
        return 0.0, p.args.environment_model.planet.T_ref, SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    _gram_warn_once(
        "outside_grid:$(surrogate.planet_name)",
        "State is outside GRAM surrogate grid. Falling back to point-to-point GRAM for this sample.",
        planet=surrogate.planet_name,
        h_m=h,
        lat_deg=rad2deg(lat),
        lon_deg=rad2deg(lon),
        grid_alt_min_m=first(surrogate.alt_nodes_m),
        grid_alt_max_m=grid_alt_max_m
    )
    return _gram_point_density(model, h, lat, lon, el_time, wind)
end

function getDensity(model::NRLMSISE00AtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}

    if config.cnf.drag_state == false && args[:keplerian] == false
        rho , T , wind_vec = density_exp(h, p)
        rho = 0.0
    elseif config.cnf.drag_state == true || args[:keplerian] == true
        jd = datetime2julian(current_time)
        atmo = SatelliteToolbox.AtmosphericModels.nrlmsise00(jd, h, lat, lon, 150, 150, 3)
        rho = atmo.total_density
        T = atmo.temperature
        wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    end

    return rho, T, wind_vec
end

function density_polyfit(h::Float64, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    # Load the polynomial coefficients for the planet
    coeffs = p.args.environment_model.planet.polyfit_coeffs
    poly_model = PolynomialFitAtmosphereModel(vec(coeffs))
    ρ, T, wind_vec = getDensity(poly_model, h, 0.0, 0.0, 0.0, false, p) # Latitude, longitude, elapsed time, and wind don't affect the density in this model
    return ρ, T, wind_vec
end
# end # module
