module GRAMSuite

using LinearAlgebra
using Serialization
using StaticArrays

export InitialTime
export GRAMAtmosphereModel, GRAMAtmosphereModelSurrogate
export point_density_state, density_state, surrogate_density_state
export precompute_gram_static_grids!, clear_gram_static_grid_cache!, clear_gram_offline_surrogate_cache!
export gram_default_surrogate_file, gram_default_point_fallback_below_m
export gram_use_global_lock, gram_static_grid_enabled, resolve_path

Base.@kwdef struct InitialTime
    year::Int32 = 2000
    month::Int16 = 1
    day::Int16 = 1
    hour::Int16 = 0
    minute::Int16 = 0
    second::Float32 = 0.0
end

struct GRAMAtmosphereModel{G, GA}
    gram::G
    gram_atmosphere::GA
    gram_root::String
    gram_data_root::String
    spice_root::String
    planet_name::String
    initial_time::InitialTime
    offline_surrogate_supported::Bool
    offline_surrogate_unsupported_reason::String
end

struct GRAMAtmosphereModelSurrogate{M}
    base_model::M
    surrogate_file::String
    point_fallback_below_m::Union{Nothing, Float64}
end

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
    wind_mode::String
    monte_carlo::String
    dust::String
    alt_nodes_m::Vector{Float64}
    lat_nodes_rad::Vector{Float64}
    lon_nodes_rad::Vector{Float64}
    rho::Array{Float64, 3}
    T::Array{Float64, 3}
    wind_e::Array{Float64, 3}
    wind_n::Array{Float64, 3}
    wind_u::Array{Float64, 3}
end

const _GRAM_WRAPPER = Ref{Any}(nothing)
const _GRAM_WRAPPER_FILE = Ref{String}("")
const _GRAM_SEED_WARNING_EMITTED = Ref(false)
const _GRAM_WIND_WARNING_EMITTED = Ref(false)
const _GRAM_NONFINITE_WIND_WARNING_EMITTED = Ref(false)
const _GRAM_LOCK_OFF_WARNING_EMITTED = Ref(false)
const _GRAM_STATIC_GRID_LOGGED = Ref(false)
const _GRAM_STATIC_GRID_CACHE = Dict{Any, Any}()
const _GRAM_STATIC_GRID_LOCK = ReentrantLock()
const _GRAM_STATIC_GRID_PREBUILD_IN_PROGRESS = Ref(false)
const _GRAM_OFFLINE_SURROGATE_CACHE = Dict{Any, Any}()
const _GRAM_OFFLINE_SURROGATE_LOCK = ReentrantLock()
const _GRAM_OFFLINE_SURROGATE_LOGGED = Ref(false)
const _GRAM_OFFLINE_SURROGATE_WARNED = Dict{String, Bool}()
const _GRAM_OFFLINE_SURROGATE_WARNED_LOCK = ReentrantLock()
const _GRAM_INTERNAL_LOCK = ReentrantLock()

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
const _GRAM_OFFLINE_SURROGATE_PROFILE = "p175_mid_all_planets"
const _GRAM_OFFLINE_SURROGATE_POINT_FALLBACK_BELOW_M_DEFAULT = Dict{String, Float64}(
    "titan" => 260_000.0
)
const _GRAM_FROZEN_PLANET_ALT_RANGE_M = Dict{String, Tuple{Float64, Float64}}(
    "earth" => (5_000.0, 1_115_000.0),
    "mars" => (0.0, 365_000.0),
    "venus" => (0.0, 460_000.0),
    "jupiter" => (0.0, 1_000_000.0),
    "titan" => (0.0, 2_500_000.0),
    "neptune" => (0.0, 4_000_000.0),
    "uranus" => (0.0, 7_000_000.0)
)

@inline _package_root() = normpath(joinpath(@__DIR__, ".."))
@inline _workspace_root() = normpath(joinpath(_package_root(), ".."))
@inline _gram_lib_extension() = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")

function resolve_path(path::AbstractString)::String
    if isempty(path)
        return ""
    end
    if isabspath(path)
        return normpath(path)
    end
    return normpath(joinpath(_workspace_root(), path))
end

@inline function _resolve_path_candidates(path::AbstractString)
    if isempty(path)
        return String[]
    end
    if isabspath(path)
        return String[normpath(path)]
    end
    return String[
        normpath(joinpath(_workspace_root(), path)),
        normpath(joinpath(_package_root(), path)),
        normpath(joinpath(pwd(), path))
    ]
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

@inline function gram_use_global_lock()::Bool
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

@inline gram_static_grid_enabled() = _parse_bool_env("SPACEAGORA_GRAM_STATIC_GRID", false)
@inline _gram_static_grid_prebuild_all_planets_enabled() = _parse_bool_env("SPACEAGORA_GRAM_STATIC_GRID_PREBUILD_ALL_PLANETS", false)
@inline _gram_static_grid_prebuild_strict() = _parse_bool_env("SPACEAGORA_GRAM_STATIC_GRID_PREBUILD_STRICT", false)

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

@inline function _coerce_initial_time(initial_time)::InitialTime
    if initial_time isa InitialTime
        return initial_time
    end
    return InitialTime(
        year=Int32(getproperty(initial_time, :year)),
        month=Int16(getproperty(initial_time, :month)),
        day=Int16(getproperty(initial_time, :day)),
        hour=Int16(getproperty(initial_time, :hour)),
        minute=Int16(getproperty(initial_time, :minute)),
        second=Float32(getproperty(initial_time, :second))
    )
end

@inline function _with_gram_lock(f::F, lock_obj) where {F}
    if gram_use_global_lock()
        if lock_obj === nothing
            return lock(_GRAM_INTERNAL_LOCK) do
                f()
            end
        end
        return lock(lock_obj) do
            f()
        end
    end
    return f()
end

@inline function Base.getproperty(model::GRAMAtmosphereModelSurrogate, name::Symbol)
    if name === :base_model || name === :surrogate_file || name === :point_fallback_below_m
        return getfield(model, name)
    end
    return getproperty(getfield(model, :base_model), name)
end

@inline function Base.propertynames(model::GRAMAtmosphereModelSurrogate, private::Bool=false)
    wrapped = propertynames(getfield(model, :base_model), private)
    return (:base_model, :surrogate_file, :point_fallback_below_m, wrapped...)
end

@inline function _gram_wind_mode()::Symbol
    raw = lowercase(strip(get(ENV, "SPACEAGORA_GRAM_WIND_MODE", "auto")))
    if raw in ("perturbed", "pert", "stochastic")
        return :perturbed
    elseif raw in ("nominal", "mean", "deterministic", "base")
        return :nominal
    elseif raw == "auto"
        return :perturbed
    end
    throw(ArgumentError(
        "Unsupported SPACEAGORA_GRAM_WIND_MODE='$raw'. Use one of: auto, nominal, perturbed."
    ))
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

@inline function _is_gram_root(path::AbstractString)::Bool
    return isdir(joinpath(path, "Build")) && isdir(joinpath(path, "Julia"))
end

function _resolve_gram_root(gram_root_directory::String, gram_directory::String)::String
    candidates = String[]

    !isempty(gram_root_directory) && append!(candidates, _resolve_path_candidates(gram_root_directory))
    env_root = strip(get(ENV, "GRAM_ROOT", ""))
    !isempty(env_root) && append!(candidates, _resolve_path_candidates(env_root))

    if !isempty(gram_directory)
        append!(candidates, _resolve_path_candidates(gram_directory))
        append!(candidates, _resolve_path_candidates(joinpath(gram_directory, "..")))
        append!(candidates, _resolve_path_candidates(joinpath(gram_directory, "..", "..")))
    end

    append!(candidates, _resolve_path_candidates("GRAM Suite 2.0"))

    for candidate in candidates
        _is_gram_root(candidate) && return candidate
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
    !isempty(spice_directory) && append!(candidates, _resolve_path_candidates(spice_directory))
    push!(candidates, normpath(joinpath(gram_root, "SPICE")))
    !isempty(gram_data_directory) && append!(candidates, _resolve_path_candidates(joinpath(gram_data_directory, "SPICE")))
    append!(candidates, _resolve_path_candidates("GRAM Suite 2.0/SPICE"))
    append!(candidates, _resolve_path_candidates("GRAM_Data/SPICE"))

    for candidate in candidates
        isdir(candidate) && return candidate
    end

    throw(ArgumentError("Unable to find a valid SPICE directory for GRAM initialization."))
end

function _resolve_gram_data_root(gram_root::String, gram_data_directory::String)::String
    candidates = String[normpath(gram_root)]
    !isempty(gram_data_directory) && append!(candidates, _resolve_path_candidates(gram_data_directory))
    append!(candidates, _resolve_path_candidates("GRAM Suite 2.0"))
    append!(candidates, _resolve_path_candidates("GRAM_Data"))

    for candidate in candidates
        if isdir(joinpath(candidate, "Earth", "data")) || isdir(joinpath(candidate, "Mars", "data"))
            return candidate
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

@inline function _gram_offline_surrogate_features_supported(
    planet_key::String;
    gram_perturbation_scales=nothing,
    mars_map_year=nothing,
    mars_mgcm_dust_levels=nothing,
    mars_dust_storm=nothing,
    mars_f107=nothing,
    mars_wind_scales=nothing,
    mars_mola_heights=nothing,
    mars_min_max=nothing
)::Tuple{Bool, String}
    unsupported = String[]

    if gram_perturbation_scales !== nothing
        push!(unsupported, "gram_perturbation_scales")
    end

    if planet_key == "mars"
        mars_map_year !== nothing && push!(unsupported, "mars_map_year")
        mars_mgcm_dust_levels !== nothing && push!(unsupported, "mars_mgcm_dust_levels")
        mars_dust_storm !== nothing && push!(unsupported, "mars_dust_storm")
        mars_f107 !== nothing && push!(unsupported, "mars_f107")
        mars_wind_scales !== nothing && push!(unsupported, "mars_wind_scales")
        mars_mola_heights !== nothing && push!(unsupported, "mars_mola_heights")
        mars_min_max !== nothing && push!(unsupported, "mars_min_max")
    end

    return isempty(unsupported), join(unsupported, ",")
end

function GRAMAtmosphereModel(;
    gram_directory::String="GRAM Suite 2.0",
    gram_data_directory::String="GRAM Suite 2.0",
    gram_root_directory::String="",
    gram_library_path::String="",
    spice_directory::String="",
    planet_name::String="earth",
    seed::Int=1001,
    initial_time=InitialTime(),
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
    it = _coerce_initial_time(initial_time)
    planet_key = lowercase(strip(planet_name))

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
            initial_time=it,
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
        gram.set_library!(resolve_path(gram_library_path))
    else
        local_lib = joinpath(gram_root, "Build", "lib", "libGRAM.$(_gram_lib_extension())")
        if isfile(local_lib)
            gram.set_library!(local_lib)
        end
    end

    spice_root = _resolve_spice_directory(gram_root, spice_directory, gram_data_directory)
    gram.initialize!(spice_root)

    gram_data_root = _resolve_gram_data_root(gram_root, gram_data_directory)
    body = _gram_body(gram, planet_key)
    data_path = _gram_data_path(gram_data_root, planet_key)
    gram_atmosphere = data_path === nothing ? gram.create_atmosphere(body) : gram.create_atmosphere(body; data_path=data_path)

    gram.set_start_time!(
        gram_atmosphere;
        year=Int(it.year),
        month=Int(it.month),
        day=Int(it.day),
        hour=Int(it.hour),
        minute=Int(it.minute),
        seconds=Float64(it.second),
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

    offline_surrogate_supported, offline_surrogate_unsupported_reason = _gram_offline_surrogate_features_supported(
        planet_key;
        gram_perturbation_scales=gram_perturbation_scales,
        mars_map_year=mars_map_year,
        mars_mgcm_dust_levels=mars_mgcm_dust_levels,
        mars_dust_storm=mars_dust_storm,
        mars_f107=mars_f107,
        mars_wind_scales=mars_wind_scales,
        mars_mola_heights=mars_mola_heights,
        mars_min_max=mars_min_max
    )

    model = GRAMAtmosphereModel(
        gram,
        gram_atmosphere,
        gram_root,
        gram_data_root,
        spice_root,
        planet_key,
        it,
        offline_surrogate_supported,
        offline_surrogate_unsupported_reason
    )

    if gram_static_grid_enabled() && _gram_static_grid_prebuild_all_planets_enabled() && !_GRAM_STATIC_GRID_PREBUILD_IN_PROGRESS[]
        wind_enabled = _parse_bool_env("SPACEAGORA_GRAM_STATIC_GRID_WIND", true)
        precompute_gram_static_grids!(model; wind=wind_enabled)
    end

    return model
end

function GRAMAtmosphereModelSurrogate(;
    surrogate_file::String="",
    point_fallback_below_m::Union{Nothing, Real}=nothing,
    kwargs...
)
    base_model = GRAMAtmosphereModel(; kwargs...)
    if !base_model.offline_surrogate_supported
        throw(ArgumentError(
            "GRAMAtmosphereModelSurrogate does not support this GRAM configuration. " *
            "Unsupported feature(s): $(base_model.offline_surrogate_unsupported_reason)"
        ))
    end

    file = isempty(strip(surrogate_file)) ?
        gram_default_surrogate_file(base_model.planet_name; gram_root=base_model.gram_root) :
        resolve_path(surrogate_file)
    point_fallback = if point_fallback_below_m === nothing
        gram_default_point_fallback_below_m(base_model.planet_name)
    else
        value = Float64(point_fallback_below_m)
        value >= 0.0 || throw(ArgumentError("point_fallback_below_m must be >= 0.0 m, got $value"))
        value
    end

    return GRAMAtmosphereModelSurrogate(base_model, file, point_fallback)
end

@inline function gram_default_point_fallback_below_m(planet::AbstractString)::Union{Nothing, Float64}
    planet_key = lowercase(strip(String(planet)))
    return get(_GRAM_OFFLINE_SURROGATE_POINT_FALLBACK_BELOW_M_DEFAULT, planet_key, nothing)
end

function gram_default_surrogate_file(
    planet::AbstractString;
    gram_root::String="",
    gram_root_directory::String="",
    gram_directory::String="GRAM Suite 2.0"
)::String
    planet_key = _gram_normalize_planet_name(planet)
    root = if !isempty(strip(gram_root))
        normpath(strip(gram_root))
    else
        _resolve_gram_root(gram_root_directory, gram_directory)
    end

    planet_dir = _gram_planet_dir_name(planet_key)
    preferred = joinpath(root, planet_dir, "$(planet_key)_surrogate.jls")
    isfile(preferred) && return preferred

    legacy = joinpath(
        root,
        "simulation",
        "GRAM",
        "static_grids",
        _GRAM_OFFLINE_SURROGATE_PROFILE,
        "surrogates",
        "$(planet_key)_surrogate.jls"
    )
    isfile(legacy) && return legacy

    return preferred
end

@inline function _gram_offline_surrogate_enabled(gram_root::String)::Bool
    mode = lowercase(strip(get(ENV, "SPACEAGORA_GRAM_OFFLINE_SURROGATE", "off")))
    if mode in ("on", "true", "1", "yes")
        return true
    elseif mode in ("off", "false", "0", "no")
        return false
    elseif mode == "auto"
        env_lib = strip(get(ENV, "GRAM_LIB", ""))
        if !isempty(env_lib)
            return !isfile(resolve_path(env_lib))
        end
        local_lib = joinpath(gram_root, "Build", "lib", "libGRAM.$(_gram_lib_extension())")
        return !isfile(local_lib)
    end
    throw(ArgumentError(
        "Unsupported SPACEAGORA_GRAM_OFFLINE_SURROGATE='$mode'. Use one of: off, on, auto."
    ))
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
    planet_key = lowercase(strip(model.planet_name))
    alt_min_default, alt_max_default = get(
        _GRAM_FROZEN_PLANET_ALT_RANGE_M,
        planet_key,
        (0.0, 2_000_000.0)
    )
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
        planet_key,
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

@inline function _gram_offline_surrogate_key(planet::String, file::String)::NamedTuple{(:planet, :file), Tuple{String, String}}
    return (planet=lowercase(strip(planet)), file=normpath(file))
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

    return GRAMOfflineSurrogate(
        payload_planet,
        file,
        lowercase(strip(String(get(dict, "wind_mode", "unknown")))),
        lowercase(strip(String(get(dict, "monte_carlo", "unknown")))),
        lowercase(strip(String(get(dict, "dust", "unknown")))),
        alt_nodes_m,
        lat_nodes_rad,
        lon_nodes_rad,
        rho,
        T,
        wind_e,
        wind_n,
        wind_u
    )
end

function _gram_offline_surrogate_get_or_load!(planet::String, file::String)::GRAMOfflineSurrogate
    key = _gram_offline_surrogate_key(planet, file)
    if !isfile(key.file)
        throw(ArgumentError("GRAM surrogate payload not found for planet='$(key.planet)': $(key.file)"))
    end

    lock(_GRAM_OFFLINE_SURROGATE_LOCK) do
        if haskey(_GRAM_OFFLINE_SURROGATE_CACHE, key)
            return _GRAM_OFFLINE_SURROGATE_CACHE[key]::GRAMOfflineSurrogate
        end
        surrogate = _gram_load_offline_surrogate(key.file, key.planet)
        _GRAM_OFFLINE_SURROGATE_CACHE[key] = surrogate
        if !_GRAM_OFFLINE_SURROGATE_LOGGED[]
            _GRAM_OFFLINE_SURROGATE_LOGGED[] = true
            @info "GRAM offline surrogate interpolation loaded."
        end
        @info "Loaded GRAM offline surrogate payload." planet=surrogate.planet_name file=surrogate.source_file
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

    rho = _gram_trilerp(
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
    return rho, Ti, SVector{3, Float64}(wE, wN, wU)
end

function _gram_density_state_native(
    model::GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
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
        wind_mode = _gram_wind_mode()
        use_perturbed_winds = wind && wind_mode == :perturbed
        wind_vec_local = if use_perturbed_winds
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
        if !all(isfinite, wind_vec_local)
            if !_GRAM_NONFINITE_WIND_WARNING_EMITTED[]
                _GRAM_NONFINITE_WIND_WARNING_EMITTED[] = true
                @warn(
                    "GRAM returned non-finite wind component(s). Replacing non-finite values with 0.0.",
                    planet=model.planet_name,
                    h_m=h,
                    lat_deg=rad2deg(lat),
                    lon_deg=rad2deg(lon),
                    elapsed_time_s=el_time,
                    wind_raw=wind_vec_local
                )
            end
            wind_vec_local = SVector{3, Float64}(
                isfinite(wind_vec_local[1]) ? wind_vec_local[1] : 0.0,
                isfinite(wind_vec_local[2]) ? wind_vec_local[2] : 0.0,
                isfinite(wind_vec_local[3]) ? wind_vec_local[3] : 0.0
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

@inline function point_density_state(
    model::GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool;
    lock_obj=nothing
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return _with_gram_lock(() -> _gram_density_state_native(model, h, lat, lon, el_time, wind), lock_obj)
end

function _gram_static_grid_build(model::GRAMAtmosphereModel, key::GRAMStaticGridKey; lock_obj=nothing)::GRAMStaticGrid
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
                rho_i, T_i, wv = point_density_state(model, h, lat, lon, key.elapsed_time_s, key.include_wind; lock_obj=lock_obj)
                rho[ia, ilat, ilon] = rho_i
                T[ia, ilat, ilon] = T_i
                wind_e[ia, ilat, ilon] = wv[1]
                wind_n[ia, ilat, ilon] = wv[2]
                wind_u[ia, ilat, ilon] = wv[3]
            end
        end
    end

    @info "Built GRAM static grid." planet=key.planet_name n_alt=key.n_alt n_lat=key.n_lat n_lon=key.n_lon elapsed_s=(time_ns() - t0_ns) * 1e-9

    return GRAMStaticGrid(key, alt_nodes, lat_nodes, lon_nodes, rho, T, wind_e, wind_n, wind_u)
end

function _gram_static_grid_get_or_build!(model::GRAMAtmosphereModel, wind::Bool; lock_obj=nothing)::GRAMStaticGrid
    key = _gram_static_grid_key(model, wind)
    lock(_GRAM_STATIC_GRID_LOCK) do
        if haskey(_GRAM_STATIC_GRID_CACHE, key)
            return _GRAM_STATIC_GRID_CACHE[key]::GRAMStaticGrid
        end
        grid = _gram_static_grid_build(model, key; lock_obj=lock_obj)
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

function clear_gram_offline_surrogate_cache!()
    lock(_GRAM_OFFLINE_SURROGATE_LOCK) do
        empty!(_GRAM_OFFLINE_SURROGATE_CACHE)
        _GRAM_OFFLINE_SURROGATE_LOGGED[] = false
    end
    lock(_GRAM_OFFLINE_SURROGATE_WARNED_LOCK) do
        empty!(_GRAM_OFFLINE_SURROGATE_WARNED)
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
    wind::Bool=true,
    lock_obj=nothing
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
                _gram_static_grid_get_or_build!(model, wind; lock_obj=lock_obj)
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

    rho = _gram_trilerp(
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

    return rho, Ti, SVector{3, Float64}(wE, wN, wU)
end

function density_state(
    model::GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool;
    lock_obj=nothing,
    vacuum_temperature::Float64=200.0
)::Tuple{Float64, Float64, SVector{3, Float64}}
    if _gram_offline_surrogate_enabled(model.gram_root)
        if model.offline_surrogate_supported
            file = gram_default_surrogate_file(model.planet_name; gram_root=model.gram_root)
            if isfile(file)
                surrogate = _gram_offline_surrogate_get_or_load!(model.planet_name, file)
                point_fallback_below_m = gram_default_point_fallback_below_m(model.planet_name)
                if point_fallback_below_m === nothing || h > point_fallback_below_m
                    s_eval = _gram_offline_surrogate_eval(surrogate, h, lat, lon)
                    if s_eval !== nothing
                        return s_eval
                    end
                    grid_alt_max_m = last(surrogate.alt_nodes_m)
                    if h > grid_alt_max_m
                        _gram_warn_once(
                            "above_ceiling_zero_rho:$(surrogate.planet_name)",
                            "State altitude is above GRAM offline surrogate ceiling. Returning rho=0 for this sample.",
                            planet=surrogate.planet_name,
                            h_m=h,
                            grid_alt_max_m=grid_alt_max_m
                        )
                        return 0.0, vacuum_temperature, SVector{3, Float64}(0.0, 0.0, 0.0)
                    end
                else
                    _gram_warn_once(
                        "point_fallback_below:$(surrogate.planet_name)",
                        "State altitude is below configured surrogate point-fallback altitude. Using point-to-point GRAM for this sample.",
                        planet=surrogate.planet_name,
                        h_m=h,
                        fallback_below_m=point_fallback_below_m
                    )
                end
            else
                _gram_warn_once(
                    "missing:$(model.planet_name):$(file)",
                    "GRAM offline surrogate payload not found; using point-to-point GRAM.",
                    planet=model.planet_name,
                    file=file
                )
            end
        else
            _gram_warn_once(
                "unsupported:$(model.planet_name):$(model.offline_surrogate_unsupported_reason)",
                "GRAM offline surrogate disabled for unsupported feature configuration; using point-to-point GRAM.",
                planet=model.planet_name,
                unsupported=model.offline_surrogate_unsupported_reason
            )
        end
    end

    if gram_static_grid_enabled()
        grid = _gram_static_grid_get_or_build!(model, wind; lock_obj=lock_obj)
        grid_alt_max_m = last(grid.alt_nodes)
        if h > grid_alt_max_m
            _gram_warn_once(
                "above_ceiling_zero_rho_static:$(grid.key.planet_name)",
                "State altitude is above GRAM static-grid ceiling. Returning rho=0 for this sample.",
                planet=grid.key.planet_name,
                h_m=h,
                grid_alt_max_m=grid_alt_max_m
            )
            return 0.0, vacuum_temperature, SVector{3, Float64}(0.0, 0.0, 0.0)
        end
        return _gram_static_grid_eval(grid, h, lat, lon)
    end

    return point_density_state(model, h, lat, lon, el_time, wind; lock_obj=lock_obj)
end

function _planet_name_from_base_model(base_model)::String
    if hasproperty(base_model, :planet_name)
        return lowercase(strip(String(getproperty(base_model, :planet_name))))
    end
    throw(ArgumentError("Base model for surrogate path must expose `planet_name`."))
end

function _surrogate_point_fallback(
    base_model,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool;
    lock_obj=nothing,
    point_density_fallback=nothing
)::Tuple{Float64, Float64, SVector{3, Float64}}
    if base_model isa GRAMAtmosphereModel
        return point_density_state(base_model, h, lat, lon, el_time, wind; lock_obj=lock_obj)
    end
    if point_density_fallback === nothing
        throw(ArgumentError("No point-density fallback is available for base model type $(typeof(base_model))."))
    end
    return point_density_fallback(base_model, h, lat, lon, el_time, wind)
end

function surrogate_density_state(
    base_model,
    surrogate_file::String,
    point_fallback_below_m::Union{Nothing, Float64},
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool;
    lock_obj=nothing,
    point_density_fallback=nothing,
    vacuum_temperature::Float64=200.0
)::Tuple{Float64, Float64, SVector{3, Float64}}
    file = resolve_path(surrogate_file)
    planet = _planet_name_from_base_model(base_model)
    surrogate = _gram_offline_surrogate_get_or_load!(planet, file)

    if point_fallback_below_m !== nothing && h <= point_fallback_below_m
        _gram_warn_once(
            "point_fallback_below:$(surrogate.planet_name)",
            "State altitude is below configured surrogate point-fallback altitude. Using point-to-point GRAM for this sample.",
            planet=surrogate.planet_name,
            h_m=h,
            fallback_below_m=point_fallback_below_m
        )
        return _surrogate_point_fallback(base_model, h, lat, lon, el_time, wind; lock_obj=lock_obj, point_density_fallback=point_density_fallback)
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
        return 0.0, vacuum_temperature, SVector{3, Float64}(0.0, 0.0, 0.0)
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
    return _surrogate_point_fallback(base_model, h, lat, lon, el_time, wind; lock_obj=lock_obj, point_density_fallback=point_density_fallback)
end

function density_state(
    model::GRAMAtmosphereModelSurrogate,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool;
    lock_obj=nothing,
    point_density_fallback=nothing,
    vacuum_temperature::Float64=200.0
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return surrogate_density_state(
        model.base_model,
        model.surrogate_file,
        model.point_fallback_below_m,
        h,
        lat,
        lon,
        el_time,
        wind;
        lock_obj=lock_obj,
        point_density_fallback=point_density_fallback,
        vacuum_temperature=vacuum_temperature
    )
end

function Base.deepcopy_internal(model::GRAMAtmosphereModel, stackdict::IdDict)
    if haskey(stackdict, model)
        return stackdict[model]
    end

    copied = Base.invokelatest(
        GRAMAtmosphereModel;
        gram_root_directory=model.gram_root,
        gram_data_directory=model.gram_data_root,
        spice_directory=model.spice_root,
        planet_name=model.planet_name,
        initial_time=model.initial_time
    )
    stackdict[model] = copied
    return copied
end

function Base.deepcopy_internal(model::GRAMAtmosphereModelSurrogate, stackdict::IdDict)
    if haskey(stackdict, model)
        return stackdict[model]
    end

    copied = GRAMAtmosphereModelSurrogate(
        deepcopy(model.base_model),
        model.surrogate_file,
        model.point_fallback_below_m
    )
    stackdict[model] = copied
    return copied
end

end # module GRAMSuite
