module SpaceAGORAGRAMSuiteExt

using SpaceAGORA
using GRAMSuite
using StaticArrays

const EM = SpaceAGORA.SimulationModel.EnvironmentModels
const GRAM_LOCK = SpaceAGORA.RuntimeServices.GRAM_LOCK

function __init__()
    EM._GRAM_USE_GLOBAL_LOCK_FN[] = () -> GRAMSuite.gram_use_global_lock()
    EM._GRAM_DEFAULT_SURROGATE_FILE_FN[] = planet -> GRAMSuite.gram_default_surrogate_file(planet)
    EM._CLEAR_GRAM_STATIC_GRID_CACHE_FN[] = () -> GRAMSuite.clear_gram_static_grid_cache!()
    EM._CLEAR_GRAM_OFFLINE_SURROGATE_CACHE_FN[] = () -> GRAMSuite.clear_gram_offline_surrogate_cache!()
end

# ---------------------------------------------------------------------------
# Constructors
# ---------------------------------------------------------------------------

function EM.GRAMAtmosphereModel(; kwargs...)
    return EM.GRAMAtmosphereModel(GRAMSuite.GRAMAtmosphereModel(; kwargs...))
end

function EM.GRAMAtmosphereModelSurrogate(;
    surrogate_file::String="",
    point_fallback_below_m::Union{Nothing, Real}=nothing,
    kwargs...
)
    base_model = EM.GRAMAtmosphereModel(; kwargs...)
    file = isempty(strip(surrogate_file)) ?
        GRAMSuite.gram_default_surrogate_file(base_model.planet_name; gram_root=base_model.gram_root) :
        GRAMSuite.resolve_path(surrogate_file)
    point_fallback = if point_fallback_below_m === nothing
        GRAMSuite.gram_default_point_fallback_below_m(base_model.planet_name)
    else
        value = Float64(point_fallback_below_m)
        value >= 0.0 || throw(ArgumentError("point_fallback_below_m must be >= 0.0 m, got $value"))
        value
    end
    return EM.GRAMAtmosphereModelSurrogate(base_model, file, point_fallback)
end

# ---------------------------------------------------------------------------
# Precomputation / cache management
# ---------------------------------------------------------------------------

function EM.precompute_gram_static_grids!(
    base_model::EM.GRAMAtmosphereModel;
    planets::Union{Nothing, AbstractVector{<:AbstractString}}=nothing,
    wind::Bool=true
)
    return GRAMSuite.precompute_gram_static_grids!(
        base_model.core;
        planets=planets,
        wind=wind,
        lock_obj=GRAM_LOCK
    )
end


# ---------------------------------------------------------------------------
# deepcopy
# ---------------------------------------------------------------------------

function Base.deepcopy_internal(model::EM.GRAMAtmosphereModel, stackdict::IdDict)
    haskey(stackdict, model) && return stackdict[model]
    copied = lock(GRAM_LOCK) do
        EM.GRAMAtmosphereModel(deepcopy(model.core))
    end
    stackdict[model] = copied
    return copied
end

function Base.deepcopy_internal(model::EM.GRAMAtmosphereModelSurrogate, stackdict::IdDict)
    haskey(stackdict, model) && return stackdict[model]
    copied = lock(GRAM_LOCK) do
        EM.GRAMAtmosphereModelSurrogate(
            deepcopy(model.base_model),
            model.surrogate_file,
            model.point_fallback_below_m
        )
    end
    stackdict[model] = copied
    return copied
end

# ---------------------------------------------------------------------------
# _gram_core_density_state — used by the isolated-pool callback path
# ---------------------------------------------------------------------------

function EM._gram_core_density_state(
    core::GRAMSuite.GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    lock_obj,
    vacuum_temperature::Float64
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return GRAMSuite.density_state(
        core, h, lat, lon, el_time, wind;
        lock_obj=lock_obj,
        vacuum_temperature=vacuum_temperature
    )
end

# ---------------------------------------------------------------------------
# _gram_point_density
# ---------------------------------------------------------------------------

@inline function EM._gram_point_density(
    model::EM.GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return GRAMSuite.point_density_state(model.core, h, lat, lon, el_time, wind; lock_obj=GRAM_LOCK)
end

# ---------------------------------------------------------------------------
# getDensity
# ---------------------------------------------------------------------------

function EM.getDensity(
    model::EM.GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p::params
)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0

    if h > 2000.0e3
        rho = 0.0
        T = p.args.environment_model.planet.T_ref
        wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif !drag_state && !p.args.mission_configuration.keplerian
        rho, T, wind_vec = EM.density_polyfit(h, p)
    else
        rho, T, wind_vec = GRAMSuite.density_state(
            model.core,
            h,
            lat,
            lon,
            el_time,
            wind;
            lock_obj=GRAM_LOCK,
            vacuum_temperature=p.args.environment_model.planet.T_ref
        )
    end

    return rho, T, wind_vec
end

function EM.getDensity(
    model::EM.GRAMAtmosphereModelSurrogate,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p::params
)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0

    if h > 2000.0e3
        return 0.0, p.args.environment_model.planet.T_ref, SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif !drag_state && !p.args.mission_configuration.keplerian
        return EM.density_polyfit(h, p)
    end

    base_model = model.base_model isa EM.GRAMAtmosphereModel ? model.base_model.core : model.base_model
    point_fallback = model.base_model isa EM.GRAMAtmosphereModel ? nothing :
        (m, h_i, lat_i, lon_i, t_i, w_i) -> EM._gram_point_density(m, h_i, lat_i, lon_i, t_i, w_i)

    println("GRAM density altitude = $(h) m ($(h / 1e3) km)")
    return GRAMSuite.surrogate_density_state(
        base_model,
        model.surrogate_file,
        model.point_fallback_below_m,
        h,
        lat,
        lon,
        el_time,
        wind;
        lock_obj=GRAM_LOCK,
        point_density_fallback=point_fallback,
        vacuum_temperature=p.args.environment_model.planet.T_ref
    )
end

end # module SpaceAGORAGRAMSuiteExt
