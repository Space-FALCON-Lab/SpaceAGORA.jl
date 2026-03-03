# include("../utils/Reference_system.jl")

using SatelliteToolbox
using StaticArrays
using LinearAlgebra
using GRAMSuite
using ..SimulationModel: GRAM_LOCK

struct NoAtmosphereModel <: AbstractDensityModel
end

struct ExponentialAtmosphereModel <: AbstractDensityModel
    ρ_ref::Float64
    h_ref::Float64
    H::Float64
end

"""
SpaceAGORA wrapper around `GRAMSuite.GRAMAtmosphereModel` that preserves
`AbstractDensityModel` dispatch compatibility.
"""
struct GRAMAtmosphereModel <: AbstractDensityModel
    core::GRAMSuite.GRAMAtmosphereModel
end

"""
Surrogate model wrapper kept shape-compatible with the previous implementation,
including the ability to provide a non-GRAM base model for custom point fallback.
"""
struct GRAMAtmosphereModelSurrogate{M} <: AbstractDensityModel
    base_model::M
    surrogate_file::String
    point_fallback_below_m::Union{Nothing, Float64}
end

struct PolynomialFitAtmosphereModel <: AbstractDensityModel
    polyfit_coeffs::Vector{Float64}
end

struct NRLMSISE00AtmosphereModel <: AbstractDensityModel
end

@inline function Base.getproperty(model::GRAMAtmosphereModel, name::Symbol)
    if name === :core
        return getfield(model, :core)
    end
    return getproperty(getfield(model, :core), name)
end

@inline function Base.propertynames(model::GRAMAtmosphereModel, private::Bool=false)
    wrapped = propertynames(getfield(model, :core), private)
    return (:core, wrapped...)
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

@inline _gram_use_global_lock() = GRAMSuite.gram_use_global_lock()

@inline function _gram_default_surrogate_file(planet::String)::String
    return GRAMSuite.gram_default_surrogate_file(planet)
end

function GRAMAtmosphereModel(; kwargs...)
    return GRAMAtmosphereModel(GRAMSuite.GRAMAtmosphereModel(; kwargs...))
end

function GRAMAtmosphereModelSurrogate(;
    surrogate_file::String="",
    point_fallback_below_m::Union{Nothing, Real}=nothing,
    kwargs...
)
    base_model = GRAMAtmosphereModel(; kwargs...)
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

    return GRAMAtmosphereModelSurrogate(base_model, file, point_fallback)
end

function precompute_gram_static_grids!(
    base_model::GRAMAtmosphereModel;
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

clear_gram_static_grid_cache!() = GRAMSuite.clear_gram_static_grid_cache!()
clear_gram_offline_surrogate_cache!() = GRAMSuite.clear_gram_offline_surrogate_cache!()

function Base.deepcopy_internal(model::GRAMAtmosphereModel, stackdict::IdDict)
    if haskey(stackdict, model)
        return stackdict[model]
    end

    copied = GRAMAtmosphereModel(deepcopy(model.core))
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

function interp(a, b, x)
    if abs(b - a) > 20.0
        if b <= 360.0 && b >= 350.0
            b = 360.0 - b
        elseif a <= 360.0 && a >= 350.0
            a = 360.0 - a
        end
    end

    return x * (b - a) + a
end

function temperature_linear(h, p)
    return p.T_ref
end

function getDensity(model::NoAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    rho = 0.0
    T = p.args.environment_model.planet.T_ref
    wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    return rho, T, wind_vec
end

function getDensity(model::ExponentialAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}
    rho = model.ρ_ref * exp((model.h_ref - h) / model.H)
    T = 200.0
    wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    return rho, T, wind_vec
end

function getDensity(model::PolynomialFitAtmosphereModel, h::Union{Float64, Vector{Float64}}, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    if !(h isa Float64)
        polyfit = model.polyfit_coeffs
        power = zeros(length(polyfit), length(h))
        h_km = h .* 1e-3
        for i in eachindex(polyfit)
            power[i, :] = h_km .^ (length(polyfit) - i)
        end

        exponent = zeros(length(h_km))
        @inbounds for j in eachindex(h_km)
            exponent[j] = sum(polyfit .* power[:, j])
        end

        rho = exp.(exponent)
        T = p.args.environment_model.planet.T_ref
        wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
        return rho, T, wind_vec
    end

    polyfit = model.polyfit_coeffs
    power = MVector{length(polyfit)}(zeros(length(polyfit)))
    h_km = h * 1e-3
    @inbounds for i in eachindex(polyfit)
        power[i] = h_km ^ (length(polyfit) - i)
    end
    exponent = sum(polyfit .* power)
    rho = exp(exponent)
    T = p.args.environment_model.planet.T_ref
    wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    return rho, T, wind_vec
end

function _gram_point_density(
    model,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    throw(MethodError(_gram_point_density, (model, h, lat, lon, el_time, wind)))
end

@inline function _gram_point_density(
    model::GRAMAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return GRAMSuite.point_density_state(model.core, h, lat, lon, el_time, wind; lock_obj=GRAM_LOCK)
end

@inline function _gram_point_density(
    model::GRAMAtmosphereModelSurrogate,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return _gram_point_density(model.base_model, h, lat, lon, el_time, wind)
end

function getDensity(model::GRAMAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0

    if h > 2000.0e3
        rho = 0.0
        T = p.args.environment_model.planet.T_ref
        wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif !drag_state && !p.args.mission_configuration.keplerian
        rho, T, wind_vec = density_polyfit(h, p)
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

function getDensity(model::GRAMAtmosphereModelSurrogate, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    EI = p.args.environment_model.EI * 1e3
    drag_state = h - EI <= 0.0

    if h > 2000.0e3
        return 0.0, p.args.environment_model.planet.T_ref, SVector{3, Float64}(0.0, 0.0, 0.0)
    elseif !drag_state && !p.args.mission_configuration.keplerian
        return density_polyfit(h, p)
    end

    base_model = model.base_model isa GRAMAtmosphereModel ? model.base_model.core : model.base_model
    point_fallback = model.base_model isa GRAMAtmosphereModel ? nothing :
        (m, h_i, lat_i, lon_i, t_i, w_i) -> _gram_point_density(m, h_i, lat_i, lon_i, t_i, w_i)

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

function getDensity(model::NRLMSISE00AtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}
    if config.cnf.drag_state == false && args[:keplerian] == false
        rho, T, wind_vec = density_exp(h, p)
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
    coeffs = p.args.environment_model.planet.polyfit_coeffs
    poly_model = PolynomialFitAtmosphereModel(vec(coeffs))
    rho, T, wind_vec = getDensity(poly_model, h, 0.0, 0.0, 0.0, false, p)
    return rho, T, wind_vec
end
