# include(joinpath(@__DIR__, "..", "..", "core", "interfaces", "reference_system.jl"))

using SatelliteToolbox
using StaticArrays
using LinearAlgebra
using Dates

"""
Vacuum baseline density model used by the default no-GRAM quickstart path.
"""
struct NoAtmosphereModel <: AbstractDensityModel
end

"""
    density_vanishes_above_entry_interface(model) -> Bool

Whether the split/implicit solver partition may treat this density model as
exactly zero above the configured entry interface (`EnvironmentModel.EI`).

Only models guaranteed to produce no density above EI may return `true`; the
implicit-atmosphere RHS uses this to skip aerodynamic evaluation entirely on
coast arcs. For every other model EI is a step-size/tolerance boundary, not a
force gate: aero must stay engaged at all altitudes (an orbit that never
crosses EI downward would otherwise silently fly drag-free), and genuinely
negligible densities are handled by the wrench's own `rho <= eps` short-circuit.
"""
@inline density_vanishes_above_entry_interface(::AbstractDensityModel)::Bool = false
@inline density_vanishes_above_entry_interface(::NoAtmosphereModel)::Bool = true

"""
Single-scale-height analytic atmosphere with zero winds and a constant
temperature placeholder.

`valid_min_altitude_m` and `valid_max_altitude_m` are advisory bounds that
document the altitude band the fit was intended to cover. Density evaluation
still extrapolates the same exponential outside that band.
"""
struct ExponentialAtmosphereModel <: AbstractDensityModel
    ρ_ref::Float64
    h_ref::Float64
    H::Float64
    temperature_k::Float64
    valid_min_altitude_m::Float64
    valid_max_altitude_m::Float64
end

function ExponentialAtmosphereModel(
    ρ_ref::Real,
    h_ref::Real,
    H::Real;
    temperature_k::Real=200.0,
    valid_min_altitude_m::Real=h_ref,
    valid_max_altitude_m::Real=h_ref + 5.0 * H
)
    H_f = Float64(H)
    H_f > 0.0 || throw(ArgumentError("ExponentialAtmosphereModel scale height H must be > 0.0 m, got $(Float64(H))."))
    valid_min_f = Float64(valid_min_altitude_m)
    valid_max_f = Float64(valid_max_altitude_m)
    valid_min_f <= valid_max_f || throw(ArgumentError(
        "ExponentialAtmosphereModel valid_min_altitude_m must be <= valid_max_altitude_m, got $valid_min_f > $valid_max_f."
    ))
    return ExponentialAtmosphereModel(
        Float64(ρ_ref),
        Float64(h_ref),
        H_f,
        Float64(temperature_k),
        valid_min_f,
        valid_max_f
    )
end

@inline function ExponentialAtmosphereModel(planet::AbstractPlanet; kwargs...)
    return ExponentialAtmosphereModel(Float64(planet.ρ_ref), Float64(planet.h_ref), Float64(planet.H); kwargs...)
end

"""
Piecewise multi-layer exponential atmosphere with zero winds and a constant
temperature placeholder.

`h_breaks_m` defines contiguous layer boundaries in meters. For `N` layers,
provide `N + 1` strictly increasing breakpoints, `N` reference densities, `N`
reference altitudes, and `N` scale heights. `valid_min_altitude_m` and
`valid_max_altitude_m` document the intended altitude band; density evaluation
uses the nearest layer to extrapolate outside that band.
"""
struct PiecewiseExponentialAtmosphereModel <: AbstractDensityModel
    h_breaks_m::Vector{Float64}
    ρ_refs::Vector{Float64}
    h_refs::Vector{Float64}
    Hs::Vector{Float64}
    temperature_k::Float64
    valid_min_altitude_m::Float64
    valid_max_altitude_m::Float64
end

function PiecewiseExponentialAtmosphereModel(
    h_breaks_m::AbstractVector{<:Real},
    ρ_refs::AbstractVector{<:Real},
    Hs::AbstractVector{<:Real};
    h_refs::Union{Nothing, AbstractVector{<:Real}}=nothing,
    temperature_k::Real=200.0,
    valid_min_altitude_m::Union{Nothing, Real}=nothing,
    valid_max_altitude_m::Union{Nothing, Real}=nothing
)
    length(h_breaks_m) >= 2 || throw(ArgumentError(
        "PiecewiseExponentialAtmosphereModel requires at least two breakpoints for one layer."
    ))

    h_breaks = Float64.(collect(h_breaks_m))
    n_layers = length(h_breaks) - 1
    h_refs_vec = h_refs === nothing ? h_breaks[1:n_layers] : collect(h_refs)
    length(ρ_refs) == n_layers || throw(ArgumentError(
        "PiecewiseExponentialAtmosphereModel requires $(n_layers) density references, got $(length(ρ_refs))."
    ))
    length(h_refs_vec) == n_layers || throw(ArgumentError(
        "PiecewiseExponentialAtmosphereModel requires $(n_layers) reference altitudes, got $(length(h_refs_vec))."
    ))
    length(Hs) == n_layers || throw(ArgumentError(
        "PiecewiseExponentialAtmosphereModel requires $(n_layers) scale heights, got $(length(Hs))."
    ))

    @inbounds for i in 1:n_layers
        h_breaks[i] < h_breaks[i + 1] || throw(ArgumentError(
            "PiecewiseExponentialAtmosphereModel breakpoints must be strictly increasing; found $(h_breaks[i]) >= $(h_breaks[i + 1])."
        ))
    end

    Hs_f = Float64.(collect(Hs))
    @inbounds for i in eachindex(Hs_f)
        Hs_f[i] > 0.0 || throw(ArgumentError(
            "PiecewiseExponentialAtmosphereModel scale heights must be > 0.0 m; layer $i has $(Hs_f[i])."
        ))
    end

    valid_min_f = Float64(something(valid_min_altitude_m, h_breaks[1]))
    valid_max_f = Float64(something(valid_max_altitude_m, h_breaks[end]))
    valid_min_f <= valid_max_f || throw(ArgumentError(
        "PiecewiseExponentialAtmosphereModel valid_min_altitude_m must be <= valid_max_altitude_m, got $valid_min_f > $valid_max_f."
    ))

    return PiecewiseExponentialAtmosphereModel(
        h_breaks,
        Float64.(collect(ρ_refs)),
        Float64.(collect(h_refs_vec)),
        Hs_f,
        Float64(temperature_k),
        valid_min_f,
        valid_max_f
    )
end

"""
SpaceAGORA wrapper around a GRAMSuite atmosphere model that preserves
`AbstractDensityModel` dispatch compatibility. The `core` field holds the
underlying GRAMSuite model object; constructors are provided by the
`SpaceAGORAGRAMSuiteExt` package extension when GRAMSuite is loaded.

`instance_lock` serializes native GRAM calls against this wrapper instance when
`SPACEAGORA_GRAM_LOCK_SCOPE=model` is active (see [`_gram_lock_scope`](@ref));
each construction (including `deepcopy`) gets a fresh lock.
"""
struct GRAMAtmosphereModel <: AbstractDensityModel
    core
    instance_lock::ReentrantLock
end

GRAMAtmosphereModel(core) = GRAMAtmosphereModel(core, ReentrantLock())

@kwdef struct ConstantDensityModel <: AbstractDensityModel
    density_kg_m3::Float64
    temperature_k::Float64
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

const _POLYFIT_DEFAULT_VALID_MIN_ALTITUDE_M = 50.0e3
const _POLYFIT_DEFAULT_VALID_MAX_ALTITUDE_M = 2_000.0e3
const _POLYFIT_LOG_DENSITY_MIN = log(nextfloat(0.0))
const _POLYFIT_LOG_DENSITY_MAX = log(floatmax(Float64))

"""
Analytic atmosphere model that fits `log(rho)` as a polynomial in altitude
measured in kilometers.

`valid_min_altitude_m` and `valid_max_altitude_m` document the fit domain and
are enforced as a clamp during evaluation so the polynomial does not blow up
under extrapolation. The resulting log-density is also bounded to the finite
`exp` range of `Float64`.
"""
struct PolynomialFitAtmosphereModel <: AbstractDensityModel
    polyfit_coeffs::Vector{Float64}
    valid_min_altitude_m::Float64
    valid_max_altitude_m::Float64
end

function PolynomialFitAtmosphereModel(
    polyfit_coeffs::AbstractVector{<:Real};
    valid_min_altitude_m::Real=_POLYFIT_DEFAULT_VALID_MIN_ALTITUDE_M,
    valid_max_altitude_m::Real=_POLYFIT_DEFAULT_VALID_MAX_ALTITUDE_M
)
    valid_min_f = Float64(valid_min_altitude_m)
    valid_max_f = Float64(valid_max_altitude_m)
    valid_min_f <= valid_max_f || throw(ArgumentError(
        "PolynomialFitAtmosphereModel valid_min_altitude_m must be <= valid_max_altitude_m, got $valid_min_f > $valid_max_f."
    ))
    return PolynomialFitAtmosphereModel(Float64.(collect(polyfit_coeffs)), valid_min_f, valid_max_f)
end

@inline function _planet_polyfit_valid_min_altitude_m(planet::AbstractPlanet)::Float64
    if hasproperty(planet, :polyfit_valid_min_altitude_m)
        return Float64(getproperty(planet, :polyfit_valid_min_altitude_m))
    end
    return _POLYFIT_DEFAULT_VALID_MIN_ALTITUDE_M
end

@inline function _planet_polyfit_valid_max_altitude_m(planet::AbstractPlanet)::Float64
    if hasproperty(planet, :polyfit_valid_max_altitude_m)
        return Float64(getproperty(planet, :polyfit_valid_max_altitude_m))
    end
    return _POLYFIT_DEFAULT_VALID_MAX_ALTITUDE_M
end

@inline function PolynomialFitAtmosphereModel(planet::AbstractPlanet; kwargs...)
    return PolynomialFitAtmosphereModel(
        collect(getproperty(planet, :polyfit_coeffs));
        valid_min_altitude_m=_planet_polyfit_valid_min_altitude_m(planet),
        valid_max_altitude_m=_planet_polyfit_valid_max_altitude_m(planet),
        kwargs...
    )
end

const _NRLMSISE00_DEFAULT_F107A = 150.0
const _NRLMSISE00_DEFAULT_F107 = 150.0
const _NRLMSISE00_DEFAULT_AP = 4.0
const _NRLMSISE00_DEFAULT_VALID_MIN_ALTITUDE_M = 0.0
const _NRLMSISE00_DEFAULT_VALID_MAX_ALTITUDE_M = 1_000.0e3
const _NRLMSISE00_LOW_ALTITUDE_DEFAULT_MAX_M = 80.0e3
const _NRLMSISE00_SPACE_INDICES_LOCK = ReentrantLock()
const _NRLMSISE00_SPACE_INDICES_READY = Ref(false)

"""
Built-in NRLMSISE-00 index provider backed by CelesTrak F10.7 and Ap data via
`SpaceIndices`.

This provider is intended to be installed through
`NRLMSISE00AtmosphereModel(use_space_indices=true)`. It uses:

- `F10.7a`: adjusted centered 81-day average flux
- `F10.7`: adjusted previous-day daily flux
- `Ap`: full 7-slot NRLMSISE history vector built from CelesTrak 3-hour bins

Below `80 km`, it returns the standard NRLMSISE fallback values `150 / 150 / 4`
without touching the dataset.
"""
mutable struct NRLMSISE00SpaceIndicesProvider
    force_download::Bool
end

"""
    init_nrlmsise_space_indices!(; force_download=false)

Initialize the CelesTrak space-weather dataset used by
`NRLMSISE00AtmosphereModel(use_space_indices=true)`.

This is optional because the dataset will also be initialized lazily on the
first atmosphere evaluation that needs it. Call this function before long runs
or Monte Carlo campaigns if you want any download/refresh work to happen before
the solver starts.
"""
function init_nrlmsise_space_indices!(; force_download::Bool=false)
    lock(_NRLMSISE00_SPACE_INDICES_LOCK) do
        if force_download || !_NRLMSISE00_SPACE_INDICES_READY[]
            SatelliteToolbox.AtmosphericModels.SpaceIndices.init(
                SatelliteToolbox.AtmosphericModels.SpaceIndices.Celestrak;
                force_download=force_download
            )
            _NRLMSISE00_SPACE_INDICES_READY[] = true
        end
    end
    return nothing
end

"""
NRLMSISE-00 atmosphere model with either fixed geophysical indices or a callable
index provider.

`f107a`, `f107`, and `ap` set the explicit solar and geomagnetic inputs. If
`index_provider` is not `nothing`, it is called at evaluation time and must
return either `(f107a, f107, ap)` or a named tuple with `f107a`, `f107`, and
`ap`. The provider may accept either `(instant)` or `(instant, h, lat, lon)`.

Set `use_space_indices=true` to install the built-in CelesTrak-backed provider
for F10.7 and Ap data. That path will lazily initialize the `SpaceIndices`
dataset on first use, or you can prewarm it with `init_nrlmsise_space_indices!`.

`valid_min_altitude_m` and `valid_max_altitude_m` document the standard model
validity range. The standard NRLMSISE-00 validity band is approximately
`0 m` to `1000 km`, but evaluation is not artificially clamped to that range.
"""
struct NRLMSISE00AtmosphereModel{A, P} <: AbstractDensityModel
    f107a::Float64
    f107::Float64
    ap::A
    index_provider::P
    include_anomalous_oxygen::Bool
    valid_min_altitude_m::Float64
    valid_max_altitude_m::Float64
end

@inline function _nrlmsise_ap_value(ap::Real)::Float64
    value = Float64(ap)
    isfinite(value) && value >= 0.0 || throw(ArgumentError("NRLMSISE00AtmosphereModel ap must be finite and >= 0.0, got $value."))
    return value
end

function _nrlmsise_ap_value(ap::AbstractVector{<:Real})
    length(ap) == 7 || throw(ArgumentError("NRLMSISE00AtmosphereModel vector ap input must have length 7, got $(length(ap))."))
    values = ntuple(i -> Float64(ap[i]), 7)
    all(x -> isfinite(x) && x >= 0.0, values) || throw(ArgumentError(
        "NRLMSISE00AtmosphereModel vector ap inputs must all be finite and >= 0.0."
    ))
    return SVector{7, Float64}(values)
end

@inline function _nrlmsise_flux_value(value::Real, name::AbstractString)::Float64
    value_f = Float64(value)
    isfinite(value_f) && value_f >= 0.0 || throw(ArgumentError(
        "NRLMSISE00AtmosphereModel $name must be finite and >= 0.0, got $value_f."
    ))
    return value_f
end

function NRLMSISE00AtmosphereModel(;
    f107a::Real=_NRLMSISE00_DEFAULT_F107A,
    f107::Real=_NRLMSISE00_DEFAULT_F107,
    ap::Union{Real, AbstractVector{<:Real}}=_NRLMSISE00_DEFAULT_AP,
    index_provider=nothing,
    use_space_indices::Bool=false,
    space_indices_force_download::Bool=false,
    include_anomalous_oxygen::Bool=true,
    valid_min_altitude_m::Real=_NRLMSISE00_DEFAULT_VALID_MIN_ALTITUDE_M,
    valid_max_altitude_m::Real=_NRLMSISE00_DEFAULT_VALID_MAX_ALTITUDE_M
)
    if use_space_indices && index_provider !== nothing
        throw(ArgumentError(
            "NRLMSISE00AtmosphereModel `use_space_indices=true` cannot be combined with a custom `index_provider`."
        ))
    end
    if !use_space_indices && space_indices_force_download
        throw(ArgumentError(
            "NRLMSISE00AtmosphereModel `space_indices_force_download=true` requires `use_space_indices=true`."
        ))
    end
    valid_min_f = Float64(valid_min_altitude_m)
    valid_max_f = Float64(valid_max_altitude_m)
    valid_min_f <= valid_max_f || throw(ArgumentError(
        "NRLMSISE00AtmosphereModel valid_min_altitude_m must be <= valid_max_altitude_m, got $valid_min_f > $valid_max_f."
    ))
    return NRLMSISE00AtmosphereModel(
        _nrlmsise_flux_value(f107a, "f107a"),
        _nrlmsise_flux_value(f107, "f107"),
        _nrlmsise_ap_value(ap),
        use_space_indices ? NRLMSISE00SpaceIndicesProvider(space_indices_force_download) : index_provider,
        include_anomalous_oxygen,
        valid_min_f,
        valid_max_f
    )
end

@inline function Base.getproperty(model::GRAMAtmosphereModel, name::Symbol)
    if name === :core || name === :instance_lock
        return getfield(model, name)
    end
    return getproperty(getfield(model, :core), name)
end

@inline function Base.propertynames(model::GRAMAtmosphereModel, private::Bool=false)
    wrapped = propertynames(getfield(model, :core), private)
    return (:core, :instance_lock, wrapped...)
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

"""
    _gram_lock_scope() -> Symbol

Which lock serializes native GRAM density calls: `:global` (default) contends on
the single `RuntimeServices.GRAM_LOCK` process-wide, `:model`
(`SPACEAGORA_GRAM_LOCK_SCOPE=model`) contends only on the wrapper's own
`instance_lock`, so distinct [`GRAMAtmosphereModel`](@ref) instances evaluate
concurrently. `:model` relies on the same safety premise as the isolated-pool
batch path (`SPACEAGORA_GRAM_ISOLATED_POOL`): independent GRAM model instances
may be called concurrently as long as each single instance is serialized.
Threaded Monte Carlo campaigns whose samples build their own GRAM models need
`:model` to scale; with `:global` all samples serialize on one lock.
"""
@inline function _gram_lock_scope()::Symbol
    raw = lowercase(strip(get(ENV, "SPACEAGORA_GRAM_LOCK_SCOPE", "global")))
    if raw == "" || raw == "global"
        return :global
    elseif raw == "model" || raw == "per_model" || raw == "per-model" || raw == "instance"
        return :model
    end
    throw(ArgumentError(
        "Unsupported SPACEAGORA_GRAM_LOCK_SCOPE='$raw'. Use one of: global, model."
    ))
end

# These Refs are set by SpaceAGORAGRAMSuiteExt when GRAMSuite is loaded.
const _GRAM_USE_GLOBAL_LOCK_FN = Ref{Function}(() -> false)
const _GRAM_DEFAULT_SURROGATE_FILE_FN = Ref{Function}(_ -> "")

@inline _gram_use_global_lock()::Bool = _GRAM_USE_GLOBAL_LOCK_FN[]()

@inline function _gram_default_surrogate_file(planet::String)::String
    return _GRAM_DEFAULT_SURROGATE_FILE_FN[](planet)
end

function _gram_not_loaded_error(fn_name::String)
    error(
        "GRAMSuite.jl is not loaded. " *
        "Add GRAMSuite as a dependency and load it before calling `$(fn_name)`."
    )
end

function precompute_gram_static_grids!(::AbstractDensityModel; kwargs...)
    _gram_not_loaded_error("precompute_gram_static_grids!")
end

# Zero-arg cache-clear functions use Callable Refs so the extension can wire them in
# without triggering method-overwrite warnings.
const _CLEAR_GRAM_STATIC_GRID_CACHE_FN = Ref{Function}(
    () -> _gram_not_loaded_error("clear_gram_static_grid_cache!")
)
const _CLEAR_GRAM_OFFLINE_SURROGATE_CACHE_FN = Ref{Function}(
    () -> _gram_not_loaded_error("clear_gram_offline_surrogate_cache!")
)

clear_gram_static_grid_cache!() = _CLEAR_GRAM_STATIC_GRID_CACHE_FN[]()
clear_gram_offline_surrogate_cache!() = _CLEAR_GRAM_OFFLINE_SURROGATE_CACHE_FN[]()

function _gram_core_density_state(
    _core,
    _h::Float64,
    _lat::Float64,
    _lon::Float64,
    _el_time::Float64,
    _wind::Bool,
    _lock_obj,
    _vacuum_temperature::Float64
)::Tuple{Float64, Float64, SVector{3, Float64}}
    _gram_not_loaded_error("GRAMAtmosphereModel density evaluation")
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

@inline function _exponential_density(ρ_ref::Float64, h_ref::Float64, H::Float64, h::Float64)::Float64
    return ρ_ref * exp((h_ref - h) / H)
end

@inline function _piecewise_layer_index(model::PiecewiseExponentialAtmosphereModel, h::Float64)::Int
    idx = searchsortedlast(model.h_breaks_m, h)
    return clamp(idx, 1, length(model.ρ_refs))
end

@inline function _polyfit_eval_altitude_km(model::PolynomialFitAtmosphereModel, h_m::Float64)::Float64
    return clamp(h_m, model.valid_min_altitude_m, model.valid_max_altitude_m) * 1e-3
end

@inline function _polyfit_log_density(model::PolynomialFitAtmosphereModel, h_m::Float64)::Float64
    coeffs = model.polyfit_coeffs
    isempty(coeffs) && return 0.0
    h_km = _polyfit_eval_altitude_km(model, h_m)
    exponent = coeffs[1]
    @inbounds for j in 2:length(coeffs)
        exponent = muladd(exponent, h_km, coeffs[j])
    end
    return clamp(exponent, _POLYFIT_LOG_DENSITY_MIN, _POLYFIT_LOG_DENSITY_MAX)
end

@inline function _polyfit_density(model::PolynomialFitAtmosphereModel, h_m::Float64)::Float64
    return exp(_polyfit_log_density(model, h_m))
end

@inline function _nrlmsise_eval_datetime(initial_time, el_time::Float64)::DateTime
    base = DateTime(
        Int(initial_time.year),
        Int(initial_time.month),
        Int(initial_time.day),
        Int(initial_time.hour),
        Int(initial_time.minute),
        0
    ) + Millisecond(round(Int, 1000 * Float64(initial_time.second)))
    return base + Millisecond(round(Int, 1000 * el_time))
end

@inline function _nrlmsise_space_indices_lookup(index::Val, instant::DateTime)
    return SatelliteToolbox.AtmosphericModels.SpaceIndices.space_index(index, instant)
end

@inline function _nrlmsise_space_indices_f107(lookup, instant::DateTime)::Float64
    return _nrlmsise_flux_value(lookup(Val(:F10adj), instant - Day(1)), "f107")
end

@inline function _nrlmsise_space_indices_f107a(lookup, instant::DateTime)::Float64
    return _nrlmsise_flux_value(lookup(Val(:F10adj_avg_center81), instant), "f107a")
end

@inline function _nrlmsise_ap_slot_index(instant::DateTime)::Int
    return fld(Dates.hour(instant), 3) + 1
end

function _nrlmsise_ap_bins(values)::NTuple{8, Float64}
    length(values) == 8 || throw(ArgumentError(
        "NRLMSISE00SpaceIndicesProvider Ap lookup must return 8 three-hour values, got $(length(values))."
    ))
    bins = ntuple(i -> Float64(values[i]), 8)
    all(x -> isfinite(x) && x >= 0.0, bins) || throw(ArgumentError(
        "NRLMSISE00SpaceIndicesProvider Ap lookup must return finite nonnegative values."
    ))
    return bins
end

@inline function _nrlmsise_space_indices_ap_bin(lookup, instant::DateTime)::Float64
    bins = _nrlmsise_ap_bins(lookup(Val(:Ap), instant))
    return bins[_nrlmsise_ap_slot_index(instant)]
end

function _nrlmsise_space_indices_ap_vector(lookup, instant::DateTime)::SVector{7, Float64}
    daily = _nrlmsise_ap_value(lookup(Val(:Ap_daily), instant))
    current = _nrlmsise_space_indices_ap_bin(lookup, instant)
    prev3 = _nrlmsise_space_indices_ap_bin(lookup, instant - Hour(3))
    prev6 = _nrlmsise_space_indices_ap_bin(lookup, instant - Hour(6))
    prev9 = _nrlmsise_space_indices_ap_bin(lookup, instant - Hour(9))
    avg12_33 = sum(_nrlmsise_space_indices_ap_bin(lookup, instant - Hour(hours)) for hours in 12:3:33) / 8.0
    avg36_57 = sum(_nrlmsise_space_indices_ap_bin(lookup, instant - Hour(hours)) for hours in 36:3:57) / 8.0
    return SVector{7, Float64}(daily, current, prev3, prev6, prev9, avg12_33, avg36_57)
end

function _nrlmsise_space_indices_indices(lookup, instant::DateTime)
    return (
        f107a=_nrlmsise_space_indices_f107a(lookup, instant),
        f107=_nrlmsise_space_indices_f107(lookup, instant),
        ap=_nrlmsise_space_indices_ap_vector(lookup, instant),
    )
end

@inline function _nrlmsise_provider_indices(result::Tuple)
    length(result) == 3 || throw(ArgumentError(
        "NRLMSISE00AtmosphereModel index_provider tuple result must have length 3: (f107a, f107, ap)."
    ))
    return (
        f107a=_nrlmsise_flux_value(result[1], "f107a"),
        f107=_nrlmsise_flux_value(result[2], "f107"),
        ap=_nrlmsise_ap_value(result[3]),
    )
end

function _nrlmsise_provider_indices(result::NamedTuple)
    names = propertynames(result)
    f107a_key = :f107a in names ? :f107a : (:F10A in names ? :F10A : nothing)
    f107_key = :f107 in names ? :f107 : (:F10 in names ? :F10 : nothing)
    ap_key = :ap in names ? :ap : (:Ap in names ? :Ap : nothing)
    isnothing(f107a_key) && throw(ArgumentError("NRLMSISE00AtmosphereModel index_provider named tuple result must contain `f107a` or `F10A`."))
    isnothing(f107_key) && throw(ArgumentError("NRLMSISE00AtmosphereModel index_provider named tuple result must contain `f107` or `F10`."))
    isnothing(ap_key) && throw(ArgumentError("NRLMSISE00AtmosphereModel index_provider named tuple result must contain `ap` or `Ap`."))
    return (
        f107a=_nrlmsise_flux_value(getproperty(result, f107a_key), "f107a"),
        f107=_nrlmsise_flux_value(getproperty(result, f107_key), "f107"),
        ap=_nrlmsise_ap_value(getproperty(result, ap_key)),
    )
end

@inline function _nrlmsise_provider_indices(result)
    throw(ArgumentError(
        "NRLMSISE00AtmosphereModel index_provider must return a named tuple or a 3-tuple: (f107a, f107, ap)."
    ))
end

function (provider::NRLMSISE00SpaceIndicesProvider)(
    instant::DateTime,
    h::Float64,
    lat::Float64,
    lon::Float64
)
    if h < _NRLMSISE00_LOW_ALTITUDE_DEFAULT_MAX_M
        return (
            f107a=_NRLMSISE00_DEFAULT_F107A,
            f107=_NRLMSISE00_DEFAULT_F107,
            ap=_NRLMSISE00_DEFAULT_AP,
        )
    end
    init_nrlmsise_space_indices!(; force_download=provider.force_download)
    provider.force_download = false
    return _nrlmsise_space_indices_indices(_nrlmsise_space_indices_lookup, instant)
end

function _nrlmsise_resolved_indices(
    model::NRLMSISE00AtmosphereModel,
    instant::DateTime,
    h::Float64,
    lat::Float64,
    lon::Float64
)
    provider = model.index_provider
    if provider === nothing
        return (f107a=model.f107a, f107=model.f107, ap=model.ap)
    end
    result = if applicable(provider, instant, h, lat, lon)
        provider(instant, h, lat, lon)
    elseif applicable(provider, instant)
        provider(instant)
    else
        throw(ArgumentError(
            "NRLMSISE00AtmosphereModel index_provider must accept `(instant)` or `(instant, h, lat, lon)`."
        ))
    end
    return _nrlmsise_provider_indices(result)
end

@inline function _nrlmsise_density_state(
    model::NRLMSISE00AtmosphereModel,
    instant::DateTime,
    h::Float64,
    lat::Float64,
    lon::Float64
)::Tuple{Float64, Float64, SVector{3, Float64}}
    indices = _nrlmsise_resolved_indices(model, instant, h, lat, lon)
    atmo = SatelliteToolbox.AtmosphericModels.nrlmsise00(
        instant,
        h,
        lat,
        lon,
        indices.f107a,
        indices.f107,
        indices.ap;
        include_anomalous_oxygen=model.include_anomalous_oxygen
    )
    return atmo.total_density, atmo.temperature, SVector{3, Float64}(0.0, 0.0, 0.0)
end

"""
TabulatedFlightAtmosphereModel — density from flight-measured per-pass
profiles (e.g. accelerometer-derived), for labeled replication/certification
diagnostics. Pure data: nearest pass in elapsed time, leg by side of that
pass's periapsis, log-linear interpolation in altitude, exponential tails
beyond each profile's coverage. `sigma_scale` shifts every point by that many
measurement standard errors (envelope runs). Temperature is derived from the
local log-density slope (scale height); winds are zero — the measurement
carries neither. Serves as the digital-twin regression sentinel and the
certification reference; predictive-atmosphere results come from the GRAM
configurations and are reported separately.
"""
struct TabulatedFlightAtmosphereModel <: AbstractDensityModel
    pass_peri_el_s::Vector{Float64}            # sorted periapsis elapsed times
    # per pass: inbound/outbound profiles as (alt_m ascending, log_rho, sigma_log)
    prof_alt_m::Vector{NTuple{2, Vector{Float64}}}
    prof_logrho::Vector{NTuple{2, Vector{Float64}}}
    prof_siglog::Vector{NTuple{2, Vector{Float64}}}
    sigma_scale::Float64
    g_ref_mps2::Float64
    gas_constant::Float64
end

# Tail scale heights are clamped to the physical Mars thermosphere range:
# noisy edge bins can invert (top of a profile is noise-dominated), and an
# unclamped slope would extend a near-constant density blanket along the whole
# orbit — integrated over an 18 h period that phantom drag exceeds the real
# periapsis pass severalfold.
const _TAB_FLIGHT_H_MIN_M = 2000.0
const _TAB_FLIGHT_H_MAX_M = 12000.0

@inline function _tab_flight_interp(alts::Vector{Float64}, logs::Vector{Float64}, sigs::Vector{Float64}, h::Float64, sigma_scale::Float64)::Tuple{Float64, Float64}
    n = length(alts)
    if h <= alts[1]
        H = n >= 2 ? (alts[2] - alts[1]) / max(logs[1] - logs[2], 1e-9) : 8000.0
        H = clamp(H, _TAB_FLIGHT_H_MIN_M, _TAB_FLIGHT_H_MAX_M)
        return (logs[1] + (alts[1] - h) / H + sigma_scale * sigs[1]), H
    elseif h >= alts[n]
        H = n >= 2 ? (alts[n] - alts[n-1]) / max(logs[n-1] - logs[n], 1e-9) : 8000.0
        H = clamp(H, _TAB_FLIGHT_H_MIN_M, _TAB_FLIGHT_H_MAX_M)
        return (logs[n] - (h - alts[n]) / H + sigma_scale * sigs[n]), H
    end
    j = searchsortedlast(alts, h)
    t = (h - alts[j]) / (alts[j+1] - alts[j])
    lr = logs[j] + t * (logs[j+1] - logs[j])
    sg = sigs[j] + t * (sigs[j+1] - sigs[j])
    H = (alts[j+1] - alts[j]) / max(logs[j] - logs[j+1], 1e-9)
    return (lr + sigma_scale * sg), H
end

function getDensity(model::TabulatedFlightAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}
    # Adaptive trial steps can probe a non-finite state (deep-impact blowup
    # before the termination callback fires, seen on +1 sigma envelope runs);
    # a NaN altitude passes both tail branches of _tab_flight_interp and
    # indexes past the profile end. Vacuum lets the solver reject the step.
    isfinite(h) && isfinite(el_time) || return 0.0, 150.0, SVector{3, Float64}(0.0, 0.0, 0.0)
    ts = model.pass_peri_el_s
    j = clamp(searchsortedlast(ts, el_time), 1, length(ts))
    if j < length(ts) && abs(ts[j+1] - el_time) < abs(el_time - ts[j])
        j += 1
    end
    leg = el_time < ts[j] ? 1 : 2
    alts = model.prof_alt_m[j][leg]
    if isempty(alts)
        leg = leg == 1 ? 2 : 1
        alts = model.prof_alt_m[j][leg]
    end
    isempty(alts) && return 0.0, 150.0, SVector{3, Float64}(0.0, 0.0, 0.0)
    logrho, H_m = _tab_flight_interp(alts, model.prof_logrho[j][leg], model.prof_siglog[j][leg], h, model.sigma_scale)
    rho = exp(logrho)
    T = clamp(H_m * model.g_ref_mps2 / model.gas_constant, 80.0, 400.0)
    return rho, T, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline function getDensity(model::TabulatedFlightAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    return getDensity(model, h, lat, lon, el_time, wind)
end

"""
TimeTabulatedAtmosphereModel — density as a pure function of scenario elapsed
time, log-linearly interpolated from a sorted (t_s, rho) table and held
constant beyond the table's ends. The orbital-arc counterpart of
[`TabulatedFlightAtmosphereModel`](@ref) (which is keyed by altitude within
periapsis passes): use it to replay a flight-inferred along-track density
history, an assimilated product sampled along the trajectory, or any other
externally supplied rho(t), e.g. as the drag channel of a digital-twin run.
`scale` multiplies every density; `temperature_k` is reported as the constant
gas temperature; winds are zero — the table carries neither.
"""
struct TimeTabulatedAtmosphereModel <: AbstractDensityModel
    t_el_s::Vector{Float64}     # sorted elapsed times from the scenario epoch
    log_rho::Vector{Float64}    # log density at each node [log(kg/m^3)]
    scale::Float64
    temperature_k::Float64
end

function TimeTabulatedAtmosphereModel(
    t_el_s::AbstractVector{<:Real},
    rho_kgm3::AbstractVector{<:Real};
    scale::Real=1.0,
    temperature_k::Real=900.0
)
    n = length(t_el_s)
    n == length(rho_kgm3) || throw(ArgumentError(
        "TimeTabulatedAtmosphereModel time and density vectors must have equal length, got $n vs $(length(rho_kgm3))."
    ))
    n >= 2 || throw(ArgumentError("TimeTabulatedAtmosphereModel needs at least 2 table nodes, got $n."))
    t = Float64.(t_el_s)
    issorted(t) || throw(ArgumentError("TimeTabulatedAtmosphereModel times must be sorted ascending."))
    all(isfinite, t) || throw(ArgumentError("TimeTabulatedAtmosphereModel times must be finite."))
    rho = Float64.(rho_kgm3)
    all(x -> isfinite(x) && x > 0.0, rho) || throw(ArgumentError(
        "TimeTabulatedAtmosphereModel densities must be finite and > 0."
    ))
    scale_f = Float64(scale)
    isfinite(scale_f) && scale_f > 0.0 || throw(ArgumentError(
        "TimeTabulatedAtmosphereModel scale must be finite and > 0, got $scale_f."
    ))
    temp_f = Float64(temperature_k)
    isfinite(temp_f) && temp_f > 0.0 || throw(ArgumentError(
        "TimeTabulatedAtmosphereModel temperature_k must be finite and > 0, got $temp_f."
    ))
    return TimeTabulatedAtmosphereModel(t, log.(rho), scale_f, temp_f)
end

function getDensity(model::TimeTabulatedAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}
    ts = model.t_el_s
    logs = model.log_rho
    n = length(ts)
    local logrho
    if el_time <= ts[1]
        logrho = logs[1]
    elseif el_time >= ts[n]
        logrho = logs[n]
    else
        j = searchsortedlast(ts, el_time)
        f = (el_time - ts[j]) / (ts[j+1] - ts[j])
        logrho = logs[j] + f * (logs[j+1] - logs[j])
    end
    return model.scale * exp(logrho), model.temperature_k, SVector{3, Float64}(0.0, 0.0, 0.0)
end

@inline function getDensity(model::TimeTabulatedAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    return getDensity(model, h, lat, lon, el_time, wind)
end

function getDensity(model::NoAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    rho = 0.0
    T = p.args.environment_model.planet.T_ref
    wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    return rho, T, wind_vec
end

function getDensity(model::ExponentialAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}
    rho = _exponential_density(model.ρ_ref, model.h_ref, model.H, h)
    T = model.temperature_k
    wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    return rho, T, wind_vec
end

@inline function getDensity(model::ExponentialAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    return getDensity(model, h, lat, lon, el_time, wind)
end

function getDensity(model::PiecewiseExponentialAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}
    idx = _piecewise_layer_index(model, h)
    rho = _exponential_density(model.ρ_refs[idx], model.h_refs[idx], model.Hs[idx], h)
    T = model.temperature_k
    wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    return rho, T, wind_vec
end

@inline function getDensity(model::PiecewiseExponentialAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    return getDensity(model, h, lat, lon, el_time, wind)
end

@inline function _batch_elapsed_time(el_time::Float64, idx::Int)::Float64
    return el_time
end

@inline function _batch_elapsed_time(el_time::AbstractVector{<:Real}, idx::Int)::Float64
    return Float64(el_time[idx])
end

@inline function _validate_density_batch_lengths(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}}
)::Int
    n = length(hs)
    length(rhos) == n || throw(ArgumentError("getDensityBatch!: rhos length $(length(rhos)) must match hs length $n."))
    length(Ts) == n || throw(ArgumentError("getDensityBatch!: Ts length $(length(Ts)) must match hs length $n."))
    length(winds) == n || throw(ArgumentError("getDensityBatch!: winds length $(length(winds)) must match hs length $n."))
    length(lats) == n || throw(ArgumentError("getDensityBatch!: lats length $(length(lats)) must match hs length $n."))
    length(lons) == n || throw(ArgumentError("getDensityBatch!: lons length $(length(lons)) must match hs length $n."))
    if el_time isa AbstractVector{<:Real}
        length(el_time) == n || throw(ArgumentError("getDensityBatch!: el_time vector length $(length(el_time)) must match hs length $n."))
    end
    return n
end

@inline function _density_scalar_for_batch(
    model,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return getDensity(model, h, lat, lon, el_time, wind, p)
end

@inline function _density_scalar_for_batch(
    model::ExponentialAtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return getDensity(model, h, lat, lon, el_time, wind)
end

@inline function _density_scalar_for_batch(
    model::NRLMSISE00AtmosphereModel,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool,
    p
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return getDensity(model, h, lat, lon, el_time, wind, p)
end

function getDensityBatch!(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    model::NoAtmosphereModel,
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}},
    wind::Bool,
    p::params
)::Nothing where params
    n = _validate_density_batch_lengths(rhos, Ts, winds, hs, lats, lons, el_time)
    T_ref = p.args.environment_model.planet.T_ref
    zero_wind = SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for i in 1:n
        rhos[i] = 0.0
        Ts[i] = T_ref
        winds[i] = zero_wind
    end
    return nothing
end

function getDensityBatch!(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    model::ExponentialAtmosphereModel,
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}},
    wind::Bool,
    p::params
)::Nothing where params
    n = _validate_density_batch_lengths(rhos, Ts, winds, hs, lats, lons, el_time)
    zero_wind = SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for i in 1:n
        h = Float64(hs[i])
        rhos[i] = _exponential_density(model.ρ_ref, model.h_ref, model.H, h)
        Ts[i] = model.temperature_k
        winds[i] = zero_wind
    end
    return nothing
end

function getDensityBatch!(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    model::PiecewiseExponentialAtmosphereModel,
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}},
    wind::Bool,
    p::params
)::Nothing where params
    n = _validate_density_batch_lengths(rhos, Ts, winds, hs, lats, lons, el_time)
    zero_wind = SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for i in 1:n
        h = Float64(hs[i])
        idx = _piecewise_layer_index(model, h)
        rhos[i] = _exponential_density(model.ρ_refs[idx], model.h_refs[idx], model.Hs[idx], h)
        Ts[i] = model.temperature_k
        winds[i] = zero_wind
    end
    return nothing
end

function getDensityBatch!(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    model::PolynomialFitAtmosphereModel,
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}},
    wind::Bool,
    p::params
)::Nothing where params
    n = _validate_density_batch_lengths(rhos, Ts, winds, hs, lats, lons, el_time)
    T_ref = p.args.environment_model.planet.T_ref
    zero_wind = SVector{3, Float64}(0.0, 0.0, 0.0)
    @inbounds for i in 1:n
        rhos[i] = _polyfit_density(model, Float64(hs[i]))
        Ts[i] = T_ref
        winds[i] = zero_wind
    end
    return nothing
end

function getDensityBatch!(
    rhos::AbstractVector{Float64},
    Ts::AbstractVector{Float64},
    winds::AbstractVector{SVector{3, Float64}},
    model::AbstractDensityModel,
    hs::AbstractVector{<:Real},
    lats::AbstractVector{<:Real},
    lons::AbstractVector{<:Real},
    el_time::Union{Float64, AbstractVector{<:Real}},
    wind::Bool,
    p::params
)::Nothing where params
    n = _validate_density_batch_lengths(rhos, Ts, winds, hs, lats, lons, el_time)
    @inbounds for i in 1:n
        h = Float64(hs[i])
        lat = Float64(lats[i])
        lon = Float64(lons[i])
        t_i = _batch_elapsed_time(el_time, i)
        rho_i, T_i, wind_i = _density_scalar_for_batch(model, h, lat, lon, t_i, wind, p)
        rhos[i] = rho_i
        Ts[i] = T_i
        winds[i] = wind_i
    end
    return nothing
end

function getDensity(model::PolynomialFitAtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    rho = _polyfit_density(model, h)
    T = p.args.environment_model.planet.T_ref
    wind_vec = SVector{3, Float64}(0.0, 0.0, 0.0)
    return rho, T, wind_vec
end

function getDensity(model::PolynomialFitAtmosphereModel, h::AbstractVector{<:Real}, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params) where params
    rho = Vector{Float64}(undef, length(h))
    @inbounds for i in eachindex(h)
        rho[i] = _polyfit_density(model, Float64(h[i]))
    end
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
    model::GRAMAtmosphereModelSurrogate,
    h::Float64,
    lat::Float64,
    lon::Float64,
    el_time::Float64,
    wind::Bool
)::Tuple{Float64, Float64, SVector{3, Float64}}
    return _gram_point_density(model.base_model, h, lat, lon, el_time, wind)
end

function getDensity(model::ConstantDensityModel, altitude_m::Float64, latitude_deg::Float64, longitude_deg::Float64, et::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    return model.density_kg_m3, model.temperature_k, SVector{3, Float64}(0.0, 0.0, 0.0)
end

# NRLMSISE-00 is calendar-dependent (space indices, solar position, seasonal
# terms), so evaluating it without a scenario epoch is a silent physics error:
# the historical fallback anchored el_time to the J2000 epoch, which put every
# caller on this path decades away from its mission dates. Callers must use the
# 7-arg method, whose epoch comes from `p.args.initial_time`.
function getDensity(model::NRLMSISE00AtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool)::Tuple{Float64, Float64, SVector{3, Float64}}
    throw(ArgumentError(
        "getDensity(::NRLMSISE00AtmosphereModel, h, lat, lon, el_time, wind) has no scenario epoch: " *
        "el_time would be anchored to the J2000 reference epoch and NRLMSISE-00 would evaluate " *
        "space indices for the wrong dates. Call the 7-argument method " *
        "getDensity(model, h, lat, lon, el_time, wind, p), which maps el_time from p.args.initial_time."
    ))
end

function getDensity(model::NRLMSISE00AtmosphereModel, h::Float64, lat::Float64, lon::Float64, el_time::Float64, wind::Bool, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    instant = _nrlmsise_eval_datetime(p.args.initial_time, el_time)
    return _nrlmsise_density_state(model, instant, h, lat, lon)
end

function density_polyfit(h::Float64, p::params)::Tuple{Float64, Float64, SVector{3, Float64}} where params
    poly_model = PolynomialFitAtmosphereModel(p.args.environment_model.planet)
    rho, T, wind_vec = getDensity(poly_model, h, 0.0, 0.0, 0.0, false, p)
    return rho, T, wind_vec
end
