module SpaceAGORAGRAMSuiteExt

using SpaceAGORA
using GRAMSuite
using StaticArrays
using SPICE
using Serialization

const EM = SpaceAGORA.SimulationModel.EnvironmentModels
# Spelled out rather than routed through `RS` below: ci_architecture_contract_gate
# asserts this exact reference textually, because the one thing that must never
# drift here is which lock this is. Same reason the gate exists at all.
const GRAM_LOCK = SpaceAGORA.RuntimeServices.GRAM_LOCK
const RS = SpaceAGORA.RuntimeServices

# Attributed views of the one shared native lock. Same critical section, same
# mutual exclusion; the site only decides which occupancy counter the time lands
# in. See RuntimeServices' native-lock section for why occupancy is attributed
# by call site rather than by density-model family.
@inline _tl(site::Symbol) = RS.tracked_lock(site)

# Lock used to serialize native GRAM calls for this model instance. With the
# default global scope every call in the process contends on GRAM_LOCK; with
# SPACEAGORA_GRAM_LOCK_SCOPE=model only calls on the same wrapper instance
# serialize, so per-sample/per-worker model copies evaluate concurrently
# (the same instance-isolation premise as the isolated-pool batch path).
#
# Return type is left unannotated because the two branches now differ: the
# global branch hands back a tracked view of the shared lock, the per-model
# branch the instance's own `ReentrantLock`. GRAMSuite's `lock_obj` keyword is
# untyped and reaches it through `lock(f, l)`, so both work; the small Union
# splits at the call site.
#
# Per-model occupancy is deliberately NOT recorded: those calls do not contend
# on the shared lock at all, which is the entire point of that mode, so their
# absence from the shared counters is the correct reading. What is lost is the
# per-instance duty cycle, which needs its own counters if that mode is ever
# measured rather than merely offered.
@inline function _gram_call_lock(model::EM.GRAMAtmosphereModel)
    return EM._gram_lock_scope() === :model ? model.instance_lock : _tl(:gram_density)
end

function __init__()
    EM._GRAM_USE_GLOBAL_LOCK_FN[] = () -> GRAMSuite.gram_use_global_lock()
    EM._GRAM_DEFAULT_SURROGATE_FILE_FN[] = planet -> GRAMSuite.gram_default_surrogate_file(planet)
    EM._CLEAR_GRAM_STATIC_GRID_CACHE_FN[] = () -> GRAMSuite.clear_gram_static_grid_cache!()
    EM._CLEAR_GRAM_OFFLINE_SURROGATE_CACHE_FN[] = () -> GRAMSuite.clear_gram_offline_surrogate_cache!()
    GRAMSuite._GRAM_EPHEMERIS_STATE_FN[] = _gram_spice_ephemeris_state
    # See the comment on RuntimeServices.GRAM_LOCK/SPICE_LOCK: libGRAM.dylib's
    # statically-linked CSPICE globally exports the same internal symbol names
    # as SpaceAGORA's own SPICE.jl bindings, so GRAM model construction (the
    # only unlocked-by-default native-call path inside GRAMSuite.jl) must
    # serialize against the same process-wide lock as every other CSPICE-
    # touching call, not just against other GRAM construction calls.
    GRAMSuite._GRAM_DEFAULT_LOCK_HOOK[] = GRAM_LOCK
end

# ---------------------------------------------------------------------------
# GRAM ephemeris-state bypass (see the long comment above
# GRAMSuite._gram_apply_user_ephemeris_state! for why this exists: GRAM's
# vendored native library has its own private, isolated CSPICE instance whose
# default kernels don't work for non-Earth bodies). Reimplements GRAM's own
# common/source/Ephemeris.cpp formulas using SpaceAGORA's working SPICE.jl
# bindings. Units/conventions matched to EphemerisStateC's fields: solarTime
# in hours, longitudeSun/subsolarLatitude/subsolarLongitude/solarZenithAngle
# in degrees, orbitalRadius in AU, oneWayLightTime in minutes, secondsPerSol
# in seconds.
# ---------------------------------------------------------------------------

const _GRAM_AU_KM = 149_597_870.7
const _GRAM_EARTH_NAIF_ID = 399

# From NSSDCA Planetary Fact Sheet -- literal constants copied from
# Ephemeris::updateSecondsPerSol() in the vendored C++ (they are treated as
# fixed per-body values there too, not recomputed from ephemeris).
const _GRAM_SECONDS_PER_SOL = Dict(
    "VENUS" => 1.00872e7,
    "EARTH" => 86400.00,
    "MARS" => 88774.92,
    "JUPITER" => 35733.24,
    "URANUS" => 62064.0,
    "NEPTUNE" => 57996.0,
    "SATURN" => 38361.6,
    "TITAN" => 1377648.0,
)

function _gram_utc_string(initial_time)::String
    return string(
        Int(initial_time.year), "-", Int(initial_time.month), "-", Int(initial_time.day), " ",
        Int(initial_time.hour), ":", Int(initial_time.minute), ":", Float64(initial_time.second),
        " UTC"
    )
end

# ---------------------------------------------------------------------------
# Epoch-level ephemeris cache.
#
# Six of the eight EphemerisStateC fields -- longitudeSun, subsolarLatitude,
# subsolarLongitude, orbitalRadius, oneWayLightTime, secondsPerSol -- depend
# only on (body, epoch). They do NOT depend on the query lat/lon. Only
# solarTime and solarZenithAngle do, and both are closed-form arithmetic on the
# subsolar point.
#
# Every satellite in a constellation is evaluated at the same `el_time` within
# one right-hand-side call, so without a cache the five SPICE calls below
# (utc2et, spkpos, lspcn, ltime, subslr) run once per satellite per stage and
# serialize on the SPICE lock. Caching the epoch-level result collapses that to
# once per epoch -- the same "compute shared environment once, amortize across
# the constellation" mechanism already used for the N-body/SRP prepass.
#
# The cache holds ONE entry per body (the most recent epoch), not a growing
# map: `el_time` advances every step, so an unbounded Dict keyed on epoch would
# grow without limit over a mission. One entry is the right size because the
# access pattern is "all satellites at epoch t, then all satellites at t'".
#
# Set SPACEAGORA_GRAM_EPHEMERIS_CACHE=off to bypass (for A/B checking); the
# cached and uncached paths compute bit-identical values.
#
# Locking: the cache is guarded by GRAM_LOCK, which is `const GRAM_LOCK =
# SPICE_LOCK` (runtime_services.jl:15) -- the same object that serializes every
# other CSPICE call in the process. That matters because the miss path calls
# SPICE.jl (utc2et/spkpos/lspcn/ltime/subslr), and CSPICE is not thread-safe:
# guarding with a private lock would serialize misses against each other but
# not against SpaceAGORA's own SRP/N-body SPICE calls on another thread. Using
# the one shared lock also removes any lock-ordering hazard, since it is
# reentrant and is already held on entry under the default global lock scope.
# ---------------------------------------------------------------------------

const _GRAM_EPHEM_CACHE = Dict{String, Tuple{Float64, NTuple{6, Float64}}}()
const _GRAM_ET0_CACHE = Dict{String, Float64}()
const _GRAM_NAIF_ID_CACHE = Dict{String, Int}()

@inline function _gram_ephemeris_cache_enabled()::Bool
    raw = lowercase(strip(get(ENV, "SPACEAGORA_GRAM_EPHEMERIS_CACHE", "on")))
    raw in ("on", "1", "true", "yes") && return true
    raw in ("off", "0", "false", "no") && return false
    throw(ArgumentError(
        "Unsupported SPACEAGORA_GRAM_EPHEMERIS_CACHE='$raw'. Use one of: on, off."
    ))
end

"""
    clear_gram_ephemeris_cache!()

Drop the cached epoch-level ephemeris, start-epoch and NAIF-id lookups. Only
needed when SPICE kernels are re-furnished mid-process.
"""
function clear_gram_ephemeris_cache!()
    lock(_tl(:gram_setup)) do
        empty!(_GRAM_EPHEM_CACHE)
        empty!(_GRAM_ET0_CACHE)
        empty!(_GRAM_NAIF_ID_CACHE)
    end
    return nothing
end

# Start epoch and NAIF id are fixed for the whole run; utc2et parses a string
# and bodn2c does a name lookup, so neither belongs in a per-call path.
@inline function _gram_et0(utc::String)::Float64
    lock(_tl(:spice_body)) do
        get!(_GRAM_ET0_CACHE, utc) do
            SPICE.utc2et(utc)
        end
    end
end

@inline function _gram_naif_id(naif_name::String)::Int
    lock(_tl(:spice_body)) do
        get!(_GRAM_NAIF_ID_CACHE, naif_name) do
            Int(SPICE.bodn2c(naif_name))
        end
    end
end

# The epoch-only half: longitudeSun, subsolarLat, subsolarLon, orbitalRadius,
# oneWayLightTime, secondsPerSol.
function _gram_body_ephemeris_uncached(naif_name::String, et::Float64)::NTuple{6, Float64}
    frame = "IAU_" * naif_name

    pos_sun, _ = SPICE.spkpos(naif_name, et, "J2000", "NONE", "SUN")
    orbital_radius_au = sqrt(sum(abs2, pos_sun)) / _GRAM_AU_KM

    longitude_sun_deg = mod(rad2deg(SPICE.lspcn(naif_name, et, "NONE")), 360.0)

    _, howlng = SPICE.ltime(et, _gram_naif_id(naif_name), "->", _GRAM_EARTH_NAIF_ID)
    one_way_light_time_min = howlng / 60.0

    spoint, _, _ = SPICE.subslr("NEAR POINT/ELLIPSOID", naif_name, et, frame, "NONE", naif_name)
    _, subsolar_lon, subsolar_lat = SPICE.reclat(spoint)
    subsolar_lon_deg = mod(rad2deg(subsolar_lon), 360.0)
    subsolar_lat_deg = rad2deg(subsolar_lat)

    return (
        longitude_sun_deg, subsolar_lat_deg, subsolar_lon_deg,
        orbital_radius_au, one_way_light_time_min, _GRAM_SECONDS_PER_SOL[naif_name]
    )
end

function _gram_body_ephemeris(naif_name::String, et::Float64)::NTuple{6, Float64}
    _gram_ephemeris_cache_enabled() || return _gram_body_ephemeris_uncached(naif_name, et)
    lock(_tl(:spice_ephemeris)) do
        hit = get(_GRAM_EPHEM_CACHE, naif_name, nothing)
        # Exact epoch match only: the value is a nonlinear function of et, so
        # interpolating or tolerating a nearby epoch would silently change the
        # atmosphere. Every satellite in one RHS call shares et exactly.
        if hit !== nothing && hit[1] === et
            return hit[2]
        end
        value = _gram_body_ephemeris_uncached(naif_name, et)
        _GRAM_EPHEM_CACHE[naif_name] = (et, value)
        return value
    end
end

function _gram_spice_ephemeris_state(
    planet_name::String,
    initial_time,
    el_time::Float64,
    lat_deg::Float64,
    lon_deg::Float64,
)
    naif_name = uppercase(planet_name)
    haskey(_GRAM_SECONDS_PER_SOL, naif_name) || return nothing

    et = _gram_et0(_gram_utc_string(initial_time)) + el_time

    longitude_sun_deg, subsolar_lat_deg, subsolar_lon_deg,
        orbital_radius_au, one_way_light_time_min, seconds_per_sol =
            _gram_body_ephemeris(naif_name, et)

    # Local solar time & solar zenith angle from the subsolar point vs. the
    # query lat/lon (standard spherical-geometry relations; matches GRAM's
    # own updateLocalSolarTime()/updateSolarZenithAngle()). These are the only
    # per-satellite terms, and they are pure arithmetic -- no SPICE.
    hour_angle_deg = mod(lon_deg - subsolar_lon_deg + 180.0, 360.0) - 180.0
    solar_time_hr = mod(12.0 + hour_angle_deg / 15.0, 24.0)

    lat_r, sublat_r = deg2rad(lat_deg), deg2rad(subsolar_lat_deg)
    dlon_r = deg2rad(lon_deg - subsolar_lon_deg)
    cos_zenith = sin(lat_r) * sin(sublat_r) + cos(lat_r) * cos(sublat_r) * cos(dlon_r)
    solar_zenith_deg = rad2deg(acos(clamp(cos_zenith, -1.0, 1.0)))

    return (
        solar_time_hr, longitude_sun_deg, subsolar_lat_deg, subsolar_lon_deg,
        orbital_radius_au, one_way_light_time_min, solar_zenith_deg, seconds_per_sol
    )
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
        lock_obj=_tl(:gram_setup)
    )
end


# ---------------------------------------------------------------------------
# deepcopy
# ---------------------------------------------------------------------------

function Base.deepcopy_internal(model::EM.GRAMAtmosphereModel, stackdict::IdDict)
    haskey(stackdict, model) && return stackdict[model]
    copied = lock(_tl(:gram_setup)) do
        EM.GRAMAtmosphereModel(deepcopy(model.core))
    end
    stackdict[model] = copied
    return copied
end

function Base.deepcopy_internal(model::EM.GRAMAtmosphereModelSurrogate, stackdict::IdDict)
    haskey(stackdict, model) && return stackdict[model]
    copied = lock(_tl(:gram_setup)) do
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
# Custom Serialization -- same rationale as GRAMSuite.jl's own
# GRAMAtmosphereModel/GRAMAtmosphereModelSurrogate methods: `core` wraps a live
# native handle, and `instance_lock` is a ReentrantLock that must never be
# serialized as-is (Task/condition-variable state has no meaning on another
# process, mirroring why deepcopy_internal above never copies it either).
# Serialize only `core` (which recurses into GRAMSuite's own serialize method)
# and reconstruct with a fresh instance_lock on the receiving side -- the same
# "fresh lock per construction" contract the single-arg constructor already
# documents. This lets a SimulationConfiguration carrying one of these models
# cross a Distributed process boundary (e.g. remotecall) transparently.
# ---------------------------------------------------------------------------

function Serialization.serialize(s::Serialization.AbstractSerializer, model::EM.GRAMAtmosphereModel)
    Serialization.writetag(s.io, Serialization.OBJECT_TAG)
    Serialization.serialize(s, EM.GRAMAtmosphereModel)
    Serialization.serialize(s, model.core)
    return nothing
end

function Serialization.deserialize(s::Serialization.AbstractSerializer, ::Type{EM.GRAMAtmosphereModel})
    core = Serialization.deserialize(s)
    return EM.GRAMAtmosphereModel(core)
end

function Serialization.serialize(s::Serialization.AbstractSerializer, model::EM.GRAMAtmosphereModelSurrogate)
    Serialization.writetag(s.io, Serialization.OBJECT_TAG)
    Serialization.serialize(s, EM.GRAMAtmosphereModelSurrogate)
    Serialization.serialize(s, model.base_model)
    Serialization.serialize(s, model.surrogate_file)
    Serialization.serialize(s, model.point_fallback_below_m)
    return nothing
end

function Serialization.deserialize(s::Serialization.AbstractSerializer, ::Type{EM.GRAMAtmosphereModelSurrogate})
    base_model = Serialization.deserialize(s)
    surrogate_file = Serialization.deserialize(s)
    point_fallback_below_m = Serialization.deserialize(s)
    return EM.GRAMAtmosphereModelSurrogate(base_model, surrogate_file, point_fallback_below_m)
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
    h_gram = max(h, -30.0)
    return GRAMSuite.point_density_state(model.core, h_gram, lat, lon, el_time, wind; lock_obj=_gram_call_lock(model))
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
        h_gram = max(h, -30.0)
        rho, T, wind_vec = GRAMSuite.density_state(
            model.core,
            h_gram,
            lat,
            lon,
            el_time,
            wind;
            lock_obj=_gram_call_lock(model),
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
    lock_obj = model.base_model isa EM.GRAMAtmosphereModel ? _gram_call_lock(model.base_model) : _tl(:gram_density)
    h_gram = max(h, -30.0)

    # println("GRAM density altitude = $(h) m ($(h / 1e3) km)")
    return GRAMSuite.surrogate_density_state(
        base_model,
        model.surrogate_file,
        model.point_fallback_below_m,
        h_gram,
        lat,
        lon,
        el_time,
        wind;
        lock_obj=lock_obj,
        point_density_fallback=point_fallback,
        vacuum_temperature=p.args.environment_model.planet.T_ref
    )
end

end # module SpaceAGORAGRAMSuiteExt
