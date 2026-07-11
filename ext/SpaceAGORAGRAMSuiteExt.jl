module SpaceAGORAGRAMSuiteExt

using SpaceAGORA
using GRAMSuite
using StaticArrays
using SPICE

const EM = SpaceAGORA.SimulationModel.EnvironmentModels
const GRAM_LOCK = SpaceAGORA.RuntimeServices.GRAM_LOCK

# Lock used to serialize native GRAM calls for this model instance. With the
# default global scope every call in the process contends on GRAM_LOCK; with
# SPACEAGORA_GRAM_LOCK_SCOPE=model only calls on the same wrapper instance
# serialize, so per-sample/per-worker model copies evaluate concurrently
# (the same instance-isolation premise as the isolated-pool batch path).
@inline function _gram_call_lock(model::EM.GRAMAtmosphereModel)::ReentrantLock
    return EM._gram_lock_scope() === :model ? model.instance_lock : GRAM_LOCK
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

function _gram_spice_ephemeris_state(
    planet_name::String,
    initial_time,
    el_time::Float64,
    lat_deg::Float64,
    lon_deg::Float64,
)
    naif_name = uppercase(planet_name)
    haskey(_GRAM_SECONDS_PER_SOL, naif_name) || return nothing

    et = SPICE.utc2et(_gram_utc_string(initial_time)) + el_time
    frame = "IAU_" * naif_name

    pos_sun, _ = SPICE.spkpos(naif_name, et, "J2000", "NONE", "SUN")
    orbital_radius_au = sqrt(sum(abs2, pos_sun)) / _GRAM_AU_KM

    longitude_sun_deg = mod(rad2deg(SPICE.lspcn(naif_name, et, "NONE")), 360.0)

    _, howlng = SPICE.ltime(et, SPICE.bodn2c(naif_name), "->", _GRAM_EARTH_NAIF_ID)
    one_way_light_time_min = howlng / 60.0

    spoint, _, _ = SPICE.subslr("NEAR POINT/ELLIPSOID", naif_name, et, frame, "NONE", naif_name)
    _, subsolar_lon, subsolar_lat = SPICE.reclat(spoint)
    subsolar_lon_deg = mod(rad2deg(subsolar_lon), 360.0)
    subsolar_lat_deg = rad2deg(subsolar_lat)

    # Local solar time & solar zenith angle from the subsolar point vs. the
    # query lat/lon (standard spherical-geometry relations; matches GRAM's
    # own updateLocalSolarTime()/updateSolarZenithAngle()).
    hour_angle_deg = mod(lon_deg - subsolar_lon_deg + 180.0, 360.0) - 180.0
    solar_time_hr = mod(12.0 + hour_angle_deg / 15.0, 24.0)

    lat_r, sublat_r = deg2rad(lat_deg), deg2rad(subsolar_lat_deg)
    dlon_r = deg2rad(lon_deg - subsolar_lon_deg)
    cos_zenith = sin(lat_r) * sin(sublat_r) + cos(lat_r) * cos(sublat_r) * cos(dlon_r)
    solar_zenith_deg = rad2deg(acos(clamp(cos_zenith, -1.0, 1.0)))

    seconds_per_sol = _GRAM_SECONDS_PER_SOL[naif_name]

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
    lock_obj = model.base_model isa EM.GRAMAtmosphereModel ? _gram_call_lock(model.base_model) : GRAM_LOCK
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
