module GRAM

export BODY_NO, BODY_VENUS, BODY_EARTH, BODY_MARS, BODY_JUPITER, BODY_SATURN, BODY_URANUS, BODY_NEPTUNE, BODY_TITAN
export GAS_ARGON, GAS_HELIUM, GAS_HYDROGEN, GAS_DIHYDROGEN, GAS_NITROGEN, GAS_DINITROGEN, GAS_OXYGEN, GAS_DIOXYGEN
export GAS_METHANE, GAS_CARBON_MONOXIDE, GAS_CARBON_DIOXIDE, GAS_OZONE, GAS_NITROUS_OXIDE, GAS_WATER
export Atmosphere, PositionC, DynamicsStateC, DensityStateC, WindsStateC, GasesStateC, ConstituentGasC, PerturbationStateC, EphemerisStateC, TrajectoryPointC
export MarsStateC, DailyDynamicsStateC
export set_library!, initialize!, create_atmosphere, copy_atmosphere, set_start_time!, set_delta!, set_position!, update!, get_start_time, get_position, get_dynamics_state, get_density_state, get_winds_state, get_gases_state, get_constituent_gas, get_ephemeris_state, get_perturbation_state
export enable_console_output!, try_get_spice_path, try_get_data_paths, set_spice_lsk!, set_spice_pck!, set_spice_kernel!, load_spice_file!
export set_seed!, set_min_relative_step_size!, set_perturbation_scales!, set_perturbation_action!, add_auxiliary_atmosphere!, set_auxiliary_values!, set_ephemeris_state!, set_ephemeris_fast_mode!, set_subsolar_update_time!
export generate_trajectory, generate_monte_carlo_trajectories
export set_map_year!, set_mgcm_dust_levels!, set_dust_storm!, set_dust_density!, set_f107!, set_wind_scales!, set_mola_heights!, set_min_max!
export set_planetary_radii!, set_height_offset_model!, set_height_above_surface!, set_perturbation_wavelength_scale!, set_exospheric_temperature!, set_wave_defaults!, set_wave_file!
export set_areoid_radius_callback!, set_topographic_height_callback!, set_callback_data!
export get_mars_state, get_daily_dynamics_state, get_mars_gases_state, get_version_string
export mars_try_get_spice_path, mars_try_get_data_paths, mars_set_spice_lsk!, mars_set_spice_pck!, mars_set_spice_kernel!, mars_initialize!, mars_load_spice_file!, create_mars_atmosphere, copy_mars_atmosphere, mars_update!, mars_get_error_message
export get_error_message, close!

const BODY_NO = 0
const BODY_VENUS = 1
const BODY_EARTH = 2
const BODY_MARS = 3
const BODY_JUPITER = 4
const BODY_SATURN = 5
const BODY_URANUS = 6
const BODY_NEPTUNE = 7
const BODY_TITAN = 8

const GAS_ARGON = 0
const GAS_HELIUM = 1
const GAS_HYDROGEN = 2
const GAS_DIHYDROGEN = 3
const GAS_NITROGEN = 4
const GAS_DINITROGEN = 5
const GAS_OXYGEN = 6
const GAS_DIOXYGEN = 7
const GAS_METHANE = 8
const GAS_CARBON_MONOXIDE = 9
const GAS_CARBON_DIOXIDE = 10
const GAS_OZONE = 11
const GAS_NITROUS_OXIDE = 12
const GAS_WATER = 13

struct GramTimeC
    year::Cint
    month::Cint
    day::Cint
    hour::Cint
    minute::Cint
    seconds::Cdouble
    scale::Cint
    frame::Cint
end

struct PositionC
    height::Cdouble
    latitude::Cdouble
    longitude::Cdouble
    elapsedTime::Cdouble
    isPlanetoCentric::Cint
    totalRadius::Cdouble
    latitudeRadius::Cdouble
    surfaceHeight::Cdouble
    gravity::Cdouble
end

struct DynamicsStateC
    temperature::Cdouble
    pressure::Cdouble
    density::Cdouble
    pressureScaleHeight::Cdouble
    densityScaleHeight::Cdouble
    speedOfSound::Cdouble
    referenceDensity::Cdouble
    referenceTemperature::Cdouble
    referencePressure::Cdouble
    sigmaLevel::Cdouble
    pressureAtSurface::Cdouble
    pressureAltitude::Cdouble
end

struct DensityStateC
    density::Cdouble
    lowDensity::Cdouble
    highDensity::Cdouble
    densityPerturbation::Cdouble
    perturbedDensity::Cdouble
    densityStandardDeviation::Cdouble
    perturbedSpeedOfSound::Cdouble
    relativeStepSize::Cdouble
    referenceDensity::Cdouble
    densityDeviation::Cdouble
    lowDensityDeviation::Cdouble
    highDensityDeviation::Cdouble
    perturbedDensityDeviation::Cdouble
    updateStatus::Cint
end

struct WindsStateC
    ewWind::Cdouble
    nsWind::Cdouble
    verticalWind::Cdouble
    ewWindPerturbation::Cdouble
    nsWindPerturbation::Cdouble
    verticalWindPerturbation::Cdouble
    perturbedEWWind::Cdouble
    perturbedNSWind::Cdouble
    perturbedVerticalWind::Cdouble
    ewStandardDeviation::Cdouble
    nsStandardDeviation::Cdouble
    verticalStandardDeviation::Cdouble
end

struct GasesStateC
    totalNumberDensity::Cdouble
    averageMolecularWeight::Cdouble
    compressibilityFactor::Cdouble
    specificGasConstant::Cdouble
    specificHeatRatio::Cdouble
end

struct ConstituentGasC
    moleFraction::Cdouble
    massFraction::Cdouble
    numberDensity::Cdouble
    averageMolecularWeight::Cdouble
    specificHeatCapacity::Cdouble
end

struct PerturbationStateC
    densityRandomNumber::Cdouble
    ewWindRandomNumber::Cdouble
    nsWindRandomNumber::Cdouble
    verticalWindRandomNumber::Cdouble
end

struct EphemerisStateC
    solarTime::Cdouble
    longitudeSun::Cdouble
    subsolarLatitude::Cdouble
    subsolarLongitude::Cdouble
    orbitalRadius::Cdouble
    oneWayLightTime::Cdouble
    solarZenithAngle::Cdouble
    secondsPerSol::Cdouble
end

struct TrajectoryPointC
    position::PositionC
    dynamics::DynamicsStateC
    winds::WindsStateC
    ephemeris::EphemerisStateC
end

struct MarsStateC
    planetoGraphicHeight::Cdouble
    planetoGraphicLatitude::Cdouble
    referenceHeight::Cdouble
    referenceRadius::Cdouble
    groundTemperature::Cdouble
    thermosphereBaseHeight::Cdouble
    thermosphereBaseTemperature::Cdouble
    exosphericTemperature::Cdouble
    f1PeakHeight::Cdouble
    albedo::Cdouble
    heightOffset::Cdouble
    localHeightOffset::Cdouble
    dustOpticalDepth::Cdouble
    dustColumnArealDensity::Cdouble
    dustMixingRatio::Cdouble
    dustMassDensity::Cdouble
    dustNumberDensity::Cdouble
    iceIsPresent::Cint
end

struct DailyDynamicsStateC
    temperatureDaily::Cdouble
    pressureDaily::Cdouble
    densityDaily::Cdouble
    ewWindDaily::Cdouble
    nsWindDaily::Cdouble
    densityMin::Cdouble
    densityMax::Cdouble
    temperatureMin::Cdouble
    temperatureMax::Cdouble
end

mutable struct Atmosphere
    ptr::Ptr{Cvoid}
    body::Cint
    function Atmosphere(ptr::Ptr{Cvoid}, body::Integer = BODY_NO)
        ptr == C_NULL && error("createAtmosphere_C returned NULL.")
        obj = new(ptr, Cint(body))
        finalizer(obj) do x
            if x.ptr != C_NULL
                ccall((:deleteAtmosphere_C, _libgram()), Cvoid, (Ptr{Cvoid},), x.ptr)
                x.ptr = C_NULL
            end
        end
        return obj
    end
end

const _LIBGRAM = Ref{String}("")
const _SPICE_PATH = Ref{String}("")
const _REGISTERED_BODIES = Set{Int}()

function _default_library_path()
    ext = Sys.iswindows() ? "dll" : (Sys.isapple() ? "dylib" : "so")
    return normpath(joinpath(@__DIR__, "..", "Build", "lib", "libGRAM.$ext"))
end

function _libgram()
    if !isempty(_LIBGRAM[])
        return _LIBGRAM[]
    end
    if haskey(ENV, "GRAM_LIB") && !isempty(ENV["GRAM_LIB"])
        set_library!(ENV["GRAM_LIB"])
        return _LIBGRAM[]
    end
    candidate = _default_library_path()
    if isfile(candidate)
        _LIBGRAM[] = candidate
        return candidate
    end
    error("GRAM shared library not found. Call set_library!(...) or set GRAM_LIB.")
end

function set_library!(path::AbstractString)
    abs = abspath(path)
    isfile(abs) || error("GRAM shared library does not exist: $abs")
    _LIBGRAM[] = abs
    return _LIBGRAM[]
end

function initialize!(spice_path::AbstractString)
    _SPICE_PATH[] = String(spice_path)
    empty!(_REGISTERED_BODIES)
    ccall((:initialize_C, _libgram()), Cvoid, (Cstring,), spice_path)
    return nothing
end

function enable_console_output!(flag::Bool = true)
    ccall((:enableConsoleOutput_C, _libgram()), Cvoid, (Cint,), flag ? Cint(1) : Cint(0))
    return nothing
end

function try_get_spice_path(buffer_size::Integer = 2048)
    n = Int(max(buffer_size, 2))
    buf = Vector{UInt8}(undef, n)
    ccall((:tryGetSpicePath_C, _libgram()), Cvoid, (Ptr{UInt8}, Cint), buf, Cint(n))
    i = findfirst(==(0x00), buf)
    if i === nothing
        return String(buf)
    end
    i == 1 && return ""
    return String(buf[1:(i - 1)])
end

function try_get_data_paths(buffer_size::Integer = 2048)
    n = Int(max(buffer_size, 2))
    spice_buf = Vector{UInt8}(undef, n)
    data_buf = Vector{UInt8}(undef, n)
    ccall((:tryGetDataPaths_C, _libgram()), Cvoid, (Ptr{UInt8}, Ptr{UInt8}, Cint), spice_buf, data_buf, Cint(n))
    spice_i = findfirst(==(0x00), spice_buf)
    data_i = findfirst(==(0x00), data_buf)
    spice = spice_i === nothing ? String(spice_buf) : (spice_i == 1 ? "" : String(spice_buf[1:(spice_i - 1)]))
    data = data_i === nothing ? String(data_buf) : (data_i == 1 ? "" : String(data_buf[1:(data_i - 1)]))
    return spice, data
end

function set_spice_lsk!(lsk::AbstractString)
    ccall((:setSpiceLsk_C, _libgram()), Cvoid, (Cstring,), lsk)
    return nothing
end

function set_spice_pck!(pck::AbstractString)
    ccall((:setSpicePck_C, _libgram()), Cvoid, (Cstring,), pck)
    return nothing
end

function set_spice_kernel!(body::Integer, bsp::AbstractString)
    ccall((:setSpiceKernel_C, _libgram()), Cvoid, (Cint, Cstring), Cint(body), bsp)
    return nothing
end

function load_spice_file!(file_name::AbstractString)
    ccall((:loadSpiceFile_C, _libgram()), Cvoid, (Cstring,), file_name)
    return nothing
end

mars_try_get_spice_path(args...; kwargs...) = try_get_spice_path(args...; kwargs...)
mars_try_get_data_paths(args...; kwargs...) = try_get_data_paths(args...; kwargs...)
mars_set_spice_lsk!(args...; kwargs...) = set_spice_lsk!(args...; kwargs...)
mars_set_spice_pck!(args...; kwargs...) = set_spice_pck!(args...; kwargs...)

function mars_set_spice_kernel!(bsp::AbstractString)
    set_spice_kernel!(BODY_MARS, bsp)
    return nothing
end

function mars_initialize!(spice_path::AbstractString)
    ccall((:initialize_M, _libgram()), Cvoid, (Cstring,), spice_path)
    push!(_REGISTERED_BODIES, BODY_MARS)
    return nothing
end

function mars_load_spice_file!(file_name::AbstractString)
    ccall((:loadSpiceFile_M, _libgram()), Cvoid, (Cstring,), file_name)
    return nothing
end

function _initialize_body!(body::Integer)
    body == BODY_NO && throw(ArgumentError("BODY_NO is not a creatable atmosphere model."))
    body == BODY_SATURN && throw(ArgumentError("SATURN is not available in this GRAM Suite build."))
    body in _REGISTERED_BODIES && return
    isempty(_SPICE_PATH[]) && error("Call initialize!(spice_path) before create_atmosphere.")
    if body == BODY_VENUS
        ccall((:initialize_V, _libgram()), Cvoid, (Cstring,), _SPICE_PATH[])
    elseif body == BODY_EARTH
        ccall((:initialize_E, _libgram()), Cvoid, (Cstring,), _SPICE_PATH[])
    elseif body == BODY_MARS
        ccall((:initialize_M, _libgram()), Cvoid, (Cstring,), _SPICE_PATH[])
    elseif body == BODY_JUPITER
        ccall((:initialize_J, _libgram()), Cvoid, (Cstring,), _SPICE_PATH[])
    elseif body == BODY_NEPTUNE
        ccall((:initialize_N, _libgram()), Cvoid, (Cstring,), _SPICE_PATH[])
    elseif body == BODY_TITAN
        ccall((:initialize_T, _libgram()), Cvoid, (Cstring,), _SPICE_PATH[])
    elseif body == BODY_URANUS
        ccall((:initialize_U, _libgram()), Cvoid, (Cstring,), _SPICE_PATH[])
    else
        throw(ArgumentError("Unsupported GRAM body: $body"))
    end
    push!(_REGISTERED_BODIES, Int(body))
end

function create_atmosphere(body::Integer; data_path::Union{Nothing, AbstractString} = nothing)
    b = Int(body)
    _initialize_body!(b)
    data_arg = data_path === nothing ? Ptr{UInt8}(C_NULL) : String(data_path)
    ptr = if b == BODY_EARTH
        ccall((:createAtmosphere_E, _libgram()), Ptr{Cvoid}, (Cstring,), data_arg)
    elseif b == BODY_MARS
        ccall((:createAtmosphere_M, _libgram()), Ptr{Cvoid}, (Cstring,), data_arg)
    elseif b == BODY_VENUS
        ccall((:createAtmosphere_V, _libgram()), Ptr{Cvoid}, ())
    elseif b == BODY_JUPITER
        ccall((:createAtmosphere_J, _libgram()), Ptr{Cvoid}, ())
    elseif b == BODY_NEPTUNE
        ccall((:createAtmosphere_N, _libgram()), Ptr{Cvoid}, ())
    elseif b == BODY_TITAN
        ccall((:createAtmosphere_T, _libgram()), Ptr{Cvoid}, ())
    elseif b == BODY_URANUS
        ccall((:createAtmosphere_U, _libgram()), Ptr{Cvoid}, ())
    else
        ccall((:createAtmosphere_C, _libgram()), Ptr{Cvoid}, (Cint,), Cint(b))
    end
    ptr == C_NULL && error("Failed to create atmosphere for body $b: $(get_error_message())")
    return Atmosphere(ptr, b)
end

function create_mars_atmosphere(; data_path::Union{Nothing, AbstractString} = nothing)
    return create_atmosphere(BODY_MARS; data_path = data_path)
end

function copy_atmosphere(atmos::Atmosphere)
    ptr = ccall((:copyAtmosphere_C, _libgram()), Ptr{Cvoid}, (Ptr{Cvoid},), atmos.ptr)
    ptr == C_NULL && error("copyAtmosphere_C failed: $(get_error_message())")
    return Atmosphere(ptr, Int(atmos.body))
end

function copy_mars_atmosphere(atmos::Atmosphere)
    _require_mars!(atmos)
    ptr = ccall((:copyAtmosphere_M, _libgram()), Ptr{Cvoid}, (Ptr{Cvoid},), atmos.ptr)
    ptr == C_NULL && error("copyAtmosphere_M failed: $(get_error_message())")
    return Atmosphere(ptr, BODY_MARS)
end

function close!(atmos::Atmosphere)
    if atmos.ptr != C_NULL
        ccall((:deleteAtmosphere_C, _libgram()), Cvoid, (Ptr{Cvoid},), atmos.ptr)
        atmos.ptr = C_NULL
    end
    return nothing
end

function set_start_time!(atmos::Atmosphere; year::Integer, month::Integer, day::Integer,
    hour::Integer = 0, minute::Integer = 0, seconds::Real = 0.0, scale::Integer = 1, frame::Integer = 1)
    t = Ref(GramTimeC(Cint(year), Cint(month), Cint(day), Cint(hour), Cint(minute), Cdouble(seconds), Cint(scale), Cint(frame)))
    ccall((:setStartTime_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GramTimeC}), atmos.ptr, t)
    return nothing
end

function set_start_time!(atmos::Atmosphere, t::GramTimeC)
    ccall((:setStartTime_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GramTimeC}), atmos.ptr, Ref(t))
    return nothing
end

function set_position!(atmos::Atmosphere; height::Real, latitude::Real, longitude::Real,
    elapsed_time::Real = 0.0, is_planetocentric::Bool = true)
    p = Ref(PositionC(
        Cdouble(height), Cdouble(latitude), Cdouble(longitude), Cdouble(elapsed_time),
        is_planetocentric ? Cint(1) : Cint(0), Cdouble(0), Cdouble(0), Cdouble(0), Cdouble(0)))
    ccall((:setPosition_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PositionC}), atmos.ptr, p)
    return nothing
end

function set_position!(atmos::Atmosphere, p::PositionC)
    ccall((:setPosition_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PositionC}), atmos.ptr, Ref(p))
    return nothing
end

function set_delta!(atmos::Atmosphere; height::Real = 0.0, latitude::Real = 0.0, longitude::Real = 0.0,
    elapsed_time::Real = 0.0, is_planetocentric::Bool = true)
    p = Ref(PositionC(
        Cdouble(height), Cdouble(latitude), Cdouble(longitude), Cdouble(elapsed_time),
        is_planetocentric ? Cint(1) : Cint(0), Cdouble(0), Cdouble(0), Cdouble(0), Cdouble(0)))
    ccall((:setDelta_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PositionC}), atmos.ptr, p)
    return nothing
end

function set_delta!(atmos::Atmosphere, p::PositionC)
    ccall((:setDelta_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PositionC}), atmos.ptr, Ref(p))
    return nothing
end

function set_seed!(atmos::Atmosphere, seed::Integer)
    ccall((:setSeed_C, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, Cint(seed))
    return nothing
end

function set_min_relative_step_size!(atmos::Atmosphere, min_relative_step_size::Real)
    ccall((:setMinRelativeStepSize_C, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble), atmos.ptr, Cdouble(min_relative_step_size))
    return nothing
end

function set_perturbation_scales!(
    atmos::Atmosphere;
    density_scale::Real = 1.0,
    ew_wind_scale::Real = 1.0,
    ns_wind_scale::Real = 1.0,
    vertical_wind_scale::Real = 1.0
)
    ccall(
        (:setPerturbationScales_C, _libgram()),
        Cvoid,
        (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble),
        atmos.ptr,
        Cdouble(density_scale),
        Cdouble(ew_wind_scale),
        Cdouble(ns_wind_scale),
        Cdouble(vertical_wind_scale)
    )
    return nothing
end

function set_perturbation_action!(atmos::Atmosphere, update::Bool)
    ccall((:setPerturbationAction_C, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, update ? Cint(1) : Cint(0))
    return nothing
end

function add_auxiliary_atmosphere!(
    atmos::Atmosphere;
    file_name::AbstractString,
    inner_radius::Real,
    outer_radius::Real,
    is_east_longitude_positive::Bool = true
)
    ccall(
        (:addAuxiliaryAtmosphere_C, _libgram()),
        Cvoid,
        (Ptr{Cvoid}, Cstring, Cdouble, Cdouble, Cint),
        atmos.ptr,
        file_name,
        Cdouble(inner_radius),
        Cdouble(outer_radius),
        is_east_longitude_positive ? Cint(1) : Cint(0)
    )
    return nothing
end

function set_auxiliary_values!(
    atmos::Atmosphere;
    dens::Real,
    pres::Real,
    temp::Real,
    ew::Real,
    ns::Real
)
    ccall(
        (:setAuxiliaryValues_C, _libgram()),
        Cvoid,
        (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble),
        atmos.ptr,
        Cdouble(dens),
        Cdouble(pres),
        Cdouble(temp),
        Cdouble(ew),
        Cdouble(ns)
    )
    return nothing
end

function set_ephemeris_state!(atmos::Atmosphere, state::EphemerisStateC)
    ccall((:setEphemerisState_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{EphemerisStateC}), atmos.ptr, Ref(state))
    return nothing
end

function set_ephemeris_fast_mode!(atmos::Atmosphere, flag::Bool = true)
    ccall((:setEphemerisFastModeOn_C, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, flag ? Cint(1) : Cint(0))
    return nothing
end

function set_subsolar_update_time!(atmos::Atmosphere, utime::Real)
    ccall((:setSubsolarUpdateTime_C, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble), atmos.ptr, Cdouble(utime))
    return nothing
end

function update!(atmos::Atmosphere)
    return Int(ccall((:update_C, _libgram()), Cint, (Ptr{Cvoid},), atmos.ptr))
end

function mars_update!(atmos::Atmosphere)
    _require_mars!(atmos)
    return Int(ccall((:update_M, _libgram()), Cint, (Ptr{Cvoid},), atmos.ptr))
end

function get_start_time(atmos::Atmosphere)
    t = Ref(GramTimeC(0, 0, 0, 0, 0, 0.0, 1, 1))
    ccall((:getStartTime_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GramTimeC}), atmos.ptr, t)
    return t[]
end

function get_position(atmos::Atmosphere)
    p = Ref(PositionC(0, 0, 0, 0, 1, 0, 0, 0, 0))
    ccall((:getPosition_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PositionC}), atmos.ptr, p)
    return p[]
end

function get_dynamics_state(atmos::Atmosphere)
    s = Ref(DynamicsStateC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getDynamicsState_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{DynamicsStateC}), atmos.ptr, s)
    return s[]
end

function get_density_state(atmos::Atmosphere)
    s = Ref(DensityStateC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getDensityState_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{DensityStateC}), atmos.ptr, s)
    return s[]
end

function get_winds_state(atmos::Atmosphere)
    w = Ref(WindsStateC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getWindsState_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{WindsStateC}), atmos.ptr, w)
    return w[]
end

function get_gases_state(atmos::Atmosphere)
    s = Ref(GasesStateC(0, 0, 0, 0, 0))
    ccall((:getGasesState_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GasesStateC}), atmos.ptr, s)
    return s[]
end

function get_constituent_gas(atmos::Atmosphere, gas_type::Integer)
    g = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    ccall((:getConstituentGas_C, _libgram()), Cvoid, (Ptr{Cvoid}, Cint, Ref{ConstituentGasC}), atmos.ptr, Cint(gas_type), g)
    return g[]
end

function get_ephemeris_state(atmos::Atmosphere)
    s = Ref(EphemerisStateC(0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getEphemerisState_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{EphemerisStateC}), atmos.ptr, s)
    return s[]
end

function get_perturbation_state(atmos::Atmosphere)
    s = Ref(PerturbationStateC(0, 0, 0, 0))
    ccall((:getPerturbationState_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PerturbationStateC}), atmos.ptr, s)
    return s[]
end

@inline function _position_c(height::Real, latitude::Real, longitude::Real, elapsed_time::Real, is_planetocentric::Bool)
    return PositionC(
        Cdouble(height), Cdouble(latitude), Cdouble(longitude), Cdouble(elapsed_time),
        is_planetocentric ? Cint(1) : Cint(0),
        Cdouble(0), Cdouble(0), Cdouble(0), Cdouble(0)
    )
end

function generate_trajectory(
    atmos::Atmosphere,
    initial::PositionC,
    delta::PositionC;
    n_points::Integer,
    update_initial_perturbations::Bool = true
)
    n = Int(n_points)
    n >= 1 || throw(ArgumentError("n_points must be >= 1."))

    output = Vector{TrajectoryPointC}(undef, n)
    err = ccall(
        (:generateTrajectory_C, _libgram()),
        Cint,
        (Ptr{Cvoid}, Ref{PositionC}, Ref{PositionC}, Cint, Cint, Ptr{TrajectoryPointC}),
        atmos.ptr,
        Ref(initial),
        Ref(delta),
        Cint(n),
        update_initial_perturbations ? Cint(1) : Cint(0),
        output
    )
    err == 0 || error("GRAM trajectory generation failed (code=$err): $(get_error_message())")
    return output
end

function generate_trajectory(
    atmos::Atmosphere;
    initial_height::Real,
    initial_latitude::Real,
    initial_longitude::Real,
    initial_elapsed_time::Real = 0.0,
    delta_height::Real = 0.0,
    delta_latitude::Real = 0.0,
    delta_longitude::Real = 0.0,
    delta_elapsed_time::Real = 0.0,
    n_points::Integer,
    is_planetocentric::Bool = true,
    update_initial_perturbations::Bool = true
)
    initial = _position_c(initial_height, initial_latitude, initial_longitude, initial_elapsed_time, is_planetocentric)
    delta = _position_c(delta_height, delta_latitude, delta_longitude, delta_elapsed_time, true)
    return generate_trajectory(
        atmos,
        initial,
        delta;
        n_points = n_points,
        update_initial_perturbations = update_initial_perturbations
    )
end

function generate_monte_carlo_trajectories(
    atmos::Atmosphere,
    initial::PositionC,
    delta::PositionC;
    n_points::Integer,
    seeds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    n_runs::Integer = 1,
    initial_seed::Integer = 1001,
    seed_step::Integer = 1,
    update_initial_perturbations::Bool = true
)
    seed_list = if seeds === nothing
        n = Int(n_runs)
        n >= 1 || throw(ArgumentError("n_runs must be >= 1."))
        [initial_seed + (i - 1) * seed_step for i in 1:n]
    else
        collected = Int.(collect(seeds))
        isempty(collected) && throw(ArgumentError("seeds must not be empty."))
        collected
    end

    tracks = Vector{Vector{TrajectoryPointC}}(undef, length(seed_list))
    for i in eachindex(seed_list)
        set_seed!(atmos, seed_list[i])
        tracks[i] = generate_trajectory(
            atmos,
            initial,
            delta;
            n_points = n_points,
            update_initial_perturbations = update_initial_perturbations
        )
    end
    return tracks
end

function generate_monte_carlo_trajectories(
    atmos::Atmosphere;
    initial_height::Real,
    initial_latitude::Real,
    initial_longitude::Real,
    initial_elapsed_time::Real = 0.0,
    delta_height::Real = 0.0,
    delta_latitude::Real = 0.0,
    delta_longitude::Real = 0.0,
    delta_elapsed_time::Real = 0.0,
    n_points::Integer,
    seeds::Union{Nothing, AbstractVector{<:Integer}} = nothing,
    is_planetocentric::Bool = true,
    update_initial_perturbations::Bool = true,
    n_runs::Integer = 1,
    initial_seed::Integer = 1001,
    seed_step::Integer = 1
)
    initial = _position_c(initial_height, initial_latitude, initial_longitude, initial_elapsed_time, is_planetocentric)
    delta = _position_c(delta_height, delta_latitude, delta_longitude, delta_elapsed_time, true)
    return generate_monte_carlo_trajectories(
        atmos,
        initial,
        delta;
        n_points = n_points,
        seeds = seeds,
        n_runs = n_runs,
        initial_seed = initial_seed,
        seed_step = seed_step,
        update_initial_perturbations = update_initial_perturbations
    )
end

@inline function _require_mars!(atmos::Atmosphere)
    Int(atmos.body) == BODY_MARS || throw(ArgumentError("This function is Mars-specific and requires a Mars atmosphere handle."))
    return nothing
end

function set_map_year!(atmos::Atmosphere, year::Integer)
    _require_mars!(atmos)
    ccall((:setMapYear_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, Cint(year))
    return nothing
end

function set_planetary_radii!(atmos::Atmosphere; equatorial::Real, polar::Real)
    _require_mars!(atmos)
    ccall((:setPlanetaryRadii_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble, Cdouble), atmos.ptr, Cdouble(equatorial), Cdouble(polar))
    return nothing
end

function set_height_offset_model!(atmos::Atmosphere; model::Integer, height_offset::Real)
    _require_mars!(atmos)
    ccall((:setHeightOffsetModel_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cint, Cdouble), atmos.ptr, Cint(model), Cdouble(height_offset))
    return nothing
end

function set_height_above_surface!(atmos::Atmosphere, height_above_surface::Real)
    _require_mars!(atmos)
    ccall((:setHeightAboveSurface_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble), atmos.ptr, Cdouble(height_above_surface))
    return nothing
end

function set_mgcm_dust_levels!(atmos::Atmosphere; constant_level::Real, min_level::Real = 0.0, max_level::Real = 0.0)
    _require_mars!(atmos)
    ccall(
        (:setMGCMDustLevels_M, _libgram()),
        Cvoid,
        (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble),
        atmos.ptr,
        Cdouble(constant_level),
        Cdouble(min_level),
        Cdouble(max_level)
    )
    return nothing
end

function set_dust_storm!(
    atmos::Atmosphere;
    longitude_sun::Real,
    duration::Real,
    intensity::Real,
    max_radius::Real,
    latitude::Real,
    longitude::Real
)
    _require_mars!(atmos)
    ccall(
        (:setDustStorm_M, _libgram()),
        Cvoid,
        (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble),
        atmos.ptr,
        Cdouble(longitude_sun),
        Cdouble(duration),
        Cdouble(intensity),
        Cdouble(max_radius),
        Cdouble(latitude),
        Cdouble(longitude)
    )
    return nothing
end

function set_dust_density!(atmos::Atmosphere; nu::Real, diameter::Real, density::Real)
    _require_mars!(atmos)
    ccall((:setDustDensity_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble), atmos.ptr, Cdouble(nu), Cdouble(diameter), Cdouble(density))
    return nothing
end

function set_f107!(atmos::Atmosphere, f107::Real)
    _require_mars!(atmos)
    ccall((:setF107_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble), atmos.ptr, Cdouble(f107))
    return nothing
end

function set_perturbation_wavelength_scale!(atmos::Atmosphere, scale::Real)
    _require_mars!(atmos)
    ccall((:setPerturbationWaveLengthScale_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble), atmos.ptr, Cdouble(scale))
    return nothing
end

function set_wind_scales!(atmos::Atmosphere; mean_winds::Real = 1.0, boundary_layer_winds::Real = 1.0)
    _require_mars!(atmos)
    ccall((:setWindScales_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble, Cdouble), atmos.ptr, Cdouble(mean_winds), Cdouble(boundary_layer_winds))
    return nothing
end

function set_mola_heights!(atmos::Atmosphere, enabled::Bool)
    _require_mars!(atmos)
    ccall((:setMOLAHeights_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, enabled ? Cint(1) : Cint(0))
    return nothing
end

function set_min_max!(atmos::Atmosphere, minmax::Integer)
    _require_mars!(atmos)
    ccall((:setMinMax_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, Cint(minmax))
    return nothing
end

function set_exospheric_temperature!(atmos::Atmosphere; offset::Real, factor::Real)
    _require_mars!(atmos)
    ccall((:setExosphericTemperature_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble, Cdouble), atmos.ptr, Cdouble(offset), Cdouble(factor))
    return nothing
end

function set_wave_defaults!(
    atmos::Atmosphere;
    date::Real,
    scale::Real,
    mean::Real,
    a1::Real,
    p1::Real,
    r1::Real,
    a2::Real,
    p2::Real,
    r2::Real,
    a3::Real,
    p3::Real,
    r3::Real
)
    _require_mars!(atmos)
    ccall(
        (:setWaveDefaults_M, _libgram()),
        Cvoid,
        (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble),
        atmos.ptr,
        Cdouble(date),
        Cdouble(scale),
        Cdouble(mean),
        Cdouble(a1),
        Cdouble(p1),
        Cdouble(r1),
        Cdouble(a2),
        Cdouble(p2),
        Cdouble(r2),
        Cdouble(a3),
        Cdouble(p3),
        Cdouble(r3)
    )
    return nothing
end

function set_wave_file!(atmos::Atmosphere, wave_file::AbstractString)
    _require_mars!(atmos)
    ccall((:setWaveFile_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cstring), atmos.ptr, wave_file)
    return nothing
end

function set_areoid_radius_callback!(atmos::Atmosphere, callback::Ptr{Cvoid})
    _require_mars!(atmos)
    ccall((:setAreoidRadiusCallback_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), atmos.ptr, callback)
    return nothing
end

function set_topographic_height_callback!(atmos::Atmosphere, callback::Ptr{Cvoid})
    _require_mars!(atmos)
    ccall((:setTopographicHeightCallback_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), atmos.ptr, callback)
    return nothing
end

function set_callback_data!(atmos::Atmosphere, data::Ptr{Cvoid})
    _require_mars!(atmos)
    ccall((:setCallbackData_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), atmos.ptr, data)
    return nothing
end

function get_version_string(atmos::Atmosphere; buffer_size::Integer = 2048)
    _require_mars!(atmos)
    n = Int(max(buffer_size, 2))
    buf = Vector{UInt8}(undef, n)
    out_len = ccall((:getVersionString_M, _libgram()), Csize_t, (Ptr{Cvoid}, Ptr{UInt8}, Csize_t), atmos.ptr, buf, Csize_t(n))
    out_len_int = Int(out_len)
    if out_len_int <= 0
        return ""
    end
    last_idx = min(out_len_int, n)
    return String(buf[1:last_idx])
end

function get_daily_dynamics_state(atmos::Atmosphere)
    _require_mars!(atmos)
    s = Ref(DailyDynamicsStateC(0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getDailyDynamicsState_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{DailyDynamicsStateC}), atmos.ptr, s)
    return s[]
end

function get_mars_state(atmos::Atmosphere)
    _require_mars!(atmos)
    s = Ref(MarsStateC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getMarsState_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{MarsStateC}), atmos.ptr, s)
    return s[]
end

function get_mars_gases_state(atmos::Atmosphere)
    _require_mars!(atmos)
    state = Ref(GasesStateC(0, 0, 0, 0, 0))
    argon = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    carbon_dioxide = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    carbon_monoxide = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    dihydrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    dinitrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    dioxygen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    helium = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    hydrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    oxygen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    water = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    ccall(
        (:getGasesState_M, _libgram()),
        Cvoid,
        (
            Ptr{Cvoid},
            Ref{GasesStateC},
            Ref{ConstituentGasC},
            Ref{ConstituentGasC},
            Ref{ConstituentGasC},
            Ref{ConstituentGasC},
            Ref{ConstituentGasC},
            Ref{ConstituentGasC},
            Ref{ConstituentGasC},
            Ref{ConstituentGasC},
            Ref{ConstituentGasC},
            Ref{ConstituentGasC}
        ),
        atmos.ptr,
        state,
        argon,
        carbon_dioxide,
        carbon_monoxide,
        dihydrogen,
        dinitrogen,
        dioxygen,
        helium,
        hydrogen,
        oxygen,
        water
    )
    return (
        state = state[],
        argon = argon[],
        carbon_dioxide = carbon_dioxide[],
        carbon_monoxide = carbon_monoxide[],
        dihydrogen = dihydrogen[],
        dinitrogen = dinitrogen[],
        dioxygen = dioxygen[],
        helium = helium[],
        hydrogen = hydrogen[],
        oxygen = oxygen[],
        water = water[]
    )
end

function mars_get_error_message(buffer_size::Integer = 2048)
    buf = Vector{UInt8}(undef, buffer_size)
    n = ccall((:getErrorMessage_M, _libgram()), Csize_t, (Ptr{UInt8}, Csize_t), buf, Csize_t(buffer_size))
    return String(buf[1:Int(n)])
end

function get_error_message(buffer_size::Integer = 2048)
    buf = Vector{UInt8}(undef, buffer_size)
    n = ccall((:getErrorMessage_C, _libgram()), Csize_t, (Ptr{UInt8}, Csize_t), buf, Csize_t(buffer_size))
    return String(buf[1:Int(n)])
end

end # module GRAM
