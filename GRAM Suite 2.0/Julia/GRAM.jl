module GRAM

export BODY_NO, BODY_VENUS, BODY_EARTH, BODY_MARS, BODY_JUPITER, BODY_SATURN, BODY_URANUS, BODY_NEPTUNE, BODY_TITAN
export Atmosphere, PositionC, DynamicsStateC, WindsStateC
export set_library!, initialize!, create_atmosphere, set_start_time!, set_position!, update!, get_position, get_dynamics_state, get_winds_state
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

mutable struct Atmosphere
    ptr::Ptr{Cvoid}
    function Atmosphere(ptr::Ptr{Cvoid})
        ptr == C_NULL && error("createAtmosphere_C returned NULL.")
        obj = new(ptr)
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
    return Atmosphere(ptr)
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

function set_position!(atmos::Atmosphere; height::Real, latitude::Real, longitude::Real,
    elapsed_time::Real = 0.0, is_planetocentric::Bool = true)
    p = Ref(PositionC(
        Cdouble(height), Cdouble(latitude), Cdouble(longitude), Cdouble(elapsed_time),
        is_planetocentric ? Cint(1) : Cint(0), Cdouble(0), Cdouble(0), Cdouble(0), Cdouble(0)))
    ccall((:setPosition_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PositionC}), atmos.ptr, p)
    return nothing
end

function update!(atmos::Atmosphere)
    return Int(ccall((:update_C, _libgram()), Cint, (Ptr{Cvoid},), atmos.ptr))
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

function get_winds_state(atmos::Atmosphere)
    w = Ref(WindsStateC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getWindsState_C, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{WindsStateC}), atmos.ptr, w)
    return w[]
end

function get_error_message(buffer_size::Integer = 2048)
    buf = Vector{UInt8}(undef, buffer_size)
    n = ccall((:getErrorMessage_C, _libgram()), Csize_t, (Ptr{UInt8}, Csize_t), buf, Csize_t(buffer_size))
    return String(buf[1:Int(n)])
end

end # module GRAM
