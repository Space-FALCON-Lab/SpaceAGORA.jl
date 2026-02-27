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

function mars_try_get_spice_path(buffer_size::Integer = 2048)
    n = Int(max(buffer_size, 2))
    buf = Vector{UInt8}(undef, n)
    ccall((:tryGetSpicePath_M, _libgram()), Cvoid, (Ptr{UInt8}, Cint), buf, Cint(n))
    i = findfirst(==(0x00), buf)
    if i === nothing
        return String(buf)
    end
    i == 1 && return ""
    return String(buf[1:(i - 1)])
end

function mars_try_get_data_paths(buffer_size::Integer = 2048)
    n = Int(max(buffer_size, 2))
    spice_buf = Vector{UInt8}(undef, n)
    data_buf = Vector{UInt8}(undef, n)
    ccall((:tryGetDataPaths_M, _libgram()), Cvoid, (Ptr{UInt8}, Ptr{UInt8}, Cint), spice_buf, data_buf, Cint(n))
    spice_i = findfirst(==(0x00), spice_buf)
    data_i = findfirst(==(0x00), data_buf)
    spice = spice_i === nothing ? String(spice_buf) : (spice_i == 1 ? "" : String(spice_buf[1:(spice_i - 1)]))
    data = data_i === nothing ? String(data_buf) : (data_i == 1 ? "" : String(data_buf[1:(data_i - 1)]))
    return spice, data
end

function mars_set_spice_lsk!(lsk::AbstractString)
    ccall((:setSpiceLsk_M, _libgram()), Cvoid, (Cstring,), lsk)
    return nothing
end

function mars_set_spice_pck!(pck::AbstractString)
    ccall((:setSpicePck_M, _libgram()), Cvoid, (Cstring,), pck)
    return nothing
end

function mars_set_spice_kernel!(bsp::AbstractString)
    ccall((:setSpiceKernel_M, _libgram()), Cvoid, (Cstring,), bsp)
    return nothing
end

function mars_initialize!(spice_path::AbstractString)
    _SPICE_PATH[] = String(spice_path)
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
        if Int(atmos.body) == BODY_MARS
            ccall((:deleteAtmosphere_M, _libgram()), Cvoid, (Ptr{Cvoid},), atmos.ptr)
        else
            ccall((:deleteAtmosphere_C, _libgram()), Cvoid, (Ptr{Cvoid},), atmos.ptr)
        end
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

@inline function _require_body!(atmos::Atmosphere, body::Integer, name::AbstractString)
    Int(atmos.body) == body || throw(ArgumentError("This function is $(name)-specific and requires a $(name) atmosphere handle."))
    return nothing
end

@inline function _require_mars!(atmos::Atmosphere)
    _require_body!(atmos, BODY_MARS, "Mars")
    return nothing
end

function get_error_message(buffer_size::Integer = 2048)
    buf = Vector{UInt8}(undef, buffer_size)
    n = ccall((:getErrorMessage_C, _libgram()), Csize_t, (Ptr{UInt8}, Csize_t), buf, Csize_t(buffer_size))
    return String(buf[1:Int(n)])
end
