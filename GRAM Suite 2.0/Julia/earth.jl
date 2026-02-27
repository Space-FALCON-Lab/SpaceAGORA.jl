@inline function _cstring_until_nul(buf::Vector{UInt8})
    i = findfirst(==(0x00), buf)
    if i === nothing
        return String(buf)
    end
    i == 1 && return ""
    return String(buf[1:(i - 1)])
end

function earth_try_get_spice_path(buffer_size::Integer = 2048)
    n = Int(max(buffer_size, 2))
    buf = Vector{UInt8}(undef, n)
    ccall((:tryGetSpicePath_E, _libgram()), Cvoid, (Ptr{UInt8}, Cint), buf, Cint(n))
    return _cstring_until_nul(buf)
end

function earth_try_get_data_paths(buffer_size::Integer = 2048)
    n = Int(max(buffer_size, 2))
    spice_buf = Vector{UInt8}(undef, n)
    data_buf = Vector{UInt8}(undef, n)
    ccall((:tryGetDataPaths_E, _libgram()), Cvoid, (Ptr{UInt8}, Ptr{UInt8}, Cint), spice_buf, data_buf, Cint(n))
    return _cstring_until_nul(spice_buf), _cstring_until_nul(data_buf)
end

function earth_set_spice_lsk!(lsk::AbstractString)
    ccall((:setSpiceLsk_E, _libgram()), Cvoid, (Cstring,), lsk)
    return nothing
end

function earth_set_spice_pck!(pck::AbstractString)
    ccall((:setSpicePck_E, _libgram()), Cvoid, (Cstring,), pck)
    return nothing
end

function earth_set_spice_kernel!(bsp::AbstractString)
    ccall((:setSpiceKernel_E, _libgram()), Cvoid, (Cstring,), bsp)
    return nothing
end

function earth_initialize!(spice_path::AbstractString)
    _SPICE_PATH[] = String(spice_path)
    ccall((:initialize_E, _libgram()), Cvoid, (Cstring,), spice_path)
    push!(_REGISTERED_BODIES, BODY_EARTH)
    return nothing
end

function earth_load_spice_file!(file_name::AbstractString)
    ccall((:loadSpiceFile_E, _libgram()), Cvoid, (Cstring,), file_name)
    return nothing
end

function create_earth_atmosphere(; data_path::Union{Nothing, AbstractString} = nothing)
    return create_atmosphere(BODY_EARTH; data_path = data_path)
end

function copy_earth_atmosphere(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return copy_atmosphere(atmos)
end

function earth_close!(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    if atmos.ptr != C_NULL
        ccall((:deleteAtmosphere_E, _libgram()), Cvoid, (Ptr{Cvoid},), atmos.ptr)
        atmos.ptr = C_NULL
    end
    return nothing
end

function earth_set_thermosphere_model!(atmos::Atmosphere, model::Integer)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setThermosphereModel_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, Cint(model))
    return nothing
end

function earth_set_surface_roughness!(atmos::Atmosphere, z0::Real)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setSurfaceRoughness_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble), atmos.ptr, Cdouble(z0))
    return nothing
end

function earth_set_atmos_path!(atmos::Atmosphere, path::AbstractString)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setAtmosPath_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cstring), atmos.ptr, path)
    return nothing
end

function earth_set_merra2_path!(atmos::Atmosphere, path::AbstractString)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setMERRA2Path_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cstring), atmos.ptr, path)
    return nothing
end

function earth_set_ncep_path!(atmos::Atmosphere, path::AbstractString)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setNCEPPath_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cstring), atmos.ptr, path)
    return nothing
end

function earth_set_merra2_parameters!(atmos::Atmosphere; hour::Integer, lat_min::Real, lat_max::Real, lon_min::Real, lon_max::Real)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setMERRA2Parameters_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cint, Cdouble, Cdouble, Cdouble, Cdouble), atmos.ptr, Cint(hour), Cdouble(lat_min), Cdouble(lat_max), Cdouble(lon_min), Cdouble(lon_max))
    return nothing
end

function earth_set_use_ncep!(atmos::Atmosphere, use_ncep::Bool)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setUseNCEP_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, use_ncep ? Cint(1) : Cint(0))
    return nothing
end

function earth_set_ncep_parameters!(atmos::Atmosphere; year::Integer, hour::Integer)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setNCEPParameters_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cint, Cint), atmos.ptr, Cint(year), Cint(hour))
    return nothing
end

function earth_set_rra_path!(atmos::Atmosphere, path::AbstractString)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setRRAPath_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cstring), atmos.ptr, path)
    return nothing
end

function earth_set_rra_site_list!(atmos::Atmosphere, file_name::AbstractString)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setRRASiteList_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cstring), atmos.ptr, file_name)
    return nothing
end

function earth_set_rra_parameters!(atmos::Atmosphere; year::Integer, inner_radius::Real, outer_radius::Real)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setRRAParameters_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cint, Cdouble, Cdouble), atmos.ptr, Cint(year), Cdouble(inner_radius), Cdouble(outer_radius))
    return nothing
end

function earth_set_use_rra!(atmos::Atmosphere, enabled::Bool)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setUseRRA_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, enabled ? Cint(1) : Cint(0))
    return nothing
end

function earth_set_solar_parameters!(atmos::Atmosphere; daily_f10::Real, mean_f10::Real, ap::Real)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setSolarParameters_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble), atmos.ptr, Cdouble(daily_f10), Cdouble(mean_f10), Cdouble(ap))
    return nothing
end

function earth_set_jb2008_parameters!(atmos::Atmosphere; daily_s10::Real, mean_s10::Real, daily_xm10::Real, mean_xm10::Real, daily_y10::Real, mean_y10::Real, dstdtc::Real)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setJB2008Parameters_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), atmos.ptr, Cdouble(daily_s10), Cdouble(mean_s10), Cdouble(daily_xm10), Cdouble(mean_xm10), Cdouble(daily_y10), Cdouble(mean_y10), Cdouble(dstdtc))
    return nothing
end

function earth_set_initial_perturbations!(atmos::Atmosphere; density::Real, temperature::Real, ew_wind::Real, ns_wind::Real, vertical_wind::Real)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setInitialPerturbations_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble), atmos.ptr, Cdouble(density), Cdouble(temperature), Cdouble(ew_wind), Cdouble(ns_wind), Cdouble(vertical_wind))
    return nothing
end

function earth_set_patchiness!(atmos::Atmosphere, enabled::Bool)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setPatchiness_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, enabled ? Cint(1) : Cint(0))
    return nothing
end

function earth_set_start_time!(atmos::Atmosphere; kwargs...)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_start_time!(atmos; kwargs...)
end

function earth_set_start_time!(atmos::Atmosphere, t::GramTimeC)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_start_time!(atmos, t)
end

function earth_set_position!(atmos::Atmosphere; kwargs...)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_position!(atmos; kwargs...)
end

function earth_set_position!(atmos::Atmosphere, p::PositionC)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_position!(atmos, p)
end

function earth_set_delta!(atmos::Atmosphere; kwargs...)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_delta!(atmos; kwargs...)
end

function earth_set_delta!(atmos::Atmosphere, p::PositionC)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_delta!(atmos, p)
end

function earth_set_seed!(atmos::Atmosphere, seed::Integer)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_seed!(atmos, seed)
end

function earth_set_perturbation_scales!(atmos::Atmosphere; random_scale::Real = 1.0, horizontal_wind_scale::Real = 1.0, vertical_wind_scale::Real = 1.0)
    _require_body!(atmos, BODY_EARTH, "Earth")
    ccall((:setPerturbationScales_E, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble), atmos.ptr, Cdouble(random_scale), Cdouble(horizontal_wind_scale), Cdouble(vertical_wind_scale))
    return nothing
end

function earth_add_auxiliary_atmosphere!(atmos::Atmosphere; kwargs...)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return add_auxiliary_atmosphere!(atmos; kwargs...)
end

function earth_set_auxiliary_values!(atmos::Atmosphere; kwargs...)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_auxiliary_values!(atmos; kwargs...)
end

function earth_set_perturbation_action!(atmos::Atmosphere, update::Bool)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_perturbation_action!(atmos, update)
end

function earth_set_ephemeris_state!(atmos::Atmosphere, state::EphemerisStateC)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_ephemeris_state!(atmos, state)
end

function earth_set_ephemeris_fast_mode!(atmos::Atmosphere, flag::Bool)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_ephemeris_fast_mode!(atmos, flag)
end

function earth_set_subsolar_update_time!(atmos::Atmosphere, utime::Real)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return set_subsolar_update_time!(atmos, utime)
end

function earth_update!(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return Int(ccall((:update_E, _libgram()), Cint, (Ptr{Cvoid},), atmos.ptr))
end

function earth_get_error_message(buffer_size::Integer = 2048)
    buf = Vector{UInt8}(undef, buffer_size)
    n = ccall((:getErrorMessage_E, _libgram()), Csize_t, (Ptr{UInt8}, Csize_t), buf, Csize_t(buffer_size))
    return String(buf[1:Int(n)])
end

function earth_get_start_time(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return get_start_time(atmos)
end

function earth_get_version_string(atmos::Atmosphere; buffer_size::Integer = 2048)
    _require_body!(atmos, BODY_EARTH, "Earth")
    n = Int(max(buffer_size, 2))
    buf = Vector{UInt8}(undef, n)
    out_len = ccall((:getVersionString_E, _libgram()), Csize_t, (Ptr{Cvoid}, Ptr{UInt8}, Csize_t), atmos.ptr, buf, Csize_t(n))
    out_len_int = Int(out_len)
    out_len_int <= 0 && return ""
    return String(buf[1:min(out_len_int, n)])
end

function earth_get_position(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return get_position(atmos)
end

function earth_get_dynamics_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return get_dynamics_state(atmos)
end

function earth_get_density_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return get_density_state(atmos)
end

function earth_get_winds_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return get_winds_state(atmos)
end

function earth_get_ephemeris_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return get_ephemeris_state(atmos)
end

function earth_get_perturbation_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    return get_perturbation_state(atmos)
end

function earth_get_gases_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    state = Ref(GasesStateC(0, 0, 0, 0, 0))
    argon = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    carbon_dioxide = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    carbon_monoxide = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    dinitrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    dioxygen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    helium = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    hydrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    methane = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    nitrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    nitrous_oxide = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    oxygen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    ozone = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    water = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    ccall(
        (:getGasesState_E, _libgram()),
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
        dinitrogen,
        dioxygen,
        helium,
        hydrogen,
        methane,
        nitrogen,
        nitrous_oxide,
        oxygen,
        ozone,
        water
    )
    return (
        state = state[],
        argon = argon[],
        carbon_dioxide = carbon_dioxide[],
        carbon_monoxide = carbon_monoxide[],
        dinitrogen = dinitrogen[],
        dioxygen = dioxygen[],
        helium = helium[],
        hydrogen = hydrogen[],
        methane = methane[],
        nitrogen = nitrogen[],
        nitrous_oxide = nitrous_oxide[],
        oxygen = oxygen[],
        ozone = ozone[],
        water = water[]
    )
end

function earth_get_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    s = Ref(EarthStateC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, ntuple(_ -> Cchar(0), 6), 0, 0, 0, 0))
    ccall((:getEarthState_E, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{EarthStateC}), atmos.ptr, s)
    return s[]
end

function earth_get_perts(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    s = Ref(EarthPertsC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getEarthPerts_E, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{EarthPertsC}), atmos.ptr, s)
    return s[]
end

function earth_get_surface(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    s = Ref(EarthSurfaceC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getEarthSurface_E, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{EarthSurfaceC}), atmos.ptr, s)
    return s[]
end

function earth_get_boundary_layer(atmos::Atmosphere)
    _require_body!(atmos, BODY_EARTH, "Earth")
    s = Ref(EarthBoundaryLayerC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getEarthBoundaryLayer_E, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{EarthBoundaryLayerC}), atmos.ptr, s)
    return s[]
end
