function mars_close!(atmos::Atmosphere)
    _require_mars!(atmos)
    if atmos.ptr != C_NULL
        ccall((:deleteAtmosphere_M, _libgram()), Cvoid, (Ptr{Cvoid},), atmos.ptr)
        atmos.ptr = C_NULL
    end
    return nothing
end

function mars_set_start_time!(atmos::Atmosphere; year::Integer, month::Integer, day::Integer,
    hour::Integer = 0, minute::Integer = 0, seconds::Real = 0.0, scale::Integer = 1, frame::Integer = 1)
    _require_mars!(atmos)
    t = Ref(GramTimeC(Cint(year), Cint(month), Cint(day), Cint(hour), Cint(minute), Cdouble(seconds), Cint(scale), Cint(frame)))
    ccall((:setStartTime_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GramTimeC}), atmos.ptr, t)
    return nothing
end

function mars_set_start_time!(atmos::Atmosphere, t::GramTimeC)
    _require_mars!(atmos)
    ccall((:setStartTime_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GramTimeC}), atmos.ptr, Ref(t))
    return nothing
end

function mars_set_position!(atmos::Atmosphere; height::Real, latitude::Real, longitude::Real,
    elapsed_time::Real = 0.0, is_planetocentric::Bool = true)
    _require_mars!(atmos)
    p = Ref(PositionC(
        Cdouble(height), Cdouble(latitude), Cdouble(longitude), Cdouble(elapsed_time),
        is_planetocentric ? Cint(1) : Cint(0), Cdouble(0), Cdouble(0), Cdouble(0), Cdouble(0)))
    ccall((:setPosition_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PositionC}), atmos.ptr, p)
    return nothing
end

function mars_set_position!(atmos::Atmosphere, p::PositionC)
    _require_mars!(atmos)
    ccall((:setPosition_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PositionC}), atmos.ptr, Ref(p))
    return nothing
end

function mars_set_delta!(atmos::Atmosphere; height::Real = 0.0, latitude::Real = 0.0, longitude::Real = 0.0,
    elapsed_time::Real = 0.0, is_planetocentric::Bool = true)
    _require_mars!(atmos)
    p = Ref(PositionC(
        Cdouble(height), Cdouble(latitude), Cdouble(longitude), Cdouble(elapsed_time),
        is_planetocentric ? Cint(1) : Cint(0), Cdouble(0), Cdouble(0), Cdouble(0), Cdouble(0)))
    ccall((:setDelta_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PositionC}), atmos.ptr, p)
    return nothing
end

function mars_set_delta!(atmos::Atmosphere, p::PositionC)
    _require_mars!(atmos)
    ccall((:setDelta_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PositionC}), atmos.ptr, Ref(p))
    return nothing
end

function mars_set_seed!(atmos::Atmosphere, seed::Integer)
    _require_mars!(atmos)
    ccall((:setSeed_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, Cint(seed))
    return nothing
end

function mars_set_min_relative_step_size!(atmos::Atmosphere, min_relative_step_size::Real)
    _require_mars!(atmos)
    ccall((:setMinRelativeStepSize_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble), atmos.ptr, Cdouble(min_relative_step_size))
    return nothing
end

function mars_set_perturbation_scales!(
    atmos::Atmosphere;
    density::Real = 1.0,
    ew_wind::Real = 1.0,
    ns_wind::Real = 1.0,
    vertical_wind::Real = 1.0
)
    _require_mars!(atmos)
    ccall(
        (:setPerturbationScales_M, _libgram()),
        Cvoid,
        (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble),
        atmos.ptr,
        Cdouble(density),
        Cdouble(ew_wind),
        Cdouble(ns_wind),
        Cdouble(vertical_wind)
    )
    return nothing
end

function mars_add_auxiliary_atmosphere!(
    atmos::Atmosphere;
    file_name::AbstractString,
    inner_radius::Real,
    outer_radius::Real,
    is_east_longitude_positive::Bool = true
)
    _require_mars!(atmos)
    ccall(
        (:addAuxiliaryAtmosphere_M, _libgram()),
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

function mars_set_auxiliary_values!(
    atmos::Atmosphere;
    density::Real,
    pressure::Real,
    temperature::Real,
    ew_wind::Real,
    ns_wind::Real
)
    _require_mars!(atmos)
    ccall(
        (:setAuxiliaryValues_M, _libgram()),
        Cvoid,
        (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble, Cdouble, Cdouble),
        atmos.ptr,
        Cdouble(density),
        Cdouble(pressure),
        Cdouble(temperature),
        Cdouble(ew_wind),
        Cdouble(ns_wind)
    )
    return nothing
end

function mars_set_perturbation_action!(atmos::Atmosphere, update::Bool)
    _require_mars!(atmos)
    ccall((:setPerturbationAction_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, update ? Cint(1) : Cint(0))
    return nothing
end

function mars_set_ephemeris_state!(atmos::Atmosphere, state::EphemerisStateC)
    _require_mars!(atmos)
    ccall((:setEphemerisState_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{EphemerisStateC}), atmos.ptr, Ref(state))
    return nothing
end

function mars_set_ephemeris_fast_mode!(atmos::Atmosphere, flag::Bool)
    _require_mars!(atmos)
    ccall((:setEphemerisFastModeOn_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, flag ? Cint(1) : Cint(0))
    return nothing
end

function mars_set_subsolar_update_time!(atmos::Atmosphere, utime::Real)
    _require_mars!(atmos)
    ccall((:setSubsolarUpdateTime_M, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble), atmos.ptr, Cdouble(utime))
    return nothing
end

function mars_get_start_time(atmos::Atmosphere)
    _require_mars!(atmos)
    t = Ref(GramTimeC(0, 0, 0, 0, 0, 0.0, 1, 1))
    ccall((:getStartTime_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GramTimeC}), atmos.ptr, t)
    return t[]
end

function mars_get_position(atmos::Atmosphere)
    _require_mars!(atmos)
    p = Ref(PositionC(0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getPosition_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PositionC}), atmos.ptr, p)
    return p[]
end

function mars_get_dynamics_state(atmos::Atmosphere)
    _require_mars!(atmos)
    s = Ref(DynamicsStateC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getDynamicsState_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{DynamicsStateC}), atmos.ptr, s)
    return s[]
end

function mars_get_density_state(atmos::Atmosphere)
    _require_mars!(atmos)
    s = Ref(DensityStateC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getDensityState_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{DensityStateC}), atmos.ptr, s)
    return s[]
end

function mars_get_winds_state(atmos::Atmosphere)
    _require_mars!(atmos)
    w = Ref(WindsStateC(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getWindsState_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{WindsStateC}), atmos.ptr, w)
    return w[]
end

function mars_get_ephemeris_state(atmos::Atmosphere)
    _require_mars!(atmos)
    s = Ref(EphemerisStateC(0, 0, 0, 0, 0, 0, 0, 0))
    ccall((:getEphemerisState_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{EphemerisStateC}), atmos.ptr, s)
    return s[]
end

function mars_get_perturbation_state(atmos::Atmosphere)
    _require_mars!(atmos)
    s = Ref(PerturbationStateC(0, 0, 0, 0))
    ccall((:getPerturbationState_M, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{PerturbationStateC}), atmos.ptr, s)
    return s[]
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
