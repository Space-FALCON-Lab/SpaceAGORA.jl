for (planet, suffix, body, label) in (
    (:venus, :V, BODY_VENUS, "Venus"),
    (:jupiter, :J, BODY_JUPITER, "Jupiter"),
    (:uranus, :U, BODY_URANUS, "Uranus"),
    (:neptune, :N, BODY_NEPTUNE, "Neptune"),
    (:titan, :T, BODY_TITAN, "Titan")
)
    pstr = String(planet)
    try_sym = Symbol("tryGetSpicePath_", suffix)
    set_lsk_sym = Symbol("setSpiceLsk_", suffix)
    set_pck_sym = Symbol("setSpicePck_", suffix)
    set_kernel_sym = Symbol("setSpiceKernel_", suffix)
    init_sym = Symbol("initialize_", suffix)
    load_spice_sym = Symbol("loadSpiceFile_", suffix)
    del_sym = Symbol("deleteAtmosphere_", suffix)
    set_pert_sym = Symbol("setPerturbationScales_", suffix)
    update_sym = Symbol("update_", suffix)
    err_sym = Symbol("getErrorMessage_", suffix)
    version_sym = Symbol("getVersionString_", suffix)

    @eval begin
        function $(Symbol(pstr * "_try_get_spice_path"))(buffer_size::Integer = 2048)
            n = Int(max(buffer_size, 2))
            buf = Vector{UInt8}(undef, n)
            ccall(($(QuoteNode(try_sym)), _libgram()), Cvoid, (Ptr{UInt8}, Cint), buf, Cint(n))
            return _cstring_until_nul(buf)
        end

        function $(Symbol(pstr * "_set_spice_lsk!"))(lsk::AbstractString)
            ccall(($(QuoteNode(set_lsk_sym)), _libgram()), Cvoid, (Cstring,), lsk)
            return nothing
        end

        function $(Symbol(pstr * "_set_spice_pck!"))(pck::AbstractString)
            ccall(($(QuoteNode(set_pck_sym)), _libgram()), Cvoid, (Cstring,), pck)
            return nothing
        end

        function $(Symbol(pstr * "_set_spice_kernel!"))(bsp::AbstractString)
            ccall(($(QuoteNode(set_kernel_sym)), _libgram()), Cvoid, (Cstring,), bsp)
            return nothing
        end

        function $(Symbol(pstr * "_initialize!"))(spice_path::AbstractString)
            _SPICE_PATH[] = String(spice_path)
            ccall(($(QuoteNode(init_sym)), _libgram()), Cvoid, (Cstring,), spice_path)
            push!(_REGISTERED_BODIES, $body)
            return nothing
        end

        function $(Symbol(pstr * "_load_spice_file!"))(file_name::AbstractString)
            ccall(($(QuoteNode(load_spice_sym)), _libgram()), Cvoid, (Cstring,), file_name)
            return nothing
        end

        function $(Symbol("create_" * pstr * "_atmosphere"))()
            return create_atmosphere($body)
        end

        function $(Symbol("copy_" * pstr * "_atmosphere"))(atmos::Atmosphere)
            _require_body!(atmos, $body, $label)
            return copy_atmosphere(atmos)
        end

        function $(Symbol(pstr * "_close!"))(atmos::Atmosphere)
            _require_body!(atmos, $body, $label)
            if atmos.ptr != C_NULL
                ccall(($(QuoteNode(del_sym)), _libgram()), Cvoid, (Ptr{Cvoid},), atmos.ptr)
                atmos.ptr = C_NULL
            end
            return nothing
        end

        function $(Symbol(pstr * "_set_start_time!"))(atmos::Atmosphere; kwargs...)
            _require_body!(atmos, $body, $label)
            return set_start_time!(atmos; kwargs...)
        end

        function $(Symbol(pstr * "_set_start_time!"))(atmos::Atmosphere, t::GramTimeC)
            _require_body!(atmos, $body, $label)
            return set_start_time!(atmos, t)
        end

        function $(Symbol(pstr * "_set_position!"))(atmos::Atmosphere; kwargs...)
            _require_body!(atmos, $body, $label)
            return set_position!(atmos; kwargs...)
        end

        function $(Symbol(pstr * "_set_position!"))(atmos::Atmosphere, p::PositionC)
            _require_body!(atmos, $body, $label)
            return set_position!(atmos, p)
        end

        function $(Symbol(pstr * "_set_delta!"))(atmos::Atmosphere; kwargs...)
            _require_body!(atmos, $body, $label)
            return set_delta!(atmos; kwargs...)
        end

        function $(Symbol(pstr * "_set_delta!"))(atmos::Atmosphere, p::PositionC)
            _require_body!(atmos, $body, $label)
            return set_delta!(atmos, p)
        end

        function $(Symbol(pstr * "_set_seed!"))(atmos::Atmosphere, seed::Integer)
            _require_body!(atmos, $body, $label)
            return set_seed!(atmos, seed)
        end

        function $(Symbol(pstr * "_set_min_relative_step_size!"))(atmos::Atmosphere, min_relative_step_size::Real)
            _require_body!(atmos, $body, $label)
            return set_min_relative_step_size!(atmos, min_relative_step_size)
        end

        function $(Symbol(pstr * "_set_perturbation_scales!"))(atmos::Atmosphere; density::Real = 1.0, ew_wind::Real = 1.0, ns_wind::Real = 1.0)
            _require_body!(atmos, $body, $label)
            ccall(($(QuoteNode(set_pert_sym)), _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble, Cdouble, Cdouble), atmos.ptr, Cdouble(density), Cdouble(ew_wind), Cdouble(ns_wind))
            return nothing
        end

        function $(Symbol(pstr * "_add_auxiliary_atmosphere!"))(atmos::Atmosphere; kwargs...)
            _require_body!(atmos, $body, $label)
            return add_auxiliary_atmosphere!(atmos; kwargs...)
        end

        function $(Symbol(pstr * "_set_auxiliary_values!"))(atmos::Atmosphere; kwargs...)
            _require_body!(atmos, $body, $label)
            return set_auxiliary_values!(atmos; kwargs...)
        end

        function $(Symbol(pstr * "_set_perturbation_action!"))(atmos::Atmosphere, update::Bool)
            _require_body!(atmos, $body, $label)
            return set_perturbation_action!(atmos, update)
        end

        function $(Symbol(pstr * "_set_ephemeris_state!"))(atmos::Atmosphere, state::EphemerisStateC)
            _require_body!(atmos, $body, $label)
            return set_ephemeris_state!(atmos, state)
        end

        function $(Symbol(pstr * "_set_ephemeris_fast_mode!"))(atmos::Atmosphere, flag::Bool)
            _require_body!(atmos, $body, $label)
            return set_ephemeris_fast_mode!(atmos, flag)
        end

        function $(Symbol(pstr * "_set_subsolar_update_time!"))(atmos::Atmosphere, utime::Real)
            _require_body!(atmos, $body, $label)
            return set_subsolar_update_time!(atmos, utime)
        end

        function $(Symbol(pstr * "_update!"))(atmos::Atmosphere)
            _require_body!(atmos, $body, $label)
            return Int(ccall(($(QuoteNode(update_sym)), _libgram()), Cint, (Ptr{Cvoid},), atmos.ptr))
        end

        function $(Symbol(pstr * "_get_error_message"))(buffer_size::Integer = 2048)
            buf = Vector{UInt8}(undef, buffer_size)
            n = ccall(($(QuoteNode(err_sym)), _libgram()), Csize_t, (Ptr{UInt8}, Csize_t), buf, Csize_t(buffer_size))
            return String(buf[1:Int(n)])
        end

        function $(Symbol(pstr * "_get_start_time"))(atmos::Atmosphere)
            _require_body!(atmos, $body, $label)
            return get_start_time(atmos)
        end

        function $(Symbol(pstr * "_get_version_string"))(atmos::Atmosphere; buffer_size::Integer = 2048)
            _require_body!(atmos, $body, $label)
            n = Int(max(buffer_size, 2))
            buf = Vector{UInt8}(undef, n)
            out_len = ccall(($(QuoteNode(version_sym)), _libgram()), Csize_t, (Ptr{Cvoid}, Ptr{UInt8}, Csize_t), atmos.ptr, buf, Csize_t(n))
            out_len_int = Int(out_len)
            out_len_int <= 0 && return ""
            return String(buf[1:min(out_len_int, n)])
        end

        function $(Symbol(pstr * "_get_position"))(atmos::Atmosphere)
            _require_body!(atmos, $body, $label)
            return get_position(atmos)
        end

        function $(Symbol(pstr * "_get_dynamics_state"))(atmos::Atmosphere)
            _require_body!(atmos, $body, $label)
            return get_dynamics_state(atmos)
        end

        function $(Symbol(pstr * "_get_density_state"))(atmos::Atmosphere)
            _require_body!(atmos, $body, $label)
            return get_density_state(atmos)
        end

        function $(Symbol(pstr * "_get_winds_state"))(atmos::Atmosphere)
            _require_body!(atmos, $body, $label)
            return get_winds_state(atmos)
        end

        function $(Symbol(pstr * "_get_ephemeris_state"))(atmos::Atmosphere)
            _require_body!(atmos, $body, $label)
            return get_ephemeris_state(atmos)
        end

        function $(Symbol(pstr * "_get_perturbation_state"))(atmos::Atmosphere)
            _require_body!(atmos, $body, $label)
            return get_perturbation_state(atmos)
        end
    end
end

function venus_get_gases_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_VENUS, "Venus")
    state = Ref(GasesStateC(0, 0, 0, 0, 0))
    hydrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    helium = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    oxygen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    nitrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    dinitrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    carbon_monoxide = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    carbon_dioxide = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    ccall(
        (:getGasesState_V, _libgram()),
        Cvoid,
        (Ptr{Cvoid}, Ref{GasesStateC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}),
        atmos.ptr, state, hydrogen, helium, oxygen, nitrogen, dinitrogen, carbon_monoxide, carbon_dioxide
    )
    return (
        state = state[],
        hydrogen = hydrogen[],
        helium = helium[],
        oxygen = oxygen[],
        nitrogen = nitrogen[],
        dinitrogen = dinitrogen[],
        carbon_monoxide = carbon_monoxide[],
        carbon_dioxide = carbon_dioxide[]
    )
end

function jupiter_get_gases_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_JUPITER, "Jupiter")
    return (state = get_gases_state(atmos),)
end

function uranus_get_gases_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_URANUS, "Uranus")
    state = Ref(GasesStateC(0, 0, 0, 0, 0))
    dihydrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    methane = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    helium = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    ccall((:getGasesState_U, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GasesStateC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}), atmos.ptr, state, dihydrogen, methane, helium)
    return (state = state[], dihydrogen = dihydrogen[], methane = methane[], helium = helium[])
end

function neptune_get_gases_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_NEPTUNE, "Neptune")
    state = Ref(GasesStateC(0, 0, 0, 0, 0))
    dihydrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    methane = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    helium = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    dinitrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    ccall((:getGasesState_N, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GasesStateC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}), atmos.ptr, state, dihydrogen, methane, helium, dinitrogen)
    return (state = state[], dihydrogen = dihydrogen[], methane = methane[], helium = helium[], dinitrogen = dinitrogen[])
end

function titan_get_gases_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_TITAN, "Titan")
    state = Ref(GasesStateC(0, 0, 0, 0, 0))
    argon = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    methane = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    dinitrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    ccall((:getGasesState_T, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GasesStateC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}), atmos.ptr, state, argon, methane, dinitrogen)
    return (state = state[], argon = argon[], methane = methane[], dinitrogen = dinitrogen[])
end

function neptune_set_dinitrogen_mole_fraction!(atmos::Atmosphere, n2mf::Real)
    _require_body!(atmos, BODY_NEPTUNE, "Neptune")
    ccall((:setDinitrogenMoleFraction_N, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble), atmos.ptr, Cdouble(n2mf))
    return nothing
end

function neptune_set_min_max_factor!(atmos::Atmosphere; factor::Real, compute::Bool = true)
    _require_body!(atmos, BODY_NEPTUNE, "Neptune")
    ccall((:setMinMaxFactor_N, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble, Cint), atmos.ptr, Cdouble(factor), compute ? Cint(1) : Cint(0))
    return nothing
end

function neptune_get_min_max_factor(atmos::Atmosphere)
    _require_body!(atmos, BODY_NEPTUNE, "Neptune")
    f = Ref{Cdouble}(0.0)
    ccall((:getMinMaxFactor_N, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{Cdouble}), atmos.ptr, f)
    return f[]
end

function titan_set_min_max_factor!(atmos::Atmosphere; factor::Real, compute::Bool = true)
    _require_body!(atmos, BODY_TITAN, "Titan")
    ccall((:setMinMaxFactor_T, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble, Cint), atmos.ptr, Cdouble(factor), compute ? Cint(1) : Cint(0))
    return nothing
end

function titan_get_min_max_factor(atmos::Atmosphere)
    _require_body!(atmos, BODY_TITAN, "Titan")
    f = Ref{Cdouble}(0.0)
    ccall((:getMinMaxFactor_T, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{Cdouble}), atmos.ptr, f)
    return f[]
end

function titan_set_methane_mole_fraction!(atmos::Atmosphere, mmf::Real)
    _require_body!(atmos, BODY_TITAN, "Titan")
    ccall((:setMethaneMoleFraction_T, _libgram()), Cvoid, (Ptr{Cvoid}, Cdouble), atmos.ptr, Cdouble(mmf))
    return nothing
end

function titan_set_model_type!(atmos::Atmosphere, model_type::Integer)
    _require_body!(atmos, BODY_TITAN, "Titan")
    ccall((:setModelType_T, _libgram()), Cvoid, (Ptr{Cvoid}, Cint), atmos.ptr, Cint(model_type))
    return nothing
end
