_define_outer_planet_bindings(:titan, :T, BODY_TITAN, "Titan")

function titan_get_gases_state(atmos::Atmosphere)
    _require_body!(atmos, BODY_TITAN, "Titan")
    state = Ref(GasesStateC(0, 0, 0, 0, 0))
    argon = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    methane = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    dinitrogen = Ref(ConstituentGasC(0, 0, 0, 0, 0))
    ccall((:getGasesState_T, _libgram()), Cvoid, (Ptr{Cvoid}, Ref{GasesStateC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}, Ref{ConstituentGasC}), atmos.ptr, state, argon, methane, dinitrogen)
    return (state = state[], argon = argon[], methane = methane[], dinitrogen = dinitrogen[])
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
