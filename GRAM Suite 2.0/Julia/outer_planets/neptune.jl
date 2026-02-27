_define_outer_planet_bindings(:neptune, :N, BODY_NEPTUNE, "Neptune")

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
